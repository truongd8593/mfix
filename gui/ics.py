# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import QLabel, QLineEdit, QPushButton, QGridLayout
from qtpy.QtGui import QPicture

UserRole = QtCore.Qt.UserRole

from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit

from tools.general import (set_item_noedit, set_item_enabled,
                           get_selected_row,
                           widget_iter)

# We don't need extended JSON here
from json import JSONDecoder, JSONEncoder

def safe_float(val):
    try:
        return float(val)
    except ValueError:
        return 0.0

FLUID_TAB = 0
SOLIDS_TAB_DUMMY_L = 1
SOLIDS_TAB = 2
SOLIDS_TAB_DUMMY_R = 3
SCALAR_TAB = 4

class ICS(object):
    # Initial Conditions Task Pane Window: This section allows a user to define the initial conditions
    # for the described model. This section relies on regions named in the Regions section.

    def init_ics(self):
        ui = self.ui
        ics = ui.initial_conditions

        self.ics = {} # key: index.  value: data dictionary for initial cond
        self.ics_current_indices = [] # List of IC indices
        self.ics_current_regions = [] # And the names of the regions which define them
        self.ics_region_dict = None

        # The top of the task pane is where users define/select IC regions
        # (see handle_ics_region_selection)
        #
        #Icons to add/remove/duplicate regions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #the region to apply the initial condition.
        #Implementation Idea: Allow the user to select multiple regions in the popup window so that they can
        #define one IC region over multiple regions. Managing the MFIX IC indices on the back end could
        #get messy, and if so, scrap this idea.
        # (Note- Suggestion implemented)
        # TODO: implement 'duplicate' (what does this do?)
        ics.toolbutton_add.clicked.connect(self.ics_show_regions_popup)
        ics.toolbutton_delete.clicked.connect(self.ics_delete_regions)
        ics.toolbutton_delete.setEnabled(False) # Need a selection

        ics.tablewidget_regions.itemSelectionChanged.connect(self.handle_ics_region_selection)

        self.ics_current_tab = 0 # #? "Fluid" tab.  If fluid is disabled, we will switch
        self.ics_current_solid = None
        ics.pushbutton_fluid.clicked.connect(lambda: self.ics_change_tab(FLUID_TAB,0))
        ics.pushbutton_scalar.clicked.connect(lambda: self.ics_change_tab(SCALAR_TAB,0))

        # Trim width of "Fluid" and "Scalar" buttons, like we do for
        # dynamically-created "Solid #" buttons
        for b in (ics.pushbutton_fluid, ics.pushbutton_scalar):
            w = b.fontMetrics().boundingRect(b.text()).width() + 20
            b.setMaximumWidth(w)


    def ics_set_volume_fraction_limit(self):
        # Set ic_ep_g from ic_ep_s (issues/121)
        if not self.ics_current_indices:
            return
        if not self.ics_current_solid:
            return
        IC0 = self.ics_current_indices[0]
        P = self.ics_current_solid
        key = 'ic_ep_s'
        widget = self.ui.initial_conditions.lineedit_keyword_ic_ep_s_args_IC_P

        s = sum(safe_float(self.project.get_value(key, default=0, args=[IC0, s]))
                for s in range(1, len(self.solids)+1) if s != P)

        lim = max(0, 1.0 - s)
        lim = round(lim, 10) # avoid problem with 1 - 0.9 != 0.1

        widget.min = 0.0
        widget.max = lim

    def handle_ics_volume_fraction(self, widget, val, args):
        # We may have been called before submit_change
        if widget is None:
            return
        self.project.submit_change(widget, val, args)
        if not self.ics_current_indices:
            return
        IC0 = self.ics_current_indices[0]
        if not self.ics_current_solid:
            return
        P = self.ics_current_solid
        s = sum(safe_float(self.project.get_value('ic_ep_s', default=0, args=[IC0, s]))
                for s in range(1, len(self.solids)+1))
        if s > 1.0:
            self.warning("Volume fractions sum to %s, must be <= 1.0" % s,
                         popup=True)
            return # ?
        val = round(1.0 - s, 10)
        for IC in self.ics_current_indices:
            self.update_keyword('ic_ep_g', val, args=[IC])


    def ics_show_regions_popup(self):
        #Users cannot select inapplicable regions.
        #IC regions must be volumes (not planes, points, or STLs)
        #No region can define more than one initial condition.
        ui = self.ui
        ics = ui.initial_conditions
        rp = self.regions_popup
        rp.clear()
        for (name,data) in self.ics_region_dict.items():
            shape = data.get('type', '---')
            # Assume available if unmarked
            available = data.get('available', True) and (shape == 'box')
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.ics_add_regions)
        rp.cancel.connect(self.ics_cancel_add)
        for item in (ics.tablewidget_regions,
                     ics.scrollarea_detail,
                     ics.toolbutton_add,
                     ics.toolbutton_delete):
            item.setEnabled(False)
        rp.popup(boundary=False)


    def ics_cancel_add(self):
        ui = self.ui
        ics = ui.initial_conditions

        for item in (ics.toolbutton_add,
                     ics.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ics.tablewidget_regions) is not None:
            for item in (ics.scrollarea_detail,
                         ics.toolbutton_delete):
                item.setEnabled(True)


    def ics_add_regions(self):
        # Interactively add regions to define ICs
        ui = self.ui
        ics = ui.initial_conditions
        rp = self.regions_popup
        self.ics_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        if not selections:
            return
        self.ics_add_regions_1(selections) # Indices will be assigned
        self.ics_setup_current_tab() # Update the widgets

    def ics_add_regions_1(self, selections, indices=None):
        # Used by both interactive and load-time add-region handlers
        ui = self.ui
        if self.ics_region_dict is None:
            self.ics_region_dict = ui.regions.get_region_dict()

        ics = ui.initial_conditions
        tw = ics.tablewidget_regions
        nrows = tw.rowCount()
        tw.setRowCount(nrows+1)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        item = make_item('+'.join(selections))

        if indices is None:
            indices = [None] * len(selections)
        else:
            assert len(selections) == len(indices)

        for (i, region_name) in enumerate(selections):
            idx = indices[i]
            if idx is None:
                idx = self.ics_find_index()
                indices[i] = idx
            self.ics[idx] = {'region': region_name}
            region_data = self.ics_region_dict.get(region_name)
            if region_data is None: # ?
                self.warn("no data for region %s" % region_name)
                continue
            self.ics_set_region_keys(region_name, idx, region_data)

            self.ics_region_dict[region_name]['available'] = False # Mark as in-use

        item.setData(UserRole, (tuple(indices), tuple(selections)))

        self.ics_current_regions = selections
        self.ics_current_indices = indices
        tw.setItem(nrows, 0, item)
        self.fixup_ics_table(tw)
        tw.setCurrentCell(nrows, 0) # Might as well make it selected

    def ics_find_index(self):
        n = 1
        while n in self.ics:
            n += 1
        return n

    def ics_delete_regions(self):
        tw = self.ui.initial_conditions.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())[:]
        for kw in kwlist:
            key, args = kw.key, kw.args
            # TODO use keyword_args here instead of startswith
            if key.startswith('ic_') and args and args[0] in self.ics_current_indices:
                self.unset_keyword(key, args=args)

        # TODO: fix any resulting holes in index sequence!

        for r in self.ics_current_regions:
            if r in self.ics_region_dict:
                self.ics_region_dict[r]['available'] = True


        for i in self.ics_current_indices:
            del self.ics[i]

        self.ics_current_regions = []
        self.ics_current_indices = []

        tw.removeRow(row)
        self.ics_setup_current_tab()

    def handle_ics_region_selection(self):
        ics = self.ui.initial_conditions
        table = ics.tablewidget_regions
        row = get_selected_row(table)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = table.item(row,0).data(UserRole)
        self.ics_current_indices, self.ics_current_regions = indices, regions
        enabled = (row is not None)
        ics.toolbutton_delete.setEnabled(enabled)
        ics.scrollarea_detail.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ics.scrollarea_detail):
                if isinstance(widget, QLineEdit): # Does this work for LineEdit?
                    widget.setText('')
            return

        self.ics_setup_current_tab() # reinitialize all widgets in current tab


    def fixup_ics_table(self, tw, stretch_column=0):
        hv = QtWidgets.QHeaderView
        if PYQT5:
            resize = tw.horizontalHeader().setSectionResizeMode
        else:
            resize = tw.horizontalHeader().setResizeMode
        ncols = tw.columnCount()
        for n in range(0, ncols):
            resize(n, hv.Stretch if n==stretch_column else hv.ResizeToContents)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # TODO FIXME scrollbar handling is not right - scrollbar status can change
        # outside of this function.  We need to call this everytime window geometry changes
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()

        if nrows==0:
            height = header_height+scrollbar_height
        else:
            height =  (header_height+scrollbar_height
                       + nrows*tw.rowHeight(0) + 4) # extra to avoid unneeded scrollbar
        tw.setMaximumHeight(height) # Works for tablewidget inside groupbox
        tw.setMinimumHeight(height) #? needed? should we allow scrollbar?
        tw.updateGeometry() #? needed?

    def ics_update_enabled(self):
        # If there are no solids, no scalar equations, and the fluid solver is disabled,
        # then we have no input tabs on the ICs pane, so disable it completely
        disabled = (len(self.ui.regions)==0
                    or (self.fluid_solver_disabled
                        and self.project.get_value('nscalar',default=0)==0
                        and len(self.solids)==0))

        self.find_navigation_tree_item("Initial Conditions").setDisabled(disabled)


    def ics_change_tab(self, tab, solid):
        ics = self.ui.initial_conditions
        index = (0 if tab==FLUID_TAB
                 else len(self.solids)+1 if tab==SCALAR_TAB
                 else solid)

        for i in range(ics.tab_box.columnCount()):
            item = ics.tab_box.itemAtPosition(0, i)
            if not item:
                continue
            widget = item.widget()
            if not widget:
                continue
            font = widget.font()
            font.setBold(i==index)
            widget.setFont(font)

        current_index = ics.stackedwidget.currentIndex()
        # If we're switching from solid m to solid n, we need some
        # special handling, because both tabs are really the same
        # widget.  We make a picture of the current tab, display that
        # in a dummy pane, then slide back to the solids tab
        if tab == current_index == SOLIDS_TAB:
            if solid == self.ics_current_solid:
                return # Really nothing to do

            if solid > self.ics_current_solid:
                dummy_label = ics.label_dummy_solids_L
                dummy_tab = SOLIDS_TAB_DUMMY_L
            else:
                dummy_label = ics.label_dummy_solids_R
                dummy_tab = SOLIDS_TAB_DUMMY_R

            picture = QPicture() #ics.page_solids.size())
            ics.page_solids.render(picture)
            dummy_label.setPicture(picture)
            ics.stackedwidget.setCurrentIndex(dummy_tab)

        self.ics_current_tab = tab
        self.ics_current_solid = solid

        #update tab contents
        if tab==FLUID_TAB:
            self.setup_ics_fluid_tab()
        elif tab==SOLIDS_TAB:
            self.setup_ics_solids_tab(self.ics_current_solid)
        elif tab==SCALAR_TAB:
            self.setup_ics_scalar_tab()

        # change stackedwidget contents
        self.animate_stacked_widget(
            ics.stackedwidget,
            ics.stackedwidget.currentIndex(),
            tab,
            direction='horizontal',
            line = ics.tab_underline,
            to_btn = ics.tab_box.itemAtPosition(0, index),
            btn_layout = ics.tab_box)



    def ics_check_region_in_use(self, region):
        return any(data.get('region')==region for data in self.ics.values())

    def ics_update_region(self, name, data):
        for (i,ic) in self.ics.items():
            if ic.get('region') == name:
                self.ics_set_region_keys(name, i, data)

    def ics_set_region_keys(self, name, index,  data):
        # Update the keys which define the box-shaped region the IC applies to
        for (key, val) in zip(('x_w', 'y_s', 'z_b', 'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            key = 'ic_' + key
            self.update_keyword(key, val, args=[index])

    def reset_ics(self):
        pass
        # TODO implement
        # Clear regions table, remove solids tabs, disable all inputs

    def ics_to_str(self):
        ics = self.ui.initial_conditions
        tw = ics.tablewidget_regions
        data = [tw.item(i,0).data(UserRole)
                for i in range(tw.rowCount())]
        return JSONEncoder().encode(data)

    def ics_regions_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)
        for (indices, regions) in data:
            self.ics_add_regions_1(regions, indices)


    def setup_ics(self):
        ui = self.ui
        ics = ui.initial_conditions
        # Grab a fresh copy, may have been updated
        self.ics_region_dict = ui.regions.get_region_dict()

        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.ics.items():
            region = data['region']
            if region in self.ics_region_dict:
                self.ics_region_dict[region]['available'] = False

        self.fixup_ics_table(ics.tablewidget_regions)
        row = get_selected_row(ics.tablewidget_regions)
        enabled = (row is not None)
        ics.toolbutton_delete.setEnabled(enabled)
        ics.scrollarea_detail.setEnabled(enabled)

        #Tabs group initial condition parameters for phases and additional equations.
        # Tabs are unavailable if no input is required from the user.

        #Fluid tab - Unavailable if the fluid phase was disabled.
        b = ics.pushbutton_fluid
        b.setText(self.fluid_phase_name)
        b.setEnabled(not self.fluid_solver_disabled)
        if self.fluid_solver_disabled:
            if self.ics_current_tab == 0: # Don't stay on disabled tab
                self.ics_change_tab(*(SOLIDS_TAB, 1) if self.solids else (SCALAR,0))
        font = b.font()
        font.setBold(self.ics_current_tab == 0)
        b.setFont(font)

        #Each solid phase will have its own tab. The tab name should be the name of the solid
        # (Could do this only on solid name change)
        n_cols = ics.tab_box.columnCount()
        # Clear out the old ones
        for i in range(n_cols-1, 0, -1):
            item = ics.tab_box.itemAtPosition(0, i)
            if not item:
                continue
            widget = item.widget()
            if not widget:
                continue
            if widget in (ics.pushbutton_fluid, ics.pushbutton_scalar):
                continue
            ics.tab_box.removeWidget(widget)
            widget.setParent(None)
            widget.deleteLater()
        # And make new ones
        for (i, solid_name) in enumerate(self.solids.keys(),1):
            b = QPushButton(text=solid_name)
            w = b.fontMetrics().boundingRect(solid_name).width() + 20
            b.setMaximumWidth(w)
            b.setFlat(True)
            font = b.font()
            font.setBold(self.ics_current_tab==1 and i==self.ics_current_solid)
            b.setFont(font)
            b.clicked.connect(lambda clicked, i=i: self.ics_change_tab(SOLIDS_TAB, i))
            ics.tab_box.addWidget(b, 0, i)
        # Don't stay on disabled tab TODO
        # if self.ics_current_tab == 1 and ...

        #Scalar (tab) - Tab only available if scalar equations are solved
        # Move the 'Scalar' button to the right of all solids, if needed
        b = ics.pushbutton_scalar
        font = b.font()
        font.setBold(self.ics_current_tab==SCALAR_TAB)
        b.setFont(font)
        nscalar = self.project.get_value('nscalar', default=0)
        enabled = (nscalar > 0)
        b.setEnabled(enabled)
        if len(self.solids) > 0:
            ics.tab_box.removeWidget(b)
            ics.tab_box.addWidget(b, 0, 1+len(self.solids))
        # Don't stay on a disabled tab TODO
        # if self.ics_current_tab == 2 and nscalar == 0:
        #
        self.ics_setup_current_tab()


    def ics_setup_current_tab(self):
        if self.ics_current_tab == FLUID_TAB:
            self.setup_ics_fluid_tab()
        elif self.ics_current_tab == SOLIDS_TAB:
            self.setup_ics_solids_tab(self.ics_current_solid)
        elif self.ics_current_tab == SCALAR_TAB:
            self.setup_ics_scalar_tab()


    def update_ics_fluid_mass_fraction_table(self):
        ics = self.ui.initial_conditions
        table = ics.tablewidget_fluid_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        if not (self.fluid_species and self.ics_current_indices):
            self.fixup_ics_table(table)
            table.setEnabled(False)
            return
        table.setEnabled(True)
        IC0 = self.ics_current_indices[0]
        species = self.fluid_species
        if species:
            nrows = len(species) + 1 # 'Total' row at end
        else:
            nrows = 0
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        for (row, (species,data)) in enumerate(species.items()):
            alias = data.get('alias', species) # default to species if no alias
            table.setItem(row, 0, make_item(alias))
            # mass fraction
            le = LineEdit()
            le.setdtype('dp')
            le.setValInfo(min=0.0, max=1.0)
            key = 'ic_x_g'
            le.key = key
            le.args = [self.ics_current_indices, row+1]
            val = self.project.get_value(key, args=[IC0, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)
            le.value_updated.connect(self.handle_ics_fluid_mass_fraction)
            table.setCellWidget(row, 1, le)
        if species:
            table.setItem(nrows-1, 0, make_item("Total"))
            table.setItem(nrows-1, 1, make_item(''))
            item = table.item(nrows-1, 0)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            self.update_ics_fluid_mass_fraction_total()

        self.fixup_ics_table(table)

    def handle_ics_fluid_mass_fraction(self, widget, value_dict, args):
        ics = self.ui.initial_conditions
        key = 'ic_x_g'
        val = value_dict[key]
        table = ics.tablewidget_fluid_mass_fraction
        widget.updateValue(key, val)
        if val == '':
            self.unset_keyword(key, args=args)
        else:
            self.update_keyword(key, val, args=args)
        self.update_ics_fluid_mass_fraction_total()

    def update_ics_fluid_mass_fraction_total(self):
        if not self.ics_current_indices:
            return
        if not self.fluid_species:
            return
        IC0 = self.ics_current_indices[0]
        ics = self.ui.initial_conditions
        key = 'ic_x_g'
        table = ics.tablewidget_fluid_mass_fraction
        if table.rowCount() == 0:
            return
        total = sum(float(self.project.get_value(key, default=0.0, args=[IC0,i]))
                    for i in range(1,len(self.fluid_species)+1))
        item =  table.item(table.rowCount()-1, 1)
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setText(str(total))

        # DEFAULT - last defined species has mass fraction of 1.0
        # (only enforce this if no mass fractions are set)

        if total == 0.0 and self.fluid_species:
            for IC in self.ics_current_indices:
                for i in range(1, len(self.fluid_species)):
                    self.update_keyword('ic_x_g', 0.0, args=[IC, i])
                self.update_keyword('ic_x_g', 1.0, args=[IC, len(self.fluid_species)]) # Last defined species
            self.update_ics_fluid_mass_fraction_table()

    # TODO DRY out fluid/solids code
    def update_ics_solids_mass_fraction_table(self):
        ics = self.ui.initial_conditions
        table = ics.tablewidget_solids_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        P = self.ics_current_solid
        if not (P and self.solids_species.get(P) and self.ics_current_indices):
            self.fixup_ics_table(table)
            table.setEnabled(False)
            return
        table.setEnabled(True)
        IC0 = self.ics_current_indices[0]
        species = self.solids_species[P]
        if species:
            nrows = len(species) + 1 # 'Total' row at end
        else:
            nrows = 0
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        for (row, (species,data)) in enumerate(species.items()):
            alias = data.get('alias', species) # default to species if no alias
            table.setItem(row, 0, make_item(alias))
            # mass fraction
            le = LineEdit()
            le.setdtype('dp')
            le.setValInfo(min=0.0, max=1.0)
            key = 'ic_x_s'
            le.key = key
            le.args = [self.ics_current_indices, P, row+1]
            val = self.project.get_value(key, args=[IC0, P, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)
            # register_widget ? TODO FIXME
            le.value_updated.connect(self.handle_ics_solids_mass_fraction)
            table.setCellWidget(row, 1, le)
        if species:
            table.setItem(nrows-1, 0, make_item("Total"))
            table.setItem(nrows-1, 1, make_item(''))
            item = table.item(nrows-1, 0)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            self.update_ics_solids_mass_fraction_total()

        self.fixup_ics_table(table)

    def handle_ics_solids_mass_fraction(self, widget, value_dict, args):
        ics = self.ui.initial_conditions
        key = 'ic_x_s'
        val = value_dict[key]
        table = ics.tablewidget_fluid_mass_fraction
        widget.updateValue(key, val)
        if val == '':
            self.unset_keyword(key, args=args)
        else:
            self.update_keyword(key, val, args=args)
        self.update_ics_solids_mass_fraction_total()

    def update_ics_solids_mass_fraction_total(self):
        if not self.ics_current_indices:
            return
        IC0 = self.ics_current_indices[0]
        P = self.ics_current_solid
        if P is None:
            return
        species = self.solids_species.get(P)
        if not P:
            return
        ics = self.ui.initial_conditions
        key = 'ic_x_s'
        table = ics.tablewidget_solids_mass_fraction
        if table.rowCount() == 0:
            return
        total = sum(float(self.project.get_value(key, default=0.0, args=[IC0,P,i]))
                    for i in range(1,len(species)+1))
        item = table.item(table.rowCount()-1, 1)
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setText(str(total))

        # DEFAULT - last defined species has mass fraction of 1.0
        # (only enforce this if no mass fractions are set)
        # TODO set ic_x_g = 1 - sum(ic_x_s)
        if total == 0.0 and species:
            for IC in self.ics_current_indices:
                for i in range(1, len(species)):
                    self.update_keyword('ic_x_s', 0.0, args=[IC, P, i])
                self.update_keyword('ic_x_s', 1.0, args=[IC, P, len(species)]) # Last defined species
            self.update_ics_solids_mass_fraction_table()



    def ics_extract_regions(self):
        if self.ics:
            # We assume that ic regions have been initialized correctly
            # from mfix_gui_comments.
            # TODO: verify that there is an IC region for each IC
            return

        if self.ics_region_dict is None:
            self.ics_region_dict = self.ui.regions.get_region_dict()

        # TODO: if we wanted to be fancy, we could find regions where
        # IC values matched, and merge into a new IC region.  That
        # is only needed for projects created outside the GUI (otherwise
        # we have already stored the IC regions).  Also woule be nice
        # to offer a way to split compound regions.
        for ic in self.project.ics:

            d = ic.keyword_dict
            extent = [d.get(k,None) for k in ('ic_x_w', 'ic_y_s', 'ic_z_b',
                                              'ic_x_e', 'ic_y_n', 'ic_z_t')]
            extent = [None if x is None else x.value for x in extent]
            if any (x is None for x in extent):
                self.warn("initial condition %s: invalid extents %s" %
                           (ic.ind, extent))
                continue
            for (region_name, data) in self.ics_region_dict.items():
                ext2 = data.get('from',[]) + data.get('to',[])
                if data.get('from',[]) + data.get('to',[]) == extent:
                    if data.get('available', True):
                        self.ics_add_regions_1([region_name], [ic.ind])
                        break
            else:
                self.warn("initial condition %s: could not match defined region %s" %
                          (ic.ind, extent))


    def setup_ics_fluid_tab(self):
        #Fluid (tab)
        ics = self.ui.initial_conditions

        if self.fluid_solver_disabled:
            # we shouldn't be on this tab!
            return

        tw = ics.tablewidget_fluid_mass_fraction
        enabled = (self.fluid_species is not None)
        if not enabled:
            tw.clearContents()
            tw.setRowCount(0)
        self.fixup_ics_table(tw)
        tw.setEnabled(enabled)

        if not self.ics_current_indices:
            # Nothing selected.  What can we do? (Clear out all lineedits?)
            return

        IC0 = self.ics_current_indices[0]
        # Note - value may not be consistent acros grouped regions
        #  For now we're going to assume that it is, and just check
        #  first subregion of IC group

        #  If we can make this code generic enough perhaps someday it can
        # be autogenerated from SRS doc
        def get_widget(key):
            return getattr(ics, 'lineedit_keyword_%s_args_IC' % key)

        def setup_key_widget(key, default=None, enabled=True):
            for name in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_IC'):
                item = getattr(ics, name%key, None)
                if item:
                    item.setEnabled(enabled)
            if not enabled:
                get_widget(key).setText('')
                return

            val = self.project.get_value(key, args=[IC0])
            if val is None:
                val = default
            if val is not None:
                for IC in self.ics_current_indices:
                    self.update_keyword(key, val, args=[IC])
                get_widget(key).updateValue(key, val, args=[IC0])

        #Define volume fraction
        # Specification always available
        # Sets keyword IC_EP_G(#)
        # DEFAULT value of 1.0
        # (terminology:  is ic_ep_g volume fraction or void fraction?)
        key = 'ic_ep_g'
        default = 1.0
        setup_key_widget(key, default)
        get_widget(key).setReadOnly(True)
        get_widget(key).setEnabled(False)
        # Issues/121 ; make non-editable

        #Define temperature
        # Specification always available
        # Input required for any of the following
        # Fluid density model: Ideal Gas Law
        # Fluid viscosity model: Sutherland's Law
        # Energy equations are solved
        # Sets keyword IC_T_G(#)
        # DEFAULT value of 293.15
        #TODO: how do we enforce "required" inputs?
        key = 'ic_t_g'
        default = 293.15
        setup_key_widget(key, default)

        #Define pressure (optional)  TODO "optional" in label?
        # Specification always available
        # Sets keyword IC_P_g(#)
        # DEFAULT - no input
        key = 'ic_p_g'
        default = None
        setup_key_widget(key, default)

        #Define velocity components (required)
        # Specification always available
        # Sets keywords IC_U_G(#), IC_V_G(#), IC_W_G(#)
        # DEFAULT value of 0.0
        default = 0.0
        for key in ('ic_u_g', 'ic_v_g', 'ic_w_g'):
            setup_key_widget(key, default)

        #Select species and set mass fractions (table format)
        # Specification always available
        # Input required for species equations TODO: implement 'required'
        # Drop down menu of fluid species
        # DEFAULT - last defined species has mass fraction of 1.0  # implemented in update_ics_fluid_mass_fraction_table
        # Error check: mass fractions must sum to one   # we show total but don't warn if != 1.0
        self.update_ics_fluid_mass_fraction_table()

        key = 'turbulence_model'
        turbulence_model = self.project.get_value(key)
        enabled = (turbulence_model is not None)
        ics.groupbox_turbulence.setEnabled(enabled)

        #Turbulence: Define mixing length model length scale
        # Specification only available with Mixing Length turbulence model
        # Sets keyword IC_L_SCALE(#)
        # DEFAULT value of 1.0
        key = 'ic_l_scale'
        default = 1.0
        enabled = (turbulence_model == 'MIXING_LENGTH')
        setup_key_widget(key, default, enabled)

        #Turbulence: Define k-ε turbulent kinetic energy
        # Specification only available with K-Epsilon turbulence model
        # Sets keyword IC_K_TURB_G(#)
        # DEFAULT value of 0.0
        key = 'ic_k_turb_g'
        default = 0.0
        enabled = (turbulence_model == 'K_EPSILON')
        setup_key_widget(key, default, enabled)

        #Turbulence: Define k-ε turbulent dissipation
        # Specification only available with K-Epsilon turbulence model
        # Sets keywords IC_E_TURB_G(#)
        # DEFAULT value of 0.0
        key = 'ic_e_turb_g'
        default = 0.0
        enabled = (turbulence_model == 'K_EPSILON')
        setup_key_widget(key, default, enabled)

        energy_eq = self.project.get_value('energy_eq', default=True)
        enabled = bool(energy_eq)
        ics.groupbox_fluid_advanced.setEnabled(enabled)

        #Advanced: Define radiation coefficient
        # Specification only available when solving energy equations
        # Sets keyword IC_GAMA_RG(#)
        # DEFAULT value of 0.0
        key = 'ic_gama_rg'
        default = 0.0
        enabled = bool(energy_eq)
        setup_key_widget(key, default, enabled)

        #Advanced: Define radiation temperature
        # Specification only available when solving energy equations
	# Sets keyword IC_T_RG(#)
        # DEFAULT value of 293.15
        key = 'ic_t_rg'
        default =  293.15
        enabled = bool(energy_eq)
        setup_key_widget(key, default, enabled)

    def setup_ics_solids_tab(self, P):
        # Solid-# (tab) - Rename tab to user provided solids name.

        # Note, solids phases are numbered 1-N
        self.ics_current_solid = P
        if P is None: # Nothing to do
            return

        if not self.ics_current_indices: # No region selected
            # TODO clear all widgets
            return

        ics = self.ui.initial_conditions
        IC0 = self.ics_current_indices[0]

        # issues/121
        self.ics_set_volume_fraction_limit()
        widget = ics.lineedit_keyword_ic_ep_s_args_IC_P
        # Have to do this after project manager has registered widgets
        widget.value_updated.disconnect()
        widget.value_updated.connect(self.handle_ics_volume_fraction)

        # Generic!
        def get_widget(key):
            widget = getattr(ics, 'lineedit_keyword_%s_args_IC_P' % key, None)
            if widget:
                return widget
            widget = getattr(ics, 'lineedit_keyword_%s_args_IC' % key, None)
            if not widget:
                self.error('no widget for key %s' % key)
            return widget

        def setup_key_widget(key, default=None, enabled=True):
            for name in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_IC_P',
                         'lineedit_keyword_%s_args_IC'):
                item = getattr(ics, name%key, None)
                if item:
                    item.setEnabled(enabled)
            if not enabled:
                get_widget(key).setText('')
                return

            no_p_keys =  ('ic_p_star',)
            args = [IC0] if key in no_p_keys else [IC0,P]
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
            if val is not None:
                for IC in self.ics_current_indices:
                    self.update_keyword(key, val, args=[IC] if key in no_p_keys
                                        else [IC,P])
                get_widget(key).updateValue(key, val, args=args)


        #Group tab inputs by equation type (e.g., momentum, energy, species).
        # Making the grouped inputs a 'collapsible list' may make navigation easier.
        #  (Note - collaspsing not implemented)

        #Define volume fraction (required)
        # Specification always available
        # Sets keyword IC_EP_S(#,#)
        # DEFAULT value of (0.0, 1 - SUM)
        # - spec changed - default 0, set ip_ep_g to 1-sum
        key = 'ic_ep_s'
        default = 0.0
        setup_key_widget(key, default)

        #Define temperature
        # Specification always available
        # Input required when solving energy equations
        # Sets keyword IC_T_S(#,#)
        # DEFAULT value of 293.15
        energy_eq = self.project.get_value('energy_eq', default=True)
        key = 'ic_t_s'
        default = 293.15
        enabled = bool(energy_eq)
        setup_key_widget(key, default, enabled)

        #Define velocity components (required)
        # Specification always available
        # Sets keywords IC_U_S(#,#), IC_V_S(#,#), IC_W_S(#,#)
        # DEFAULT value of 0.0
        for key in 'ic_u_s', 'ic_v_s', 'ic_w_s':
            setup_key_widget(key, default=0.0)

        #Define pressure (optional) # TODO 'optional' in label?
        # Specification only available for SOLIDS_MODEL(#)='TFM'
        # Sets keyword IC_P_STAR(#)
        # DEFAULT of 0.0
        # TODO: If a value is entered, then if you edit it on another solids tab,
        # you get a warning that it applies to all solids in the IC region.

        solids_model = self.project.get_value('solids_model', args=[P])
        enabled = (solids_model=='TFM')
        key = 'ic_p_star'
        default = 0.0
        setup_key_widget(key, default, enabled)

        #Define granular temperature
        # Specification only available for SOLIDS_MODEL(#)='TFM' and non-algebraic
        # formulation viscous stress model (see continuous solids model section) or for
        # SOLIDS_MODEL(#)=DEM' or SOLIDS_MODEL(#)='PIC'
        # Sets keyword IC_THETA_M(#,#)
        # DEFAULT value of 0.0
        # TODO: use (m^2/sec^2) for solids granular energy units.
        # Some of the KT_TYPES also include a mass unit (kg) -
        # but the default model (Lun) will have units of m^s/sec^2
        solids_model = self.project.get_value('solids_model', args=[P])
        kt_type = self.project.get_value('kt_type')
        enabled = ( (solids_model=='TFM' and kt_type != 'ALGEBRAIC')
                    or solids_model=='DEM'
                    or solids_model=='PIC')
        key = 'ic_theta_m'
        default = 0.0
        setup_key_widget(key, default, enabled)

        #Define particles per parcel
        # Specification only available for SOLIDS_MODEL(#)='PIC'
        # Sets keyword IC_PIC_CONST_STATWT(#,#)
        # DEFAULT value of 10.0
        enabled = (solids_model=='PIC')
        key = 'ic_pic_const_statwt'
        default = 10.0
        setup_key_widget(key, default, enabled)

        #Select species and set mass fractions (table format)
        # Specification always available
        # Input required for species equations
        # Drop down menu of solids species
        # DEFAULT - last defined species has mass fraction of 1.0
        # Error check: mass fractions must sum to one
        self.update_ics_solids_mass_fraction_table()

        enabled = (solids_model=='DEM' or bool(energy_eq))
        ics.groupbox_solids_advanced.setEnabled(enabled)

        #Advanced: Option to enable fitting DES particles to region
        # Option only available for DEM solids
        # Sets keyword: IC_DES_FIT_TO_REGION
        # Disabled [DEFAULT]
        enabled = (solids_model=='DEM')
        item = ics.checkbox_keyword_ic_des_fit_to_region_args_IC
        item.setEnabled(enabled)
        key = 'ic_des_fit_to_region'
        default = False
        val = self.project.get_value(key, args=[IC0])
        if val is None:
            val = default
            if enabled:
                for ic in self.ics_current_indices:
                    self.update_keyword(key, val, args=ic)
        item.setChecked(val)
        # TODO: popup dialog to warn that this apples to all phases

        #Advanced: Define radiation coefficient
        # Specification only available when solving energy equations
        # Sets keyword IC_GAMA_S(#,#)   # emended to IC_GAMA_RS
        # DEFAULT value of 0.0
        enabled = bool(energy_eq)
        key = 'ic_gama_rs'
        default = 0.0
        setup_key_widget(key, default, enabled)

        #Advanced: Define radiation temperature
        # Specification only available when solving energy equations
        # Sets keyword IC_T_RS(#,#)
        # DEFAULT value of 293.15
        enabled = bool(energy_eq)
        key = 'ic_t_rs'
        default = 293.15
        setup_key_widget(key, default, enabled)


    def setup_ics_scalar_tab(self):
        #Note that this tab should only be available if scalar equations are being solved.
        #Furthermore, the number of scalars requiring input comes from the number of
        #scalar equations specified by the user.

        # Note - there's a bunch of widget registration/unregistration/creation/deletion
        # here that could be avoided - we really only need to do this when nscalar changes

        if not self.ics_current_indices:
            return # No selection
        IC0 = self.ics_current_indices[0]

        ics = self.ui.initial_conditions
        nscalar = self.project.get_value('nscalar', default=0)

        page =  ics.page_scalar
        layout = page.layout()

        spacer = None
        for i in range(layout.rowCount()-1, -1, -1):
            for j in (1,0):
                item = layout.itemAtPosition(i,j)
                if not item:
                    continue
                widget = item.widget()
                if not widget:
                    spacer = item
                    continue
                if isinstance(widget, LineEdit):
                    self.project.unregister_widget(widget)
                widget.setParent(None)
                widget.deleteLater()

        if spacer:
            layout.removeItem(spacer)

        #Define initial scalar value
        #Sets keyword IC_SCALAR(#,#)
        #DEFAULT value of 0.0
        key = 'ic_scalar'
        for i in range(nscalar):
            label = QLabel('Scalar %s' % (i+1))
            layout.addWidget(label, i, 0)
            le = LineEdit()
            le.key = key
            le.args = ['IC', (i+1)]
            self.project.register_widget(le, [key], ['IC', (i+1)])
            le.setdtype('dp')
            le.default_value = 0.0

            val = self.project.get_value(key, args=[IC0, (i+1)])
            if val is None:
                val = 0.0
                for ic in self.ics_current_indices:
                    self.update_keyword(key, val, args=[ic, (i+1)])
            le.setText(str(val))
            layout.addWidget(le, i, 1)

        if spacer:
            layout.addItem(spacer, i+1, 0)
