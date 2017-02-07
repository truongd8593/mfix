# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import QLabel, QPushButton, QWidget
from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4

UserRole = QtCore.Qt.UserRole

from mfixgui.widgets.base import LineEdit

from mfixgui.tools.general import (set_item_noedit, get_selected_row,
                                   widget_iter, safe_float)

from mfixgui.tools.keyword_args import mkargs

# We don't need extended JSON here
from json import JSONDecoder, JSONEncoder

FLUID_TAB = 0
SOLIDS_TAB_DUMMY_L = 1
SOLIDS_TAB = 2
SOLIDS_TAB_DUMMY_R = 3
SCALAR_TAB = 4

class ICS(object):
    # Initial Conditions Task Pane Window: This section allows a user to define the initial conditions
    # for the described model. This section relies on regions named in the Regions section.

    def init_ics(self):
        ui = self.ui.initial_conditions

        self.ics = {} # key: index.  value: data dictionary for initial cond
        self.ics_current_indices = [] # List of IC indices
        self.ics_current_regions = [] # And the names of the regions which define them
        self.ics_region_dict = None
        self.ics_saved_solids_names = []

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
        ui.toolbutton_add.clicked.connect(self.ics_show_regions_popup)
        ui.toolbutton_delete.clicked.connect(self.ics_delete_regions)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        ui.tablewidget_regions.itemSelectionChanged.connect(self.handle_ics_region_selection)

        self.ics_current_tab = FLUID_TAB # If fluid is disabled, we will switch
        self.ics_current_solid = self.P = None
        ui.pushbutton_fluid.pressed.connect(lambda: self.ics_change_tab(FLUID_TAB,None))
        ui.pushbutton_scalar.pressed.connect(lambda: self.ics_change_tab(SCALAR_TAB,None))

        # Trim width of "Fluid" and "Scalar" buttons, like we do for
        # dynamically-created "Solid #" buttons
        for b in (ui.pushbutton_fluid, ui.pushbutton_scalar):
            w = b.fontMetrics().boundingRect(b.text()).width() + 20
            b.setMaximumWidth(w)

        # Not auto-registered with project manager
        widget = ui.lineedit_ic_ep_s_args_IC_P
        key = 'ic_ep_s'
        widget.key = key
        widget.args = ['IC', 'P']
        self.add_tooltip(widget, key)
        widget.dtype = float
        widget.value_updated.connect(self.handle_ics_volume_fraction)

        widget = ui.lineedit_ic_p_star_args_IC
        widget.key = 'ic_p_star'
        widget.args = ['IC']
        widget.dtype = float
        widget.value_updated.connect(self.handle_ic_p_star)





    def handle_ic_p_star(self, widget, data, args):
        if not self.ics_current_indices:
            return
        IC0 = self.ics_current_indices[0]
        key = 'ic_p_star'
        prev_val = self.project.get_value(key, args=[IC0])
        new_key, new_val = data.popitem()

        ui = self.ui.initial_conditions
        if len(self.solids) > 1:
            resp=self.message(text="Pressure setting applies to all solids phases\nAre you sure?",
                              buttons=['yes','no'],
                              default = 'no')
            if resp != 'yes': # Reject update, set lineedit back to previous value
                widget.updateValue(self.project.get_value(key, prev_val))
                return
        for IC in self.ics_current_indices:
            self.update_keyword(key, new_val, args=[IC])


    def ics_set_volume_fraction_limit(self):
        if not self.ics_current_indices:
            return
        if not self.ics_current_solid:
            return
        IC0 = self.ics_current_indices[0]
        P = self.ics_current_solid
        ui = self.ui.initial_conditions
        key = 'ic_ep_s'

        s = sum(safe_float(self.project.get_value(key, default=0, args=[IC0, s]))
                for s in range(1, len(self.solids)+1) if s != P)

        lim = max(0, 1.0 - s)
        lim = round(lim, 10) # avoid problem with 1 - 0.9 != 0.1

        widget = ui.lineedit_ic_ep_s_args_IC_P
        widget.min = 0.0
        widget.max = lim


    def handle_ics_volume_fraction(self, widget, val, args):
        if not self.ics_current_indices:
            return
        IC0 = self.ics_current_indices[0]
        if not self.ics_current_solid:
            return
        P = self.ics_current_solid
        ui = self.ui.initial_conditions
        key = 'ic_ep_s'
        self.project.submit_change(widget, val, args)

        s = sum(safe_float(self.project.get_value(key, default=0, args=[IC0, s]))
                for s in range(1, len(self.solids)+1))
        if s > 1.0:
            self.warning("Volume fractions sum to %s, must be <= 1.0" % s,
                         popup=True)
            return # ?
        # Set ic_ep_g from ic_ep_s (issues/121)
        val = round(1.0 - s, 10)
        for IC in self.ics_current_indices:
            self.update_keyword('ic_ep_g', val, args=[IC])


    def ics_show_regions_popup(self):
        #Users cannot select inapplicable regions.
        # IC regions must be volumes or planes (not points or STLs)
        #  Volumes are always valid IC regions
        #  XY-Planes are valid IC regions for 2D simulations (NO_K=.TRUE.)
        #  XZ- and YZ Planes are never valid IC regions
        #No region can define more than one initial condition.
        ui = self.ui.initial_conditions
        rp = self.regions_popup
        rp.clear()
        no_k = self.project.get_value('no_k', default=False)

        for (name,data) in self.ics_region_dict.items():
            shape = data.get('type', '---')
            # Assume available if unmarked
            available = (data.get('available', True)
                         #and not self.check_region_in_use(name) # allow region sharing
                         and (shape == 'box') or (no_k and shape=='XY-plane'))

            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.ics_add_regions)
        rp.cancel.connect(self.ics_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.bottom_frame,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup('initial conditions')


    def ics_cancel_add(self):
        ui = self.ui.initial_conditions

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.bottom_frame,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def ics_add_regions(self):
        # Interactively add regions to define ICs
        ui = self.ui.initial_conditions
        rp = self.regions_popup
        self.ics_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        if not selections:
            return
        self.ics_add_regions_1(selections, autoselect=True)
        self.ics_setup_current_tab() # Update the widgets


    def ics_add_regions_1(self, selections, indices=None, autoselect=False):
        # Used by both interactive and load-time add-region handlers
        ui = self.ui.initial_conditions

        if self.ics_region_dict is None:
            self.ics_region_dict = self.ui.regions.get_region_dict()

        tw = ui.tablewidget_regions
        nrows = tw.rowCount()
        tw.setRowCount(nrows+1)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        item = make_item('+'.join(selections))

        if indices is None: # interactive
            indices = [None] * len(selections)
        else: # loading file
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
        tw.setItem(nrows, 0, item)

        self.fixup_ics_table(tw)

        if autoselect:
            tw.setCurrentCell(nrows, 0)


    def ics_find_index(self):
        n = 1
        while n in self.ics:
            n += 1
        return n


    def ics_delete_regions(self):
        ui = self.ui.initial_conditions
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            if key.startswith('ic_') and args and args[0] in self.ics_current_indices:
                self.unset_keyword(key, args=args)

        for r in self.ics_current_regions:
            if r in self.ics_region_dict:
                self.ics_region_dict[r]['available'] = True

        for i in self.ics_current_indices:
            del self.ics[i]

        self.ics_current_regions = []
        self.ics_current_indices = []

        tw.removeRow(row)
        self.fixup_ics_table(tw)
        self.ics_setup_current_tab()
        self.update_nav_tree()


    def ics_delete_solids_phase(self, phase_index):
        """adjust ics_current_solid when solids phase deleted"""
        if (self.ics_current_solid is not None and
            self.ics_current_solid >= phase_index):
            self.ics_current_solid -= 1
            if self.ics_current_solid == 0:
                self.ics_current_solid = None


    def handle_ics_region_selection(self):
        ui = self.ui.initial_conditions
        table = ui.tablewidget_regions
        row = get_selected_row(table)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = table.item(row,0).data(UserRole)
        self.ics_current_indices, self.ics_current_regions = indices, regions
        enabled = (row is not None)
        for item in (ui.toolbutton_delete,
                     ui.bottom_frame):
            item.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.bottom_frame):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.ics_setup_current_tab() # reinitialize all widgets in current tab
        # Scroll to top
        ui.scrollarea_detail.ensureVisible(0, 0)


    def fixup_ics_table(self, tw, stretch_column=0):
        ui = self.ui.initial_conditions
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

        # Note - scrollbar status can change outside of this function.
        # Do we need to call this everytime window geometry changes?
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()
        if nrows==0:
            row_height = 0
            height = header_height+scrollbar_height
        else:
            row_height = tw.rowHeight(0)
            height =  (header_height+scrollbar_height
                       + nrows*row_height + 4) # extra to avoid unneeded scrollbar

        if tw == ui.tablewidget_regions: # main table, adjust top splitter
            ui.top_frame.setMaximumHeight(height+24)
            ui.top_frame.setMinimumHeight(header_height+24+row_height*min(nrows,5))
            ui.top_frame.updateGeometry()
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
        else: # mass fraction tables
            if tw == ui.tablewidget_solids_mass_fraction:
                ui.groupbox_solids_composition.setMaximumHeight(height+30)
            elif tw == ui.tablewidget_fluid_mass_fraction:
                ui.groupbox_fluid_composition.setMaximumHeight(height+30)
            tw.setMaximumHeight(height) # Works for tablewidget inside groupbox
            tw.setMinimumHeight(height) #? needed? should we allow scrollbar?
        tw.updateGeometry() #? needed?


    def ics_update_enabled(self):
        if self.ics:
            # Never disable if there are ICs defined
            disabled = False
        else:
            # If there are no solids, no scalar equations, and the fluid solver is disabled,
            # then we have no input tabs on the ICs pane, so disable it completely
            regions = self.ui.regions.get_region_dict()
            nregions = sum(1 for (name, r) in regions.items()
                           if not self.check_region_in_use(name)
                           and r.get('type')=='box')
            disabled = (nregions==0
                        or (self.fluid_solver_disabled
                            and self.project.get_value('nscalar',default=0)==0
                            and len(self.solids)==0))
        self.find_navigation_tree_item("Initial Conditions").setDisabled(disabled)


    def ics_change_tab(self, tab, solid):
        ui = self.ui.initial_conditions
        index = (0 if tab==FLUID_TAB
                 else len(self.solids)+1 if tab==SCALAR_TAB
                 else solid)

        for i in range(ui.tab_layout.columnCount()):
            item = ui.tab_layout.itemAtPosition(0, i)
            if not item:
                continue
            widget = item.widget()
            if not widget:
                continue
            font = widget.font()
            font.setBold(i==index)
            widget.setFont(font)

        current_index = ui.stackedwidget.currentIndex()
        # If we're switching from solid m to solid n, we need some
        # special handling, because both tabs are really the same
        # widget.  We make a picture of the current tab, display that
        # in a dummy pane, then slide back to the solids tab
        if tab == current_index == SOLIDS_TAB:
            if solid == self.ics_current_solid:
                return # Really nothing to do

            if solid > (self.ics_current_solid or 0):
                dummy_label = ui.label_dummy_solids_L
                dummy_tab = SOLIDS_TAB_DUMMY_L
            else:
                dummy_label = ui.label_dummy_solids_R
                dummy_tab = SOLIDS_TAB_DUMMY_R

            pixmap = QPixmap(ui.page_solids.size())
            pixmap.fill() # fill bg with white
            ui.page_solids.render(pixmap, flags=QWidget.DrawChildren)  #avoid rendering bg
            dummy_label.setPixmap(pixmap)
            ui.stackedwidget.setCurrentIndex(dummy_tab)

        self.ics_current_tab = tab
        self.ics_current_solid = self.P = solid if tab==SOLIDS_TAB else None
        self.ics_setup_current_tab()

        # change stackedwidget contents
        self.animate_stacked_widget(
            ui.stackedwidget,
            ui.stackedwidget.currentIndex(),
            tab,
            direction='horizontal',
            line = ui.tab_underline,
            to_btn = ui.tab_layout.itemAtPosition(0, index),
            btn_layout = ui.tab_layout)
        # Scroll to top
        ui.scrollarea_detail.ensureVisible(0, 0)


    def ics_check_region_in_use(self, name):
        return any(data.get('region')==name for data in self.ics.values())


    def ics_update_region(self, name, data):
        for (i,ic) in self.ics.items():
            if ic.get('region') == name:
                self.ics_set_region_keys(name, i, data)


    def ics_set_region_keys(self, name, idx, data):
        # Update the keys which define the region the IC applies to
        no_k = self.project.get_value('no_k')
        for (key, val) in zip(('x_w', 'y_s', 'z_b',
                               'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            # ic_z_t and ic_z_b keywords should not be added when no_k=True
            if no_k and key in ('z_t', 'z_b'):
                continue
            self.update_keyword('ic_'+key, val, args=[idx])


    def ics_change_region_name(self, old_name, new_name):
        ui = self.ui.initial_conditions
        for (key, val) in self.ics.items():
            if val.get('region') == old_name:
                self.ics[key]['region'] = new_name
                tw = ui.tablewidget_regions
                for i in range(tw.rowCount()):
                    data = tw.item(i,0).data(UserRole)
                    indices, names = data
                    if key in indices:
                        item = tw.item(i,0)
                        names = [new_name if n==old_name else n for n in names]
                        item.setData(UserRole, (indices, names))
                        item.setText('+'.join(names))
                        break
                break


    def reset_ics(self):
        self.ics.clear()
        self.ics_current_indices = []
        self.ics_current_regions = []
        self.ics_region_dict = None
        self.ics_current_solid = self.P = None
        ui = self.ui.initial_conditions
        ui.tablewidget_regions.clearContents()
        ui.tablewidget_regions.setRowCount(0)
        # anything else to do here?


    def ic_regions_to_str(self):
        ui = self.ui.initial_conditions
        tw = ui.tablewidget_regions
        data = [tw.item(i,0).data(UserRole)
                for i in range(tw.rowCount())]
        return JSONEncoder().encode(data)


    def ics_regions_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)
        for (indices, regions) in data:
            self.ics_add_regions_1(regions, indices, autoselect=False)


    def setup_ics(self):
        ui = self.ui.initial_conditions

        # Grab a fresh copy, may have been updated
        self.ics_region_dict = self.ui.regions.get_region_dict()

        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.ics.items():
            region = data['region']
            if region in self.ics_region_dict:
                self.ics_region_dict[region]['available'] = False

        self.fixup_ics_table(ui.tablewidget_regions)
        row = get_selected_row(ui.tablewidget_regions)
        # Autoselect if only 1 row
        if row is None and ui.tablewidget_regions.rowCount() == 1:
            row = 0
            ui.tablewidget_regions.setCurrentCell(row, 0)
        enabled = (row is not None)
        for item in (ui.toolbutton_delete,
                     ui.bottom_frame):
            item.setEnabled(enabled)

        #Tabs group initial condition parameters for phases and additional equations.
        # Tabs are unavailable if no input is required from the user.

        #Fluid tab - Unavailable if the fluid phase was disabled.
        setup_done = False
        b = ui.pushbutton_fluid
        b.setText(self.fluid_phase_name)
        b.setEnabled(not self.fluid_solver_disabled)
        if self.fluid_solver_disabled:
            if self.ics_current_tab == FLUID_TAB: # Don't stay on disabled tab
                self.ics_change_tab(*self.ics_find_valid_tab())
                setup_done = True
        font = b.font()
        font.setBold(self.ics_current_tab == FLUID_TAB)
        b.setFont(font)

        #Each solid phase will have its own tab. The tab name should be the name of the solid
        solids_names = list(self.solids.keys())
        if self.ics_saved_solids_names != solids_names:
            # Clear out the old ones
            n_cols = ui.tab_layout.columnCount()
            for i in range(n_cols-1, 0, -1):
                item = ui.tab_layout.itemAtPosition(0, i)
                if not item:
                    continue
                widget = item.widget()
                if not widget:
                    continue
                if widget in (ui.pushbutton_fluid, ui.pushbutton_scalar):
                    continue
                ui.tab_layout.removeWidget(widget)
                widget.setParent(None)
                widget.deleteLater()
            # And make new ones
            for (i, solid_name) in enumerate(solids_names, 1):
                b = QPushButton(text=solid_name)
                w = b.fontMetrics().boundingRect(solid_name).width() + 20
                b.setMaximumWidth(w)
                b.setFlat(True)
                font = b.font()
                font.setBold(self.ics_current_tab==SOLIDS_TAB and i==self.ics_current_solid)
                b.setFont(font)
                b.pressed.connect(lambda i=i: self.ics_change_tab(SOLIDS_TAB, i))
                ui.tab_layout.addWidget(b, 0, i)

        # Don't stay on a disabled tab
        if self.ics_current_tab == SOLIDS_TAB and not self.solids:
            self.ics_change_tab(*self.ics_find_valid_tab())
            setup_done = True


        #Scalar (tab) - Tab only available if scalar equations are solved
        # Move the 'Scalar' button to the right of all solids, if needed
        b = ui.pushbutton_scalar
        font = b.font()
        font.setBold(self.ics_current_tab==SCALAR_TAB)
        b.setFont(font)
        nscalar = self.project.get_value('nscalar', default=0)
        enabled = (nscalar > 0)
        b.setEnabled(enabled)
        if len(self.solids) != len(self.ics_saved_solids_names):
            ui.tab_layout.removeWidget(b)
            ui.tab_layout.addWidget(b, 0, 1+len(self.solids))

        # Don't stay on a disabled tab
        if self.ics_current_tab == SCALAR_TAB and nscalar == 0:
            self.ics_change_tab(*self.ics_find_valid_tab())
            setup_done = True

        self.ics_saved_solids_names = solids_names
        self.P = self.ics_current_solid

        if not setup_done:
            self.ics_setup_current_tab()

        # make sure underline is in the right place, as # of solids may
        # have changed (lifted from animate_stacked_widget, which we
        # don't want to call here)
        tab = self.ics_current_tab
        line_to = (0 if tab==FLUID_TAB
                   else len(self.solids)+1 if tab==SCALAR_TAB
                   #else len(self.solids)+2 if tab==CYCLIC_TAB
                   else self.ics_current_solid)
        line = ui.tab_underline
        btn_layout = ui.tab_layout
        if line_to is not None:
            btn_layout.addItem(btn_layout.takeAt(
                btn_layout.indexOf(line)), 1, line_to)


    def ics_find_valid_tab(self): # Don't stay on a disabled tab
        if not self.fluid_solver_disabled:
            return (FLUID_TAB, None)
        elif self.solids:
            return (SOLIDS_TAB, 0)
        elif self.project.get_value('nscalar', default=0) != 0:
            return (SCALAR_TAB, None)
        else:
            self.error("Initial condition:  all tabs disabled!")
            return (FLUID_TAB, None) # What else to do?


    def ics_setup_current_tab(self):
        if self.ics_current_tab == FLUID_TAB:
            self.setup_ics_fluid_tab()
        elif self.ics_current_tab == SOLIDS_TAB:
            self.setup_ics_solids_tab(self.ics_current_solid)
        elif self.ics_current_tab == SCALAR_TAB:
            self.setup_ics_scalar_tab()


    def update_ics_fluid_mass_fraction_table(self):
        ui = self.ui.initial_conditions
        table = ui.tablewidget_fluid_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        if not (self.fluid_species and self.ics_current_indices):
            self.fixup_ics_table(table)
            ui.groupbox_fluid_composition.setEnabled(False)
            return
        ui.groupbox_fluid_composition.setEnabled(True)
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
            le.setValInfo(min=0.0, max=1.0) # TODO adjust max dynamically
            key = 'ic_x_g'
            le.key = key
            le.args = [self.ics_current_indices, row+1]
            self.add_tooltip(le, key)
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
        ui = self.ui.initial_conditions
        key = 'ic_x_g'
        val = value_dict[key]
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
        ui = self.ui.initial_conditions
        key = 'ic_x_g'
        table = ui.tablewidget_fluid_mass_fraction
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


    def update_ics_solids_mass_fraction_table(self):
        ui = self.ui.initial_conditions
        table = ui.tablewidget_solids_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        P = self.ics_current_solid
        if not (P and self.solids_species.get(P) and self.ics_current_indices):
            self.fixup_ics_table(table)
            table.setEnabled(False)
            ui.groupbox_solids_composition.setEnabled(False)
            return
        ui.groupbox_solids_composition.setEnabled(True)
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
            self.add_tooltip(le, key)
            val = self.project.get_value(key, args=[IC0, P, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)
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
        ui = self.ui.initial_conditions
        key = 'ic_x_s'
        val = value_dict[key]
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
        ui = self.ui.initial_conditions
        key = 'ic_x_s'
        table = ui.tablewidget_solids_mass_fraction
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
        if total == 0.0 and species:
            for IC in self.ics_current_indices:
                for i in range(1, len(species)):
                    self.update_keyword('ic_x_s', 0.0, args=[IC, P, i])
                self.update_keyword('ic_x_s', 1.0, args=[IC, P, len(species)]) # Last defined species
            self.update_ics_solids_mass_fraction_table()


    def ics_extract_regions(self):
        if self.ics:
            # We assume that IC regions have been initialized correctly
            # from mfix_gui_comments.
            # TODO: verify that there is an IC region for each IC
            return

        if self.ics_region_dict is None:
            self.ics_region_dict = self.ui.regions.get_region_dict()

        # TODO: if we wanted to be fancy, we could find regions where
        # IC values matched, and merge into a new IC region.  That
        # is only needed for projects created outside the GUI (otherwise
        # we have already stored the IC regions).  Also would be nice
        # to offer a way to split compound regions.
        for ic in self.project.ics:

            d = ic.keyword_dict
            extent = [d.get('ic_'+k,None) for k in ('x_w', 'y_s', 'z_b',
                                                    'x_e', 'y_n', 'z_t')]
            extent = [0 if x is None else x.value for x in extent]
            #if any (x is None for x in extent):
            #    self.warn("initial condition %s: invalid extents %s" %
            #               (ic.ind, extent))
            #    continue
            for (region_name, data) in self.ics_region_dict.items():
                ext2 = [0 if x is None else x for x in
                        (data.get('from',[]) + data.get('to',[]))]
                if ext2 == extent:
                    if data.get('available', True):
                        self.ics_add_regions_1([region_name], indices=[ic.ind], autoselect=False)
                        break
            else:
                self.warn("initial condition %s: could not match defined region %s" %
                          (ic.ind, extent))


    def setup_ics_fluid_tab(self):
        #Fluid (tab)
        if self.fluid_solver_disabled:
            # we shouldn't be on this tab!
            return
        ui = self.ui.initial_conditions
        tw = ui.tablewidget_fluid_mass_fraction
        enabled = bool(self.fluid_species)
        if not enabled:
            tw.clearContents()
            tw.setRowCount(0)
        self.fixup_ics_table(tw)
        tw.setEnabled(enabled)

        if not self.ics_current_indices:
            # Nothing selected.  What can we do? (Clear out all lineedits?)
            return

        IC0 = self.ics_current_indices[0]
        # Note - value may not be consistent across grouped regions
        #  For now we're going to assume that it is, and just check
        #  first subregion of IC group

        #  If we can make this code generic enough perhaps someday it can
        # be autogenerated from SRS doc

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_IC',
                        'lineedit_%s_args_IC'):
                widget = getattr(ui, pat % key)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True):
            for pat in ('label_%s', 'label_%s_units',
                        'lineedit_keyword_%s_args_IC'):
                name = pat%key
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            if not enabled:
                get_widget(key).setText('') #?
                return
            args = mkargs(key, ic=IC0)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for IC in self.ics_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, ic=IC))
            get_widget(key).updateValue(key, val, args=args)

        #Define volume fraction
        # Specification always available
        # Sets keyword IC_EP_G(#)
        # DEFAULT value of 1.0
        # (terminology:  is ic_ep_g volume fraction or void fraction?)
        key = 'ic_ep_g'
        default = 1.0
        setup_key_widget(key, default)
        get_widget(key).setReadOnly(True)
        get_widget(key).setEnabled(False) # better way to indicate read-only?
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
        ui.groupbox_turbulence.setEnabled(enabled)

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
        ui.groupbox_fluid_advanced.setEnabled(enabled)

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
        self.ics_current_solid = self.P = P
        if P is None: # Nothing to do
            return

        if not self.ics_current_indices: # No region selected
            # TODO clear all widgets (?)
            return

        ui = self.ui.initial_conditions
        IC0 = self.ics_current_indices[0]

        # issues/121
        self.ics_set_volume_fraction_limit()

        # Generic!
        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_IC_P',
                        'lineedit_keyword_%s_args_IC',
                        'lineedit_%s_args_IC_P',
                        'lineedit_%s_args_IC'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_IC_P',
                         'lineedit_keyword_%s_args_IC',
                         'lineedit_%s_args_IC_P',
                         'lineedit_%s_args_IC'):
                name = pat%key
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            if not enabled:
                get_widget(key).setText('') #?
                return
            args = mkargs(key, ic=IC0, phase=P)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for IC in self.ics_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, ic=IC, phase=P))
            get_widget(key).updateValue(key, val, args=args)

        #Group tab inputs by equation type (e.g., momentum, energy, species).
        # Making the grouped inputs a 'collapsible list' may make navigation easier.
        #  (Note - collaspsing not implemented)

        #Define volume fraction (required)
        # Specification always available
        # Sets keyword IC_EP_S(#,#)
        # DEFAULT value of 0.0
        key = 'ic_ep_s'
        default = 0.0
        # Some input decks may or may not contain IC_EP_S keyword:
        #  Volume fraction is specified using the solids bulk density
        #    IC_EP_S(#,#) == IC_ROP_S(#,#) / IC_ROs(#)
        #    Solids density IC_ROs is determined by the solids density model
        # IC_ROs(#) = RO_S0(#)
        # IC_ROs(#) = -- FINISH LATER -
        # (above section of spec NOT IMPLEMENTED)

        #  Volume fraction may be inferred from IC_EP_G
        #    IC_EP_S(#,#) = 1.0 - IC_EP_G(#)
        #    Only valid for one solids phase (MMAX=1)
        # (note, this is handled in project_manager.load_project_file, see issues/142)
        #ic_ep_s = self.project.get_value(key, args=[IC0, P])
        #ic_ep_g = self.project.get_value('ic_ep_g', args=[IC0])
        #if ic_ep_s is None and ic_ep_g is not None and len(self.solids)==1:
        #    for IC in self.ics_current_indices:
        #        self.update_keyword(key, 1.0-ic_ep_g, args=[IC, P])

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

        #Define pressure (optional)
        # Specification only available for SOLIDS_MODEL(#)='TFM'
        # Sets keyword IC_P_STAR(#)
        # DEFAULT of 0.0
        # Common to all phases - Warn user if changed.
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
        # Use (m^2/sec^2) for solids granular energy units.
        # # Note:
        # # Some of the KT_TYPES also include a mass unit (kg) -
        # # but the default model (Lun) will have units of m^s/sec^2
        solids_model = self.project.get_value('solids_model', args=[P])
        kt_type = self.project.get_value('kt_type', default='ALGEBRAIC')
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
        # Sets keyword IC_X_S
        # DEFAULT - last defined species has mass fraction of 1.0
        # Error check: mass fractions must sum to one
        self.update_ics_solids_mass_fraction_table()

        enabled = (solids_model=='DEM' or bool(energy_eq))
        ui.groupbox_solids_advanced.setEnabled(enabled)

        #Advanced: Option to enable fitting DES particles to region
        # Option only available for DEM solids
        # Sets keyword: IC_DES_FIT_TO_REGION
        # Disabled [DEFAULT]
        enabled = (solids_model=='DEM')
        item = ui.checkbox_keyword_ic_des_fit_to_region_args_IC
        item.setEnabled(enabled)
        key = 'ic_des_fit_to_region'
        default = False
        val = self.project.get_value(key, args=[IC0])
        if val is None:
            val = default
            if enabled:
                for IC in self.ics_current_indices:
                    self.update_keyword(key, val, args=[IC])
        item.setChecked(val)
        # TODO: popup dialog to warn that this apples to all phases

        #Advanced: Define radiation coefficient
        # Specification only available when solving energy equations
        # Sets keyword IC_GAMA_RS(#,#)
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

        if not self.ics_current_indices:
            return # No selection
        IC0 = self.ics_current_indices[0]

        ui = self.ui.initial_conditions
        nscalar = self.project.get_value('nscalar', default=0)
        old_nscalar = getattr(ui, 'nscalar', None)
        ui.nscalar = nscalar

        key = 'ic_scalar'
        if nscalar == old_nscalar:
            # What do we have to do?  Just make the lineedits reflect current vals
            for i in range(1, nscalar+1):
                le = getattr(ui, "lineedit_ic_scalar_%s" % i, None)
                val = self.project.get_value(key, args=[IC0, i])
                if le:
                    le.updateValue(key, val)
            return

        page =  ui.page_scalar
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
            row = 0
            for i in range(1, nscalar+1):
                label = QLabel('Scalar %s' % i)
                self.add_tooltip(label, key)
                layout.addWidget(label, row, 0)
                le = LineEdit()
                le.key = key
                le.args = ['IC', i]
                self.add_tooltip(le, key)
                le.dtype = float
                le.default_value = 0.0
                self.project.register_widget(le, [key], ['IC', i])
                setattr(ui, 'lineedit_ic_scalar_%s'%i, le)

                self.add_tooltip(le, key)
                val = self.project.get_value(key, args=[IC0, i])
                if val is None:
                    val = 0.0
                    for IC in self.ics_current_indices:
                        self.update_keyword(key, val, args=[IC, i])
                le.updateValue(key, val)
                layout.addWidget(le, row, 1)
                row += 1

            if spacer:
                layout.addItem(spacer, row, 0)
