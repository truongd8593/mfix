# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import QLabel, QLineEdit, QPushButton, QGridLayout, QWidget
from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4

UserRole = QtCore.Qt.UserRole

from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit

from project import Equation # used for default phi_w

from tools.general import (set_item_noedit, set_item_enabled,
                           get_selected_row, widget_iter)

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

# Move to "constants"?
BC_TYPES = ['MI', 'PO', 'NSW', 'FSW', 'PSW', 'PI', 'MO']
BC_NAMES = ['Mass Inflow', 'Pressure Outflow', 'No Slip Wall',
            'Free Slip Wall', 'Partial Slip Wall',
            'Pressure Inflow', 'Mass Outflow']

(MASS_INFLOW, PRESSURE_OUTFLOW,
 NO_SLIP_WALL, FREE_SLIP_WALL, PARTIAL_SLIP_WALL,
 PRESSURE_INFLOW, MASS_OUTFLOW) = range(7)

(NO_FLUX, SPECIFIED_TEMPERATURE, SPECIFIED_FLUX, CONVECTIVE_FLUX) = range(4)
SPECIFIED_MASS_FRACTION = 1

DEFAULT_BC_TYPE = NO_SLIP_WALL

class BCS(object):
    #Boundary Conditions Task Pane Window: This section allows a user to define the boundary
    #conditions for the described model. This section relies on regions named in the Regions section.

    def init_bcs(self):
        ui = self.ui.boundary_conditions

        self.bcs = {} # key: index.  value: data dictionary for boundary cond
        self.bcs_current_indices = [] # List of BC indices
        self.bcs_current_regions = [] # And the names of the regions which define them
        self.bcs_region_dict = None

        # The top of the task pane is where users define/select BC regions
        # (see handle_bcs_region_selection)
        #
        #Icons to add/remove/duplicate boundary conditions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #a region to apply the boundary condition.
        ui.toolbutton_add.clicked.connect(self.bcs_show_regions_popup)
        ui.toolbutton_delete.clicked.connect(self.bcs_delete_regions)
        # TODO implement 'duplicate' (what does this do?)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        ui.tablewidget_regions.itemSelectionChanged.connect(self.handle_bcs_region_selection)

        self.bcs_current_tab = FLUID_TAB # #  If fluid is disabled, we will switch
        self.bcs_current_solid = self.P = None
        ui.pushbutton_fluid.pressed.connect(lambda: self.bcs_change_tab(FLUID_TAB,None))
        ui.pushbutton_scalar.pressed.connect(lambda: self.bcs_change_tab(SCALAR_TAB,None))

        # Trim width of "Fluid" and "Scalar" buttons, like we do for
        # dynamically-created "Solid #" buttons
        for b in (ui.pushbutton_fluid, ui.pushbutton_scalar):
            w = b.fontMetrics().boundingRect(b.text()).width() + 20
            b.setMaximumWidth(w)

        # Checkbox callbacks
        ui.checkbox_bc_jj_ps_args_BC.dtype = int
        ui.checkbox_bc_jj_ps_args_BC.clicked.connect(self.set_bc_jj_ps)

        # Combobox callbacks
        for name in ('fluid_energy_eq', 'fluid_species_eq',
                     'bc_jj_ps',
                     'solids_energy_eq', 'solids_granular_energy_eq'): # no solids_species_eq?
            item = getattr(ui, 'combobox_%s_type' % name)
            setter = getattr(self, 'set_bcs_%s_type' % name)
            item.currentIndexChanged.connect(setter)


    def bcs_show_regions_popup(self):
        # Users cannot select inapplicable regions.
        # BC regions must be planes or STLs (not volumes or points)
        # No region can define more than one boundary condition.
        ui = self.ui.boundary_conditions
        rp = self.regions_popup
        rp.clear()
        for (name,data) in self.bcs_region_dict.items():
            shape = data.get('type', '---')
            # Assume available if unmarked
            available = data.get('available', True) and (shape=='STL' or 'plane' in shape)
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.bcs_add_regions)
        rp.cancel.connect(self.bcs_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.scrollarea_detail,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup(boundary=True)


    def bcs_cancel_add(self):
        ui = self.ui.boundary_conditions

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.scrollarea_detail,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def bcs_add_regions(self):
        #Select boundary type
        # Selection is required
        # Available selections:
        #  Mass Inflow
        #    Plane regions set keyword BC_TYPE(#) to 'MI'
        #    STL regions set keyword BC_TYPE(#) to 'CG_MI'
        #  Pressure Outflow
        #    Plane regions set keyword BC_TYPE(#) to 'PO'
        #    STL regions set keyword BC_TYPE(#) to 'CG_PO'
        #  No Slip Wall
        #    Plane regions set keyword BC_TYPE(#) to 'NSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_NSW'
        #  Free Slip Wall
        #    Plane regions set keyword BC_TYPE(#) to 'FSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_FSW'
        #  Partial Slip Wall
        #    Plane regions set keyword BC_TYPE(#) to 'PSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_PSW'
        #  Pressure Inflow
        #    Plane regions set keyword BC_TYPE(#) to 'PI'
        #    Not available for STL regions
        # Mass Outflow
        #    Plane regions set keyword BC_TYPE(#) to 'MO'
        #    STL regions set keyword BC_TYPE(#) to 'CG_MO'
        # Specification always available
        # DEFAULT - No slip wall
        # Error check: mass fractions must sum to one

        # Interactively add regions to define BCs
        ui = self.ui.boundary_conditions
        rp = self.regions_popup
        self.bcs_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        bc_type = rp.combobox_boundary_type.currentIndex()
        if not selections:
            return
        self.bcs_add_regions_1(selections, bc_type) # Indices will be assigned
        self.bcs_setup_current_tab() # Update the widgets


    def bcs_add_regions_1(self, selections,
                          bc_type=DEFAULT_BC_TYPE, indices=None):
        # Used by both interactive and load-time add-region handlers
        if self.bcs_region_dict is None:
            self.bcs_region_dict = self.ui.regions.get_region_dict()

        ui = self.ui.boundary_conditions
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
            autoselect = True
        else: # loading file
            assert len(selections) == len(indices)
            autoselect = False

        for (i, region_name) in enumerate(selections):
            idx = indices[i]
            if idx is None:
                idx = self.bcs_find_index()
                indices[i] = idx
            self.bcs[idx] = {'region': region_name}

            region_data = self.bcs_region_dict.get(region_name)
            if region_data is None: # ?
                self.warn("no data for region %s" % region_name)
                continue
            self.bcs_set_region_keys(region_name, idx, region_data, bc_type)
            self.bcs_region_dict[region_name]['available'] = False # Mark as in-use

        item.setData(UserRole, (tuple(indices), tuple(selections)))
        tw.setItem(nrows, 0, item)

        item = make_item(BC_NAMES[bc_type])
        tw.setItem(nrows, 1, item)

        self.fixup_bcs_table(tw)

        if autoselect:
            tw.setCurrentCell(nrows, 0)


    def bcs_find_index(self):
        n = 1
        while n in self.bcs:
            n += 1
        return n


    def bcs_delete_regions(self):
        ui = self.ui.boundary_conditions
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            # TODO use keyword_args here instead of startswith
            if key.startswith('bc_') and args and args[0] in self.bcs_current_indices:
                self.unset_keyword(key, args=args)

        # TODO: fix any resulting holes in index sequence!

        for r in self.bcs_current_regions:
            if r in self.bcs_region_dict:
                self.bcs_region_dict[r]['available'] = True

        for i in self.bcs_current_indices:
            del self.bcs[i]

        self.bcs_current_regions = []
        self.bcs_current_indices = []

        tw.removeRow(row)
        self.fixup_bcs_table(tw)
        self.bcs_setup_current_tab()


    def handle_bcs_region_selection(self):
        ui = self.ui.boundary_conditions
        table = ui.tablewidget_regions
        row = get_selected_row(table)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = table.item(row,0).data(UserRole)
        self.bcs_current_indices = indices
        self.bcs_current_regions = regions
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.scrollarea_detail.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.scrollarea_detail):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.bcs_setup_current_tab() # reinitialize all widgets in current tab


    def fixup_bcs_table(self, tw, stretch_column=0):
        # TODO fix and unify all the fixup_*_table functions
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

        ui = self.ui.boundary_conditions
        ui.top_frame.setMaximumHeight(height+(40 if nrows==0 else 32))
        ui.top_frame.updateGeometry()

    def bcs_update_enabled(self):
        # If there are no solids, no scalar equations, and the fluid solver is disabled,
        regions = self.ui.regions.get_region_dict()
        nregions = len([r for r in regions.values()
                       if r.get('type')=='STL' or 'plane' in r.get('type')])
        disabled = (nregions==0
                    or (self.fluid_solver_disabled
                        and self.project.get_value('nscalar',default=0)==0
                        and len(self.solids)==0))
        self.find_navigation_tree_item("Boundary Conditions").setDisabled(disabled)


    def bcs_change_tab(self, tab, solid):
        ui = self.ui.boundary_conditions
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
            if solid == self.bcs_current_solid:
                return # Really nothing to do

            if solid > self.bcs_current_solid:
                dummy_label = ui.label_dummy_solids_L
                dummy_tab = SOLIDS_TAB_DUMMY_L
            else:
                dummy_label = ui.label_dummy_solids_R
                dummy_tab = SOLIDS_TAB_DUMMY_R

            pixmap = QPixmap(ui.page_solids.size())
            pixmap.fill() #fill bg with white
            ui.page_solids.render(pixmap, flags=QWidget.DrawChildren) # avoid rendering bg
            dummy_label.setPixmap(pixmap)
            ui.stackedwidget.setCurrentIndex(dummy_tab)

        self.bcs_current_tab = tab
        self.bcs_current_solid = self.P = solid if tab==SOLIDS_TAB else None
        self.bcs_setup_current_tab()

        # change stackedwidget contents
        self.animate_stacked_widget(
            ui.stackedwidget,
            ui.stackedwidget.currentIndex(),
            tab,
            direction='horizontal',
            line = ui.tab_underline,
            to_btn = ui.tab_layout.itemAtPosition(0, index),
            btn_layout = ui.tab_layout)


    def bcs_check_region_in_use(self, name):
        # Should we allow any change of region type?  eg. xy plane -> xz plane?
        #  Probably not
        return any(data.get('region')==name for data in self.bcs.values())


    def bcs_update_region(self, name, data):
        for (i,bc) in self.bcs.items():
            if bc.get('region') == name:
                self.bcs_set_region_keys(name, i, data)


    def bcs_change_region_name(self, old_name, new_name):
        ui = self.ui.boundary_conditions
        for (key, val) in self.bcs.items():
            if val.get('region') == old_name:
                self.bcs[key]['region'] = new_name
                tw = ui.tablewidget_regions
                for i in range(tw.rowCount()):
                    data = tw.item(i,0).data(UserRole)
                    indices, names = data
                    if old_name in names:
                        item = tw.item(i,0)
                        names = (new_name if n==old_name else n for n in names)
                        item.setData(UserRole, (indices, names))
                        item.setText('+'.join(names))
                        break
                break


    def bcs_set_region_keys(self, name, idx, data, bc_type=None):
        # Update the keys which define the box-shaped region the BC applies to
        if bc_type is not None:
            val = "%s%s" %('CG_' if data.get('type') == 'STL' else '',
                           BC_TYPES[bc_type])
            if val== 'CG_PI': # Shouldn't happen!
                self.error("Invalid bc_type %s" % val)
                return
            self.update_keyword('bc_type', val, args=[idx])

        for (key, val) in zip(('x_w', 'y_s', 'z_b', 'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            key = 'bc_' + key
            self.update_keyword(key, val, args=[idx])


    def reset_bcs(self):
        self.bcs.clear()
        self.bcs_current_indices = []
        self.bcs_current_regions = []
        self.bcs_region_dict = None
        ui = self.ui.boundary_conditions
        ui.tablewidget_regions.clearContents()
        ui.tablewidget_regions.setRowCount(0)
        # anything else to do here?


    def bcs_to_str(self):
        ui = self.ui.boundary_conditions
        tw = ui.tablewidget_regions
        data = [tw.item(i,0).data(UserRole) for i in range(tw.rowCount())]
        return JSONEncoder().encode(data)


    def bcs_regions_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)
        for (indices, regions) in data:
            # bc_type keyword should be set already when we call this
            self.bcs_add_regions_1(regions, indices)


    def setup_bcs(self):
        # Grab a fresh copy, may have been updated
        self.bcs_region_dict = self.ui.regions.get_region_dict()
        ui = self.ui.boundary_conditions
        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.bcs.items():
            region = data['region']
            if region in self.bcs_region_dict:
                self.bcs_region_dict[region]['available'] = False

        self.fixup_bcs_table(ui.tablewidget_regions)
        row = get_selected_row(ui.tablewidget_regions)
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.scrollarea_detail.setEnabled(enabled)
        # Tabs group boundary condition parameters for phases and additional equations. Tabs are
        # unavailable if no input is required from the user.
        #
        #Fluid tab - Unavailable if the fluid phase was disabled.
        b = ui.pushbutton_fluid
        b.setText(self.fluid_phase_name)
        b.setEnabled(not self.fluid_solver_disabled)
        if self.fluid_solver_disabled:
            if self.bcs_current_tab == 0: # Don't stay on disabled tab
                self.bcs_change_tab(*(SOLIDS_TAB, 1) if self.solids else (SCALAR,None))
        font = b.font()
        font.setBold(self.bcs_current_tab == 0)
        b.setFont(font)

        #  Each solid phase will have its own tab. The tab name should be the name of the solid
        # (Could do this only on solid name change)
        n_cols = ui.tab_layout.columnCount()
        # Clear out the old ones
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
        for (i, solid_name) in enumerate(self.solids.keys(),1):
            b = QPushButton(text=solid_name)
            w = b.fontMetrics().boundingRect(solid_name).width() + 20
            b.setMaximumWidth(w)
            b.setFlat(True)
            font = b.font()
            font.setBold(self.bcs_current_tab==SOLIDS_TAB and i==self.bcs_current_solid)
            b.setFont(font)
            b.pressed.connect(lambda i=i: self.bcs_change_tab(SOLIDS_TAB, i))
            ui.tab_layout.addWidget(b, 0, i)
        # Don't stay on disabled tab TODO
        # if self.bcs_current_tab == 1 and ...

        #Scalar (tab) - Tab only available if scalar equations are solved
        # Move the 'Scalar' button to the right of all solids, if needed
        b = ui.pushbutton_scalar
        font = b.font()
        font.setBold(self.bcs_current_tab==SCALAR_TAB)
        b.setFont(font)
        nscalar = self.project.get_value('nscalar', default=0)
        enabled = (nscalar > 0)
        b.setEnabled(enabled)
        if len(self.solids) > 0:
            ui.tab_layout.removeWidget(b)
            ui.tab_layout.addWidget(b, 0, 1+len(self.solids))
        # Don't stay on a disabled tab TODO
        # if self.bcs_current_tab == 2 and nscalar == 0:
        #
        self.P = self.bcs_current_solid
        self.bcs_setup_current_tab()


    def bcs_setup_current_tab(self):
        if self.bcs_current_tab == FLUID_TAB:
            self.setup_bcs_fluid_tab()
        elif self.bcs_current_tab == SOLIDS_TAB:
            self.setup_bcs_solids_tab(self.bcs_current_solid)
        elif self.bcs_current_tab == SCALAR_TAB:
            self.setup_bcs_scalar_tab()


    def bcs_set_volume_fraction_limit(self):
        pass

    def handle_bcs_volume_fraction(self, widget, val, args):
        pass

    def update_bcs_fluid_mass_fraction_table(self):
        pass

    def handle_bcs_fluid_mass_fraction(self, widget, value_dict, args):
        pass
    def update_bcs_fluid_mass_fraction_total(self):
        pass
    def update_bcs_solids_mass_fraction_table(self):
        pass

    def handle_bcs_solids_mass_fraction(self, widget, value_dict, args):
        pass
    def update_bcs_solids_mass_fraction_total(self):
        pass

    def bcs_extract_regions(self):
        if self.bcs:
            # We assume that bc regions have been initialized correctly
            # from mfix_gui_comments.
            # TODO: verify that there is an BC region for each BC
            return

        if self.bcs_region_dict is None:
            self.bcs_region_dict = self.ui.regions.get_region_dict()

        # TODO: if we wanted to be fancy, we could find regions where
        # BC values matched, and merge into a new BC region.  That
        # is only needed for projects created outside the GUI (otherwise
        # we have already stored the BC regions).  Also would be nice
        # to offer a way to split compound regions.
        for bc in self.project.bcs:
            d = bc.keyword_dict
            extent = [d.get(k, None) for k in ('bc_x_w', 'bc_y_s', 'bc_z_b',
                                              'bc_x_e', 'bc_y_n', 'bc_z_t')]
            # should we distinguish 0 from unset?  in the region_dict we get
            #  from the regions_widget, the values are 0, while in the project,
            #  keywords are simply unset (value None) rather than set to 0

            extent = [0.0 if x is None else x.value for x in extent]

            bc_type = d.get('bc_type')
            if bc_type is None:
                self.error("No type for boundary condition %s" % bc.ind)
                continue
            bc_type = bc_type.value

            is_stl = (bc_type.startswith('CG_'))
            if is_stl:
                bc_type = bc_type[3:]

            # Check dimensionality?
            #if any (x is None for x in extent):
            #    self.warn("boundary condition %s: invalid extents %s" %
            #               (bc.ind, extent))
            #    continue
            for (region_name, data) in self.bcs_region_dict.items():

                ext2 = data.get('from',[]) + data.get('to',[])

                # TODO this only works for a single STL
                if (is_stl and data.get('type')=='STL'
                    or not is_stl and ext2==extent):

                    if data.get('available', True):
                        if bc_type is None:
                            self.warn("no bc_type for region %s" % bc.ind)
                        if bc_type not in BC_TYPES:
                            self.warn("invalid bc_type %s for region %s" % (bc_type, bc.ind))
                        else:
                            self.bcs_add_regions_1([region_name], BC_TYPES.index(bc_type), [bc.ind])
                            break
            else:
                self.warn("boundary condition %s: could not match defined region %s" %
                          (bc.ind, extent))


    def set_bcs_fluid_energy_eq_type(self, eq_type):
        if not self.bcs_current_indices:
            return
        #Select energy equation boundary type:
        # Available selections:
        #  No-Flux (adiabatic) [DEFAULT]
        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_T_G(#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_T_G(#) to 0.0
            c = 0.0
            #    Sets keyword BC_TW_G(#) to UNDEFINED
            tw = None
        #  Specified Temperature
        elif eq_type == SPECIFIED_TEMPERATURE:
            #    Sets keyword BC_HW_T_G(#) to UNDEFINED
            hw = None
            #    Sets keyword BC_C_T_G(#) to 0.0
            c = 0.0
            #    Requires BC_TW_G(#)
            tw = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_T_G(#) to 0.0
            hw = 0.0
            #    Requires BC_C_T_G(#)
            c = True
            #    Sets keyword BC_TW_G(#) to UNDEFINED
            tw = None
        elif eq_type == CONVECTIVE_FLUX:
            #    Requires BC_HW_T_G(#)
            hw = True
            #    Sets keyword BC_C_T_G(#) to 0.0
            c = 0.0
            #    Requires BC_TW_G(#)
            tw = True
        else:
            self.error("Invalid fluid energy_eq type %s" % eq_type)
            return

        for BC in self.bcs_current_indices:
            self.bcs[BC]['fluid_energy_eq_type'] = eq_type
            for (key, val) in (('bc_hw_t_g', hw), ('bc_c_t_g', c), ('bc_tw_g', tw)):
                if val is True:
                    pass
                else:
                    self.update_keyword(key, val, args=[BC])

        self.setup_bcs_fluid_tab()


    def set_bcs_fluid_species_eq_type(self, eq_type):
        if not self.bcs_current_indices:
            return
        #Select species equation boundary type:
        # Selection only available when solving species equations
        # Available selections:
        #  No-Flux [DEFAULT]
        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_X_G(#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_X_G(#) to 0.0
            c = 0
            #    Sets keyword BC_XW_G(#) to UNDEFINED
            xw = None
        #  Specified Mass Fraction
        elif eq_type == SPECIFIED_MASS_FRACTION:
            #    Sets keyword BC_XW_T_G(#) to UNDEFINED # should be BC_H_X_G ?
            hw = None
            #    Sets keyword BC_C_X_G(#) to 0.0
            c = 0.0
            #    Requires BC_XW_G(#)
            xw = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_X_G(#) to 0.0
            hw = 0.0
            #    Requires BC_C_X_G(#)
            c = True
            #    Sets keyword BC_XW_G(#) to UNDEFINED
            xw = None
        #  Convective Flux
        elif eq_type == CONVECTIVE_FLUX:
            #    Requires BC_HW_X_G(#)
            hw = True
            #    Sets keyword BC_C_X_G(#) to 0.0
            c = 0.0
            #    Requires BC_XW_G(#)
            xw = True
        else:
            self.error("Invalid fluid species_eq type %s" % eq_type)
            return

        for BC in self.bcs_current_indices:
            self.bcs[BC]['fluid_species_eq_type'] = eq_type
            for (key, val) in (('bc_hw_x_g', hw), ('bc_c_x_g', c), ('bc_xw_g', xw)):
                if val is True:
                    pass # 'required'
                else:
                    self.update_keyword(key, val, args=[BC])

        self.setup_bcs_fluid_tab()


    def set_bc_jj_ps(self, val):
        if not self.bcs_current_indices:
            return
        key = 'bc_jj_ps'
        for BC in self.bcs_current_indices:
            self.update_keyword(key, val, args=BC)
        self.setup_bcs_solids_tab(self.bcs_current_solid)


    def set_bcs_bc_jj_ps_type(self, val):
        # Available selections:
        #  Default Jackson-Johnson BC [DEFAULT]
        #    Sets keyword BC_JJ_M to .FALSE.
        #    Sets keyword JENKINS to .FALSE.
        #  Variable specularity coefficient
        #    Sets keyword BC_JJ_M to .TRUE.
        #    Sets keyword JENKINS to .FALSE.
        #  Jenkins small frictional boundary
        #    Sets keyword BC_JJ_M to .FALSE.
        #    Sets keyword JENKINS to .TRUE.
        if val < 0 or val > 2:
            self.warn("Invalid bc_jj_ps_type % val")
            return
        jenkins, bc_jj_m = int(val/2), val%2
        self.update_keyword('bc_jj_m', bool(bc_jj_m))
        self.update_keyword('jenkins', bool(jenkins))
        self.setup_bcs_solids_tab(self.bcs_current_solid)


    def set_bcs_solids_energy_eq_type(self, eq_type):
        if not self.bcs_current_indices:
            return
        P = self.bcs_current_solid
        if P is None:
            return
        # Available selections:
        #  No-Flux (adiabatic) [DEFAULT]
        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_T_S(#,#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_T_S(#,#) to 0.0
            c = 0.0
            #    Sets keyword BC_TW_S(#,#) to UNDEFINED
            tw = 0.0
        #  Specified Temperature
        elif eq_type == SPECIFIED_TEMPERATURE:
            #    Sets keyword BC_HW_T_S(#,#) to UNDEFINED
            hw = None
            #    Sets keyword BC_C_T_S(#,#) to 0.0
            c = 0.0
            #    Requires BC_TW_S(#,#)
            tw = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_T_S(#,#) to 0.0
            hw = 0.0
            #    Requires BC_C_T_S(#)
            c = True
            #    Sets keyword BC_TW_S(#,#) to UNDEFINED
            tw = None
        #  Convective Flux
        elif eq_type == CONVECTIVE_FLUX:
            #    Requires BC_HW_T_S(#,#)
            hw = True
            #    Sets keyword BC_C_T_S(#,#) to 0.0
            c = 0.0
            #    Requires BC_TW_S(#,#)
            tw = True
        else:
            self.error("Invalid solid %s energy_eq type %s" % (P, eq_type))
            return

        for BC in self.bcs_current_indices:
            # FIXME this won't survive index remapping when solids are deleted!
            self.bcs[BC]['solid_%s_energy_eq_type'%P] = eq_type
            for (key, val) in (('bc_hw_t_s', hw), ('bc_c_t_s', c), ('bc_tw_s', tw)):
                if val is True:
                    pass # 'required'
                else:
                    self.update_keyword(key, val, args=[BC,P])

        self.setup_bcs_solids_tab(self.bcs_current_solid)


    def set_bcs_solids_granular_energy_eq_type(self, eq_type):
        if not self.bcs_current_indices:
            return
        P = self.bcs_current_solid
        if P is None:
            return
        # Available selections:
        #  No-Flux [DEFAULT]

        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_THETA_M(#,#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_THETA_M (#,#) to 0.0
            c = 0.0
            #    Sets keyword BC_THETAW_M(#,#) to UNDEFINED
            theta = None
        #  Specified Temperature
        elif eq_type == SPECIFIED_TEMPERATURE:
            #    Sets keyword BC_HW_THETA_M(#,#) to UNDEFINED
            hw = None
            #    Sets keyword BC_C_THETA_M(#,#) to 0.0
            c = 0.0
            #    Requires BC_THETAW_M(#,#)
            theta = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_THETA_M(#,#) to 0.0
            hw = 0.0
            #    Requires BC_C_THETA_M(#)
            c = True
            #    Sets keyword BC_THETAW_M(#,#) to UNDEFINED
            theta = None
        else:
            self.error("Invalid solid %s granualar energy_eq type %s" % (P, eq_type))
            return

        for BC in self.bcs_current_indices:
            self.bcs[BC]['solid_%s_granular_energy_eq_type'%P] = eq_type
            for (key, val) in (('bc_hw_theta_m', hw), ('bc_c_theta_m', c), ('bc_thetaw_m', theta)):
                if val is True:
                    pass # 'required'
                else:
                    self.update_keyword(key, val, args=[BC,P])

        self.setup_bcs_solids_tab(self.bcs_current_solid)


    def setup_bcs_fluid_tab(self):
        #Fluid (tab)
        if self.fluid_solver_disabled:
            # we shouldn't be on this tab!
            return
        ui = self.ui.boundary_conditions
        # Subtask Pane Tab for Wall type (NSW, FSW, PSW, CG_NSW, CG_FSW, CG_PSW) Boundary Condition Regions
        #
        if not self.bcs_current_indices:
            return # Nothing selected.  What can we do? (Clear out all lineedits?)

        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        # as generic as possible - see comments in ics.py

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC',
                        'lineedit_%s_args_BC'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for name in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC'):
                item = getattr(ui, name%(key+suffix), None)
                if item:
                    item.setEnabled(enabled)

            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default

            for BC in self.bcs_current_indices:
                self.update_keyword(key, val, args=[BC])
            get_widget(key+suffix).updateValue(key, val, args=[BC0])

        #    Define transfer coefficient
        # Specification only available with PSW
        # Sets keyword BC_HW_G(#)
        # DEFAULT value of 0.0
        #
        #    Define Wall U-velocity
        # Specification only available with PSW
        # Sets keyword BC_UW_G(#)
        # DEFAULT value of 0.0
        #
        #    Define Wall V-velocity
        # Specification only available with PSW
        # Sets keyword BC_VW_G(#)
        # DEFAULT value of 0.0
        #
        #    Define Wall W-velocity
        # Specification only available with PSW
        # Sets keyword BC_WW_G(#)
        # DEFAULT value of 0.0
        enabled = (bc_type=='PSW')
        default = 0.0 if enabled else None
        for c in 'huvw':
            key = 'bc_%sw_g' % c
            setup_key_widget(key, default, enabled)
        # Enable/disable entire groupbox
        ui.groupbox_fluid_momentum_eq.setEnabled(enabled) # Hide it???

        #Select energy equation boundary type:
        # Selection only available when solving energy equations
        energy_eq = self.project.get_value('energy_eq', default=True)
        enabled = bool(energy_eq)
        ui.groupbox_fluid_energy_eq.setEnabled(enabled)
        eq_type = None
        if enabled:
            eq_type = self.bcs[BC0].get('fluid_energy_eq_type')
            if eq_type is None: # Attempt to infer from keywords
                hw = self.project.get_value('bc_hw_t_g', args=[BC0])
                c = self.project.get_value('bc_c_t_g', args=[BC0])
                tw = self.project.get_value('bc_tw_g', args=[BC0])
                if hw==0.0 and c==0.0 and tw is None:
                    eq_type = NO_FLUX
                elif hw is None and c==0.0 and tw is not None:
                    eq_type = SPECIFIED_TEMPERATURE
                elif hw==0.0 and c!=0.0 and tw is None:
                    eq_type = SPECIFIED_FLUX
                elif hw is not None and c==0.0 and tw is not None:
                    eq_type = CONVECTIVE_FLUX
                else:
                    self.error("Cannot determine type for fluid energy boundary equation %s" % BC0)
                    eq_type = NO_FLUX # Default
                    self.set_bcs_fluid_energy_eq_type(eq_type)

            if eq_type is not None:
                ui.combobox_fluid_energy_eq_type.setCurrentIndex(eq_type)

        if energy_eq:
            #Define wall temperature
            # Specification only available with 'Specified Temperature' BC type
            # Sets keyword BC_TW_G(#)
            # DEFAULT value of 293.15
            enabled = (eq_type==SPECIFIED_TEMPERATURE)

            key = 'bc_tw_g'
            default = 293.15 if enabled else None
            self.print_internal("%s %s %s" % (key, default, enabled))
            setup_key_widget(key, default, enabled)
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_tw_g_2_args_BC.setText('')
            else:
                ui.lineedit_keyword_bc_tw_g_args_BC.setText('')
            #Define constant flux
            # Specification only available with 'Specified Flux' BC type
            # Sets keyword BC_C_T_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_FLUX)
            key = 'bc_c_t_g'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define transfer coefficient
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_HW_T_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_hw_t_g'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define free stream temperature
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_TW_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_tw_g'
            default = 0.0 if enabled else None
            setup_key_widget(key, default, enabled, suffix='_2')
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_tw_g_args_BC.setText('')
            else:
                ui.lineedit_keyword_bc_tw_g_2_args_BC.setText('')

        #Select species equation boundary type:
        # Selection only available when solving species equations
        species_eq = self.project.get_value('species_eq', default=True, args=[0])
        enabled = bool(species_eq)
        ui.groupbox_fluid_species_eq.setEnabled(enabled)
        eq_type = None
        if enabled:
            eq_type = self.bcs[BC0].get('fluid_species_eq_type')
            if eq_type is None: # Attempt to infer from keywords
                hw = self.project.get_value('bc_hw_x_g', args=[BC0])
                c = self.project.get_value('bc_c_x_g', args=[BC0])
                xw = self.project.get_value('bc_xw_g', args=[BC0])
                if hw==0.0 and c==0.0 and xw is None:
                    eq_type = NO_FLUX
                elif hw is None and c==0.0 and xw is not None:
                    eq_type = SPECIFIED_MASS_FRACTION
                elif hw==0.0 and c!=0.0 and xw is None:
                    eq_type = SPECIFIED_FLUX
                elif hw is not None and c==0.0 and xw is not None:
                    eq_type = CONVECTIVE_FLUX
                else:
                    self.error("Cannot determine type for fluid species boundary equation %s" % BC0)
                    eq_type = NO_FLUX # Default
                    self.set_bcs_fluid_species_eq_type(eq_type)

            if eq_type is not None:
                ui.combobox_fluid_species_eq_type.setCurrentIndex(eq_type)

        if species_eq:
            #Define wall mass fraction
            # Specification only available with 'Specified Mass Fraction' BC type
            # Sets keyword BC_XW_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_MASS_FRACTION)
            key = 'bc_xw_g'
            default = 0.0 if enabled else None
            setup_key_widget(key, default, enabled)
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_xw_g_2_args_BC.setText('')
            else:
                ui.lineedit_keyword_bc_xw_g_args_BC.setText('')

            #Define constant flux
            # Specification only available with 'Specified Flux' BC type
            # Sets keyword BC_C_X_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_FLUX)
            key = 'bc_c_x_g'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define transfer coefficient
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_HW_X_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_hw_x_g'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define free stream mass fraction
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_XW_G(#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_xw_g'
            default = 0.0 if enabled else None
            setup_key_widget(key, default, enabled, suffix='_2')
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_xw_g_args_BC.setText('')
            else:
                ui.lineedit_keyword_bc_xw_g_2_args_BC.setText('')


    def setup_bcs_solids_tab(self, P):
        #Solids-# (tab) - (Replace with phase name defined by the user)
        # Note, solids phases are numbered 1-N
        self.bcs_current_solid = self.P = P
        if P is None: # Nothing to do
            return

        if not self.bcs_current_indices: # No region selected
            # TODO clear all widgets (?)
            return

        ui = self.ui.boundary_conditions
        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC_P',
                        'lineedit_keyword_%s_args_BC',
                        'lineedit_%s_args_BC_P',
                        'lineedit_%s_args_BC',
                        'lineedit_keyword_%s',
                        'lineedit_%s'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for name in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC_P',
                         'lineedit_keyword_%s_args_BC',
                         'lineedit_%s_args_BC_P',
                         'lineedit_%s_args_BC',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                item = getattr(ui, name%(key+suffix), None)
                if item:
                    item.setEnabled(enabled)

            no_p_keys =  () # use keyword_args db here
            args = [BC0] if key in no_p_keys else [BC0,P]
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
            for BC in self.ics_current_indices:
                self.update_keyword(key, val, args=[BC] if key in no_p_keys
                                    else [BC,P])
            get_widget(key+suffix).updateValue(key, val, args=args)


        #    Enable Jackson-Johnson partial slip boundary
        # Disabled (0.0) for CARTESIAN_GRID = .TRUE.
        # Disabled (0.0) for KT_TYPE = 'ALGEBRAIC'
        # Disabled (0.0) for KT_TYPE = 'GHD_2007'
        # Sets keyword BC_JJ_PS(#)
        # DEFAULT value of 1.0 when not disabled
        #
        # Note, according to doc it's an int, not float
        cartesian_grid = bool(self.project.get_value('cartesian_grid', default=False))
        kt_type = self.project.get_value('kt_type', default='ALGEBRAIC')
        enabled = not cartesian_grid and kt_type not in ('ALGEBRAIC', 'GHD_2007')
        key = 'bc_jj_ps'
        default = 1 if enabled else 0
        widget = ui.checkbox_bc_jj_ps_args_BC
        val = self.project.get_value(key, default, args=BC0) if enabled else 0
        if val and not enabled:
            val = 0
            for BC in self.bcs_current_indices:
                self.update_keyword('bc_jj_ps', 0, args=[BC])
        if enabled: # Default value of 1.0 when not disabled
            for BC in self.bcs_current_indices:
                self.update_keyword('bc_jj_ps', val, args=[BC])
        widget.setChecked(val)

        widget.setEnabled(enabled)
        # TODO should this also depend on energy_eq ?
        ui.groupbox_solids_momentum_eq.setEnabled(enabled or bc_type=='PSW')

        #    Select type of Jackson and Johnson BC:
        # Selection only available BC_JJ_PS(#) = 1.0
        bc_jj_ps = int(self.project.get_value('bc_jj_ps', default=0, args=BC0))
        enabled = (bc_jj_ps==1)
        for widget in (ui.label_bc_jj_ps_type, ui.combobox_bc_jj_ps_type):
            widget.setEnabled(enabled)

        # Available selections:
        #  Default Jackson-Johnson BC [DEFAULT]
        #    Sets keyword BC_JJ_M to .FALSE.
        #    Sets keyword JENKINS to .FALSE.
        #  Variable specularity coefficient
        #    Sets keyword BC_JJ_M to .TRUE.
        #    Sets keyword JENKINS to .FALSE.
        #  Jenkins small frictional boundary
        #    Sets keyword BC_JJ_M to .FALSE.
        #    Sets keyword JENKINS to .TRUE.
        bc_jj_m = self.project.get_value('bc_jj_m', default=False)
        jenkins = self.project.get_value('jenkins', default=False)
        val = 2*jenkins + bc_jj_m
        if val > 2:
            self.warn("boundary condition %s: invalid combination of bc_jj_m=%s, jenkins=%s" %
                      (BC0, bc_jj_m, jenkins))
        else:
            widget.setCurrentIndex(val)

        #Define restitution coefficient
        # Specification only available with BC_JJ_PS(#) = 1.0
        # Sets keyword E_W
        # DEFAULT value of 1.0
        # Required when available # TODO implement 'required'
        enabled = (bc_jj_ps==1)
        key = 'e_w'
        default = 1.0
        setup_key_widget(key, default, enabled)

        #Define specularity coefficient
        # Specification only available with BC_JJ_PS(#)=1.0 and JENKINS=.FALSE.
        # Sets keyword PHIP
        # DEFAULT value of 0.6
        # Required when available
        jenkins = self.project.get_value('jenkins', default=False)
        enabled = bc_jj_ps==1 and jenkins==False
        key = 'phip'
        default = 0.6
        setup_key_widget(key, default, enabled)

        #Define specularity coefficient at zero slip
        # Specification only available with BC_JJ_PS(#)=1.0 and BC_JJ_M=.TRUE.
        # Sets keyword PHIP0
        # DEFAULT -blank
        # Optional when available
        bc_jj_m = self.project.get_value('bc_jj_m')
        enabled = bc_jj_ps==1 and bc_jj_m
        key = 'phip0'
        default = None
        setup_key_widget(key, default, enabled)

        #Define angle of internal friction
        # Specification only available with BC_JJ_PS(#)=1.0 and (JENKINS=.TRUE. FRICTION_MODEL=SRIVASTAVA) # ??? or ?
        friction_model = self.project.get_value('friction_model')
        enabled = (bc_jj_ps==1) and (jenkins or friction_model=='SRIVASTAVA') #? correct reading of srs?
        # DEFAULT value of 11.31 = atan(0.2) * (180/pi)
        default = Equation('atan(0.2) * (180/pi)') # 11.31
        # Sets keyword PHI_W # ?
        key = 'phi_w'
        # Required when available
        setup_key_widget(key, default, enabled)

        #Define transfer coefficient
        # Specification only available with PSW
        # Sets keyword BC_HW_S(#,#)
        # DEFAULT value of 0.0
        enabled = (bc_type=='PSW')
        key = 'bc_hw_s'
        default = 0.0 if enabled else None
        setup_key_widget(key, default, enabled)

        #Define Wall U-velocity
        # Specification only available with PSW or BC_JJ_PS(#) = 1.0
        # Sets keyword BC_UW_S(#,#)
        # DEFAULT value of 0.0
        #
        #Define Wall V-velocity
        # Specification only available with PSW or BC_JJ_PS(#) = 1.0
        # Sets keyword BC_VW_S(#,#)
        # DEFAULT value of 0.0
        #
        #Define Wall W-velocity
        # Specification only available with PSW or BC_JJ_PS(#) = 1.0
        # Sets keyword BC_WW_S(#,#)
        # DEFAULT value of 0.0
        enabled = (bc_type=='PSW') or (bc_jj_ps==1)
        default = 0 if enabled else None
        for c in 'uvw':
            key = 'bc_%sw_s' % c
            setup_key_widget(key, default, enabled)

        #Select energy equation boundary type:
        # Selection only available when solving energy equations
        energy_eq = self.project.get_value('energy_eq', default=True)
        enabled = bool(energy_eq)
        ui.groupbox_solids_energy_eq.setEnabled(enabled)
        eq_type = None
        if enabled:
            eq_type = self.bcs[BC0].get('solid_%s_energy_eq_type' % P)
            if eq_type is None: # Attempt to infer from keywords
                hw = self.project.get_value('bc_hw_t_s', args=[BC0,P])
                c =  self.project.get_value('bc_c_t_s', args=[BC0,P])
                tw = self.project.get_value('bc_tw_s', args=[BC0,P])
                # Available selections:
                #  No-Flux (adiabatic) [DEFAULT]
                #    Sets keyword BC_HW_T_S(#,#) to 0.0
                #    Sets keyword BC_C_T_S(#,#) to 0.0
                #    Sets keyword BC_TW_S(#,#) to UNDEFINED
                if hw==0.0 and c==0.0 and tw==None:
                    eq_type = NO_FLUX
                #  Specified Temperature
                #    Sets keyword BC_HW_T_S(#,#) to UNDEFINED
                #    Sets keyword BC_C_T_S(#,#) to 0.0
                #    Requires BC_TW_S(#,#)
                elif hw is None and c==0.0 and tw is not None:
                    eq_type = SPECIFIED_TEMPERATURE
                #  Specified Flux
                #    Sets keyword BC_HW_T_S(#,#) to 0.0
                #    Requires BC_C_T_S(#)
                #    Sets keyword BC_TW_S(#,#) to UNDEFINED
                elif hw==0.0 and c is not None and tw is None:
                    eq_type = SPECIFIED_FLUX
                #  Convective Flux
                #    Requires BC_HW_T_S(#,#)
                #    Sets keyword BC_C_T_S(#,#) to 0.0
                #    Requires BC_TW_S(#,#)
                elif hw is not None and c==0.0 and tw is None:
                    eq_type = CONVECTIVE_FLUX
                else:
                    self.error("Cannot determine type for solid %s energy boundary equation %s" % (P,BC0))
                    eq_type = NO_FLUX # Default
                    self.set_bcs_solids_energy_eq_type(eq_type)

            if eq_type is not None:
                ui.combobox_solids_energy_eq_type.setCurrentIndex(eq_type)

        if energy_eq:
            #Define wall temperature
            # Specification only available with 'Specified Temperature' BC type
            # Sets keyword BC_TW_S(#,#)
            # DEFAULT value of 293.15
            enabled = (eq_type==SPECIFIED_TEMPERATURE)
            key = 'bc_tw_s'
            default = 293.15 if enabled else None
            setup_key_widget(key, default, enabled)
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_tw_s_args_BC_P.setText('')
            else:
                ui.lineedit_keyword_bc_tw_s_args_BC_P.setText('')

            #Define constant flux
            # Specification only available with 'Specified Flux' BC type
            # Sets keyword BC_C_T_S(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_FLUX)
            key = 'bc_c_t_s'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define transfer coefficient
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_HW_T_S(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_hw_t_s'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define free stream temperature
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_TW_S(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_tw_s'
            default = 0.0 if enabled else None
            setup_key_widget(key, default, enabled, suffix='_2')
            # Hack to prevent dup. display
            if enabled:
                ui.lineedit_keyword_bc_tw_s_args_BC_P.setText('')
            else:
                ui.lineedit_keyword_bc_tw_s_2_args_BC_P.setText('')

        # FIXME should this also depend on 'granular_energy' keyword ?
        #Select granular energy equation boundary type:
        # Selection only available with BC_JJ_PS(#)=0.0 and KT_TYPE /= 'ALGEBRAIC'
        enabled = (bc_jj_ps==0) and (kt_type != 'ALGEBRAIC')
        ui.groupbox_solids_granular_energy_eq.setEnabled(enabled)
        eq_type = None
        if enabled:
            eq_type = self.bcs[BC0].get('solid_%d_granular_energy_eq_type'%P)
            if eq_type is None:
                hw = self.project.get_value('bc_hw_theta_m', args=[BC0,P])
                c =  self.project.get_value('bc_c_theta_m', args=[BC0,P])
                theta =  self.project.get_value('bc_thetaw_m', args=[BC0,P])
                #  No-Flux [DEFAULT]
                #    Sets keyword BC_HW_THETA_M(#,#) to 0.0
                #    Sets keyword BC_C_THETA_M (#,#) to 0.0
                #    Sets keyword BC_THETAW_M(#,#) to UNDEFINED
                if (hw==0.0) and (c==0.0) and (theta is None):
                    eq_type = NO_FLUX
                #  Specified Temperature
                #    Sets keyword BC_HW_THETA_M(#,#) to UNDEFINED
                #    Sets keyword BC_C_THETA_M(#,#) to 0.0
                #    Requires BC_THETAW_M(#,#)
                elif (hw is None) and (c==0.0) and theta is not None:
                    eq_type == SPECIFIED_TEMPERATURE
                #  Specified Flux
                #    Sets keyword BC_HW_THETA_M(#,#) to 0.0
                #    Requires BC_C_THETA_M(#)
                #    Sets keyword BC_THETAW_M(#,#) to UNDEFINED
                elif (hw==0.0) and (c is not None) and (theta is None):
                    eq_type = SPECIFIED_FLUX
                else:
                    self.error("Cannot determine type for solid %s granular energy boundary equation %s" % (P,BC0))
                    eq_type = NO_FLUX # Default
                    self.set_bcs_solids_granular_energy_eq_type(eq_type)

            if eq_type is not None:
                ui.combobox_solids_granular_energy_eq_type.setCurrentIndex(eq_type)

        if enabled:
            #Define granular temperature
            # Specification only available with 'Specified Temperature' BC type
            # Sets keyword BC_THETAW_M(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_TEMPERATURE)
            key = 'bc_thetaw_m'
            default = 0.0
            setup_key_widget(key, default, enabled)

            #Define constant flux
            # Specification only available with 'Specified Flux' BC type
            # Sets keyword BC_C_THETA_M(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_FLUX)
            key = 'bc_c_theta_m'


        #When solving solids species equations:
        #    Set keyword BC_HW_X_S(#,#,#) to 0.0
        #    Set keyword BC_C_X_S(#,#,#) to 0.0
        #    Set keyword BC_XW_S(#,#,#) to UNDEFINED
        #
        # TODO:  Where do we set this?
        #  1:  when toggling solids species eq
        #  2:  when adding a new species
        #  3:  when adding a new BC

    def setup_bcs_scalar_tab(self):
        pass

"""
Mockup of Task pane for specifying the Solid-# properties for WALL boundary condition regions.

Scalar (tab) - Tab only available if scalar equations are solved
#    Select scalar boundary type:
# Available selections:
#  No-Flux [DEFAULT]
#    Sets keyword BC_HW_SCALAR(#,#) to 0.0
#    Sets keyword BC_C_SCALAR(#,#) to 0.0
#    Sets keyword BC_SCALARW(#,#) to UNDEFINED
#  Specified Temperature
#    Sets keyword BC_HW_T_S(#,#) to UNDEFINED
#    Sets keyword BC_C_SCALAR (#,#) to 0.0
#    Requires BC_SCALARW (#,#)
#  Specified Flux
#    Sets keyword BC_HW_T_S(#,#) to 0.0
#    Requires BC_C_SCALAR (#)
#    Sets keyword BC_SCALARW (#,#) to UNDEFINED
#  Convective Flux
#    Requires BC_HW_T_S(#,#)
#    Sets keyword BC_C_SCALAR (#,#) to 0.0
#    Requires BC_SCALARW (#,#)
#    Define wall temperature
# Specification only available with 'Specified Temperature' BC type
# Sets keyword BC_SCALARW (#,#)
# DEFAULT value of 0.0
#    Define constant flux
# Specification only available with 'Specified Flux' BC type
# Sets keyword BC_C_SCALAR (#,#)
# DEFAULT value of 0.0
#    Define transfer coefficient
# Specification only available with 'Convective Flux' BC type
# Sets keyword BC_HW_SCALAR(#,#)
# DEFAULT value of 0.0
#    Define free stream temperature

# Specification only available with 'Convective Flux' BC type
# Sets keyword BC_SCALARW (#,#)
# DEFAULT value of 0.0

Mockup of Task pane for specifying the Scalar properties for WALL boundary condition regions.
Note that this tab should only be available if scalar equations are being solved. Furthermore, the
number of scalars requiring input (here 2) comes from the number of scalar equations specified by
the user.
Subtask Pane Tab for INFLOW type (MI, PI, CG_MI) Boundary Condition Regions
Fluid (tab)
#    Define volume fraction
# Specification always available
# Sets keyword BC_EP_G(#)
#  DEFAULT value of 1.0 for MI and CG_MI; leave [UNDEFINED] for PI
#  Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0
#  Error Check: For PI - either all are defined and sum to 1, or all are undefined
#    Define inflow properties: Mass inflow specification changes based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
# For BC_TYPE='MI' and XZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Y-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:

Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='MI' and YZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    X-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='MI' and XY-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Z-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='CG_MI' or 'PI'
#  Specify all velocity components:
#    Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0

Define temperature
# Specification always available
# Input required for any of the following
#  Fluid density model: Ideal Gas Law
#  Fluid viscosity model: Sutherland's Law
#  Energy equations are solved
# Sets keyword BC_T_G(#)
# DEFAULT value of 293.15

Define pressure
# Specification always available
# Input required when combining ideal gas law and specified mass inflow rate
# Input required for BC_TYPE = PI
# Sets keyword BC_P_G(#)
# DEFAULT 101.325d3

Select species and set mass fractions (table format)
# Specification always available
# Input required for species equations
# Drop down menu of fluid species
# Sets keyword BC_X_G(#,#)
# DEFAULT - last defined species has mass fraction of 1.0
# Error check: mass fractions must sum to one

Turbulence: Define k-ε turbulent kinetic energy
# Specification only available with K-Epsilon turbulence model
# Sets keyword BC_K_TURB_G(#)
# DEFAULT value of 0.0

Turbulence: Define k-ε turbulent dissipation
# Specification only available with K-Epsilon turbulence model
# Sets keywords BC_E_TURB_G(#)
# DEFAULT value of 0.0

Mockup of Task pane for specifying the Fluid properties for MI, boundary condition regions.

Solid-# (tab) - Rename tab to user provided solids name.

Define volume fraction
# Specification always available
# Sets keyword BC_EP_S(#,#)
#  DEFAULT value of 1.0 - (sum of previous tabs) for MI and CG_MI; leave [UNDEFINED] for PI
#  Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0
#  Error Check: For PI - either all are defined and sum to 1, or all are undefined
Define inflow properties: Mass inflow specification changes based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
# For BC_TYPE='MI' and XZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Y-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_S(#,#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_G(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='MI' and YZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    X-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_U_S(#,#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_S#, (#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define Y-Axial Velocity
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_S(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='MI' and XY-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Z-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_W_S(#,#)
# DEFAULT value of 0.0
Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_S(#,#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_S(#,#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='CG_MI' or 'PI'
#  Specify all velocity components:
#    Define X-Axial Velocity
# Sets keyword BC_U_S(#,#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0

Define temperature
# Specification always available
# Input required when energy equations are solved
# Sets keyword BC_T_S(#,#)
# DEFAULT value of 293.15
Select species and set mass fractions (table format)
# Specification always available
# Input required for species equations
# Drop down menu of fluid species
# Sets keyword BC_X_S(#,#,#)
# DEFAULT - last defined species has mass fraction of 1.0
# Error check: mass fractions must sum to one

Mockup of Task pane for specifying the Solids properties for MI, CG_MI, and PI boundary
condition regions.

Scalar (tab) - Tab only available if scalar equations are solved
#    Define initial scalar value
# Sets keyword BC_SCALAR(#,#)
# DEFAULT value of 0.0
Subtask Pane Tab for PRESSURE OUTFLOW type (PO) Boundary Condition Regions
Fluid (tab)
#    Define pressure
# Specification always available
# Input required
# Sets keyword BC_P_G(#)
# DEFAULT 101.325d3

The remaining inputs are "optional." They do not have default values, because MFIX will calculate
appropriate values if they are unspecified and 'backflow' occurs at the outlet.
#    Define volume fraction
# Specification always available
# Sets keyword BC_EP_G(#)
# No DEFAULT value
# Error Check: If any volume fraction for the BC region is specified, then all volume fractions for the BC region must be specified and must sum to one.
#    Define temperature
# Specification always available
# NO DEFAULT value
# Sets keyword BC_T_G(#)
#    Select species and set mass fractions (table format)
# Specification always available
# NO DEFAULT value
# Sets keyword BC_X_G(#,#)
# Error check: if specified, mass fractions must sum to one_

Solids-# (tab)

All inputs are optional. They do not have default values, because MFIX will calculate appropriate values if they are unspecified and 'backflow' occurs at the outlet.
#    Define volume fraction
# Specification always available
# Sets keyword BC_EP_S(#,#)
# No DEFAULT value
# Error Check: If any volume fraction for the BC region is specified, then all volume fractions for the BC region must be specified and must sum to one.
#    Define temperature
# Specification always available
# NO DEFAULT value
# Sets keyword BC_T_S(#,#)
#    Select species and set mass fractions (table format)
# Specification always available
# NO DEFAULT value
# Sets keyword BC_X_S(#,#,#)
# Error check: if specified, mass fractions must sum to one_
Scalar (tab) - Tab only available if scalar equations are solved
All inputs are optional. They do not have default values, because MFIX will calculate appropriate  values if they are unspecified and 'backflow' occurs at the outlet.

Define scalar value
# Sets keyword BC_SCALAR(#,#)
# NO DEFAULT value

Subtask Pane Tab for MASS OUTFLOW type (MO) Boundary Condition Regions
Fluid (tab)
#    Define inflow properties: Mass inflow specification changes based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
# For BC_TYPE='MO' and XZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Y-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='MO' and YZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    X-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='MO' and XY-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Z-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_W_G(#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_G(#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_G(#)
# DEFAULT value of 0.0
Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
# For BC_TYPE='CG_MO'
#  Specify all velocity components:
#    Define X-Axial Velocity
# Sets keyword BC_U_G(#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_V_G(#)
# DEFAULT value of 0.0

Define duration to average outflow rate.
# Specification always available
# Input required
# Sets keyword BC_DT_0(#)
# DEFAULT value of 0.1
# Error Check: Value should be positive (non-negative)
# BC_DT_0 specification should persist across the gas and solids tabs. If the user sets
it in the gas phase tab, but then changes it under a solids tab, a warning messing
indicating that this value is 'constant' across all phases should be given.

The remaining inputs are only required when either the mass or the volumetric flowrates are
specified. They are not required if the velocities are given for the outlet.

Define volume fraction
# Specification only available with mass or volumetric flowrates.
# Input required
# Sets keyword BC_EP_G(#)
# DEFAULT value 1.0
# Error Check: If any volume fraction for the BC region is specified, then all volume
fractions for the BC region must be specified and must sum to one.

Define temperature
# Specification only available with mass or volumetric flowrates and R_G0 is UNDEFINED
# DEFAULT value 293.15
# Sets keyword BC_T_G(#)
Select species and set mass fractions (table format)
# Specification only available with mass or volumetric flowrates and R_G0 is UNDEFINED
# DEFAULT value 1.0 of last defined species
# Sets keyword BC_X_G(#,#)
# Error check: if specified, mass fractions must sum to one

Solids-# (tab)

Define inflow properties: Mass inflow specification changes based on the BC_TYPE and
Region orientation (e.g., XZ-Plane)
# For BC_TYPE='MO' and XZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Y-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_V_S(#,#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_ S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_ S(#,#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_ S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_ S(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='MO' and YZ-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    X-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_U_ S(#,#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_ S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_ S(#,#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define Y-Axial Velocity
# Sets keyword BC_V_ S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_W_ S(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='MO' and XY-Plane region
#  Select mass inflow specification type:
#    Available selections:
#    Z-Axial Velocity (m/s) [DEFAULT]
# Sets keyword BC_W_ S(#,#)
# DEFAULT value of 0.0
#    Volumetric Flowrate (m3/s)
# Sets keyword BC_VOLFLOW_ S(#,#)
# DEFAULT value of 0.0
#    Mass Flowrate (kg/s)
# Sets keyword BC_MASSFLOW_ S(#,#)
# DEFAULT value of 0.0
#  Define Tangential Velocities:
#    Define X-Axial Velocity
# Sets keyword BC_U_ S(#,#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_ S(#,#)
# DEFAULT value of 0.0
# For BC_TYPE='CG_MO'
#  Specify all velocity components:
#    Define X-Axial Velocity
# Sets keyword BC_U_ S(#,#)
# DEFAULT value of 0.0
#    Define Y-Axial Velocity
# Sets keyword BC_V_ S(#,#)
# DEFAULT value of 0.0
#    Define Z-Axial Velocity
# Sets keyword BC_V_ S(#,#)
# DEFAULT value of 0.0

Define duration to average outflow rate.
# Specification always available
# Input required
# Sets keyword BC_DT_0(#)
# DEFAULT value of 0.1
# Error Check: Value should be positive (non-negative)
BC_DT_0 specification should persist across the gas and solids tabs. If the user sets it in the gas
phase tab, but then changes it under a solids tab, a warning messing indicating that this value is
'constant' across all phases should be given.
"""
