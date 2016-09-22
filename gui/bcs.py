# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QComboBox, QGridLayout, QGroupBox, QHBoxLayout,
                            QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget)


from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4

UserRole = QtCore.Qt.UserRole

from constants import *
from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit, ComboBox

from project import Equation, FloatExp

from tools.general import (set_item_noedit, set_item_enabled,
                           get_combobox_item, get_selected_row,
                           widget_iter)

from tools.keyword_args import mkargs

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

PAGE_WALL, PAGE_INFLOW, PAGE_PO, PAGE_MO = range(4) #page indices in stackedwidget

SPECIFIED_MASS_FRACTION = 1


xmap = {'X':'u', 'Y':'v', 'Z': 'w'}

class BCS(object):
    #Boundary Conditions Task Pane Window: This section allows a user to define the boundary
    #conditions for the described model. This section relies on regions named in the Regions section.

    def init_bcs(self):
        ui = self.ui.boundary_conditions

        self.bcs = {} # key: index.  value: data dictionary for boundary cond
        self.bcs_current_indices = [] # List of BC indices
        self.bcs_current_regions = [] # And the names of the regions which define them
        self.bcs_region_dict = None
        self.bcs_fluid_species_boxes = OrderedDict() # In the 'wall' tab each species needs a groupbox

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
        ui.checkbox_bc_jj_ps_args_BC.key = 'bc_jj_ps'
        ui.checkbox_bc_jj_ps_args_BC.args = ['BC']
        self.add_tooltip(ui.checkbox_bc_jj_ps_args_BC, 'bc_jj_ps')

        # Combobox callbacks
        for name in ('fluid_energy_eq', #'fluid_species_eq',
                     'bc_jj_ps',
                     'solids_energy_eq', 'solids_granular_energy_eq'): # no solids_species_eq?
            item = getattr(ui, 'combobox_%s_type' % name)
            setter = getattr(self, 'set_bcs_%s_type' % name)
            item.currentIndexChanged.connect(setter)

        # Tangential velocities are not auto-registered
        for suffix in ('1', '2', '1_3', '2_3'):
            widget = getattr(ui, 'lineedit_fluid_tangential_velocity_'+suffix)
            widget.value_updated.connect(self.project.submit_change)
        for suffix in ('1', '2', '1_4', '2_4'):
            widget = getattr(ui, 'lineedit_solids_tangential_velocity_'+suffix)
            widget.value_updated.connect(self.project.submit_change)


        # Inflow/outflow pane has some special inputs
        for phase_type in 'fluid', 'solids': # type
            for flow_direction in 'inflow', 'outflow': #direction
                #le = ui.lineedit_fluid_inflow
                le = getattr(ui, 'lineedit_%s_%s' % (phase_type, flow_direction), None)
                le.value_updated.connect(self.bcs_handle_flow_input)
                le.key = 'Unset' # diagnostic
                le.args = ['BC'] if phase_type=='fluid' else ['BC', 'P']
                le.dtype = float
                # Tooltip added dynamically
                cb = getattr(ui, 'combobox_%s_%s_type' % (phase_type, flow_direction))
                cb.value_updated.connect(self.bcs_handle_flow_type)
                cb.key = '%s_%s' % (phase_type, flow_direction)
                cb.help_text = 'Select %s %s specification type' % (phase_type, flow_direction)
                cb.setToolTip(cb.help_text)

        le = ui.lineedit_solids_inflow
        cb = ui.combobox_solids_inflow_type
        le.value_updated.connect(self.bcs_handle_flow_input)
        le.key = 'Unset' # diagnostic
        le.args = ['BC', 'P']
        le.dtype = float
        # Tooltip added dynamically
        cb.value_updated.connect(self.bcs_handle_flow_type)
        cb.key = 'solids_inflow'

        # Not auto-registered with project manager
        key = 'bc_ep_s'
        for widget in (ui.lineedit_bc_ep_s_args_BC_P, ui.lineedit_bc_ep_s_2_args_BC_P):
            widget.value_updated.connect(self.handle_bcs_volume_fraction)
            widget.key = key
            widget.args = ['BC', 'P']
            widget.dtype = float
            self.add_tooltip(widget, key)

        key = 'bc_dt_0'
        for widget in (ui.lineedit_bc_dt_0, ui.lineedit_bc_dt_0_4):
            widget.value_updated.connect(self.handle_bc_dt_0)
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)
            # Error Check: Value should be positive (non-negative)
            widget.min = 0


    def bcs_set_volume_fraction_limit(self):
        # Set bc_ep_g from bc_ep_s, like we do for ICs
        if not self.bcs_current_indices:
            return
        if not self.bcs_current_solid:
            return
        BC0 = self.bcs_current_indices[0]
        P = self.bcs_current_solid
        ui = self.ui.boundary_conditions
        key = 'bc_ep_s'

        s = sum(safe_float(self.project.get_value(key, default=0, args=[BC0, p]))
                for p in range(1, len(self.solids)+1) if p != P)

        lim = max(0, 1.0 - s)
        lim = round(lim, 10) # avoid problem with 1 - 0.9 != 0.1

        for widget in (ui.lineedit_bc_ep_s_args_BC_P, ui.lineedit_bc_ep_s_2_args_BC_P):
            widget.min = 0.0
            widget.max = lim


    def handle_bcs_volume_fraction(self, widget, val, args):
        if not self.bcs_current_indices:
            return
        BC0 = self.bcs_current_indices[0]
        if not self.bcs_current_solid:
            return
        P = self.bcs_current_solid
        ui = self.ui.boundary_conditions
        key = 'bc_ep_s'
        self.project.submit_change(widget, val, args)

        s = sum(safe_float(self.project.get_value(key, default=0, args=[BC0, s]))
                for s in range(1, len(self.solids)+1))
        if s > 1.0:
            self.warning("Volume fractions sum to %s, must be <= 1.0" % s,
                         popup=True)
            return # ?
        val = round(1.0 - s, 10)
        for BC in self.bcs_current_indices:
            self.update_keyword('bc_ep_g', val, args=[BC])


    def bcs_show_regions_popup(self):
        # Users cannot select inapplicable regions.
        # BC regions must be planes, volumes, or STLs (not volumes or points)
        # No region can define more than one boundary condition.
        ui = self.ui.boundary_conditions
        rp = self.regions_popup
        rp.clear()
        for (name,data) in self.bcs_region_dict.items():
            shape = data.get('type', '---')
            available = (data.get('available', True)
                         and not self.check_region_in_use(name)
                         and (shape in ('STL', 'box')
                              or 'plane' in shape))
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.bcs_add_regions)
        rp.cancel.connect(self.bcs_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.detail_pane,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup('boundary conditions')


    def bcs_cancel_add(self):
        ui = self.ui.boundary_conditions

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.detail_pane,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def bcs_add_regions(self):
        #Select boundary type
        # Selection is required
        # Available selections:
        #  Mass Inflow
        #    Plane regions set keyword BC_TYPE(#) to 'MI'
        #    STL regions set keyword BC_TYPE(#) to 'CG_MI'
        #    Not available for volume regions
        #  Pressure Outflow
        #    Plane regions set keyword BC_TYPE(#) to 'PO'
        #    STL regions set keyword BC_TYPE(#) to 'CG_PO'
        #    Not available for volume regions
        #  No Slip Wall
        #    Volume and plane regions set keyword BC_TYPE(#) to 'NSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_NSW'
        #  Free Slip Wall
        #    Volume and plane regions set keyword BC_TYPE(#) to 'FSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_FSW'
        #  Partial Slip Wall
        #    Volume and plane regions set keyword BC_TYPE(#) to 'PSW'
        #    STL regions set keyword BC_TYPE(#) to 'CG_PSW'
        #  Pressure Inflow
        #    Plane regions set keyword BC_TYPE(#) to 'PI'
        #    Not available for volume regions
        #    Not available for STL regions
        # Mass Outflow
        #    Plane regions set keyword BC_TYPE(#) to 'MO'
        #    STL regions set keyword BC_TYPE(#) to 'CG_MO'
        #    Not available for volume regions

        # Specification always available
        # DEFAULT - No slip wall
        # Error check: mass fractions must sum to one
        # (selection logic implemented in regions_popup.py)

        # Interactively add regions to define BCs
        ui = self.ui.boundary_conditions
        rp = self.regions_popup
        self.bcs_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        bc_type = rp.combobox.currentIndex()
        if not selections:
            return
        self.bcs_add_regions_1(selections, bc_type=bc_type, indices=None) # Indices will be assigned
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
        self.update_nav_tree()

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
        ui.detail_pane.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.detail_pane):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.bcs_setup_current_tab() # reinitialize all widgets in current tab


    def fixup_bcs_table(self, tw, stretch_column=0):
        ui = self.ui.boundary_conditions
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

        if tw == ui.tablewidget_regions:
            ui.top_frame.setMaximumHeight(height+40)
            ui.top_frame.setMinimumHeight(header_height+40)
            ui.top_frame.updateGeometry()
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height)
        else:
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(height)
        tw.updateGeometry() #? needed?


    def bcs_update_enabled(self):
        if self.bcs:
            # Never disable if there are BCs defined
            disabled = False
        else:
            # If there are no solids, no scalar equations, and the fluid solver is disabled,
            # then we have no input tabs on the BCs pane, so disable it completely
            regions = self.ui.regions.get_region_dict()
            nregions = sum(1 for (name, r) in regions.items()
                           if not self.check_region_in_use(name)
                           and (r.get('type')=='STL' or 'plane' in r.get('type')))
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
        # Update the keys which define the region the BC applies to
        if bc_type is not None:
            val = BC_TYPES[bc_type]
            if data.get('type') == 'STL':
                val = 'CG_' + val
            if val== 'CG_PI': # Shouldn't happen!
                self.error("Invalid bc_type %s" % val)
                return
            self.update_keyword('bc_type', val, args=[idx])

        no_k = self.project.get_value('no_k')

        for (key, val) in zip(('x_w', 'y_s', 'z_b',
                               'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            # bc_z_t and bc_z_b keywords should not be added when no_k=True
            if no_k and key in ('z_t', 'z_b'):
                continue
            self.update_keyword('bc_'+key, val, args=[idx])


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
            self.bcs_add_regions_1(regions, bc_type=None, indices=indices)


    def setup_bcs(self):
        ui = self.ui.boundary_conditions

        # Grab a fresh copy, may have been updated
        self.bcs_region_dict = self.ui.regions.get_region_dict()

        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.bcs.items():
            region = data['region']
            if region in self.bcs_region_dict:
                self.bcs_region_dict[region]['available'] = False

        self.fixup_bcs_table(ui.tablewidget_regions)
        row = get_selected_row(ui.tablewidget_regions)
        # Autoselect if only 1 row
        if row is None and ui.tablewidget_regions.rowCount() == 1:
            row = 0
            ui.tablewidget_regions.setCurrentCell(row, 0)
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.detail_pane.setEnabled(enabled)

        # Tabs group boundary condition parameters for phases and additional equations. Tabs are
        # unavailable if no input is required from the user.
        #
        #Fluid tab - Unavailable if the fluid phase was disabled.
        b = ui.pushbutton_fluid
        b.setText(self.fluid_phase_name)
        b.setEnabled(not self.fluid_solver_disabled)
        if self.fluid_solver_disabled:
            if self.bcs_current_tab == 0: # Don't stay on disabled tab
                self.bcs_change_tab(*(SOLIDS_TAB, 1) if self.solids else (SCALAR_TAB,None)) # what if nscalar==0?
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
            extent = [d.get('bc_'+k, None) for k in ('x_w', 'y_s', 'z_b',
                                                     'x_e', 'y_n', 'z_t')]
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

                ext2 = [0 if x is None else x for x in
                        data.get('from',[]) + data.get('to',[])]

                # TODO this only works for a single STL
                if (is_stl and data.get('type')=='STL'
                    or not is_stl and ext2==extent):

                    if data.get('available', True):
                        if bc_type is None:
                            self.warn("no bc_type for region %s" % bc.ind)
                        if bc_type not in BC_TYPES:
                            self.warn("invalid bc_type %s for region %s" % (bc_type, bc.ind))
                        else:
                            self.bcs_add_regions_1([region_name], bc_type=BC_TYPES.index(bc_type), indices=[bc.ind])
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


    def set_bcs_fluid_species_eq_type(self, eq_type, species_index):
        if not self.bcs_current_indices:
            return
        #Select species equation boundary type:
        # Selection only available when solving species equations
        # Available selections:
        #  No-Flux [DEFAULT]
        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_X_G(#,#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_X_G(#,#) to 0.0
            c = 0
            #    Sets keyword BC_XW_G(#,#) to UNDEFINED
            xw = None
        #  Specified Mass Fraction
        elif eq_type == SPECIFIED_MASS_FRACTION:
            #    Sets keyword BC_XW_T_G(#,#) to UNDEFINED # should be BC_H_X_G ?
            hw = None
            #    Sets keyword BC_C_X_G(#,#) to 0.0
            c = 0.0
            #    Requires BC_XW_G(#,#)
            xw = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_X_G(#,#) to 0.0
            hw = 0.0
            #    Requires BC_C_X_G(#,#)
            c = True
            #    Sets keyword BC_XW_G(#,#) to UNDEFINED
            xw = None
        #  Convective Flux
        elif eq_type == CONVECTIVE_FLUX:
            #,#    Requires BC_HW_X_G(#,#)
            hw = True
            #    Sets keyword BC_C_X_G(#,#) to 0.0
            c = 0.0
            #    Requires BC_XW_G(#,#)
            xw = True
        else:
            self.error("Invalid fluid species_eq type %s" % eq_type)
            return

        for BC in self.bcs_current_indices:
            #FIXME will not survive species deletion
            self.bcs[BC]['fluid_species_%s_eq_type'%species_index] = eq_type
            for (key, val) in (('bc_hw_x_g', hw), ('bc_c_x_g', c), ('bc_xw_g', xw)):
                if val is True:
                    pass # 'required'
                else:
                    self.update_keyword(key, val, args=[BC, species_index])

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
            self.error("Invalid solid energy_eq type %s" % eq_type)
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
            self.error("Invalid granualar energy_eq type %s" % eq_type)
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
        if not self.bcs_current_indices:
            return # Disable inputs?
        BC0 = self.bcs_current_indices[0]
        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
        if bc_type.endswith('W'):
            self.setup_bcs_fluid_wall_tab()
        elif bc_type.endswith('I'):
            self.setup_bcs_fluid_inflow_tab()
        elif bc_type.endswith('PO'):
            self.setup_bcs_fluid_po_tab()
        elif bc_type.endswith('MO'):
            self.setup_bcs_fluid_mo_tab()
        else:
            self.error("Invalid bc_type %s" % bc_type)


    def setup_bcs_fluid_wall_tab(self):
        # Subtask Pane Tab for Wall type (NSW, FSW, PSW, CG_NSW, CG_FSW, CG_PSW) Boundary Condition Regions
        ui = self.ui.boundary_conditions
        ui.page_fluid.setCurrentIndex(PAGE_WALL)

        if not self.bcs_current_indices:
            return # Nothing selected.  (Clear out all lineedits?)
        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        def get_widget(key, ui):
            for pat in ('lineedit_keyword_%s_args_BC',
                        'lineedit_keyword_%s_args_BC_species',
                        'lineedit_%s_args_BC',
                        'lineedit_%s_args_BC_species'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True,
                             species_index=None, box=None, suffix=''):
            ui = box or self.ui.boundary_conditions # widgets nested into groupboxes
            for pat in ('label_%s', 'label_%s_units',
                        'lineedit_keyword_%s_args_BC',
                        'lineedit_keyword_%s_args_BC_species'):
                name = pat%(key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0, species=species_index)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    args = mkargs(key, bc=BC, species=species_index)
                    self.update_keyword(key, val, args=args)
            get_widget(key+suffix, ui).updateValue(key, val, args=args)


        def make_fluid_species_box(title, species_index):
            box = QGroupBox(title, parent=ui.page_fluid_wall)
            box_layout = QGridLayout()
            #margins = ui.groupbox_fluid_energy_eq.getContentsMargins()
            box.setFlat(True)
            margins = (5, 10, 0, 5) #left top right bottom
            box_layout.setContentsMargins(*margins)
            box.setLayout(box_layout)
            cb = box.combobox_fluid_species_eq_type = QComboBox()
            for item in ('No-Flux', 'Specified Mass Fraction',
                         'Specified Flux', 'Convective Flux'):
                cb.addItem(item)
            cb.setCurrentIndex(NO_FLUX)
            cb.currentIndexChanged.connect(
                lambda index, species_index=species_index: self.set_bcs_fluid_species_eq_type(index, species_index))
            hbox = QHBoxLayout()
            label = QLabel('Type')
            # TODO: tooltips
            hbox.addWidget(label)
            hbox.addWidget(cb)
            hbox.addStretch()
            box_layout.addLayout(hbox, 0, 0, 1, 3)
            row = 0
            for (label_text, key, units) in (('Wall mass fraction', 'bc_xw_g', None),
                                             ('Constant flux', 'bc_c_x_g', 'TBD'),
                                             ('Transfer coefficient', 'bc_hw_x_g', 'TBD'),
                                             ('Free stream mass frac.', 'bc_xw_g', None)):
                row += 1
                suffix = '_2' if label_text.startswith('Free') else ''
                label = QLabel(label_text)
                self.add_tooltip(label, key)
                box_layout.addWidget(label, row, 0, 1, 1)
                setattr(box, 'label_%s'%(key+suffix), label)
                le = LineEdit()
                le.key = key
                le.args = ['BC', species_index]
                le.dtype = float
                # TODO: LIMITS
                self.add_tooltip(le, key)
                box_layout.addWidget(le, row, 1)
                setattr(box, 'lineedit_keyword_%s_args_BC_species'%(key+suffix), le)
                self.project.register_widget(le, [key], ['BC', species_index])
                if units:
                    label = QLabel(units)
                    box_layout.addWidget(label, row, 2)
                    setattr(box, 'label_%s_units'%(key+suffix), label)
            page_layout = ui.page_fluid_wall.layout()
            page_layout.insertWidget(page_layout.count()-1, box)
            return box

        #Fluid (tab)

        #    Define transfer coefficient
        # Specification only available with PSW
        # Sets keyword BC_HW_G(#)
        # DEFAULT value of 0.0

        #    Define Wall U-velocity
        # Specification only available with PSW
        # Sets keyword BC_UW_G(#)
        # DEFAULT value of 0.0

        #    Define Wall V-velocity
        # Specification only available with PSW
        # Sets keyword BC_VW_G(#)
        # DEFAULT value of 0.0

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
                    #self.error("Cannot determine type for fluid energy boundary equation %s" % BC0)
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
        species_eq = self.project.get_value('species_eq', default=True, args=[0])#0 for fluid-phase
        enabled = bool(species_eq)

        #  Note - we need a copy of this groupbox for each fluid species
        # Let's only do this when the species have changed
        box_keys = list(self.bcs_fluid_species_boxes.keys())
        species_names = list(self.fluid_species.keys())
        species_keys = [self.project.get_value('species_alias_g', default=name, args=[i])
                           for (i, name) in enumerate(species_names, 1)]
        if box_keys != species_keys:
            n_boxes = len(box_keys)
            n_species = len(species_keys)
            if n_boxes > n_species:
                for i in range(n_boxes-1, n_species-1, -1):
                    box_key = box_keys[i]
                    box = self.bcs_fluid_species_boxes[box_key]
                    for w in widget_iter(box):
                        self.unregister_widget(w)
                        w.deleteLater()
                    box.hide()
                    box = None
                    self.bcs_fluid_species_boxes[box_key].deleteLater()

            elif n_boxes < n_species:
                for i in range(n_boxes, n_species):
                    species_key = species_keys[i]
                    box = make_fluid_species_box(species_key, i+1)
                    self.bcs_fluid_species_boxes[species_key] = box

            tmp_dict = OrderedDict()
            for (key, box) in zip(species_keys, self.bcs_fluid_species_boxes.values()):
                tmp_dict[key] = box
                box.setTitle("Species Equation: %s" % key)
            self.bcs_fluid_species_boxes = tmp_dict

        for (species_index, box) in enumerate(self.bcs_fluid_species_boxes.values(), 1):
            enabled = bool(species_eq)
            box.setEnabled(enabled)
            eq_type = None
            if enabled:
                eq_type = self.bcs[BC0].get('fluid_species_%s_eq_type'%species_index)
                if eq_type is None: # Attempt to infer from keywords
                    hw = self.project.get_value('bc_hw_x_g', args=[BC0, species_index])
                    c = self.project.get_value('bc_c_x_g', args=[BC0, species_index])
                    xw = self.project.get_value('bc_xw_g', args=[BC0, species_index])
                    if hw==0.0 and c==0.0 and xw is None:
                        eq_type = NO_FLUX
                    elif hw is None and c==0.0 and xw is not None:
                        eq_type = SPECIFIED_MASS_FRACTION
                    elif hw==0.0 and c!=0.0 and xw is None:
                        eq_type = SPECIFIED_FLUX
                    elif hw is not None and c==0.0 and xw is not None:
                        eq_type = CONVECTIVE_FLUX
                    else:
                        #self.error("Cannot determine type for fluid species boundary equation %s" % BC0)
                        eq_type = NO_FLUX # Default
                    self.set_bcs_fluid_species_eq_type(eq_type, species_index)

            if eq_type is not None:
                box.combobox_fluid_species_eq_type.setCurrentIndex(eq_type)

            if species_eq:
                #Define wall mass fraction
                # Specification only available with 'Specified Mass Fraction' BC type
                # Sets keyword BC_XW_G(#)
                # DEFAULT value of 0.0
                enabled = (eq_type==SPECIFIED_MASS_FRACTION)
                key = 'bc_xw_g'
                default = 0.0 if enabled else None
                setup_key_widget(key, default, enabled, species_index=species_index, box=box)
                # Hack to prevent dup. display
                if enabled:
                    box.lineedit_keyword_bc_xw_g_2_args_BC_species.setText('')
                else:
                    box.lineedit_keyword_bc_xw_g_args_BC_species.setText('')

                #Define constant flux
                # Specification only available with 'Specified Flux' BC type
                # Sets keyword BC_C_X_G(#)
                # DEFAULT value of 0.0
                enabled = (eq_type==SPECIFIED_FLUX)
                key = 'bc_c_x_g'
                default = 0.0
                setup_key_widget(key, default, enabled, species_index=species_index, box=box)

                #Define transfer coefficient
                # Specification only available with 'Convective Flux' BC type
                # Sets keyword BC_HW_X_G(#)
                # DEFAULT value of 0.0
                enabled = (eq_type==CONVECTIVE_FLUX)
                key = 'bc_hw_x_g'
                default = 0.0
                setup_key_widget(key, default, enabled, species_index=species_index, box=box)

                #Define free stream mass fraction
                # Specification only available with 'Convective Flux' BC type
                # Sets keyword BC_XW_G(#)
                # DEFAULT value of 0.0
                enabled = (eq_type==CONVECTIVE_FLUX)
                key = 'bc_xw_g'
                default = 0.0 if enabled else None
                setup_key_widget(key, default, enabled, species_index=species_index, box=box, suffix='_2')
                # Hack to prevent dup. display
                if enabled:
                    box.lineedit_keyword_bc_xw_g_args_BC_species.setText('')
                else:
                    box.lineedit_keyword_bc_xw_g_2_args_BC_species.setText('')


    def setup_bcs_solids_tab(self, P):
        #Solids-# (tab) - (Replace with phase name defined by the user)
        # Note, solids phases are numbered 1-N
        self.bcs_current_solid = self.P = P
        if P is None: # Nothing to do
            return

        if not self.bcs_current_indices: # No region selected
            # TODO clear all widgets (?)
            return
        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return
        if bc_type.endswith('W'):
            self.setup_bcs_solids_wall_tab(P)
        elif bc_type.endswith('I'):
            self.setup_bcs_solids_inflow_tab(P)
        elif bc_type.endswith('PO'):
            self.setup_bcs_solids_po_tab(P)
        elif bc_type.endswith('MO'):
            self.setup_bcs_solids_mo_tab(P)
        else:
            self.error("Invalid bc_type %s" % bc_type)


    def setup_bcs_solids_wall_tab(self, P):
        #Subtask Pane Tab for Wall type (NSW, FSW, PSW, CG_NSW, CG_FSW, CG_PSW) Boundary
        #Solids-# (tab) - (Replace with phase name defined by the user)
        # Note, solids phases are numbered 1-N
        ui = self.ui.boundary_conditions
        ui.page_solids.setCurrentIndex(PAGE_WALL)

        self.bcs_current_solid = self.P = P
        if P is None: # Nothing to do
            return

        if not self.bcs_current_indices: # No region selected
            # TODO clear all widgets (?)
            return

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
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC_P',
                         'lineedit_keyword_%s_args_BC',
                         'lineedit_%s_args_BC_P',
                         'lineedit_%s_args_BC',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0, phase=P)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC, phase=P))
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
        enabled = bc_jj_ps==1 and bc_jj_m==True
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
                    #self.error("Cannot determine type for solid %s energy boundary equation %s" % (P,BC0))
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
                    #self.error("Cannot determine type for solid %s granular energy boundary equation %s" % (P,BC0))
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
        #Scalar (tab) - Tab only available if scalar equations are solved
        nscalar = self.project.get_value('nscalar', default=0)
        if nscalar==0: # We shouldn't be here
            # Clear out pane?
            return
        if not self.bcs_current_indices: # nothing selected
            return
        BC0 = self.bcs_current_indices[0]
        bc_type = self.project.get_value('bc_type', args=BC0)
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return
        if bc_type.endswith('W'):
            self.setup_bcs_scalar_wall_tab()
        else:
            pass # TODO implement


    def clear_bcs_scalar_tab(self):
        ui = self.ui.boundary_conditions
        nscalar = self.project.get_value('nscalar', default=0)
        if nscalar == 0:
            # Clear contents?
            return

        page = ui.page_scalar
        layout = page.layout()

        spacer = None
        for i in range(layout.count()-1, -1, -1):
            item = layout.itemAt(i)
            if not item:
                continue
            widget = item.widget()
            if not widget:
                spacer = item
                continue
            if isinstance(widget, LineEdit):
                self.project.unregister_widget(widget)
            widget.setParent(None)
            for w2 in widget_iter(widget):
                if isinstance(w2, (LineEdit, ComboBox)):
                    self.project.unregister_widget(w2)
            widget.deleteLater()

        if spacer:
            layout.removeItem(spacer)
        return spacer

    def set_bcs_scalar_eq_type(self, eq_type, i):
        if not self.bcs_current_indices:
            return

        ui = self.ui.boundary_conditions
        cb = getattr(ui, 'combobox_scalar_%s_eq_type' % i, None)
        if cb is None:
            self.error("Invalid scalar_eq %s" % i)
            return
        # Available selections
        #  No-Flux [DEFAULT]
        if eq_type == NO_FLUX:
            #    Sets keyword BC_HW_SCALAR(#,#) to 0.0
            hw = 0.0
            #    Sets keyword BC_C_SCALAR(#,#) to 0.0
            c = 0.0
            #    Sets keyword BC_SCALARW(#,#) to UNDEFINED
            sw = None
        #  Specified Temperature
        elif eq_type == SPECIFIED_TEMPERATURE:
            #    Sets keyword BC_HW_T_S(#,#) to UNDEFINED
            #                 ^^^^ presumably BC_HW_SCALAR
            hw = None
            #    Sets keyword BC_C_SCALAR (#,#) to 0.0
            c = 0.0
            #    Requires BC_SCALARW (#,#)
            sw = True
        #  Specified Flux
        elif eq_type == SPECIFIED_FLUX:
            #    Sets keyword BC_HW_T_S(#,#) to 0.0
            #                 ^^^^ presumably BC_HW_SCALAR
            hw = 0.0
            #    Requires BC_C_SCALAR (#)
            c = True
            #    Sets keyword BC_SCALARW (#,#) to UNDEFINED
            sw = None
        #  Convective Flux
        elif eq_type == CONVECTIVE_FLUX:
            #    Requires BC_HW_T_S(#,#)
            #                 ^^^^ presumably BC_HW_SCALAR
            hw = True
            #    Sets keyword BC_C_SCALAR (#,#) to 0.0
            c = 0.0
            #    Requires BC_SCALARW (#,#)
            sw = True
        else:
            self.error("Invalid scalar energy_eq type %s" % eq_type)
            return

        for BC in self.bcs_current_indices:
            # FIXME this won't survive index remapping when solids are deleted!
            self.bcs[BC]['scalar_%s_eq_type'%i] = eq_type
            for (key, val) in (('bc_hw_scalar', hw), ('bc_c_scalar', c), ('bc_scalarw', sw)):
                if val is True:
                    pass # 'required'
                else:
                    self.update_keyword(key, val, args=[BC,i])

        self.setup_bcs_scalar_tab()


    def setup_bcs_scalar_wall_tab(self):
        # Subtask Pane Tab for Wall type (NSW, FSW, PSW, CG_NSW, CG_FSW, CG_PSW) Boundary
        #Scalar (tab) - Tab only available if scalar equations are solved
        #Note that this tab should only be available if scalar
        #equations are being solved.  Furthermore, the number of
        #scalars requiring input (here 2) comes from the number of
        #scalar equations specified by the user.
        if not self.bcs_current_indices:
            return # No selection
        BC0 = self.bcs_current_indices[0]

        ui = self.ui.boundary_conditions
        nscalar = self.project.get_value('nscalar', default=0)
        old_nscalar = getattr(ui, 'nscalar', None)
        ui.nscalar = nscalar

        if nscalar == 0:
            # Clear contents?
            return

        page = ui.page_scalar
        page_layout = page.layout()

        spacer = None
        if nscalar != old_nscalar:
            spacer = self.clear_bcs_scalar_tab()

            for i in range(1, nscalar+1):
                groupbox = QGroupBox('Scalar %s' % i)
                groupbox.setFlat(True)
                # TODO adjust margins/spacing
                page_layout.addWidget(groupbox)
                groupbox_layout = QGridLayout()
                groupbox.setLayout(groupbox_layout)
                #    Select scalar boundary type:
                hbox = QHBoxLayout()
                label = QLabel("Type")
                hbox.addWidget(label)
                cb = ComboBox()
                setattr(ui, "combobox_scalar_%s_eq_type"%i, cb)
                # Available selections:
                cb.addItem("No-Flux")
                cb.addItem("Specified Temperature")
                cb.addItem("Specified Flux")
                cb.addItem("Convective Flux")
                cb.setCurrentIndex(NO_FLUX)
                cb.currentIndexChanged.connect(lambda eq_type, i=i: self.set_bcs_scalar_eq_type(eq_type, i))
                hbox.addWidget(cb)
                hbox.addStretch()

                groupbox_layout.addLayout(hbox, 0, 0, 1, -1)

                row = 0
                for (key, label_text, suffix) in (('bc_scalarw', 'Wall value', ''),
                                          ('bc_c_scalar', 'Constant flux', ''),
                                          ('bc_hw_scalar', 'Transfer coefficient', ''),
                                          ('bc_scalarw', 'Free stream value', '_2')):
                    row += 1 # Skip over "type"
                    args = ['BC', i]
                    label = QLabel(label_text)
                    self.add_tooltip(label, key)
                    setattr(ui, 'label_%s_args_BC_%s' % (key+suffix, i), label)
                    groupbox_layout.addWidget(label, row, 0)
                    le = LineEdit()
                    le.key = key
                    le.args = ['BC', i]
                    le.dtype = float
                    self.add_tooltip(le, key)
                    le.default_value = 0.0 #?
                    groupbox_layout.addWidget(le, row, 1)
                    self.project.register_widget(le, [key], ['BC', i])
                    setattr(ui, 'lineedit_keyword_%s_args_BC_%s' % (key+suffix, i), le)

                    # Right column is the stretchy one
                    for col in (0, 1):
                        groupbox_layout.setColumnStretch(col, col)

        # Clear memoized data above current nscalar if nscalar decreased
        if old_nscalar is None:
            old_nscalar = 0
        for i in range(nscalar+1, old_nscalar+2):
            for bc in self.bcs.values():
                for k in list(bc.keys()):
                    if k.startswith('scalar_%s' % i):
                        del bc[k]


        def get_widget(key, i):
            for pat in ('lineedit_keyword_%s_args_BC_%s',
                        'lineedit_%s_args_BC_%s'):
                name = pat % (key, i)
                widget = getattr(ui, name, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, i, default=None, enabled=True, suffix=''):
            for pat in ('label_%s_args_BC_%s',
                         'lineedit_keyword_%s_args_BC_%s'):
                name = pat%(key+suffix, i)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0, scalar=i)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC, scalar=i))
            get_widget(key+suffix, i).updateValue(key, val, args=args)

        for i in range(1, nscalar+1):
            hw = self.project.get_value('bc_hw_scalar', args=[BC0,i])
            c = self.project.get_value('bc_c_scalar', args=[BC0,i])
            sw = self.project.get_value('bc_scalarw', args=[BC0,i])
            eq_type = self.bcs[BC0].get('scalar_%s_eq_type'%i)
            cb = getattr(ui, 'combobox_scalar_%s_eq_type'%i)
            if eq_type is None:
                #  No-Flux [DEFAULT]
                #    Sets keyword BC_HW_SCALAR(#,#) to 0.0
                #    Sets keyword BC_C_SCALAR(#,#) to 0.0
                #    Sets keyword BC_SCALARW(#,#) to UNDEFINED
                if (hw==0) and (c==0) and (sw is None):
                    eq_type = NO_FLUX
                #  Specified Temperature
                #    Sets keyword BC_HW_T_S(#,#) to UNDEFINED
                #                 ^^^^ presumably BC_HW_SCALAR
                #    Sets keyword BC_C_SCALAR (#,#) to 0.0
                #    Requires BC_SCALARW (#,#)
                elif (hw is None) and (c==0.0) and (sw is not None):
                    eq_type = SPECIFIED_TEMPERATURE
                #  Specified Flux
                #    Sets keyword BC_HW_T_S(#,#) to 0.0
                #                 ^^^^ presumably BC_HW_SCALAR
                #    Requires BC_C_SCALAR (#)
                #    Sets keyword BC_SCALARW (#,#) to UNDEFINED
                elif (hw==0.0) and (c is not None) and (sw is None):
                    eq_type = SPECIFIED_FLUX
                #  Convective Flux
                #    Requires BC_HW_T_S(#,#)
                #                 ^^^^ presumably BC_HW_SCALAR
                #    Sets keyword BC_C_SCALAR (#,#) to 0.0
                #    Requires BC_SCALARW (#,#)
                elif (hw is not None) and (c==0.0) and (sw is not None):
                    eq_type = CONVECTIVE_FLUX

                if eq_type is None:
                    eq_type = NO_FLUX # default
                    self.set_bcs_scalar_eq_type(eq_type, i)
                    # Message?

            cb.setCurrentIndex(eq_type)

            #    Define wall temperature
            # Specification only available with 'Specified Temperature' BC type
            # Sets keyword BC_SCALARW (#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_TEMPERATURE)
            key = 'bc_scalarw'
            default = 0.0 if enabled else None
            setup_key_widget(key, i,  default, enabled)
            # Hack to prevent dup. display
            if enabled:
                le = getattr(ui, 'lineedit_keyword_bc_scalarw_2_args_BC_%s' %i)
                le.setText('')
            else:
                le = getattr(ui, 'lineedit_keyword_bc_scalarw_args_BC_%s' %i)
                le.setText('')

            #    Define constant flux
            # Specification only available with 'Specified Flux' BC type
            # Sets keyword BC_C_SCALAR (#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==SPECIFIED_FLUX)
            key = 'bc_c_scalar'
            default = 0.0
            setup_key_widget(key, i, default, enabled)

            #    Define transfer coefficient
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_HW_SCALAR(#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_hw_scalar'
            default = 0.0
            setup_key_widget(key, i, default, enabled)

            #    Define free stream temperature
            # Specification only available with 'Convective Flux' BC type
            # Sets keyword BC_SCALARW (#,#)
            # DEFAULT value of 0.0
            enabled = (eq_type==CONVECTIVE_FLUX)
            key = 'bc_scalarw'
            default = 0.0 if enabled else None
            setup_key_widget(key, i, default, enabled, suffix='_2')
            # Hack to prevent dup. display
            if enabled:
                le = getattr(ui, 'lineedit_keyword_bc_scalarw_args_BC_%s' %i)
                le.setText('')
            else:
                le = getattr(ui, 'lineedit_keyword_bc_scalarw_2_args_BC_%s' %i)
                le.setText('')
        if spacer:
            page_layout.addItem(spacer)


    def bcs_handle_flow_input(self, widget, data, ignore_key):
        if not data:
            log.error('bcs_handle_flow_input: no data')
            return
        key, val = data.popitem()
        P = self.bcs_current_solid
        for BC in self.bcs_current_indices:
            self.update_keyword(widget.key, val, args=mkargs(key, bc=BC, phase=P))


    def handle_bc_dt_0(self, widget, data, ignore_key):
        ui = self.ui.boundary_conditions
        if not self.bcs_current_indices:
            return
        BC0 = self.bcs_current_indices[0]
        # BC_DT_0 specification should persist across the gas and solids tabs.
        # If the user sets it in the gas phase tab, but then changes it under a solids tab,
        # a warning message indicating that this value is 'constant' across all phases should be given.
        key = 'bc_dt_0'
        prev_val = self.project.get_value(key, args=[BC0])
        new_key, new_val = data.popitem()
        if new_val == prev_val:
            return

        resp = self.message(text="Setting bc_dt_0 applies to all fluid and solid phases\nAre you sure?",
                            buttons=['yes', 'no'],
                            default='no')

        if resp != 'yes':
            widget.updateValue(key, prev_val)
            return


        for BC in self.bcs_current_indices:
            self.update_keyword(key, new_val, args=[BC])


    def bcs_handle_flow_type(self, widget, data, ignore_key):
        # This handles inflow and outflow for fluid and solids phases
        #  The combobox has a (non-keyword) 'key' set to distinguish where
        #  this is being called from
        if not data:
            log.error('bcs_handle_flow_type: no data')
            return
        ui = self.ui.boundary_conditions
        if not self.bcs_current_indices:
            return
        index = widget.currentIndex()

        BC0 = self.bcs_current_indices[0]
        P = self.bcs_current_solid
        key, val = data.popitem()
        phase_type, flow_direction = widget.key.split('_')
        axis = get_combobox_item(widget, 0).text()[0]

        if phase_type == 'fluid':
            subkeys = ['bc_%s_g' % xmap[axis], 'bc_volflow_g', 'bc_massflow_g']
        elif phase_type == 'solids':
            subkeys = ['bc_%s_s' % xmap[axis], 'bc_volflow_s', 'bc_massflow_s']
        else:
            raise ValueError(phase_type)

        subkey = subkeys[index]
        le = getattr(ui, 'lineedit_%s' % widget.key)
        prev_val = le.value

        # Memoize selection in 'bcs' object.  FIXME, this will not survive
        #  deleting/remapping solids
        bc_key = ('solid_%s_%s_type' % (P, flow_direction) if phase_type=='solids'
                  else 'fluid_%s_type'%flow_direction)

        for BC in self.bcs_current_indices:
            self.bcs[BC][bc_key] = index

        #    Available selections:
        # FLUID                                # SOLIDS
        #    Y-Axial Velocity (m/s) [DEFAULT]
        # Sets keyword BC_V_G(#)               # Sets keyword BC_V_S(#,#)
        # DEFAULT value of 0.0
        #    Volumetric Flowrate (m3/s)
        # Sets keyword BC_VOLFLOW_G(#)         # Sets keyword BC_VOLFLOW_S(#,#)
        # DEFAULT value of 0.0
        #    Mass Flowrate (kg/s)
        # Sets keyword BC_MASSFLOW_G(#)        # Sets keyword BC_MASSFLOW_S(#,#)
        # DEFAULT value of 0.0

        le.key = subkey
        le.args = ['BC', 'P'] if phase_type=='solids' else ['BC']
        le.dtype = float
        self.add_tooltip(le, subkey)

        for (i, k) in enumerate(subkeys):
            if i == index:
                continue
            for BC in self.bcs_current_indices:
                self.unset_keyword(k, args=mkargs(k, bc=BC, phase=P))

        val = self.project.get_value(subkeys[index], args=mkargs(subkeys[index], bc=BC0, phase=P))
        default = val if val is not None else prev_val if prev_val is not None else 0.0
        for BC in self.bcs_current_indices:
            self.update_keyword(subkey, default, args=mkargs(subkey, bc=BC, phase=P))
        le.updateValue(subkey, default)

        if widget.key == 'fluid_inflow':
            self.setup_bcs_fluid_inflow_tab()
        elif widget.key == 'fluid_outflow':
            self.setup_bcs_fluid_mo_tab()
        if widget.key == 'solids_inflow':
            self.setup_bcs_solids_inflow_tab(P)
        elif widget.key == 'solids_outflow':
            self.setup_bcs_solids_mo_tab(P)


    def update_bcs_fluid_mass_fraction_table(self):
        ui = self.ui.boundary_conditions
        table = ui.tablewidget_fluid_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        if not (self.fluid_species and self.bcs_current_indices):
            self.fixup_bcs_table(table)
            ui.groupbox_fluid_composition.setEnabled(False)
            return
        ui.groupbox_fluid_composition.setEnabled(True)
        BC0 = self.bcs_current_indices[0]
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
            key = 'bc_x_g'
            le.key = key
            le.args = [self.bcs_current_indices, row+1]
            le.dtype = float
            self.add_tooltip(le, key)
            val = self.project.get_value(key, args=[BC0, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)
            le.value_updated.connect(self.handle_bcs_fluid_mass_fraction)
            table.setCellWidget(row, 1, le)
        if species:
            table.setItem(nrows-1, 0, make_item("Total"))
            table.setItem(nrows-1, 1, make_item(''))
            item = table.item(nrows-1, 0)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            self.update_bcs_fluid_mass_fraction_total()
        self.fixup_bcs_table(table)


    def handle_bcs_fluid_mass_fraction(self, widget, value_dict, args):
        ui = self.ui.boundary_conditions
        key = 'bc_x_g'
        val = value_dict[key]
        widget.updateValue(key, val)
        if val == '':
            self.unset_keyword(key, args=args)
        else:
            self.update_keyword(key, val, args=args)
        self.update_bcs_fluid_mass_fraction_total()


    def update_bcs_fluid_mass_fraction_total(self):
        if not self.bcs_current_indices:
            return
        if not self.fluid_species:
            return
        BC0 = self.bcs_current_indices[0]
        ui = self.ui.boundary_conditions
        key = 'bc_x_g'
        table = ui.tablewidget_fluid_mass_fraction
        if table.rowCount() == 0:
            return
        total = sum(float(self.project.get_value(key, default=0.0, args=[BC0,i]))
                    for i in range(1,len(self.fluid_species)+1))
        item =  table.item(table.rowCount()-1, 1)
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setText(str(total))

        # DEFAULT - last defined species has mass fraction of 1.0
        # (only enforce this if no mass fractions are set)

        if total == 0.0 and self.fluid_species:
            for BC in self.bcs_current_indices:
                for i in range(1, len(self.fluid_species)):
                    self.update_keyword('bc_x_g', 0.0, args=[BC, i])
                self.update_keyword('bc_x_g', 1.0, args=[BC, len(self.fluid_species)]) # Last defined species
            self.update_bcs_fluid_mass_fraction_table()


    def setup_bcs_fluid_inflow_tab(self):
        #Subtask Pane Tab for INFLOW type (MI, PI, CG_MI) Boundary Condition Regions
        #Fluid (tab)
        ui = self.ui.boundary_conditions
        ui.page_fluid.setCurrentIndex(PAGE_INFLOW)

        if not self.bcs_current_indices:
            return # Nothing selected.  (Clear out all lineedits?)
        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC',
                        'lineedit_%s_args_BC'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC'):
                name = pat%(key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC))
            get_widget(key+suffix).updateValue(key, val, args=args)

        #    Define volume fraction
        # Specification always available
        # Sets keyword BC_EP_G(#)
        #  DEFAULT value of 1.0 for MI and CG_MI; leave [UNDEFINED] for PI
        #  Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0 #(see bcs_set_volume_fraction_limit)
        # TODO  Error Check: For PI - either all are defined and sum to 1, or all are undefined
        enabled = True
        key = 'bc_ep_g'
        default = 1.0 if bc_type in ('MI', 'CG_MI') else None
        setup_key_widget(key, default, enabled)
        get_widget(key).setReadOnly(True)
        get_widget(key).setEnabled(False) # better way to indicate read-only?

        #    Define inflow properties: Mass inflow specification changes
        # based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
        #
        region_name = self.bcs[BC0].get('region')
        if not region_name:  # should not happen!
            self.error("No region defined for BC %s" % BC0)
            return
        region_info = self.bcs_region_dict.get(region_name)
        if not region_info:
            self.error("No definition for region %s" % region_name)
            return
        region_type = region_info.get('type')
        if not region_type:
            self.error("No type for region %s" % region_name)
            return

        # For BC_TYPE='MI' and XZ-Plane region
        # For BC_TYPE='MI' and YZ-Plane region
        # For BC_TYPE='MI' and XY-Plane region
        #  Select mass inflow specification type:
        if bc_type == 'MI':
            ui.stackedwidget_fluid_inflow.setCurrentIndex(0) # subpage_fluid_inflow_MI
            ui.label_fluid_tangential_velocities.show()
            if not region_type.endswith('plane'):
                self.error("Invalid type %s for region %s" % (region_type, region_name))

            tangents = [C for C in 'XYZ' if C in region_type]
            normal = [C for C in 'XYZ' if C not in tangents]
            if len(normal) != 1 or len(tangents) != 2:
                self.error("Invalid type %s for region" % (region_type, region_name))
            normal = normal[0]

            cb = ui.combobox_fluid_inflow_type
            item = get_combobox_item(cb, 0)
            item.setText('%s-axial velocity' % normal)
            inflow_type = self.bcs[BC0].get('fluid_inflow_type')
            keys = ['bc_%s_g'%xmap[normal], 'bc_volflow_g', 'bc_massflow_g']
            vals = [self.project.get_value(k, args=[BC0]) for k in keys]
            if inflow_type is None:
                vel_flow, vol_flow, mass_flow = vals
                if (vel_flow is not None) and (vol_flow is None) and (mass_flow is None):
                    inflow_type = 0
                elif (vel_flow is None) and (vol_flow is not None) and (mass_flow is None):
                    inflow_type = 1
                elif (vel_flow is None) and (vol_flow is None) and (mass_flow is not None):
                    inflow_type = 2
                else: # warn?
                    inflow_type = 0 # default
                for BC in self.bcs_current_indices:
                    self.bcs[BC]['fluid_inflow_type'] = inflow_type

            cb.setCurrentIndex(inflow_type)
            key = keys[inflow_type]
            val = vals[inflow_type]
            le = ui.lineedit_fluid_inflow
            le.key = key
            le.args = ['BC']
            le.dtype = float
            self.add_tooltip(le, key)
            le.updateValue(key , 0.0 if val is None else val)
            units = ['m/s', 'm³/s', 'kg/s'][inflow_type]
            ui.label_fluid_inflow_units.setText(units)

            #  Define Tangential Velocities:
            #Define X-Axial Velocity
            # Sets keyword BC_U_G(#)
            # DEFAULT value of 0.0
            key = 'bc_%s_g'%xmap[tangents[0]]
            label = ui.label_fluid_tangential_velocity_1
            default = 0.0
            label.setText('  %s-axial velocity' % tangents[0])
            self.add_tooltip(label, key)
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget = ui.lineedit_fluid_tangential_velocity_1
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_G(#)
            # DEFAULT value of 0.0
            key = 'bc_%s_g'%xmap[tangents[1]]
            label = ui.label_fluid_tangential_velocity_2
            label.setText('  %s-axial velocity' % tangents[1])
            self.add_tooltip(label, key)
            default = 0.0
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget = ui.lineedit_fluid_tangential_velocity_2
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            widget.updateValue(key, val)

            self.add_tooltip(widget, key)

        # For BC_TYPE='CG_MI' or 'PI'
        elif bc_type in ('CG_MI', 'PI'):
            #  Specify all velocity components:
            ui.stackedwidget_fluid_inflow.setCurrentIndex(1) # subpage_fluid_inflow_CG_MI
            ui.label_fluid_tangential_velocities.hide()
            #    Define X-Axial Velocity
            # Sets keyword BC_U_G(#)
            # DEFAULT value of 0.0
            key = 'bc_u_g'
            default = 0.0
            widget = ui.lineedit_fluid_inflow
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Y-Axial Velocity
            # Sets keyword BC_V_G(#)
            # DEFAULT value of 0.0
            key = 'bc_v_g'
            label =  ui.label_fluid_tangential_velocity_1
            label.setText('Y-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_fluid_tangential_velocity_1
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword  BC_W_G
            # DEFAULT value of 0.0
            key = 'bc_w_g'
            label =  ui.label_fluid_tangential_velocity_2
            label.setText('Z-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_fluid_tangential_velocity_2
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

        #Define temperature
        # Specification always available
        # Input required for any of the following
        #  Fluid density model: Ideal Gas Law
        #  Fluid viscosity model: Sutherland's Law
        #  Energy equations are solved
        # Sets keyword BC_T_G(#)
        # DEFAULT value of 293.15
        energy_eq = self.project.get_value('energy_eq', default=True)
        enabled = (self.fluid_density_model==OTHER
                   or self.fluid_viscosity_model==OTHER
                   or bool(energy_eq))
        key = 'bc_t_g'
        default = 293.15
        setup_key_widget(key, default, enabled)

        #Define pressure
        # Specification always available
        # Input required when combining ideal gas law and specified mass inflow rate
        # Input required for BC_TYPE = PI
        # TODO "required"
        # Sets keyword BC_P_G(#)
        # DEFAULT 101.325d3
        enabled = True
        key = 'bc_p_g'
        default = FloatExp('101.325e3')
        setup_key_widget(key, default, enabled)

        #Select species and set mass fractions (table format)
        # Specification always available
        # TODO Input required for species equations
        # Drop down menu of fluid species
        # Sets keyword BC_X_G(#,#)
        # DEFAULT - last defined species has mass fraction of 1.0
        # TODO Error check: mass fractions must sum to one
        #
        # Move back the table, in case it got stolen
        comp = ui.groupbox_fluid_composition
        parent = comp.parent()
        if parent != ui.page_fluid_inflow:
            comp.hide()
            parent.layout().removeWidget(comp)
            layout = ui.page_fluid_inflow.layout()
            layout.insertWidget(layout.count()-2, comp)
            comp.show()
        self.update_bcs_fluid_mass_fraction_table()

        #Turbulence: Define k-ε turbulent kinetic energy
        turbulence_model = self.project.get_value('turbulence_model')
        enabled = (turbulence_model == 'K_EPSILON')
        ui.groupbox_inflow_turbulence.setEnabled(enabled)

        # Specification only available with K-Epsilon turbulence model
        # Sets keyword BC_K_TURB_G(#)
        # DEFAULT value of 0.0
        default = 0.0 if enabled else None
        key = 'bc_k_turb_g'
        setup_key_widget(key, default, enabled)

        #Turbulence: Define k-ε turbulent dissipation
        # Specification only available with K-Epsilon turbulence model
        # Sets keywords BC_E_TURB_G(#)
        # DEFAULT value of 0.0
        key = 'bc_e_turb_g'
        default = 0.0 if enabled else None
        setup_key_widget(key, default, enabled)


    # DRY out fluid/solids code
    def update_bcs_solids_mass_fraction_table(self):
        ui = self.ui.boundary_conditions
        table = ui.tablewidget_solids_mass_fraction
        table.clearContents()
        table.setRowCount(0)
        P = self.bcs_current_solid
        if not (P and self.solids_species.get(P) and self.bcs_current_indices):
            self.fixup_bcs_table(table)
            table.setEnabled(False)
            ui.groupbox_solids_composition.setEnabled(False)
            return
        ui.groupbox_solids_composition.setEnabled(True)
        table.setEnabled(True)
        BC0 = self.bcs_current_indices[0]
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
            key = 'bc_x_s'
            le = LineEdit()
            le.key = key
            le.args = [self.bcs_current_indices, P, row+1]
            le.dtype = float
            le.setValInfo(min=0.0, max=1.0)
            self.add_tooltip(le, key)
            val = self.project.get_value(key, args=[BC0, P, row+1], default=None)
            if val is not None:
                le.updateValue(key, val)
            le.value_updated.connect(self.handle_bcs_solids_mass_fraction)
            table.setCellWidget(row, 1, le)
        if species:
            table.setItem(nrows-1, 0, make_item("Total"))
            table.setItem(nrows-1, 1, make_item(''))
            item = table.item(nrows-1, 0)
            font = item.font()
            font.setBold(True)
            item.setFont(font)
            self.update_bcs_solids_mass_fraction_total()

        self.fixup_bcs_table(table)


    def handle_bcs_solids_mass_fraction(self, widget, value_dict, args):
        ui = self.ui.boundary_conditions
        key = 'bc_x_s'
        val = value_dict[key]
        widget.updateValue(key, val)
        if val == '':
            self.unset_keyword(key, args=args)
        else:
            self.update_keyword(key, val, args=args)
        self.update_bcs_solids_mass_fraction_total()


    def update_bcs_solids_mass_fraction_total(self):
        if not self.bcs_current_indices:
            return
        BC0 = self.bcs_current_indices[0]
        P = self.bcs_current_solid
        if P is None:
            return
        species = self.solids_species.get(P)
        if not P:
            return
        ui = self.ui.boundary_conditions
        key = 'bc_x_s'
        table = ui.tablewidget_solids_mass_fraction
        if table.rowCount() == 0:
            return
        total = sum(float(self.project.get_value(key, default=0.0, args=[BC0,P,i]))
                    for i in range(1,len(species)+1))
        item = table.item(table.rowCount()-1, 1)
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setText(str(total))

        # DEFAULT - last defined species has mass fraction of 1.0
        # (only enforce this if no mass fractions are set)
        if total == 0.0 and species:
            for BC in self.bcs_current_indices:
                for i in range(1, len(species)):
                    self.update_keyword('bc_x_s', 0.0, args=[BC, P, i])
                self.update_keyword('bc_x_s', 1.0, args=[BC, P, len(species)]) # Last defined species
            self.update_bcs_solids_mass_fraction_table()


    def setup_bcs_solids_inflow_tab(self, P):
        #Subtask Pane Tab for INFLOW type (MI, PI, CG_MI) Boundary Condition Regions
        #Solid-# (tab) - Rename tab to user provided solids name.
        ui = self.ui.boundary_conditions
        ui.page_solids.setCurrentIndex(PAGE_INFLOW)

        self.bcs_current_solid = self.P = P
        if P is None: # Nothing to do
            return
        if not self.bcs_current_indices:
            return # Nothing selected.  (Clear out all lineedits?)

        BC0 = self.bcs_current_indices[0]
        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        self.bcs_set_volume_fraction_limit()

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
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC_P',
                         'lineedit_keyword_%s_args_BC',
                         'lineedit_%s_args_BC_P',
                         'lineedit_%s_args_BC',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)

            args = mkargs(key, bc=BC0, phase=P)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
            for BC in self.bcs_current_indices:
                self.update_keyword(key, val, args=mkargs(key, bc=BC, phase=P))
            get_widget(key+suffix).updateValue(key, val, args=args)

        #Define volume fraction
        # Specification always available
        # Sets keyword BC_EP_S(#,#)
        #  DEFAULT value of 1.0 - (sum of previous tabs) for MI and CG_MI; leave [UNDEFINED] for PI
        # TODO Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0
        # TODO Error Check: For PI - either all are defined and sum to 1, or all are undefined
        enabled = True
        key = 'bc_ep_s'
        s = sum(self.project.get_value(key, default=0, args=[BC0, p]) for p in range(1, P))
        default = (1.0 - s) if bc_type in ('MI', 'CG_MI') else None
        setup_key_widget(key, default, enabled)

        #Define inflow properties: Mass inflow specification changes
        #based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
        #
        region_name = self.bcs[BC0].get('region')
        if not region_name:  # should not happen!
            self.error("No region defined for BC %s" % BC0)
            return
        region_info = self.bcs_region_dict.get(region_name)
        if not region_info:
            self.error("No definition for region %s" % region_name)
            return
        region_type = region_info.get('type')
        if not region_type:
            self.error("No type for region %s" % region_name)
            return

        # For BC_TYPE='MI' and XZ-Plane region
        # For BC_TYPE='MI' and YZ-Plane region
        # For BC_TYPE='MI' and XY-Plane region
        #  Select mass inflow specification type:
        if bc_type == 'MI':
            ui.stackedwidget_solids_inflow.setCurrentIndex(0) # subpage_solids_inflow_MI
            ui.label_solids_tangential_velocities.show()
            if not region_type.endswith('plane'):
                self.error("Invalid type %s for region %s" % (region_type, region_name))

            tangents = [C for C in 'XYZ' if C in region_type]
            normal = [C for C in 'XYZ' if C not in tangents]
            if len(normal) != 1 or len(tangents) != 2:
                self.error("Invalid type %s for region" % (region_type, region_name))
            normal = normal[0]

            cb = ui.combobox_solids_inflow_type
            item = get_combobox_item(cb, 0)
            item.setText('%s-axial velocity' % normal)
            inflow_type = self.bcs[BC0].get('solids_%s_inflow_type'%P)
            keys = ['bc_%s_s'%xmap[normal], 'bc_volflow_s', 'bc_massflow_s']
            vals = [self.project.get_value(k, args=[BC0, P]) for k in keys]
            if inflow_type is None:
                vel_flow, vol_flow, mass_flow = vals
                if (vel_flow is not None) and (vol_flow is None) and (mass_flow is None):
                    inflow_type = 0
                elif (vel_flow is None) and (vol_flow is not None) and (mass_flow is None):
                    inflow_type = 1
                elif (vel_flow is None) and (vol_flow is None) and (mass_flow is not None):
                    inflow_type = 2
                else: # warn?
                    inflow_type = 0 # default
                for BC in self.bcs_current_indices:
                    self.bcs[BC]['solid_%s_inflow_type'%P] = inflow_type

            cb.setCurrentIndex(inflow_type)
            key = keys[inflow_type]
            val = vals[inflow_type]
            le = ui.lineedit_solids_inflow
            le.key = key
            le.args = ['BC', 'P']
            le.dtype = float
            self.add_tooltip(le, key)
            le.updateValue(key, 0.0 if val is None else val)
            units = ['m/s', 'm³/s', 'kg/s'][inflow_type]
            ui.label_solids_inflow_units.setText(units)

            #  Define Tangential Velocities:
            #    Define X-Axial Velocity
            # Sets keyword BC_U_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_%s_s'%xmap[tangents[0]]
            default = 0.0
            label = ui.label_solids_tangential_velocity_1
            label.setText('  %s-axial velocity' % tangents[0])
            self.add_tooltip(label, key)
            val = self.project.get_value(key, args=[BC0, P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget = ui.lineedit_solids_tangential_velocity_1
            widget.key = key
            widget.args = ['BC', 'P']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_G(#)
            # DEFAULT value of 0.0
            key = 'bc_%s_s'%xmap[tangents[1]]
            label = ui.label_solids_tangential_velocity_2
            label.setText('  %s-axial velocity' % tangents[1])
            self.add_tooltip(label, key)
            default = 0.0
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget = ui.lineedit_solids_tangential_velocity_2
            widget.key = key
            widget.args = ['BC', 'P']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

        # For BC_TYPE='CG_MI' or 'PI'
        elif bc_type in ('CG_MI', 'PI'):
            #  Specify all velocity components:
            ui.stackedwidget_solids_inflow.setCurrentIndex(1) # subpage_solids_inflow_CG_MI
            ui.label_solids_tangential_velocities.hide()
            #    Define X-Axial Velocity
            # Sets keyword BC_U_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_u_s'
            default = 0.0
            widget = ui.lineedit_solids_inflow
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.updateValue(key, val)
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Y-Axial Velocity
            # Sets keyword BC_V_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_v_s'
            label =  ui.label_solids_tangential_velocity_1
            label.setText('Y-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_solids_tangential_velocity_1
            val = self.project.get_value(key, args=[BC0, P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.key = key
            widget.args = ['BC', 'P']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_w_s'
            label =  ui.label_solids_tangential_velocity_2
            label.setText('Z-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_solids_tangential_velocity_2
            val = self.project.get_value(key, args=[BC0, P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.key = key
            widget.args = ['BC', 'P']
            widget.dtype = float
            self.add_tooltip(widget, key)

        #Define temperature
        # Specification always available
        # TODO Input required when energy equations are solved
        # Sets keyword BC_T_S(#,#)
        # DEFAULT value of 293.15
        enabled = True
        key = 'bc_t_s'
        default = 293.15
        setup_key_widget(key, default, enabled)

        #Select species and set mass fractions (table format)
        # Specification always available
        # Input required for species equations
        # Drop down menu of solids species
        # Sets keyword BC_X_S(#,#,#)
        # DEFAULT - last defined species has mass fraction of 1.0
        # TODO Error check: mass fractions must sum to one
        self.update_bcs_solids_mass_fraction_table()


    def setup_bcs_scalar_inflow_tab(self):
        #Subtask Pane Tab for INFLOW type (MI, PI, CG_MI) Boundary Condition Regions
        #Scalar (tab) - Tab only available if scalar equations are solved
        #    Define initial scalar value
        # Sets keyword BC_SCALAR(#,#)
        # DEFAULT value of 0.0
        # TODO implement this tab
        pass


    def setup_bcs_fluid_po_tab(self):
        #Subtask Pane Tab for PRESSURE OUTFLOW type (PO) Boundary Condition Regions
        #Fluid (tab)
        ui = self.ui.boundary_conditions
        ui.page_fluid.setCurrentIndex(PAGE_PO)

        if not self.bcs_current_indices:
            return
        BC0 = self.bcs_current_indices[0]

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC',
                        'lineedit_%s_args_BC',
                        'lineedit_keyword_%s',
                        'lineedit_%s'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC',
                         'lineedit_%s_args_BC',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC))
            get_widget(key+suffix).updateValue(key, val, args=args)

        #    Define pressure
        # Specification always available
        # TODO Input required
        # Sets keyword BC_P_G(#)
        # DEFAULT 101.325d3
        enabled = True
        key = 'bc_p_g'
        default = FloatExp('101.325e3')
        setup_key_widget(key, default, enabled, suffix='_2')

        #The remaining inputs are "optional." They do not have default values, because MFIX will calculate
        #appropriate values if they are unspecified and 'backflow' occurs at the outlet.
        #
        #    Define volume fraction
        # Specification always available
        # Sets keyword BC_EP_G(#)
        # No DEFAULT value
        # Error Check: If any volume fraction for the BC region is specified, then all volume fractions
        # for the BC region must be specified and must sum to one.
        # TODO implement the above error check.  Note that not all fractions will be set during the
        #  setup phase.
        enabled = True
        key = 'bc_ep_g'
        default = None
        setup_key_widget(key, default, enabled, suffix='_2')
        get_widget(key+'_2').setReadOnly(True)
        get_widget(key+'_2').setEnabled(False) # better way to indicate read-only?

        #    Define temperature
        # Specification always available
        # NO DEFAULT value
        # Sets keyword BC_T_G(#)
        enabled = True
        default = None
        key = 'bc_t_g'
        setup_key_widget(key, default, enabled, suffix='_2')

        #    Select species and set mass fractions (table format)
        # Specification always available
        # NO DEFAULT value
        # Sets keyword BC_X_G(#,#)
        # Error check: if specified, mass fractions must sum to one

        # We already have a mass fraction table - just steal it and move it here
        comp = ui.groupbox_fluid_composition
        parent = comp.parent()
        if parent != ui.page_fluid_po:
            comp.hide()
            parent.layout().removeWidget(comp)
            # (Should we put it inside the 'optional' box?  or underneath?)
            #ui.groupbox_fluid_po_optional.layout().addWidget(comp, 2, 0, 1, 3)
            layout = ui.page_fluid_po.layout() #Underneath - nested groupbox looks weird
            layout.insertWidget(layout.count()-1, comp)
            comp.show()
        self.update_bcs_fluid_mass_fraction_table()


    def setup_bcs_solids_po_tab(self, P):
        #Subtask Pane Tab for PRESSURE OUTFLOW type (PO) Boundary Condition Regions
        #    Solids-# (tab)
        ui = self.ui.boundary_conditions
        ui.page_solids.setCurrentIndex(PAGE_PO)

        self.bcs_current_solid = self.P = P
        if P is None:
            return
        if not self.bcs_current_indices:
            return # Nothing selected.  (Clear out all lineedits?)

        BC0 = self.bcs_current_indices[0]

        self.bcs_set_volume_fraction_limit()

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC_P',
                        'lineedit_%s_args_BC_P',
                        'lineedit_keyword_%s',
                        'lineedit_%s'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC_P',
                         'lineedit_%s_args_BC_P',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = mkargs(key, bc=BC0, phase=P)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC, phase=P))
            get_widget(key+suffix).updateValue(key, val, args=args)

        #All inputs are optional. They do not have default values, because MFIX will calculate
        #appropriate values if they are unspecified and 'backflow' occurs at the outlet.
        #  (label in UI indicates this)

        #    Define volume fraction
        # Specification always available
        # Sets keyword BC_EP_S(#,#)
        # No DEFAULT value
        # Error Check: If any volume fraction for the BC region is specified, then all
        # volume fractions for the BC region must be specified and must sum to one.
        key = 'bc_ep_s'
        default = None
        enabled= True
        setup_key_widget(key, default, enabled, suffix='_2')

        #    Define temperature
        # Specification always available
        # NO DEFAULT value
        # Sets keyword BC_T_S(#,#)
        key = 'bc_t_s'
        default = None
        enabled= True
        setup_key_widget(key, default, enabled, suffix='_2')

        #    Select species and set mass fractions (table format)
        # Specification always available
        # NO DEFAULT value
        # Sets keyword BC_X_S(#,#,#)
        # Error check: if specified, mass fractions must sum to one
        #
        # We already have a mass fraction table - just steal it and move it here
        comp = ui.groupbox_solids_composition
        parent = comp.parent()
        if parent != ui.page_solids_po:
            comp.hide()
            parent.layout().removeWidget(comp)
            # (Should we put it inside the 'optional' box?  or underneath?)
            layout = ui.page_solids_po.layout() #Underneath - nested groupbox looks weird
            layout.insertWidget(layout.count()-1, comp)
            comp.show()
        self.update_bcs_solids_mass_fraction_table()


    def setup_bcs_scalar_po_tab(self):
        #Subtask Pane Tab for PRESSURE OUTFLOW type (PO) Boundary Condition Regions
        #Scalar (tab) - Tab only available if scalar equations are solved
        #All inputs are optional. They do not have default values, because MFIX will calculate appropriate
        # values if they are unspecified and 'backflow' occurs at the outlet.
        #Define scalar value
        # Sets keyword BC_SCALAR(#,#)
        # NO DEFAULT value
        pass


    def setup_bcs_fluid_mo_tab(self):
        #Subtask Pane Tab for MASS OUTFLOW type (MO) Boundary Condition Regions
        #Fluid (tab)
        #  NB, we use the '_3' suffix in this section for consistency, even though it is not
        #  needed everywhere
        ui = self.ui.boundary_conditions
        ui.page_fluid.setCurrentIndex(PAGE_MO)

        if not self.bcs_current_indices:
            return
        BC0 = self.bcs_current_indices[0]

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC',
                        'lineedit_%s_args_BC',
                        'lineedit_keyword_%s',
                        'lineedit_%s'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC',
                         'lineedit_%s_args_BC',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)
            args = [BC0]
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            get_widget(key+suffix).updateValue(key, val, args=args)

        #    Define outflow properties: Mass outflow specification changes
        # based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
        #
        region_name = self.bcs[BC0].get('region')
        if not region_name:  # should not happen!
            self.error("No region defined for BC %s" % BC0)
            return
        region_info = self.bcs_region_dict.get(region_name)
        if not region_info:
            self.error("No definition for region %s" % region_name)
            return
        region_type = region_info.get('type')
        if not region_type:
            self.error("No type for region %s" % region_name)
            return

        # For BC_TYPE='MO' and XZ-Plane region
        # For BC_TYPE='MO' and YZ-Plane region
        # For BC_TYPE='MO' and XY-Plane region
        #  Select mass outflow specification type:
        if bc_type == 'MO':
            ui.stackedwidget_fluid_outflow.setCurrentIndex(0) # subpage_fluid_outflow_MO
            ui.label_fluid_tangential_velocities_3.show()
            if not region_type.endswith('plane'):
                self.error("Invalid type %s for region %s" % (region_type, region_name))

            tangents = [C for C in 'XYZ' if C in region_type]
            normal = [C for C in 'XYZ' if C not in tangents]
            if len(normal) != 1 or len(tangents) != 2:
                self.error("Invalid type %s for region" % (region_type, region_name))
            normal = normal[0]

            cb = ui.combobox_fluid_outflow_type
            item = get_combobox_item(cb, 0)
            item.setText('%s-axial velocity' % normal)
            outflow_type = self.bcs[BC0].get('fluid_outflow_type')
            keys = ['bc_%s_g'%xmap[normal], 'bc_volflow_g', 'bc_massflow_g']
            vals = [self.project.get_value(k, args=[BC0]) for k in keys]
            if outflow_type is None:
                vel_flow, vol_flow, mass_flow = vals
                if (vel_flow is not None) and (vol_flow is None) and (mass_flow is None):
                    outflow_type = 0
                elif (vel_flow is None) and (vol_flow is not None) and (mass_flow is None):
                    outflow_type = 1
                elif (vel_flow is None) and (vol_flow is None) and (mass_flow is not None):
                    outflow_type = 2
                else: # warn?
                    outflow_type = 0 # default
                for BC in self.bcs_current_indices:
                    self.bcs[BC]['fluid_outflow_type'] = outflow_type

            cb.setCurrentIndex(outflow_type)
            key = keys[outflow_type]
            val = vals[outflow_type]
            le = ui.lineedit_fluid_outflow
            le.key = key
            le.args = ['BC']
            le.dtype = float
            self.add_tooltip(le, key)
            le.updateValue(key, 0.0 if val is None else val)
            units = ['m/s', 'm³/s', 'kg/s'][outflow_type]
            ui.label_fluid_outflow_units.setText(units)

            #  Define Tangential Velocities:
            #    Define X-Axial Velocity
            # Sets keyword BC_U_G(#)
            # DEFAULT value of 0.0
            key = 'bc_%s_g'%xmap[tangents[0]]
            default = 0.0
            label = ui.label_fluid_tangential_velocity_1_3
            label.setText('  %s-axial velocity' % tangents[0])
            self.add_tooltip(label, key)
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget = ui.lineedit_fluid_tangential_velocity_1_3
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_G(#)
            # DEFAULT value of 0.0
            key = 'bc_%s_g'%xmap[tangents[1]]
            label = ui.label_fluid_tangential_velocity_2_3
            label.setText('  %s-axial velocity' % tangents[1])
            self.add_tooltip(label, key)
            default = 0.0
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget = ui.lineedit_fluid_tangential_velocity_2_3
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

        # For BC_TYPE='CG_MO'
        elif bc_type == 'CG_MO':
            #  Specify all velocity components:
            ui.stackedwidget_fluid_outflow.setCurrentIndex(1) # subpage_fluid_outflow_CG_MO
            ui.label_fluid_tangential_velocities_3.hide()

            #    Define X-Axial Velocity
            # Sets keyword BC_U_G(#)
            # DEFAULT value of 0.0
            key = 'bc_u_g'
            default = 0.0
            widget = ui.lineedit_fluid_outflow
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.updateValue(key, val)
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Y-Axial Velocity
            # Sets keyword BC_V_G(#)
            # DEFAULT value of 0.0
            key = 'bc_v_g'
            label =  ui.label_fluid_tangential_velocity_1_3
            label.setText('Y-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_fluid_tangential_velocity_1_3
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_V_G(#) # fixed, BC_W_G
            # DEFAULT value of 0.0
            key = 'bc_w_g'
            label =  ui.label_fluid_tangential_velocity_2_3
            label.setText('Z-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_fluid_tangential_velocity_2
            val = self.project.get_value(key, args=[BC0])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC])
            widget.key = key
            widget.args = ['BC']
            widget.dtype = float
            self.add_tooltip(widget, key)

        #Define duration to average outflow rate.
        # Specification always available
        # Input required
        # Sets keyword BC_DT_0(#)
        # DEFAULT value of 0.1
        # TODO: See comments in mfix_user_guide re: BC_DT_0 and restarting
        enabled = True
        key = 'bc_dt_0'
        default = 0.1
        setup_key_widget(key, default, enabled)

        #The remaining inputs are only required when either the mass or the volumetric flowrates are
        #specified. They are not required if the velocities are given for the outlet.

        #Define volume fraction
        # Specification only available with mass or volumetric flowrates.
        # Input required  # TODO decide if we allow specifiying fluid vol. fractions,
        #                   or set them from 1-bc_ep_s
        # Sets keyword BC_EP_G(#)
        # DEFAULT value 1.0
        # Error Check: If any volume fraction for the BC region is specified, then all volume
        #fractions for the BC region must be specified and must sum to one.
        enabled = True
        key = 'bc_ep_g'
        default = None
        setup_key_widget(key, default, enabled, suffix='_3')
        get_widget(key+'_3').setReadOnly(True)
        get_widget(key+'_3').setEnabled(False) # better way to indicate read-only?

        #Define temperature
        # Specification only available with mass or volumetric flowrates and R_G0 is UNDEFINED (RO_G0)
        # DEFAULT value 293.15
        # Sets keyword BC_T_G(#)
        ro_g0 = self.project.get_value('ro_g0')
        mass_flow = self.project.get_value('bc_massflow_g', args=[BC0])
        vol_flow = self.project.get_value('bc_volflow_g', args=[BC0])
        enabled = (ro_g0 is None) and (mass_flow is not None or vol_flow is not None)
        default = 293.15 if enabled else None
        key = 'bc_t_g'
        setup_key_widget(key, default, enabled, suffix='_3')

        #Select species and set mass fractions (table format)
        # Specification only available with mass or volumetric flowrates and R_G0 is UNDEFINED (RO_G0)
        # DEFAULT value 1.0 of last defined species
        # Sets keyword BC_X_G(#,#)
        # Error check: if specified, mass fractions must sum to one
        #
        enabled = (ro_g0 is None) and (mass_flow is not None or vol_flow is not None)
        comp = ui.groupbox_fluid_composition
        parent = comp.parent()
        if enabled:
            if parent != ui.page_fluid_mo:
                # Steal the table widget from wherever it is
                comp.hide()
                parent.layout().removeWidget(comp)
                layout = ui.page_fluid_mo.layout()
                layout.insertWidget(layout.count()-1, comp)
            comp.show()
            self.update_bcs_fluid_mass_fraction_table()
        elif parent == ui.page_fluid_mo:
            # Table widget is here, but we don't want it
            #comp.hide()
            comp.setEnabled(False)


    def setup_bcs_solids_mo_tab(self, P):
        #Subtask Pane Tab for MASS OUTFLOW type (MO) Boundary Condition Regions
        #Solids-# (tab)
        #  NB, we use the '_4' suffix in this section for consistency, even though it is not
        #  needed everywhere
        ui = self.ui.boundary_conditions
        ui.page_solids.setCurrentIndex(PAGE_MO)

        self.bcs_current_solid = self.P = P
        if P is None:
            return
        if not self.bcs_current_indices:
            return

        BC0 = self.bcs_current_indices[0]

        def get_widget(key):
            for pat in ('lineedit_keyword_%s_args_BC_P',
                        'lineedit_%s_args_BC_P',
                        'lineedit_keyword_%s',
                        'lineedit_%s'):
                widget = getattr(ui, pat % key, None)
                if widget:
                    return widget
            self.error('no widget for key %s' % key)

        def setup_key_widget(key, default=None, enabled=True, suffix=''):
            for pat in ('label_%s', 'label_%s_units',
                         'lineedit_keyword_%s_args_BC_P',
                         'lineedit_%s_args_BC_P',
                         'lineedit_keyword_%s',
                         'lineedit_%s'):
                name = pat % (key+suffix)
                item = getattr(ui, name, None)
                if item:
                    item.setEnabled(enabled)

            args = mkargs(key, bc=BC0, phase=P)
            val = self.project.get_value(key, args=args)
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=mkargs(key, bc=BC, phase=P))
            get_widget(key+suffix).updateValue(key, val, args=args)

        bc_type = self.project.get_value('bc_type', args=[BC0])
        if bc_type is None:
            self.error("bc_type not set for region %s" % BC0)
            return

        #Define outflow properties: Mass outflow specification changes
        #based on the BC_TYPE and Region orientation (e.g., XZ-Plane)
        #
        region_name = self.bcs[BC0].get('region')
        if not region_name:  # should not happen!
            self.error("No region defined for BC %s" % BC0)
            return
        region_info = self.bcs_region_dict.get(region_name)
        if not region_info:
            self.error("No definition for region %s" % region_name)
            return
        region_type = region_info.get('type')
        if not region_type:
            self.error("No type for region %s" % region_name)
            return

        # For BC_TYPE='MO' and XZ-Plane region
        # For BC_TYPE='MO' and YZ-Plane region
        # For BC_TYPE='MO' and XY-Plane region
        if bc_type == 'MO':
            ui.stackedwidget_solids_outflow.setCurrentIndex(0) # subpage_solids_outflow_MO
            ui.label_solids_tangential_velocities_4.show()
            if not region_type.endswith('plane'):
                self.error("Invalid type %s for region %s" % (region_type, region_name))

            tangents = [C for C in 'XYZ' if C in region_type]
            normal = [C for C in 'XYZ' if C not in tangents]
            if len(normal) != 1 or len(tangents) != 2:
                self.error("Invalid type %s for region" % (region_type, region_name))
            normal = normal[0]

            cb = ui.combobox_solids_outflow_type
            item = get_combobox_item(cb, 0)
            item.setText('%s-axial velocity' % normal)
            outflow_type = self.bcs[BC0].get('solid_%s_outflow_type'%P)
            keys = ['bc_%s_s'%xmap[normal], 'bc_volflow_s', 'bc_massflow_s']
            vals = [self.project.get_value(k, args=[BC0,P]) for k in keys]
            if outflow_type is None:
                vel_flow, vol_flow, mass_flow = vals
                if (vel_flow is not None) and (vol_flow is None) and (mass_flow is None):
                    outflow_type = 0
                elif (vel_flow is None) and (vol_flow is not None) and (mass_flow is None):
                    outflow_type = 1
                elif (vel_flow is None) and (vol_flow is None) and (mass_flow is not None):
                    outflow_type = 2
                else: # warn?
                    outflow_type = 0 # default
                for BC in self.bcs_current_indices:
                    self.bcs[BC]['solid_%s_outflow_type'%P] = outflow_type

            cb.setCurrentIndex(outflow_type)
            key = keys[outflow_type]
            val = vals[outflow_type]
            le = ui.lineedit_solids_outflow
            le.key = key
            le.args = ['BC', 'P']
            le.dtype = float
            self.add_tooltip(le, key)
            le.updateValue(key, 0.0 if val is None else val)
            units = ['m/s', 'm³/s', 'kg/s'][outflow_type]
            ui.label_solids_outflow_units.setText(units)

            #  Define Tangential Velocities:
            #    Define X-Axial Velocity
            # Sets keyword BC_U_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_%s_s'%xmap[tangents[0]]
            default = 0.0
            label = ui.label_solids_tangential_velocity_1_4
            label.setText('  %s-axial velocity' % tangents[0])
            self.add_tooltip(label, key)
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget = ui.lineedit_solids_tangential_velocity_1_4
            widget.key = key
            widget.args = ['BC','P']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_%s_s'%xmap[tangents[1]]
            label = ui.label_solids_tangential_velocity_2_4
            label.setText('  %s-axial velocity' % tangents[1])
            self.add_tooltip(label, key)
            default = 0.0
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget = ui.lineedit_solids_tangential_velocity_2_4
            widget.key = key
            widget.args = ['BC','P']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

        # For BC_TYPE='CG_MO'
        elif bc_type == 'CG_MO':
            #  Specify all velocity components:
            ui.stackedwidget_solids_outflow.setCurrentIndex(1) # subpage_solids_outflow_CG_MO
            ui.label_solids_tangential_velocities_4.hide()

            #    Define X-Axial Velocity
            # Sets keyword BC_U_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_u_s'
            default = 0.0
            widget = ui.lineedit_solids_outflow
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.key = key
            widget.args = ['BC','P']
            widget.dtype = float
            widget.updateValue(key, val)
            self.add_tooltip(widget, key)

            #    Define Y-Axial Velocity
            # Sets keyword BC_V_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_v_s'
            label =  ui.label_solids_tangential_velocity_1_4
            label.setText('Y-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_solids_tangential_velocity_1_4
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.key = key
            widget.args = ['BC','P']
            widget.dtype = float
            self.add_tooltip(widget, key)

            #    Define Z-Axial Velocity
            # Sets keyword BC_W_S(#,#)
            # DEFAULT value of 0.0
            key = 'bc_w_s'
            label =  ui.label_solids_tangential_velocity_2_4
            label.setText('Z-axial velocity')
            self.add_tooltip(label, key)
            default = 0.0
            widget = ui.lineedit_solids_tangential_velocity_2_4
            val = self.project.get_value(key, args=[BC0,P])
            if val is None:
                val = default
                for BC in self.bcs_current_indices:
                    self.update_keyword(key, val, args=[BC, P])
            widget.key = key
            widget.args = ['BC','P']
            widget.dtype = float
            self.add_tooltip(widget, key)

        #Define duration to average outflow rate.
        # Specification always available
        # Input required
        # Sets keyword BC_DT_0(#)
        # DEFAULT value of 0.1
        # TODO: See comments in mfix_user_guide re: BC_DT_0 and restarting
        enabled = True
        key = 'bc_dt_0'
        default = 0.1
        setup_key_widget(key, default, enabled, suffix='_4')

        # TODO: No mass fraction table?  No volume fraction?


    def setup_bcs_scalar_mo_tab(self):
        pass # ? not mentioned in SRS
