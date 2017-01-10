# -*- coding: utf-8 -*-

"""Output Task Pane Window"""

from __future__ import print_function, absolute_import, unicode_literals, division

import logging
log = logging.getLogger(__name__)

#import Qt
from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QPushButton, QWidget

from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4
UserRole = QtCore.Qt.UserRole

#local imports
from mfixgui.constants import *
from mfixgui.widgets.base import LineEdit
from mfixgui.tools import keyword_args
from mfixgui.tools.general import (widget_iter,
                                   get_selected_row,
                                   get_combobox_item,
                                   set_item_noedit,
                                   safe_float)

from mfixgui.widgets.base import (BaseWidget, LineEdit, CheckBox)

TAB_BASIC, TAB_VTK, TAB_SPX, TAB_NETCDF = range(4)

MAX_SP = 11

class Output(object):
    #Output Task Pane Window:
    #The output input is split into tabs.
    # Output sub-pane navigation

    def output_change_tab(self, tabnum, to_btn):
        ui = self.ui.output
        self.output_current_tab = tabnum
        self.animate_stacked_widget(
            ui.stackedwidget_output,
            ui.stackedwidget_output.currentIndex(),
            tabnum,
            direction='horizontal',
            line=ui.line_output,
            to_btn=to_btn,
            btn_layout=ui.gridlayout_tab_btns)
        self.setup_output_tab(tabnum)
        for btn in self.output_pushbuttons:
            btn.setChecked(btn==to_btn)
            font = btn.font()
            font.setBold(btn==to_btn)
            btn.setFont(font)


    def init_output(self):
        ui = self.ui.output

        self.output_current_tab = TAB_BASIC

        self.outputs = {} # key: index.  value: data dictionary for point source
        self.output_current_indices = [] # List of PS indices
        self.output_current_regions = [] # And the names of the regions which define them
        self.output_region_dict = None
        self.output_current_solid = self.P = None

        self.output_pushbuttons = (ui.pushbutton_basic,
                                   ui.pushbutton_vtk,
                                   ui.pushbutton_spx,
                                   ui.pushbutton_netcdf)

        # connect tab buttons
        for i, btn in enumerate(self.output_pushbuttons):
            btn.pressed.connect(lambda i=i, btn=btn: self.output_change_tab(i, btn))

        self.init_output_basic_tab()
        self.init_output_vtk_tab()
        self.init_output_spx_tab()
        self.init_output_netcdf_tab()


    def init_output_basic_tab(self):
        ui = self.ui.output

        for item in widget_iter(ui.page_basic):
            if isinstance(item, BaseWidget):
                item.post_update = self.setup_output_basic_tab

        gb = ui.groupbox_write_vtk_files
        key = 'write_vtk_files'
        gb.clicked.connect(self.output_enable_vtk)
        self.add_tooltip(gb, key)

        cb = ui.checkbox_binary_spx_files
        cb.clicked.connect(self.output_enable_spx)

        cb = ui.checkbox_netcdf_files
        cb.clicked.connect(self.output_enable_netcdf)


    def init_output_vtk_tab(self):
        ui = self.ui.output
        #VTK (tab)
        #Icons and table similar to IC/BC/PS/IS for adding VTK regions. This section requires
        #WRITE_VTK_FILES = .TRUE.
        #Icons to add/remove/duplicate regions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #a VTK region.

        ui.toolbutton_add.clicked.connect(self.output_show_regions_popup)
        ui.toolbutton_delete.clicked.connect(self.output_delete_regions)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        ui.tablewidget_regions.itemSelectionChanged.connect(self.handle_output_region_selection)

        #Select Output type
        # Selection is required
        # Available selections:
        #  Cell data
        #    Selection always available
        #    Set keyword VTK_DATA(#) to 'C'
        #  Particle data
        #    Selection only available with DEM or PIC solids
        #    Sets keyword VTK_DATA(#) to 'P'

    def handle_output_region_selection(self):
        ui = self.ui.output
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = tw.item(row,0).data(UserRole)
        self.output_current_indices, self.output_current_regions = indices, regions
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.detail_pane.setEnabled(enabled)
        if not enabled:
            # Clear
            for widget in widget_iter(ui.detail_pane):
                if isinstance(widget, LineEdit):
                    widget.setText('')
            return
        self.setup_output_vtk_tab() # reinitialize all widgets in current tab


    def fixup_output_table(self, tw, stretch_column=0):
        ui = self.ui.output
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
        ui.top_frame.setMaximumHeight(height+24)
        ui.top_frame.setMinimumHeight(header_height+24)
        ui.top_frame.updateGeometry()
        tw.setMaximumHeight(height)
        tw.setMinimumHeight(header_height)
        tw.updateGeometry()


    def output_show_regions_popup(self):
        # Users cannot select inapplicable regions.
        # VTK regions can be points, planes, or volumes (not STLs)
        # Regions can define multiple VTK regions.

        ui = self.ui.output
        rp = self.regions_popup
        rp.clear()
        for (name,data) in self.output_region_dict.items():
            shape = data.get('type', '---')
            # Assume available if unmarked
            available = (data.get('available', True)
                         and not self.check_region_in_use(name)
                         and (shape in ('point', 'box')
                              or 'plane' in shape))
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.output_add_regions)
        rp.cancel.connect(self.output_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.detail_pane,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup('VTK output')


    def output_delete_solids_phase(self, phase):
        pass # XXX TODO WRITEME


    def output_cancel_add(self):
        ui = self.ui.output

        for item in (ui.toolbutton_add,
                     ui.tablewidget_regions):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_regions) is not None:
            for item in (ui.detail_pane,
                         ui.toolbutton_delete):
                item.setEnabled(True)


    def output_add_regions(self):
        # Interactively add regions to define VTK Output
        ui = self.ui.output
        rp = self.regions_popup
        self.output_cancel_add() # Reenable input widgets
        selections = rp.get_selection_list()
        if not selections:
            return
        self.output_add_regions_1(selections, indices=None, autoselect=False)
        self.setup_output_current_tab() # Update the widgets


    def output_add_regions_1(self, selections, indices=None, autoselect=False):
        # Used by both interactive and load-time add-region handlers
        ui = self.ui.output

        if self.output_region_dict is None:
            self.output_region_dict = self.ui.regions.get_region_dict()

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
                idx = self.output_find_index()
                indices[i] = idx
            self.outputs[idx] = {'region': region_name}
            region_data = self.output_region_dict.get(region_name)
            if region_data is None: # ?
                self.warn("no data for region %s" % region_name)
                continue
            self.output_region_dict[region_name]['available'] = False # Mark as in-use
        item.setData(UserRole, (tuple(indices), tuple(selections)))
        tw.setItem(nrows, 0, item)

        #self.fixup_output_table(tw)

        if autoselect:
            tw.setCurrentCell(nrows, 0)


    def output_find_index(self):
        n = 1
        while n in self.outputs:
            n += 1
        return n


    def output_delete_regions(self):
        ui = self.ui.output
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        # Unset keywords
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            if key.startswith('vtk_') and args and args[0] in self.output_current_indices:
                self.unset_keyword(key, args=args)

        for r in self.output_current_regions:
            if r in self.output_region_dict:
                self.output_region_dict[r]['available'] = True

        for i in self.output_current_indices:
            del self.outputs[i]

        self.output_current_regions = []
        self.output_current_indices = []

        tw.removeRow(row)
        #self.fixup_output_table(tw)
        self.setup_output_vtk_tab()
        #self.update_nav_tree()



    def init_output_spx_tab(self):
        ui = self.ui.output
        gb = ui.groupbox_print_des_data # aka "Write ASCII particle data"
        gb.clicked.connect(lambda val:  self.update_keyword('print_des_data', val))
        self.add_tooltip(gb, key = 'print_des_data')

        cb = ui.combobox_des_output_type
        cb.currentIndexChanged.connect(lambda val: self.update_keyword('des_output_type', DES_OUTPUT_TYPES[val]))
        self.add_tooltip(cb, key='des_output_type')


    def init_output_netcdf_tab(self):
        ui = self.ui.output


    def output_enable_vtk(self, enabled):
        self.update_keyword('write_vtk_files', enabled)
        self.setup_output() # enable/disable gui widgets


    def output_enable_spx(self, enabled):
        ui = self.ui.output
        if enabled:
            res_dt = self.project.get_value('res_dt', default=1.0)
            default = max(1.0, res_dt)
            for i in range(1, MAX_SP+1): # Only set keys not already set
                if self.project.get_value('spx_dt', args=[i]) is None:
                    # Try to restore value from disabled lineedit
                    le = getattr(ui, 'lineedit_keyword_spx_dt_args_%s' % i)
                    val = le.value
                    if val == '' or val is None or val == 0:
                        val = default
                    self.update_keyword('spx_dt', val, args=[i])
        else:
            for i in range(1, MAX_SP+1):
                self.unset_keyword('spx_dt', args=[i])
        self.setup_output()


    def output_enable_netcdf(self, enabled):
        ui = self.ui.output


    def setup_output(self):
        ui = self.ui.output
        # Grab a fresh copy, may have been updated
        self.output_region_dict = self.ui.regions.get_region_dict()

        spx_enabled = any(self.project.get_value('spx_dt', args=[i]) is not None
                          for i in range(1,MAX_SP+1))
        ui.pushbutton_spx.setEnabled(spx_enabled)

        vtk_enabled = self.project.get_value('write_vtk_files', default=False)
        ui.pushbutton_vtk.setEnabled(vtk_enabled)

        netcdf_enabled = False
        ui.pushbutton_netcdf.setEnabled(netcdf_enabled)

        # TODO don't stay on disabled tab
        self.setup_output_current_tab()


    def setup_output_tab(self, tabnum):
        if tabnum == TAB_BASIC:
            self.setup_output_basic_tab()
        elif tabnum == TAB_VTK:
            self.setup_output_vtk_tab()
        elif tabnum == TAB_SPX:
            self.setup_output_spx_tab()
        elif tabnum == TAB_NETCDF:
            self.setup_output_netcdf_tab()
        else:
            raise ValueError(tabnum)


    def setup_output_current_tab(self):
        self.setup_output_tab(self.output_current_tab)


    def setup_output_basic_tab(self):
        #Basic (tab)
        ui = self.ui.output
        #
        #Specify Restart/checkpoint write interval
        #    Specification always available (required)
        #    Sets keyword RES_DT
        #    DEFAULT value of 1.0
        key = 'res_dt'
        res_dt = self.project.get_value(key, default=1.0)
        # 'reverse constraint', res_dt must be less than a bunch of other keys
        # (which must be greater than it)

        m = min(safe_float(self.project.get_value('spx_dt', default=1.0, args=[i]))
                for i in range(1,MAX_SP+1))
        if self.project.get_value('res_backups', default=0) > 0:
            m = min(m, safe_float(self.project.get_value('res_backup_dt', default=1.0)))
        ui.lineedit_keyword_res_dt.max = m
        #if res_dt > m:
        #    res_dt = m
        #    self.update_keyword(key, res_dt)

        #Specify the number of backup copies
        #    Specification always available
        #    Sets keyword RES_BACKUPS
        #    DEFAULT value of 0
        #    Error check: value is greater than or equal to 0

        #Specify the backup interval
        #    Specification only available if RES_BACKUPS > 0
        enabled = self.project.get_value('res_backups', default=0) > 0
        #    Sets keyword RES_BACKUP_DT
        key = 'res_backup_dt'
        #    DEFAULT value of 1.0
        default = 1.0
        #    Error check: value must be greater than or equal to RES_DT
        ui.lineedit_keyword_res_backup_dt.min = res_dt
        for item in (ui.label_res_backup_dt,
                     ui.lineedit_keyword_res_backup_dt,
                     ui.label_res_backup_dt_units):
            item.setEnabled(enabled)
            if enabled:
                val = self.project.get_value(key)
                if val is None:
                    val = default
                    if val < res_dt:
                        val = res_dt
                    self.update_keyword(key, val)

        #Enable VTK output
        #    Specification always available
        #    Sets keyword WRITE_VTK_FILES
        #    DEFAULT value of .FALSE.
        #    Enables VTK tab
        # (handled by keyword widget and post_update)
        write_vtk_files = self.project.get_value('write_vtk_files', default=False)
        ui.groupbox_write_vtk_files.setChecked(write_vtk_files)
        ui.pushbutton_vtk.setEnabled(write_vtk_files)

        #Enable time-dependent VTK files
        #    Specification only if WRITE_VTK_FILES = .TRUE.
        enabled = bool(write_vtk_files)

        #    Sets keyword TIME_DEPENDENT_FILENAME
        key = 'time_dependent_filename'
        #    DEFAULT value .TRUE.
        default = True
        cb = ui.checkbox_keyword_time_dependent_filename
        #cb.setEnabled(enabled)        #(handled by checkable groupbox)
        value = self.project.get_value(key)
        self.add_tooltip(cb, key, value=True if value is None else value)
        if enabled:
            if value is None:
                value = default
                self.update_keyword(key, value)
        else:
            self.unset_keyword(key)
        ui.pushbutton_vtk.setEnabled(enabled)

        #Specify VTK Directory
        #    Specification only if WRITE_VTK_FILES = .TRUE.
        #    Sets keyword VTU_DIR
        #    No default (empty string)
        key = 'vtu_dir'
        enabled = bool(write_vtk_files)
        value = self.project.get_value(key, default='')
        #for item in (ui.label_vtu_dir, ui.lineedit_keyword_vtu_dir):
        #    item.setEnabled(enabled)        #(handled by checkable groupbox)
        if enabled:
            if value is None or value=='':
                value = ui.lineedit_keyword_vtu_dir.value # saved value in GUI
                if value is None:
                    value = default
                self.update_keyword(key, value)
        else:
            self.unset_keyword(key)
        ui.pushbutton_vtk.setEnabled(enabled)

        #Write binary Single Precision files (SPx)
        #    No keyword association
        #    Enables SPx tab
        #    Backwards compatibility: Enabled if any SPx time values are specified
        enabled = any(self.project.get_value('spx_dt', args=[i]) is not None
                      for i in range(1,MAX_SP+1)) # Note, enabled in template! *** Jordan?
        ui.checkbox_binary_spx_files.setChecked(enabled)
        ui.pushbutton_spx.setEnabled(enabled)

        #Enable NetCDF output files
        #    Not available when SPx output is enabled
        #    No keyword association.
        #    Enables NetCDF tab



    def setup_output_spx_tab(self):
        ui = self.ui.output
        #SPx (tab)
        #Note: Technically, MFIX will now permit a user to mix-and-match the SPx output files meaning that
        #some can be written and others not. However, this is likely to break the ParaView reader.
        #Therefore, if the “Write binary SPx” checkbox is enabled, output is required for all SPx files.
        #Otherwise, all should remain unspecified to skip writing the SPx files.

        #Write interval for gas volume fraction
        #    Sets keyword SPX_DT(1)
        #    DEFAULT value of 1.0
        #    Required if SPx data is enabled.
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for gas and solids pressure
        #    Sets keyword SPX_DT(2)
        #    Required if SPx data is enabled.
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for gas velocity
        #    Sets keyword SPX_DT(3)
        #    Required if SPx data is enabled.
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for solids velocity
        #    Sets keyword SPX_DT(4)
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for solids bulk density
        #    Sets keyword SPX_DT(5)
        #    Required if SPx data is enabled.
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for gas and solids temperature
        #    Sets keyword SPX_DT(6)
        #    Required if SPx data is enabled and solving any energy equations.
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for gas and solids mass fractions
        #    Sets keyword SPX_DT(7)
        #    Required if SPx data is enabled and solving any species equations.
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for granular temperature
        #    Sets keyword SPX_DT(8)
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for user defined scalars
        #    Sets keyword SPX_DT(9)
        #    Required if SPx data is enabled and solving any user defined scalar equations
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Write interval for reaction rates
        #    Required if SPx data is enabled and NRR > 0 (see below)
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        #Number of reaction rates to write
        #    Specification always available
        #    Sets keyword NRR
        #    DEFAULT value of 0
        #    Error check: value must be greater than or equal to 0

        #Write interval for turbulence quantities
        #    Sets keyword SPX_DT(11)
        #    Required if SPx data is enabled and TURBULENCE_MODEL = “K_EPSILON”
        #    DEFAULT value of 1.0
        #    Error check: value must be greater than or equal to RES_DT

        res_dt = self.project.get_value('res_dt', default=1.0)

        key = 'spx_dt'
        default = 1.0
        for sp in range(1, MAX_SP+1):
            le = getattr(ui, 'lineedit_keyword_spx_dt_args_%s' % sp)
            le.min = res_dt
            val = self.project.get_value(key, args=[sp])
            if val is None:
                val = default
                if val < res_dt:
                    val = res_dt
                self.update_keyword(key, val, args=[sp])

        #Write ASCII particle data
        #    Selection only available if DEM or PIC solids
        #    Sets keyword PRINT_DES_DATA
        #    DEFAULT value of .TRUE.
        key = 'print_des_data'
        solids_models = set(self.project.get_value('solids_model', args=[i])
                            for i in range(1, len(self.solids)+1))
        enabled = 'DEM' in solids_models or 'PIC' in solids_models
        gb = ui.groupbox_print_des_data
        if enabled:
            checked = self.project.get_value(key, default=True)
        else:
            checked = self.project.get_value(key, default=False)
        gb.setEnabled(enabled)
        gb.setChecked(checked)


        #Specify VTP Directory
        #    Specification only if PRINT_DES_DATA = .TRUE.
        #    Sets keyword VTP_DIR
        #    No default (empty string)
        #Select particle data format
        #
        #Selection only available if DEM or PIC solids and PRINT_DES_DATA = .TRUE.
        #  (note, containing groupbox is disabled if this condition is not met)
        #Sets keyword DES_OUTPUT_TYPE
        #Available Selections
        #  ParaView - VTK/.vtp [DEFAULT]
        #  Tecplot - .dat
        cb = ui.combobox_des_output_type
        key = 'des_output_type'
        default = DES_OUTPUT_TYPES[0] # "PARAVIEW"
        val = self.project.get_value(key, default)
        if val not in DES_OUTPUT_TYPES:
            val = default
            self.update_keyword(key, val)
        cb.setCurrentIndex(DES_OUTPUT_TYPES.index(val))


    def setup_output_vtk_tab(self):
        ui = self.ui.output
        self.fixup_output_table(ui.tablewidget_regions)

        #Cell data sub-pane
        #There is a need for some hand waving here. Many mfix.dat files may use a different specification
        #for VTK input. There will need to be a way of catching the 'old' format and converting it to this
        #input style.

        #Specify filename base
        # Specification is required.
        # Sets keyword VTK_FILEBASE(#)
        # DEFAULT value of region name

        #Specify write interval
        # Specification is required
        # Sets keyword VTK_DT(#)
        # DEFAULT value of 1.0 (must write)

        #Specify region x-axis slices
        # Specification always available
        # Sets keyword VTK_NXS(#)
        # DEFAULT value of 0

        #Specify region y-axis slices
        # Specification always available
        # Sets keyword VTK_NYS(#)
        # DEFAULT value of 0

        #Specify region z-axis slices
        # Specification always available
        # Sets keyword VTK_NZS(#)
        # DEFAULT value of 0

        # Fluid Phase (tab?)
        #Enable writing gas volume fraction
        # Selection always available
        # Sets keyword VTK_EP_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas pressure
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_P_G(#)
        # DEFAULT value .FALSE.

        #Enable writing solids pressure
        # Requires TFM solids
        # Sets keyword VTK_P_STAR
        # DEFAULT value .FALSE.

        #Enable writing gas velocity vector
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_VEL_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas velocity x-component
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_U_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas velocity y-component
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_V_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas velocity z-component
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_W_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas temperature
        # Requires fluid solver (RO_G0 /= 0.0) and ENERGY_EQ = .TRUE.
        # Sets keyword VTK_T_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas species N (an entry for each defined species)
        # Requires defined gas phase species
        # Sets keyword VTK_X_G(#,N)
        # DEFAULT value .FALSE.

        #Enable writing gas temperature
        # Requires fluid solver (RO_G0 /= 0.0) and ENERGY_EQ = .TRUE.
        # Sets keyword VTK_T_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas species N
        # Requires fluid solver (RO_G0 /= 0.0) and SPECIES_EQ(0) = .TRUE.
        # Sets keyword VTK_X_G(#,N)
        # DEFAULT value .FALSE.

        #Enable writing turbulent kinetic energy
        #    Requires fluid solver (RO_G0 /= 0.0) and TURBULENCE_MODEL='K_EPSILON'
        #    Sets keyword VTK_K_TURB_G(#)
        #    DEFAULT value of .FALSE.

        #Enable writing turbulent dissipation
        #    Requires fluid solver (RO_G0 /= 0.0) and TURBULENCE_MODEL='K_EPSILON'
        #    Sets keyword VTK_E_TURB_G(#)
        #    DEFAULT value of .FALSE.

        #Enable writing reaction rates
        # Requires nRR > 0
        # Sets keyword VTK_RRATE(#) ## requires 'reaction' index
        # DEFAULT value .FALSE.

        #Solids Phase (tab?)

        #Enable writing solids velocity vector
        # Requires TMF solids
        # Sets keyword VTK_VEL_S(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids velocity x-component
        # Requires TMF solids
        # Sets keyword VTK_U_S(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids velocity y-component
        # Requires TMF solids
        # Sets keyword VTK_V_S(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids velocity z-component
        # Requires TMF solids
        # Sets keyword VTK_W_S(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids bulk density
        # Requires TMF solids
        # Sets keyword VTK_ROP_S(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids temperature
        # Requires TFM solids and ENERGY_EQ = .TRUE.
        # Sets keyword VTK_S_G(#,#)
        # DEFAULT value .FALSE.

        #Enable writing solids phase M, species N
        # Requires TFM solids and SPECIES_EQ(#) = .TRUE.
        # Sets keyword VTK_X_S(#,M,N)
        # DEFAULT value .FALSE.

        #Enable writing solids phase granular temperature
        # Requires TFM solids and KT_TYPE /= “ALGEBRAIC”
        # Sets keyword VTK_THETA_M(#,#)
        # DEFAULT value .FALSE.

        #Scalar (tab?)

        #Enable writing user defined scalar
        # Requires NSCALAR > 0
        # Sets keyword VTK_SCALAR(#, #) # requires Scalar index
        # DEFAULT value .FALSE.

        #Other (tab?)

        #Enable writing vorticity
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_VORTICITY (#)
        # Sets keyword VTK_LAMBDA_2(#)
        # DEFAULT value .FALSE.

        #Enable writing partition
        # Sets keyword VTK_PARTITION(#)
        # DEFAULT value .FALSE.

        #Enable writing boundary ID
        # Sets keyword VTK_BC_ID(#)
        # DEFAULT value .FALSE.

        #Enable writing wall distance
        # Sets keyword VTK_DWALL(#)
        # DEFAULT value .FALSE.

        #Enable writing cell index
        # Sets keyword VTK_IJK(#)
        # DEFAULT value .FALSE.

        #Particle data sub-pane
        #There is a need for some hand waving here. Many mfix.dat files may use a different specification
        #for VTK input. There will need to be a way of catching the 'old' format and converting it to this
        #input style.

        #Specify filename base
        # Specification is required.
        # Sets keyword VTK_FILEBASE(#)
        # DEFAULT value of region name

        #Specify write interval
        # Specification is required
        # Sets keyword VTK_DT(#)
        # DEFAULT value of 1.0 (must write)

        #Specify region x-axis slices
        # Specification always available
        # Sets keyword VTK_NXS(#)
        # DEFAULT value of 0

        #Specify region y-axis slices
        # Specification always available
        # Sets keyword VTK_NYS(#)
        # DEFAULT value of 0

        #Specify region z-axis slices
        # Specification always available
        # Sets keyword VTK_NZS(#)
        # DEFAULT value of 0

        #Enable writing particle diameter
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_DIAMETER(#)
        # DEFAULT value .FALSE.

        #Enable writing particle translational velocity
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_VEL(#)
        # DEFAULT value .FALSE.

        #Enable writing particle rotational velocity
        # Requires DEM or PIC solids
        # Sets keyword VTK_ANGULAR_VEL(#)
        # DEFAULT value .FALSE.

        #Enable writing particle orientation
        # Requires DEM or PIC solids
        # Sets keyword PARTICLE_ORIENTATION = .TRUE.
        # Sets keyword VTK_ORIENTATION(#)
        # DEFAULT value .FALSE.

        #Enable writing particle user variable
        # Requires DEM or PIC solids and DES_USR_VAR > 0
        # Sets keyword VTK_PART_USR_VAR(#,#)
        # DEFAULT value .FALSE.

        #Enable writing particle rotational velocity
        # Requires DEM or PIC solids
        # Sets keyword VTK_ANGULAR_VEL(#)
        # DEFAULT value .FALSE.

        #Enable writing particle temperature
        # Requires DEM or PIC solids and ENERGY_EQ=.TRUE.
        # Sets keyword VTK_PART_TEMP(#)
        # DEFAULT value .FALSE.

        #Enable writing particle species composition
        # Requires DEM or PIC solids and any SPECIES_EQ=.TRUE.
        # Sets keyword VTK_PART_X_S(#)
        # DEFAULT value .FALSE.

        #Enable writing particle MPI rank
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_RANK(#)
        # DEFAULT value .FALSE.

        #Enable writing particle global ID
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_ID(#)
        # DEFAULT value .FALSE.


    def reset_output(self):
        ui = self.ui.output
        # Set all output-related state back to default
        ui.pushbutton_vtk.setEnabled(False)
        self.output_change_tab(TAB_BASIC, ui.pushbutton_basic)
        self.outputs.clear()
        self.output_current_indices = []
        self.output_current_regions = []
        self.output_region_dict = None
        self.output_current_solid = self.P = None
        ui.tablewidget_regions.clearContents()
        ui.tablewidget_regions.setRowCount(0)
