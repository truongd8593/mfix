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

# We don't need extended JSON here
from json import JSONDecoder, JSONEncoder

#local imports
from mfixgui.constants import *
from mfixgui.tools import keyword_args
from mfixgui.tools.general import (widget_iter,
                                   get_selected_row,
                                   get_combobox_item,
                                   set_item_noedit,
                                   set_item_enabled,
                                   safe_float)

from mfixgui.widgets.base import (BaseWidget, LineEdit, CheckBox)

TAB_BASIC, TAB_VTK, TAB_SPX, TAB_NETCDF = range(4)
SUBPAGE_CELL, SUBPAGE_PARTICLE = (0, 1)

MAX_SP = 11

VTK_DATA_TYPES = ['C', 'P']

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

        self.outputs = {} # key: index.  value: data dictionary for VTK output region
        self.vtk_current_indices = [] # List of VTK output indices
        self.vtk_current_regions = [] # And the names of the regions which define them
        self.output_region_dict = None
        self.vtk_current_solid = self.P = None

        self.output_pushbuttons = (ui.pushbutton_basic,
                                   ui.pushbutton_vtk,
                                   ui.pushbutton_spx,
                                   ui.pushbutton_netcdf)
        # Dynamically created items
        self.part_usr_var_checkboxes = []
        self.part_x_s_checkboxes = []

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
        ui.checkbox_keyword_vtk_part_orientation_args_V.post_update = self.output_set_particle_orientation

        ui.stackedwidget_vtk_data.setCurrentIndex(SUBPAGE_CELL)



    def handle_output_region_selection(self):
        ui = self.ui.output
        tw = ui.tablewidget_regions
        row = get_selected_row(tw)
        if row is None:
            indices = []
            regions = []
        else:
            (indices, regions) = tw.item(row,0).data(UserRole)
        self.vtk_current_indices, self.vtk_current_regions = indices, regions
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        #ui.detail_pane.setEnabled(enabled)
        ui.bottom_frame.setEnabled(enabled) # more widgets above detail_pane
        if not enabled: # No selection, clear inputs
            for widget in widget_iter(ui.bottom_frame):
                if isinstance(widget, LineEdit):
                    widget.setText('')
                elif isinstance(widget, CheckBox):
                    widget.setChecked(False)
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
            row_height = 0
            height = header_height+scrollbar_height
        else:
            row_height = tw.rowHeight(0)
            height =  (header_height+scrollbar_height
                       + nrows*row_height + 4) # extra to avoid unneeded scrollbar (?)
        ui.top_frame.setMaximumHeight(height+24)
        ui.top_frame.setMinimumHeight(header_height+24+row_height*min(nrows,3))
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
                         and (shape in ('point', 'box')
                              or 'plane' in shape))
            row = (name, shape, available)
            rp.add_row(row)
        #Select Output type
        # Selection is required
        # Available selections:
        #  Cell data
        #    Selection always available
        #    Set keyword VTK_DATA(#) to 'C'
        #  Particle data
        #    Selection only available with DEM or PIC solids
        #    Sets keyword VTK_DATA(#) to 'P'
        solids_models = set(self.project.get_value('solids_model', args=[i])
                            for i in range(1, len(self.solids)+1))
        enabled = 'DEM' in solids_models or 'PIC' in solids_models


        rp.reset_signals()
        rp.save.connect(self.output_add_regions)
        rp.cancel.connect(self.output_cancel_add)
        for item in (ui.tablewidget_regions,
                     ui.detail_pane,
                     ui.toolbutton_add,
                     ui.toolbutton_delete):
            item.setEnabled(False)
        rp.popup('VTK output') # repopulates combobox
        set_item_enabled(get_combobox_item(rp.combobox,1), enabled)


    def output_delete_solids_phase(self, phase):
        # unset keywords associated with solid phase
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
        output_type = VTK_DATA_TYPES[rp.combobox.currentIndex()]
        self.output_add_regions_1(selections, output_type=output_type, indices=None, autoselect=True)
        self.output_setup_current_tab() # Update the widgets


    def output_add_regions_1(self, selections, output_type=None, indices=None, autoselect=False):
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
            self.output_set_region_keys(region_name, idx, region_data, output_type)
            self.output_region_dict[region_name]['available'] = False # Mark as in-use
        item.setData(UserRole, (tuple(indices), tuple(selections)))
        tw.setItem(nrows, 0, item)

        name = 'Cell data' if output_type=='C' else 'Particle data' if output_type=='P' else '???'
        item = make_item(name)
        if output_type == 'C':
            item.setToolTip('Cell data (VTU file)')
        elif output_type == 'P':
            item.setToolTip('Particle data (VTP file)')

        tw.setItem(nrows, 1, item)

        #self.fixup_output_table(tw) # avoid dup. call

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
            if key.startswith('vtk_') and args and args[0] in self.vtk_current_indices:
                self.unset_keyword(key, args=args)

        for r in self.vtk_current_regions:
            if r in self.output_region_dict:
                self.output_region_dict[r]['available'] = True

        for i in self.vtk_current_indices:
            del self.outputs[i]

        self.vtk_current_regions = []
        self.vtk_current_indices = []

        tw.removeRow(row)
        #self.fixup_output_table(tw)
        self.setup_output_vtk_tab()
        #self.update_nav_tree()


    def vtk_regions_to_str(self):
        """convert VTK output region definitions to savable form"""
        ui = self.ui.output
        tw = ui.tablewidget_regions
        data = [tw.item(i,0).data(UserRole)
                for i in range(tw.rowCount())]
        return JSONEncoder().encode(data)


    def vtk_regions_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)
        for (indices, regions) in data:
            if not indices:
                continue # should not get empty tuple
            # vtk_data (output type) keyword should be set already when we call this
            output_type = self.project.get_value('vtk_data', args=[indices[0]], default='C')
            self.output_add_regions_1(regions, output_type=output_type,
                                      indices=indices, autoselect=False)


    def vtk_extract_regions(self):
        if self.outputs:
            # We assume that output regions have been initialized correctly
            # from mfix_gui_comments.
            # TODO: verify that there is an output region for each output
            return

        if self.output_region_dict is None:
            self.output_region_dict = self.ui.regions.get_region_dict()

        # TODO: if we wanted to be fancy, we could find regions where
        # output values matched, and merge into a new output region.  That
        # is only needed for projects created outside the GUI (otherwise
        # we have already stored the output regions).  Also would be noutpute
        # to offer a way to split compound regions.
        for vtk in self.project.vtks:

            d = vtk.keyword_dict
            extent = [d.get('vtk_'+k,None) for k in ('x_w', 'y_s', 'z_b',
                                                    'x_e', 'y_n', 'z_t')]
            extent = [0 if x is None else x.value for x in extent]
            #if any (x is None for x in extent):
            #    self.warn("vtk output %s: invalid extents %s" %
            #               (vtk.ind, extent))
            #    continue
            for (region_name, data) in self.output_region_dict.items():
                ext2 = [0 if x is None else x for x in
                        (data.get('from',[]) + data.get('to',[]))]
                if ext2 == extent:
                    if data.get('available', True):
                        self.output_add_regions_1([region_name], indices=[vtk.ind], autoselect=False)
                        break

            else:
                self.warn("vtk output %s: could not match defined region %s" %
                          (vtk.ind, extent))


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

        # Mark regions which are in use (this gets reset each time we get here)
        for (i, data) in self.outputs.items():
            region = data['region']
            if region in self.output_region_dict:
                self.output_region_dict[region]['available'] = False


        spx_enabled = any(self.project.get_value('spx_dt', args=[i]) is not None
                          for i in range(1,MAX_SP+1))
        ui.pushbutton_spx.setEnabled(spx_enabled)

        vtk_enabled = self.project.get_value('write_vtk_files', default=False)
        ui.pushbutton_vtk.setEnabled(vtk_enabled)

        netcdf_enabled = False
        ui.pushbutton_netcdf.setEnabled(netcdf_enabled)

        # TODO don't stay on disabled tab
        self.output_setup_current_tab()


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


    def output_setup_current_tab(self):
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
                      for i in range(1,MAX_SP+1)) # Note, enabled in template! XXX Jordan?
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


    def output_default_vtk_filebase(self, region_names):
        # Construct default value for VTK_FILEBASE,
        # replacing possibly troublesome characters
        val = '+'.join(region_names)
        for c in ': /*?':
            val = val.replace(c, '_')
        return val


    def output_set_particle_orientation(self):
        o = any(self.project.get_value('vtk_part_orientation', args=[V])
                for V in self.outputs)
        self.update_keyword('particle_orientation', o)


    def setup_output_vtk_tab(self):
        ui = self.ui.output
        self.fixup_output_table(ui.tablewidget_regions)

        row = get_selected_row(ui.tablewidget_regions)
        # Autoselect if only 1 row
        if row is None and ui.tablewidget_regions.rowCount() == 1:
            row = 0
            ui.tablewidget_regions.setCurrentCell(row, 0)
        enabled = (row is not None)
        ui.toolbutton_delete.setEnabled(enabled)
        ui.bottom_frame.setEnabled(enabled)

        indices = self.vtk_current_indices
        if not indices:
            # Clear inputs?  should have been done in handle_selection
            return

        V0 = indices[0]
        vtk_data = self.project.get_value('vtk_data', args=[V0])

        enabled = (vtk_data=='C')
        ui.stackedwidget_vtk_data.setCurrentIndex(SUBPAGE_CELL if vtk_data=='C' else SUBPAGE_PARTICLE)

        #Cell data sub-pane
        #There is a need for some hand waving here. Many mfix.dat files may use a different specification
        #for VTK input. There will need to be a way of catching the 'old' format and converting it to this
        #input style.

        # Note, filebase through nzs are common to cell/particle

        #Specify filename base
        # Specification is required.
        # Sets keyword VTK_FILEBASE(#)
        # DEFAULT value of region name
        key = 'vtk_filebase'
        le = ui.lineedit_keyword_vtk_filebase_args_V
        val = self.project.get_value(key, args=[V0])
        if val is None: # Construct from region name
            val = self.output_default_vtk_filebase(self.vtk_current_regions)
            for V in self.vtk_current_indices:
                self.update_keyword(key, val, args=[V])
        le.updateValue(key, val)

        #Specify write interval
        # Specification is required
        # Sets keyword VTK_DT(#)
        # DEFAULT value of 1.0 (must write)
        key = 'vtk_dt'
        default = 1.0
        le = ui.lineedit_keyword_vtk_dt_args_V
        val = self.project.get_value(key, args=[V0])
        if val is None:
            val = default
            for V in self.vtk_current_indices:
                self.update_keyword(key, val, args=[V])
        le.updateValue(key, val)

        for c in 'xyz':
            #Specify region x[yz]-axis slices
            # Specification always available
            # Sets keyword VTK_NXS(#)
            # DEFAULT value of 0
            key = 'vtk_n%ss' % c
            default = 0
            le = getattr(ui, 'lineedit_keyword_vtk_n%ss_args_V'%c)
            val = self.project.get_value(key, args=[V0])
            if val is None:
                val = default
                for V in self.vtk_current_indices:
                    self.update_keyword(key, val, args=[V])
            le.updateValue(key, val)

        # Fluid Phase (tab?)
        #Enable writing gas volume fraction
        # Selection always available
        # Sets keyword VTK_EP_G(#)
        # DEFAULT value .FALSE.

        #Enable writing gas pressure
        # Requires fluid solver (RO_G0 /= 0.0)
        # Sets keyword VTK_P_G(#)
        # DEFAULT value .FALSE.

        #Enable writing solids pressure  # MOVE TO SOLIDS TAB
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


        #move to separate Reactions tab?
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

        # Note 1, filebase through nzs are common to cell/particle
        # Note 2, this groupbox is disabled completely if not DEM or PIC

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
        # Sets keyword VTK_ANGULAR_VEL(#)  # VTK_PART_ANGULAR_VEL
        # DEFAULT value .FALSE.

        #Enable writing particle orientation
        # Requires DEM or PIC solids
        # Sets keyword PARTICLE_ORIENTATION = .TRUE.
        # Sets keyword VTK_ORIENTATION(#)    ## VTK_PART_ORIENTATION
        # DEFAULT value .FALSE.

        #Enable writing particle temperature
        # Requires DEM or PIC solids and ENERGY_EQ=.TRUE.
        # Sets keyword VTK_PART_TEMP(#)
        # DEFAULT value .FALSE.

        #Enable writing particle MPI rank
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_RANK(#)
        # DEFAULT value .FALSE.

        #Enable writing particle global ID
        # Requires DEM or PIC solids
        # Sets keyword VTK_PART_ID(#)
        # DEFAULT value .FALSE.
        for key in ('vtk_part_diameter',
                    'vtk_part_vel',
                    'vtk_part_angular_vel',
                    'vtk_part_orientation',
                    'vtk_part_temp',
                    'vtk_part_rank',
                    'vtk_part_id'):
            cb = getattr(ui, 'checkbox_keyword_%s_args_V' % key)
            val = self.project.get_value(key, args=[V0], default=False)
            cb.setChecked(val)


        # disable temp if energy_eq disabled
        cb =  ui.checkbox_keyword_vtk_part_temp_args_V
        key = 'vtk_part_temp'
        enabled = self.project.get_value('energy_eq', default=True)
        cb.setEnabled(enabled)
        if not enabled:
            for V in self.vtk_current_indices:
                self.unset_keyword(key, args=[V])

        # Dynamically created GUI items - need to remove and re-add spacer
        layout = ui.groupbox_particle_data.layout()
        spacer = None
        spacer_moved = False
        # find spacer, can't do this by name for some reason
        for i in range(layout.count()-1, -1, -1):
            item = layout.itemAt(i)
            if not item:
                continue
            widget = item.widget()
            if not widget:
                spacer = item
                break

        #enable writing particle user variable
        # Requires DEM or PIC solids and DES_USR_VAR > 0 ## DES_USR_VAR_SIZE
        # Sets keyword VTK_PART_USR_VAR(#,#)
        # DEFAULT value .FALSE.

        des_usr_var_size = self.project.get_value('des_usr_var_size', default=0)
        # Remove extra widgets if number decreased
        if len(self.part_usr_var_checkboxes) > des_usr_var_size:
            for cb in self.part_usr_var_checkboxes[des_usr_var_size:]:
                layout.removeWidget(cb)
                cb.deleteLater()
            self.part_usr_var_checkboxes = self.part_usr_var_checkboxes[:des_usr_var_size]
        # If adding widgets, remove spacer first (will re-add it at end)
        if len(self.part_usr_var_checkboxes) != des_usr_var_size:
            layout.removeItem(spacer)
            spacer_moved = True
        while len(self.part_usr_var_checkboxes) < des_usr_var_size:
            n = 1+len(self.part_usr_var_checkboxes)
            cb = CheckBox("DES user scalar %s" % n)
            cb.key = 'vtk_part_usr_var'
            cb.args = ['V', n]
            cb.value_updated.connect(self.project.submit_change)
            self.part_usr_var_checkboxes.append(cb)
            self.add_tooltip(cb, key=cb.key)
            layout.addWidget(cb)

        #Enable writing particle species composition
        # Requires DEM or PIC solids and any SPECIES_EQ=.TRUE.
        # Sets keyword VTK_PART_X_S(#)
        # NOTE, VTK_PART_X_S(#,N) where N ranges from 1 to max(mmax_s)  # nmax_s
        # DEFAULT value .FALSE.
        max_n = max(self.project.get_value('nmax_s', args=[N], default=0)
                    for N in range(1, len(self.solids)+1))
        # Remove extra widgets if number decreased
        if len(self.part_x_s_checkboxes) > max_n:
            for cb in self.part_x_s_checkboxes[max_n:]:
                layout.removeWidget(cb)
                cb.deleteLater()
            self.part_x_s_checkboxes = self.part_x_s_checkboxes[:max_n]
        # If adding widgets, remove spacer first (will re-add it at end)
        if len(self.part_x_s_checkboxes) != max_n:
            if not spacer_moved:
                layout.removeItem(spacer)
                spacer_moved = True
        while len(self.part_x_s_checkboxes) < max_n:
            n = 1+len(self.part_x_s_checkboxes)
            cb = CheckBox("Species %s mass fraction" % n)
            cb.key = 'vtk_part_x_s'
            cb.args = ['V', n]
            cb.value_updated.connect(self.project.submit_change)
            self.part_x_s_checkboxes.append(cb)
            self.add_tooltip(cb, key=cb.key)
            layout.addWidget(cb)

        if spacer_moved:
            layout.addItem(spacer)



    def output_set_region_keys(self, name, idx, data, output_type=None):
        # Update the keys which define the region the PS applies to
        if output_type is not None:
            self.update_keyword('vtk_data', output_type, args=[idx])

        no_k = self.project.get_value('no_k')
        for (key, val) in zip(('x_w', 'y_s', 'z_b',
                               'x_e', 'y_n', 'z_t'),
                              data['from']+data['to']):
            # vtk_z_t and vtk_z_b keywords should not be added when no_k=True
            if no_k and key in ('z_t', 'z_b'):
                continue
            self.update_keyword('vtk_'+key, val, args=[idx])


    def output_check_region_in_use(self, name):
        return any(data.get('region')==name for data in self.outputs.values())


    def output_update_region(self, name, data):
        for (i, output) in self.outputs.items():
            if output.get('region') == name:
                self.output_set_region_keys(name, i, data)


    def output_change_region_name(self, old_name, new_name):
        ui = self.ui.output
        for (key, val) in self.outputs.items():
            if val.get('region') == old_name:
                self.outputs[key]['region'] = new_name
                tw = ui.tablewidget_regions
                for i in range(tw.rowCount()):
                    data = tw.item(i,0).data(UserRole)
                    indices, names = data
                    if key in indices:
                        item = tw.item(i,0)
                        new_names = [new_name if n==old_name else n for n in names]
                        item.setData(UserRole, (indices, new_names))
                        item.setText('+'.join(new_names))
                        # Also update vtk_filebase, if it is at the default setting
                        vtk_filebase = self.project.get_value('vtk_filebase', args=[key])
                        if vtk_filebase == self.output_default_vtk_filebase(names):
                            vtk_filebase = self.output_default_vtk_filebase(new_names)
                            for i in indices:
                                self.update_keyword('vtk_filebase', vtk_filebase, args=[i])
                        break
                break


    def reset_output(self):
        ui = self.ui.output
        # Set all output-related state back to default
        ui.pushbutton_vtk.setEnabled(False)
        self.output_change_tab(TAB_BASIC, ui.pushbutton_basic)
        self.outputs.clear()
        self.vtk_current_indices = []
        self.vtk_current_regions = []
        self.output_region_dict = None
        self.vtk_current_solid = self.P = None
        ui.tablewidget_regions.clearContents()
        ui.tablewidget_regions.setRowCount(0)

        ui.stackedwidget_vtk_data.setCurrentIndex(SUBPAGE_CELL)
