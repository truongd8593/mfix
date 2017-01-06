# -*- coding: utf-8 -*-

"""Output Task Pane Window"""

from __future__ import print_function, absolute_import, unicode_literals, division

import logging
log = logging.getLogger(__name__)

#import Qt
from qtpy import QtWidgets, PYQT5
from qtpy.QtCore import Qt

#local imports
from mfixgui.constants import *

from mfixgui.tools import keyword_args
from mfixgui.tools.general import (widget_iter, get_combobox_item, safe_float)

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

        self.output_current_tab = TAB_BASIC

        ui = self.ui.output
        self.output_pushbuttons = (ui.pushbutton_basic,
                                   ui.pushbutton_vtk,
                                   ui.pushbutton_spx,
                                   ui.pushbutton_netcdf)


        # connect tab buttons
        for i, btn in enumerate(self.output_pushbuttons):
            btn.pressed.connect(lambda i=i, btn=btn: self.output_change_tab(i, btn))

        for item in widget_iter(ui.page_basic):
            if isinstance(item, BaseWidget):
                item.post_update = self.setup_output_basic_tab

        cb = ui.checkbox_binary_spx_files
        cb.clicked.connect(self.output_enable_spx)

        gb = ui.groupbox_print_des_data
        gb.clicked.connect(lambda val:  self.update_keyword('print_des_data', val))
        self.add_tooltip(gb, key = 'print_des_data')

        cb = ui.combobox_des_output_type
        cb.currentIndexChanged.connect(lambda val: self.update_keyword('des_output_type', DES_OUTPUT_TYPES[val]))
        self.add_tooltip(cb, key='des_output_type')

        #cb = ui.checkbox_keyword_write_vtk_files
        #cb.clicked.connect(self.output_enable_vtk)
        #self.add_tooltip(cb, key, value=write_vtk_files)

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


    def output_enable_des(self, val):
        self.update_keyword('print_des_data', val)


    def setup_output(self):
        ui = self.ui.output
        spx_enabled = any(self.project.get_value('spx_dt', args=[i]) is not None
                          for i in range(1,MAX_SP+1))
        ui.pushbutton_spx.setEnabled(spx_enabled)
        self.setup_output_tab(self.output_current_tab)

        vtk_enabled = self.project.get_value('write_vtk_files', default=False)
        ui.pushbutton_vtk.setEnabled(vtk_enabled)


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
        ui.pushbutton_vtk.setEnabled(write_vtk_files)

        #Enable time-dependent VTK files
        #    Specification only if WRITE_VTK_FILES = .TRUE.
        enabled = bool(write_vtk_files)
        #    Sets keyword TIME_DEPENDENT_FILENAME
        key = 'time_dependent_filename'
        #    DEFAULT value .TRUE.
        default = True
        cb = ui.checkbox_keyword_time_dependent_filename
        cb.setEnabled(enabled)
        value = self.project.get_value(key)
        self.add_tooltip(cb, key, value=True if value is None else value)
        if enabled:
            if value is None:
                value = default
                self.update_keyword(key, value)
        else:
            self.unset_keyword(key) #?
        ui.pushbutton_vtk.setEnabled(enabled)

        #Specify VTK Directory
        #    Specification only if WRITE_VTK_FILES = .TRUE.
        #    Sets keyword VTU_DIR
        #    No default (empty string)
        key = 'vtu_dir'
        enabled = bool(write_vtk_files)
        for item in (ui.label_vtu_dir, ui.lineedit_keyword_vtu_dir):
            item.setEnabled(enabled)

        #Write binary Single Precision files (SPx)
        #    No keyword association
        #    Enables SPx tab
        #    Backwards compatibility: Enabled if any SPx time values are specified
        enabled = any(self.project.get_value('spx_dt', args=[i]) is not None
                      for i in range(1,MAX_SP+1)) # Note, enabled in template! *** ??? JORDAN
        ui.checkbox_binary_spx_files.setChecked(enabled)
        ui.pushbutton_spx.setEnabled(enabled)

        #Enable NetCDF output files
        #    Not available when SPx output is enabled
        #    No keyword association.
        #    Enables NetCDF tab



    def output_enable_vtk(self, val):
        print(val)


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
        pass


    def reset_output(self):
        ui = self.ui.output
        # Set all output-related state back to default
        ui.pushbutton_vtk.setEnabled(False)
        self.output_change_tab(TAB_BASIC, ui.pushbutton_basic)


#VTK (tab)
#
#Icons and table similar to IC/BC/PS/IS for adding VTK regions. This section requires WRITE_VTK_FILES = .TRUE.

#Icons to add/remove/duplicate regions are given at the top
#Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select a VTK region.
# Users cannot select inapplicable regions.
# VTK regions can be points, planes, or volumes (not STLs)
# Regions can define multiple VTK regions.

#Select Output type
# Selection is required
# Available selections:
#  Cell data
#    Selection always available
#    Set keyword VTK_DATA(#) to 'C'
#  Particle data
#    Selection only available with DEM or PIC solids
#    Sets keyword VTK_DATA(#) to 'P'

#Cell data sub-pane

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

#Flud Phase (tab?)
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
# Enable writing gas temperature
# Requires fluid solver (RO_G0 /= 0.0) and ENERGY_EQ = .TRUE.
# Sets keyword VTK_T_G(#)
# DEFAULT value .FALSE.
# Enable writing solids velocity vector
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
