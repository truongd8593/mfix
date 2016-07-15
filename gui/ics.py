# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5

from widgets.regions_popup import RegionsPopup
from tools.general import (set_item_noedit, set_item_enabled, get_selected_row)

"""
Initial Conditions Task Pane Window: This section allows a user to define the initial conditions
for the described model. This section relies on regions named in the Regions section.
The top of the task pane is where users define/select IC regions
"""

"""
Tabs group initial condition parameters for phases and additional equations. Tabs are unavailable if
no input is required from the user.
 Fluid tab - Unavailable if the fluid phase was disabled.

UI

Each solid phase will have its own tab. The tab name should be the name of the solid
Group tab inputs by equation type (e.g., momentum, energy, species). Making the grouped
inputs a 'collapsible list' may make navigation easier.

Fluid (tab)
 Define volume fraction
 Specification always available
 Sets keyword IC_EP_G(#)
 DEFAULT value of 1.0

Define temperature
 Specification always available
 Input required for any of the following
 Fluid density model: Ideal Gas Law
 Fluid viscosity model: Sutherland's Law
 Energy equations are solved
 Sets keyword IC_T_G(#)
 DEFAULT value of 293.15

Define pressure (optional)
 Specification always available
 Sets keyword IC_P_g(#)
 DEFAULT - no inputDefine velocity components (required)
 Specification always available
 Sets keywords IC_U_G(#), IC_V_G(#), IC_W_G(#)
 DEFAULT value of 0.0

Select species and set mass fractions (table format)
 Specification always available
 Input required for species equations
 Drop down menu of fluid species
 DEFAULT - last defined species has mass fraction of 1.0
 Error check: mass fractions must sum to one

Turbulence: Define mixing length model length scale
 Specification only available with Mixing Length turbulence model
 Sets keyword IC_L_SCALE(#)
 DEFAULT value of 1.0

Turbulence: Define k-ε turbulent kinetic energy
 Specification only available with K-Epsilon turbulence model
 Sets keyword IC_K_TURB_G(#)
 DEFAULT value of 0.0

Turbulence: Define k-ε turbulent dissipation
 Specification only available with K-Epsilon turbulence model
 Sets keywords IC_E_TURB_G(#)
 DEFAULT value of 0.0

Advanced: Define radiation coefficient
 Specification only available when solving energy equations
 Sets keyword IC_GAMA_G(#)
 DEFAULT value of 0.0
Advanced: Define radiation temperature
 Specification only available when solving energy equations
 Sets keyword IC_T_RG(#)
 DEFAULT value of 293.15


Solid-# (tab) - Rename tab to user provided solids name.
 Define volume fraction (required)
 Specification always available
 Sets keyword IC_EP_S(#,#)
 DEFAULT value of (0.0, 1 - SUM)
 Define temperature
 Specification always available
 Input required when solving energy equations
 Sets keyword IC_T_S(#,#)
 DEFAULT value of 293.15
 Define velocity components (required)
 Specification always available
 Sets keywords IC_U_S(#,#), IC_V_S(#,#), IC_W_S(#,#)
 DEFAULT value of 0.0
 Define pressure (optional)
 Specification only available for SOLIDS_MODEL(#)='TFM'
 Sets keyword IC_P_s(#,#)
 DEFAULT of 0.0
 Define granular temperature
 Specification only available for SOLIDS_MODEL(#)='TFM' and non-algebraic
 formulation viscous stress model (see continuous solids model section) or for
 SOLIDS_MODEL(#)=DEM' or SOLIDS_MODEL(#)='PIC'
 Sets keyword IC_THETA_M(#,#)
 DEFAULT value of 0.0
 Define particles per parcel
 Specification only available for SOLIDS_MODEL(#)='PIC'
 Sets keyword IC_PIC_CONST_STATWT(#,#)
 DEFAULT value of 10.0
 Select species and set mass fractions (table format)
 Specification always available
 Input required for species equations
 Drop down menu of solids species
 DEFAULT - last defined species has mass fraction of 1.0
 Error check: mass fractions must sum to one
Advanced: Option to enable fitting DES particles to region
 Option only available for DEM solids
 Sets keyword: IC_DES_FIT_TO_REGION
 Disabled [DEFAULT]
Advanced: Define radiation coefficient
 Specification only available when solving energy equations
 Sets keyword IC_GAMA_S(#,#)
 DEFAULT value of 0.0
Advanced: Define radiation temperature
 Specification only available when solving energy equations
 Sets keyword IC_T_RS(#,#)
 DEFAULT value of 293.15


Scalar (tab) - Tab only available if scalar equations are solved
 Define initial scalar value
 Sets keyword IC_SCALAR(#,#)
 DEFAULT value of 0.0
 Note that this tab should only be available if scalar equations are being solved.
 Furthermore, the number of scalars requiring input comes from the number of
 scalar equations specified by the user.

"""

class ICS(object):
    def init_ics(self):
        ui = self.ui
        ics = ui.initial_conditions
        self.ics_current_region = None

        #Icons to add/remove/duplicate regions are given at the top
        #Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
        #the region to apply the initial condition.
        #Users cannot select inapplicable regions.
        #IC regions must be volumes (not planes, points, or STLs)
        #No region can define more than one initial condition.
        #
        #Implementation Idea: Allow the user to select multiple regions in the popup window so that they can
        #define one IC region over multiple regions. Managing the MFIX IC indices on the back end could
        #get messy, and if so, scrap this idea.
        ics.toolbutton_add.clicked.connect(self.ics_add_region)
        ics.toolbutton_delete.clicked.connect(self.ics_delete_region)

        ics.tablewidget_regions.itemSelectionChanged.connect(self.handle_ics_region_selection)

    def handle_ics_region_selection(self):
        ics = self.ui.initial_conditions
        table = ics.tablewidget_regions
        row = get_selected_row(table)
        enabled = (row is not None)
        ics.toolbutton_delete.setEnabled(enabled)
        ics.scrollarea.setEnabled(enabled)

    def ics_add_region(self):
        ui = self.ui
        ics = ui.initial_conditions
        rp = self.regions_popup

        rp.clear()
        for (name,data) in self.region_dict.items():
            shape = data.get('type', '---')
            available = (shape=='box') # and not in use
            row = (name, shape, available)
            rp.add_row(row)
        rp.reset_signals()
        rp.save.connect(self.ics_update_regions)
        rp.popup("Select region(s) for initial condition")

    def ics_update_regions(self):
        ui = self.ui
        ics = ui.initial_conditions
        rp = self.regions_popup
        sel = rp.get_selection()
        if not sel:
            return
        tw = ics.tablewidget_regions
        nrows = tw.rowCount()
        tw.setRowCount(nrows+1)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        tw.setItem(nrows, 0, make_item('+'.join(sel)))
        self.fixup_table(tw)

    def fixup_table(self, tw, stretch_column=0):
        # fixme, this is getting called excessively
        # Should we just hide the entire table (including header) if no rows?
        hv = QtWidgets.QHeaderView
        if PYQT5:
            resize = tw.horizontalHeader().setSectionResizeMode
        else:
            resize = tw.horizontalHeader().setResizeMode
        ncols = tw.columnCount()
        resize(0, hv.Stretch)
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
        tw.setMinimumHeight(height) #? needed for tablewidget_des_en_input. should we allow scrollbar?
        tw.updateGeometry() #? needed?

    def ics_delete_region(self):
        pass

    def setup_ics(self):
        ui = self.ui
        self.region_dict = ui.regions.get_region_dict()
        self.fixup_table(ui.initial_conditions.tablewidget_regions)
