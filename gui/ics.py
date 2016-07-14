# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from widgets.regions_popup import RegionsPopup

"""
Initial Conditions Task Pane Window: This section allows a user to define the initial conditions
for the described model. This section relies on regions named in the Regions section.
The top of the task pane is where users define/select IC regions
"""

"""
Icons to add/remove/duplicate regions are given at the top
Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
the region to apply the initial condition.
 Users cannot select inapplicable regions.
 IC regions must be volumes (not planes, points, or STLs)
 No region can define more than one initial condition.
"""

"""
Implementation Idea: Allow the user to select multiple regions in the popup window so that they can
define one IC region over multiple regions. Managing the MFIX IC indices on the back end could
get messy, and if so, scrap this idea.

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
        self.ics_current_region = None

    def ics_add_region(self):
        rp = self
        ui = self.ui

        for (k,v) in ui.regions.get_region_dict().items():
            print(k)

    def setup_ics(self):
        ui = self.ui
        ics = self.ui.initial_conditions
        ics.toolbutton_add.clicked.connect(self.ics_add_region)
