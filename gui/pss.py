# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox)

from qtpy.QtGui import QPixmap # QPicture doesn't work with Qt4

UserRole = QtCore.Qt.UserRole

from constants import *
from widgets.regions_popup import RegionsPopup
from widgets.base import LineEdit, ComboBox

from project import Equation, FloatExp

from tools.general import (set_item_noedit, set_item_enabled,
                           get_combobox_item, get_selected_row,
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
#SCALAR_TAB = 4

class PSS(object):

    #Point Source Task Pane Window: This section allows a user to define point sources for the
    #described model. This section relies on regions named in the Regions section.

    def init_pss(self):
        ui = self.ui.point_sources

    #The top of the task pane is where users define/select PS regions



#Icons to add/remove/duplicate boundary conditions are given at the top

#Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select
#a region to apply the point source.
# Users cannot select inapplicable regions.
# PS regions can be points, planes, or volumes (not STLs)
# No region can define more than one point source.

#Tabs group point source parameters for phases. Tabs are unavailable if no input is required from the user.
#    Fluid tab - Unavailable if the fluid phase was disabled.
#    Each solid phase will have its own tab. The tab name should be the name of the solid

#Fluid (tab)
#    Define mass flowrate:
# Sets keyword PS_MASSFLOW_G(#)
# DEFAULT value of 0.0

#    Define temperature
# Specification only available when solving energy equations
# DEFAULT value 293.15
# Sets keyword PS_T_G(#)

#    Select species and set mass fractions (table format)
# Specification only available when solving species equations
# DEFAULT value 1.0 of last defined species
# Sets keyword PS_X_G(#,#)
# Error check: if specified, mass fractions must sum to one

#    Define X-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_U_G(#)
#    Define Y-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_V_G(#)
#    Define Z-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_W_G(#)

#Solids-# (tab) - At this time, only TFM solids can be defined with point sources.
#At some point in the future, this could be extended to PIC solids, but never DEM.

#Define mass flowrate:
# Select mass inflow specification type:
# Sets keyword PS_MASSFLOW_S(#,#)
# DEFAULT value of 0.0

#Define temperature
# Specification only available when solving energy equations
# DEFAULT value 293.15
# Sets keyword PS_T_S(#,#)

#Select species and set mass fractions (table format)
# Specification only available when solving species equations
# DEFAULT value 1.0 of last defined species
# Sets keyword PS_X_S(#,#,#)
# Error check: if specified, mass fractions must sum to one

#Define X-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_U_S(#,#)

#Define Y-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_V_S(#,#)

#Define Z-axial velocity:
# Specification always
# No DEFAULT value (blank)
# Sets keyword PS_W_S(#,#)

# No scalar tab for point sources
