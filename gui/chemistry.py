# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division
from collections import OrderedDict

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtWidgets import (QLabel, QLineEdit, QPushButton, QGridLayout,
                            QHBoxLayout, QWidget, QGroupBox)

#Chemistry Task Pane Window: This section allows a user to define chemical reaction input.
#Chemistry pane is disabled if any solids are specified as PIC.

#Enable the stiff chemistry solver
# Selection always available
# Sets keyword STIFF_CHEMISTRY to .TRUE.

#Chemical reaction input is handled different that keyword pair inputs. All homogeneous gas phase
#chemical reactions and all heterogeneous gas-tfm solids reactions are specified between @(RXNS)
#and @(END) reaction block. All heterogeneous gas-dem and gas-pic reactions are specified
#between @(DES_RXNS) and @(END) reaction block.

#Users use the globally unique species aliases to specify the chemical reaction equation. Each
#reaction is specified with a unique reaction identify.

#Specify the reaction identifier (Name)
# Specification always available
# DEFAULT value reactionN (for the Nth reaction)
# Reaction identifiers must be "Fortran compilable"
#  Alphanumeric combinations (no special characters excluding underscores)
#  Limited to 32 characters
#  First character must be a letter
#  No black spaces

#Specify chemical reaction reactants (table format)
# Use +/- buttons to add/remove reactants
# Column 1 - Select phase for reactant
#  Use drop down list of user defined phase names
#  Reactions are limited to homogeneous or two-phase heterogeneous reactions
#(e.g., species from three separate phases cannot be referenced by any single chemical reaction)
# Column 2 - Select reactant species
#  Use drop down list to show species in selected phase
#  A species can only appear once as a reactant in the same reaction
# Column 3 - Enter stoichiometric coefficient
#  Numerical value (integer or float)
#  Value must be non-negative
#Specify chemical reaction products (table format)
# Use +/- buttons to add/remove products
# Column 1 - Select phase for product
#  Use drop down list of user defined phase names

#Reactions are limited to homogeneous or two-phase heterogeneous reactions
#(e.g., species from three separate phases cannot be referenced by any single
#chemical reaction)

# Column 2 - Select product species
#  Use drop down list to show species in selected phase
#  A species can only appear once as a product in the same reaction

# Column 3 - Enter stoichiometric coefficient
#  Numerical value (integer or float)
#  Value must be non-negative

#Reactant/Product information is combined with stoichiometric coefficients to define the
#chemical reaction as a string.
# Sets reaction construct keyword CHEM_EQ
# Example: CHEM_EQ = "rcoeff1*reactant1 + rcoeff2*reactant2 --> pcoeff1*product1"

# Error check: Mass of the reactants equal mass of the products (within a tolerance, 1.0e-6).
# abs[(rcoeff1*MW(reactant1) + rcoeff2*MW(reactant2) - (pcoeff1*product1)] < 1.0e-6

#Enable user-defined heat of reaction
# Selection always available
# DEFAULT disabled

#Specify heat of reaction
# Only available if user-defined heat of reaction is enabled
# DEFAULT value 0.0
# Sets reaction construct keyword DH

#Specify HoR fraction assigned to -phase
# Only available if user-defined heat of reaction is enabled
# Homogeneous chemical reactions
#  Specification is not available
#  Set reaction construct keyword fracDH(#) to 1.0 where # is the phase index
# Heterogeneous chemical reactions
#  Entry for each phase referenced by the reaction
#  DEFAULT value 0.5 for both entires
#  Sets reaction construct keyword fracDH(#) for each referenced phase

# The user cannot 'save' the reaction if there are errors. After
# saving (adding?) the reaction, the reaction identifier (name) and
# chemical equation are shown in the summary box at the top. A
# chemical reaction is activated/deactivated by checking/unchecking the box. If the user 'deactivates'
# the chemical equation, the CHEM_EQ reaction construct keyword should get set to "NONE."
