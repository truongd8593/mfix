from __future__ import print_function, absolute_import, unicode_literals, division

import logging
log = logging.getLogger(__name__)

#import Qt
from qtpy import QtWidgets, PYQT5
from qtpy.QtCore import Qt

#local imports
from mfixgui.constants import *
from mfixgui.tools.general import (set_item_noedit, get_selected_row,
                                   widget_iter, format_key_with_args,
                                   get_combobox_item, set_item_enabled)


from mfixgui.tools import keyword_args


(TAB_RESIDUALS, TAB_DISCRETIZATION, TAB_LINEAR_SOLVER,
 TAB_PRECONDITIONER, TAB_ADVANCED) = range(5)

class Numerics(object):
    """Numerics Task Pane Window
    The numerics input is split into tabs to group similar inputs and reduce the amount of input needed
    on a single pane."""

    def init_numerics(self):
        ui = self.ui.numerics
        self.numerics_pushbuttons = (ui.pushbutton_residuals,
                                     ui.pushbutton_discretization,
                                     ui.pushbutton_linear_solver,
                                     ui.pushbutton_preconditioner,
                                     ui.pushbutton_advanced)
        for (i, btn) in enumerate(self.numerics_pushbuttons):
            btn.pressed.connect(lambda i=i, btn=btn: self.numerics_change_tab(i, btn))

        self.numerics_current_tab = TAB_RESIDUALS


    def numerics_update_tabs(self):
        ui = self.ui.numerics
        # Don't stay on a disabled tab!
        active_tab = ui.stackedwidget_numerics.currentIndex()


    # Numerics sub-pane navigation
    def numerics_change_tab(self, tabnum, to_btn):
        ui = self.ui.numerics
        self.numerics_current_tab = tabnum
        self.animate_stacked_widget(
            ui.stackedwidget_numerics,
            ui.stackedwidget_numerics.currentIndex(),
            tabnum,
            direction='horizontal',
            line=ui.line_numerics,
            to_btn=to_btn,
            btn_layout=ui.gridlayout_tab_btns)
        self.setup_numerics_tab(tabnum)
        for btn in self.numerics_pushbuttons:
            btn.setChecked(btn==to_btn)
            font = btn.font()
            font.setBold(btn==to_btn)
            btn.setFont(font)


    def setup_numerics(self):
        self.setup_numerics_tab(self.numerics_current_tab)


    def setup_numerics_tab(self, tabnum):

        if tabnum == TAB_RESIDUALS:
            self.numerics_setup_residuals_tab()
        elif tabnum == TAB_DISCRETIZATION:
            self.numerics_setup_discretization_tab()
        elif tabnum == TAB_LINEAR_SOLVER:
            self.numerics_setup_linear_solver_tab()
        elif tabnum == TAB_PRECONDITIONER:
            self.numerics_setup_preconditioner_tab()
        elif tabnum == TAB_ADVANCED:
            self.numerics_setup_advanced_tab()
        else:
            raise ValueError(tabnum)


    def reset_numerics(self):
        ui = self.ui.numerics
        # Set all numeric-related state back to default


    def numerics_setup_residuals_tab(self):
        """Residuals (tab)"""
        pass

        #Specify residual for continuity plus momentum equations
        #    Specification always available
        #    Sets keyword TOL_RESID
        #    DEFAULT value of 1.0e-3
        #Specify residual for energy equations
        #    Specification always available
        #    Sets keyword TOL_RESID_T
        #    DEFAULT value of 1.0e-4
        #Specify residual for species equations
        #    Specification always available
        #    Sets keyword TOL_RESID_X
        #    DEFAULT value of 1.0e-4
        #Specify residual for granular energy equations
        #    Specification always available
        #    Sets keyword TOL_RESID_TH
        #    DEFAULT value of 1.0e-4
        #Specify residual for scalar/K-Epsilon
        #    Specification always available
        #    Sets keyword TOL_RESID_SCALAR
        #    DEFAULT value of 1.0e-4

    def numerics_setup_discretization_tab(self):
        """Discretization (tab)"""

        #Select temporal discretization scheme
        #    Selection always available
        #    Available selections:
        #  Implicit Euler [DEFAULT]
        #    Sets keyword CN_ON to .FALSE.
        #  Crank-Nicholson
        #    Sets keyword CN_ON to .TRUE.
        #Specify spatial discretization and under relation factors (table format)
        #    Specification always available
        #    Column 1: List of equations
        #    Column 2: Select discretization scheme for equation #
        #  Available selections
        #    First-order upwind [DEFAULT for all equations]

        # Sets keyword DISCRETIZE(#) to 0
        #First-order upwind (dwf)
        # Sets keyword DISCRETIZE(#) to 1
        #    Superbee
        # Sets keyword DISCRETIZE(#) to 2
        #    SMART
        # Sets keyword DISCRETIZE(#) to 3
        #    ULTRA-QUICK
        # Sets keyword DISCRETIZE(#) to 4
        #    QUICKEST
        # Sets keyword DISCRETIZE(#) to 5
        #    MUSCL
        # Sets keyword DISCRETIZE(#) to 6
        #    van Leer
        # Sets keyword DISCRETIZE(#) to 7
        #    minmod
        # Sets keyword DISCRETIZE(#) to 8
        #    Central
        # Sets keyword DISCRETIZE(#) to 9
        #    Column 3: Specify under relation factors
        #  Specification always available
        #  Sets keyword UR_FAC for each equation #
        #  DEFAULTS are equation type specific
        #    1 - gas pressure: 0.8
        #    2 - volume fraction: 0.5
        #    3 - u-momentum: 0.5
        #    4 - v-momentum: 0.5
        #    5 - w-momentum: 0.5
        #    6 - energy: 1.0
        #    7 - species: 1.0
        #    8 - granular energy: 0.5
        #    9 - user-scalar/k-epsilon: 0.8
        #    10 - DES diffusion: 1.0
        #Enable deferred correction
        #    Selection only available if minval(discretize) > 0
        #    Sets keyword DEF_COR
        #    DEFAULT value .FALSE.
        #Enable chi-scheme correction
        #    Selection only available if the species equation spatial discretization is SMART or
        #MUSCL (DISCRETIZE(7) = 3 or 6)
        #    Sets keyword CHI_SCHEME
        #    DEFAULT value .FALSE.


    def numerics_setup_linear_solver_tab(self):
        """Linear Solver (tab)"""

        #Specify linear solver, number of iterations, and convergence tolerance (table format)
        #    Specification always available
        #    Column 1: List of equations
        #    Column 2: Select linear equation solver method for equation #
        #  Available selections
        #    BiCGSTAB [DEFAULT for all equations]
        # Sets keyword LEQ_METHOD(#) to 2
        #    GMRES
        #  Sets keyword LEQ_METHOD(#) to 3
        #    Column 3: Specify number of iterations
        #  Specification always available
        #  Sets keyword LEQ_IT for each equation #
        #  DEFAULTS are equation type specific
        #    1 - gas pressure: 20
        #    2 - volume fraction: 20
        #    3 - u-momentum: 5
        #    4 - v-momentum: 5
        #    5 - w-momentum: 5
        #    6 - energy: 15
        #    7 - species: 15
        #    8 - granular energy: 15
        #    9 - user-scalar/k-epsilon: 15
        #    10 - DES diffusion: 5
        #    Column 4: Specify convergence tolerance
        #  Specification always available
        #  Sets keyword LEQ_TOL
        #  DEFAULT value of 1.0E-4 for all equations
        pass

    def numerics_setup_preconditioner_tab(self):
        """Preconditioner (tab)"""

        #Specify linear solver, number of preconditioner and sweep direction (table format)
        #    Specification only available for equations using BiCGSTAB solver
        #    Column 1: List of equations
        #    Column 2: Preconditioner for equation #
        #  Available selections
        #    None
        # Sets keyword LEQ_PC(#) to 'NONE'
        #    Line Relaxation [DEFAULT for all equations]
        # Sets keyword LEQ_PC(#) to 'LINE'
        #    Diagonal Scaling
        # Sets keyword LEQ_PC(#) to 'DIAG'
        #    Column 3: Preconditioner sweep direction for equation #
        #  Selection only available for equations with LINE preconditioner
        #  Available selections
        #    'Red-black sweep' [DEFAULT for all equations]
        # Sets keyword LEQ_PC(#) to 'RSRS'
        #    All sweep
        # Sets keyword LEQ_PC(#) to 'ASAS'
        #    I-sweep
        # Sets keyword LEQ_PC(#) to 'ISIS'
        #    J-sweep
        # Sets keyword LEQ_PC(#) to 'JSJS'
        #    K-sweep
        # Sets keyword LEQ_PC(#) to 'KSKS'
        pass

    def numerics_setup_advanced_tab(self):
        """Advanced (tab)"""
        #Specify maximum inlet velocity factor
        #    Specification always available
        #    Sets keyword MAX_INLET_VEL_FAC
        #    DEFAULT value of 1.0
        #    Error check: Value greater than or equal to 1.0
        #Specify drag under relation factor
        #    Specification only available with MFIX-TFM and MFIX-Hybrid solvers
        #    Sets keyword UR_F_GS
        #    DEFAULT value of 1.0
        #    Error check: Value bounded between 0 and 1
        #Specify IA theory conductivity under relation factor
        #    Specification only available with KT_TYPE = 'IA_NONEP'
        #    Sets keyword IA_NONEP
        #    DEFAULT value of 1.0
        #    Error check: value bounded between 0 and 1
        #Enable four point, fourth order interpolation
        #    Specification always available
        #    Sets keyword FPFOI
        #    DEFAULT value of .FALSE.
        #Specify the universal limiter
        #    Specification only available if four point, forth order interpolation is enabled
        #    Sets keyword C_FAC
        #    DEFAULT value of 1.0
        #    Error check: value bounded between 0 and 1
        pass
