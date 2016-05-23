#!/usr/bin/env python
"""MFIX GUI"""

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import copy
import getopt
import logging
import os
import shutil
import signal
import sys
import time

from copy import deepcopy
from collections import OrderedDict

# Initialize logger early
SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
log = logging.getLogger(__name__)
log.debug(SCRIPT_DIRECTORY)

try:
    import requests
except ImportError:
    requests = None
    log.warn("requests module not available")

# import qt
from qtpy import QtCore, QtWidgets, QtGui, PYQT4, PYQT5
from qtpy.QtCore import (QObject, QThread, pyqtSignal, QUrl, QSettings,
                         Qt)

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16

PRECOMPILE_UI = False

if not PRECOMPILE_UI:
    try:
        from PyQt5 import uic
    except ImportError:
        from PyQt4 import uic

# Debugging hooks
def debug_trace(__name__):
    """Set a tracepoint in the Python debugger that works with Qt"""
    from qtpy.QtCore import pyqtRemoveInputHook
    from pdb import set_trace
    pyqtRemoveInputHook(__name__)
    set_trace(__name__)


# local imports
from project_manager import ProjectManager
from mfix_threads import MfixThread, MonitorThread

from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table)
from widgets.regions import RegionsWidget
from widgets.linear_equation_table import LinearEquationTable
from widgets.species_popup import SpeciesPopup


from tools.general import (make_callback, get_icon,
                           widget_iter, set_script_directory,
                           format_key_with_args, to_unicode_from_fs,
                           set_item_noedit, get_selected_row)

set_script_directory(SCRIPT_DIRECTORY)


from tools.namelistparser import buildKeywordDoc

from constants import *


# look for pyqtnode
try:
    from pyqtnode import NodeWidget
except ImportError:
    NodeWidget = None

if PRECOMPILE_UI:
    try:
        from uifiles.general import Ui_general
        from uifiles.geometry import Ui_geometry
        from uifiles.gui import Ui_MainWindow
        from uifiles.mesh import Ui_mesh
        from uifiles.model_setup import Ui_model_setup
        from uifiles.solids import Ui_solids
        from uifiles.monitors import Ui_monitors
        from uifiles.numerics import Ui_numerics
        from uifiles.output import Ui_output
        from uifiles.post_processing import Ui_post_processing
        from uifiles.run import Ui_run
        from uifiles.vtk import Ui_vtk
    except ImportError:
        print("You must compile ui files!  cd uifiles; make")
        sys.exit(1)


# --- Main Gui ---


class MfixGui(QtWidgets.QMainWindow):
    """Main window class handling all gui interactions"""

    settings = QSettings('MFIX', 'MFIX')

    def __init__(self, app, parent=None, project_file=None):
        # load settings early so get_project_file returns the right thing.
        if project_file:
            self.set_project_file(project_file)

        QtWidgets.QMainWindow.__init__(self, parent)

        # reference to qapp instance
        self.app = app
        self.use_vtk = 'MFIX_NO_VTK' not in os.environ

        self.mfix_config = None
        self.smp_enabled = False
        self.dmp_enabled = False
        self.pymfix_enabled = False

        # load ui file
        self.customWidgets = {'LineEdit':      LineEdit,
                              'CheckBox':      CheckBox,
                              'ComboBox':      ComboBox,
                              'DoubleSpinBox': DoubleSpinBox,
                              'SpinBox':       SpinBox,
                              'Table':         Table,
                              }

        if PRECOMPILE_UI:
            self.ui = Ui_MainWindow()
            self.ui.setupUi(self)

            def make_widget(cls):
                # Create an instance of a new class which is a subclass
                # of QWidget and the specified class
                class Widget(QtWidgets.QWidget, cls):
                    def __init__(self):
                        QtWidgets.QWidget.__init__(self, parent)
                        self.setupUi(self)
                return Widget()

            for cls in (Ui_general, Ui_geometry, Ui_mesh, RegionsWidget,
                        Ui_model_setup, Ui_solids, Ui_numerics, Ui_output, Ui_vtk,
                        Ui_monitors, Ui_post_processing, Ui_run):
                if cls == RegionsWidget: # not loaded from ui file
                    widget = RegionsWidget()
                    name = 'regions'
                else:
                    widget = make_widget(cls)
                    name = cls.__name__.split('_',1)[1] # part after "Ui_"
                # assign 'self.ui.general', etc
                setattr(self.ui, name, widget)
                self.ui.stackedWidgetTaskPane.addWidget(widget)

        else:  # not precompiled
            uifiles = os.path.join(SCRIPT_DIRECTORY, 'uifiles')
            self.ui = uic.loadUi(os.path.join(uifiles, 'gui.ui'))
            self.setCentralWidget(self.ui)
            assert self is not self.ui

            for name in ('general', 'geometry', 'mesh', 'regions',
                         'model_setup', 'solids', 'numerics', 'output', 'vtk',
                         'monitors', 'run'):
                if name == 'regions':  # not loaded from .ui file
                    widget = RegionsWidget()
                else:
                    widget = QtWidgets.QWidget()
                    uic.loadUi(os.path.join(uifiles, name+'.ui'), widget)
                    # assign 'self.ui.general', etc
                setattr(self.ui, name, widget)
                self.ui.stackedWidgetTaskPane.addWidget(widget)
        # end of ui loading

        self.species_popup = SpeciesPopup(QtWidgets.QDialog())
        #self.species_popup.setModal(True) # ?
        self.number_source_lines = True # No user control for this, yet

        # set title and icon
        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        # build keyword documentation from namelist docstrings
        self.keyword_doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY,
                                                        os.pardir, 'model'))

        # create project manager
        # NOTE.  it's a ProjectManager, not a Project.  But
        # ProjectManager is a subclass of Project.  Please
        # do not "fix" the code by renaming self.project to
        # self.project_manager
        self.project = ProjectManager(self, self.keyword_doc)

        # note, the set_fluid_*_model methods have a lot of repeated code
        # see 'set_fluid_mol_weight_model' to help understand this
        def make_fluid_model_setter(self, name, key):
            def setter(model):
                setattr(self, name, model) # self.fluid_<name>_model = model
                combobox = getattr(self.ui, 'combobox_' + name)
                prev_model = combobox.currentIndex()
                if model != prev_model:
                    combobox.setCurrentIndex(model)
                    return

                # Enable spinbox for constant model
                key_g0 = key + "_g0"
                key_usr = "usr_" + key + "g"
                spinbox = getattr(self.ui, 'spinbox_keyword_%s' % key_g0)
                spinbox.setEnabled(model==CONSTANT)

                if model == CONSTANT:
                    value = spinbox.value() # Possibly re-enabled gui item
                    if self.project.get_value(key_g0) != value:
                        self.set_keyword(key_g0, value) # Restore keyword value
                elif model == UDF:
                    self.unset_keyword(key_g0)
                    self.set_keyword(key_usr, True)
                else: # Ideal gas law, Sutherland, etc
                    self.unset_keyword(key_g0)
                    self.unset_keyword(key_usr)
                    # anything else to do in this case? validation?
            return setter

        # Create setters for the cases which are similar (mol. wt. handled separately)
        for (name, key) in (
                ('density', 'ro'),
                ('viscosity', 'mu'),
                ('specific_heat', 'cp'),
                ('conductivity', 'k'),
                ('diffusion', 'dif')):
            model_name = 'fluid_%s_model' % name
            setattr(self, 'set_'+model_name, make_fluid_model_setter(self, model_name, key))


        # --- data ---
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper,
                               }

        self.animation_speed = 400
        self.animating = False
        self.stack_animation = None

        # ---- parameters which do not map neatly to keywords
        self.fluid_nscalar_eq = 0
        self.solid_nscalar_eq = 0 # Infer these from phase4scalar
        # Defaults
        #self.solver = SINGLE - moved to Project
        self.fluid_density_model = CONSTANT
        self.fluid_viscosity_model = CONSTANT
        self.fluid_mol_weight_model = CONSTANT
        self.fluid_specific_heat_model = CONSTANT
        self.fluid_conductivity_model = AIR
        self.fluid_diffusion_model = AIR

        # --- icons ---
        # loop through all widgets, because I am lazy
        for widget in widget_iter(self):
            if isinstance(widget, QtWidgets.QToolButton):
                name = str(widget.objectName())
                if 'add' in name:
                    widget.setIcon(get_icon('add.png'))
                elif 'delete' in name or 'remove' in name:
                    widget.setIcon(get_icon('remove.png'))
                elif 'copy' in name:
                    widget.setIcon(get_icon('copy.png'))

        self.ui.toolbutton_new.setIcon(get_icon('newfolder.png'))
        self.ui.toolbutton_open.setIcon(get_icon('openfolder.png'))
        self.ui.toolbutton_save.setIcon(get_icon('save.png'))
        self.ui.toolbutton_save_as.setIcon(get_icon('save.png'))
        self.ui.toolbutton_export.setIcon(get_icon('export.png'))

        self.ui.toolbutton_run.setIcon(get_icon('play.png'))
        self.ui.toolbutton_reset.setIcon(get_icon('restart.png'))

        self.ui.geometry.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        self.ui.geometry.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        self.ui.geometry.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        self.ui.geometry.toolbutton_geometry_intersect.setIcon(
            get_icon('intersect.png'))
        self.ui.geometry.toolbutton_geometry_difference.setIcon(
            get_icon('difference.png'))

        # --- Connect Signals to Slots---
        # open/save/new project
        self.ui.toolbutton_open.clicked.connect(self.handle_open_action)
        self.ui.toolbutton_save.clicked.connect(self.save_project)
        self.ui.toolbutton_save_as.clicked.connect(self.handle_save_as)
        self.ui.toolbutton_export.clicked.connect(self.handle_export)

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.clicked.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        self.ui.treewidget_model_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # build/run/connect MFIX
        self.ui.run.button_run_mfix.clicked.connect(self.handle_run_stop)
        self.ui.run.button_pause_mfix.clicked.connect(self.pause_mfix)
        self.ui.run.button_reset_mfix.clicked.connect(self.remove_output_files)
        self.ui.toolbutton_run.clicked.connect(self.handle_run_stop)
        self.ui.toolbutton_reset.clicked.connect(self.remove_output_files)
        self.ui.run.spinbox_mfix_executables.activated.connect(self.handle_select_executable)

        # Print welcome message.  Do this early so it appears before any
        # other messages that may occur during this __init__
        self.print_welcome()

        # --- Threads ---
        self.run_thread = MfixThread(self)
        self.monitor_thread = MonitorThread(self)

        def make_handler(qtextbrowser):
            " make a closure to read output from external process "

            def handle_line(line, color=None): # Combine with print_internal
                " closure to read output from external process "
                log.debug(str(line).strip())
                cursor = qtextbrowser.textCursor()
                cursor.movePosition(cursor.End)
                char_format = QtGui.QTextCharFormat()
                if color:
                    char_format.setForeground(QtGui.QColor(color))
                cursor.setCharFormat(char_format)
                cursor.insertText(line)
                scrollbar = qtextbrowser.verticalScrollBar()
                scrollbar.setValue(scrollbar.maximum())

            return handle_line

        self.run_thread.line_printed.connect(self.print_internal)
        self.run_thread.update_run_options.connect(self.update_run_options)

        self.monitor_thread.sig.connect(self.update_run_options)
        self.monitor_thread.start()

        # --- setup widgets ---
        self.__setup_simple_keyword_widgets()
        self.__setup_other_widgets()  # refactor/rename - cgw

        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        if NodeWidget is not None:
            self.__setup_workflow_widget()

        # --- default ---
        self.mode_changed('modeler')
        self.change_pane('geometry')

        # some data fields, these should probably be in Project
        self.fluid_species = OrderedDict()
        self.saved_fluid_species = None
        self.solids = OrderedDict()
        self.current_solid_index = None

        # Update run options
        self.update_run_options()

    def set_keyword(self, key, value, args=None):
        """convenience function to set keyword"""
        self.project.submit_change(None, {key:value}, args)

    def update_keyword(self, key, value, args=None):
        """like set_keyword but no action if value already set"""
        if self.project.get_value(key, args=args) == value:
            return
        self.set_keyword(key, value, args)

    def unset_keyword(self, key, args=None):
        """convenience function to undefine keyword"""
        if isinstance(args, int):
            args = [args]
        elif args is None:
            args = []
        success = self.project.removeKeyword(key, args, warn=False)
        if success:
            self.print_internal("%s" % format_key_with_args(key, args),
                                font='strikeout')


    def is_project_open(self):
        """returns False if and only if no project is loaded and main window invisible"""
        return self.ui.stackedwidget_mode.isVisible()

    def update_run_options(self):
        """Updates list of of mfix executables and sets run dialog options"""

        if not self.is_project_open():
            return

        running = (self.run_thread.mfixproc is not None)
        not_running = not running

        res_file_exists = bool(self.monitor_thread.get_res())

        self.handle_select_executable()

        current_selection = self.ui.run.spinbox_mfix_executables.currentText()
        self.ui.run.spinbox_mfix_executables.clear()
        output = self.monitor_thread.get_executables()
        for executable in output:
            self.ui.run.spinbox_mfix_executables.addItem(executable)
        if current_selection in self.monitor_thread.executables:
            self.ui.run.spinbox_mfix_executables.setEditText(current_selection)
        self.handle_select_executable()

        mfix_available = bool(output)

        self.ui.run.spinbox_mfix_executables.setEnabled(not_running)
        self.ui.run.button_pause_mfix.setEnabled(running and self.pymfix_enabled)
        self.ui.run.openmp_threads.setEnabled(not_running and self.smp_enabled)
        self.ui.run.spinbox_keyword_nodesi.setEnabled(not_running and self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesj.setEnabled(not_running and self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesk.setEnabled(not_running and self.dmp_enabled)

        self.ui.run.spinbox_mfix_executables.setVisible(mfix_available)
        self.ui.run.spinbox_mfix_executables_warning.setVisible(not mfix_available)

        self.ui.run.use_spx_checkbox.setEnabled(res_file_exists)
        self.ui.run.use_spx_checkbox.setChecked(res_file_exists)
        self.ui.toolbutton_reset.setEnabled(res_file_exists and not_running)
        self.ui.run.button_reset_mfix.setEnabled(res_file_exists and not_running)

        if res_file_exists:
            if running:
                self.ui.run.button_run_mfix.setText("Stop")
                self.ui.toolbutton_run.setIcon(get_icon('stop.png'))
                self.ui.toolbutton_run.setText("Stop")
            else:
                self.ui.run.button_run_mfix.setText("Resume")
                self.ui.toolbutton_run.setIcon(get_icon('play.png'))
                self.ui.toolbutton_run.setText("Resume")
        else:
            if running:
                self.ui.run.button_run_mfix.setText("Stop")
                self.ui.toolbutton_run.setIcon(get_icon('stop.png'))
                self.ui.toolbutton_run.setText("Stop")
            else:
                self.ui.run.button_run_mfix.setText("Run")
                self.ui.toolbutton_run.setIcon(get_icon('play.png'))
                self.ui.toolbutton_run.setText("Run")


    def print_welcome(self):
        self.print_internal("Welcome to MFIX - https://mfix.netl.doe.gov",
                            color='blue')
        self.print_internal("MFIX-GUI version %s" % self.get_version(),
                            color='blue')

    def get_version(self):
        return "0.2x" # placeholder

    def closeEvent(self, event):

        log.debug('Closing GUI')

        # cleanup threads
        self.run_thread.quit()
        self.monitor_thread.quit()
        if hasattr(self, 'update_residuals_thread'):
            self.update_residuals_thread.quit()

        event.accept()

    def set_navigation_item_state(self, item_name, state):
        tree = self.ui.treewidget_model_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        items = tree.findItems(item_name, flags, 0)

        assert len(items)==1, "%s menu items matching %s"%("Multiple" if len(items)>1 else "No",
                                                           item_name)
        for item in items:
            item.setDisabled(not state)


    # Top-level "Model Setup"
    def set_solver(self, solver):
        """handler for "Solver" combobox in Model Setup"""
        self.project.solver = solver
        if solver is None: #
            log.warn("set_solver called with solver=None")
            return

        ui = self.ui
        model_setup = ui.model_setup
        cb = model_setup.combobox_solver
        if cb.currentIndex != solver:
            cb.setCurrentIndex(solver)

        solver_name = model_setup.combobox_solver.currentText()
        self.print_internal("set solver to %s" % solver_name)

        item_names =  ("Solids", "Continuum Solids Model",
                       "Discrete Element Model", "Particle in Cell Model")

        item_states = {SINGLE: (False, False, False, False),
                       TFM: (True, True, False, False),
                       DEM: (True, False, True, False),
                       PIC: (True, False, False, True),
                       HYBRID: (True, True, True, False)}

        #for item_name, item_state in zip(item_names, item_states[solver]):
        #    self.set_navigation_item_state(item_name, item_state)

        # Options which require TFM, DEM, or PIC
        enabled = solver in (TFM, DEM, PIC)
        interphase = model_setup.groupbox_interphase
        interphase.setEnabled(enabled)

        # TFM only
        # use a groupbox here, instead of accessing combobox + label?
        enabled = (solver == TFM)
        model_setup.combobox_subgrid_model.setEnabled(enabled)
        model_setup.label_subgrid_model.setEnabled(enabled)
        model_setup.groupbox_subgrid_params.setEnabled(enabled and
                                                       self.subgrid_model > 0)

        ui.checkbox_enable_fluid_scalar_eq.setEnabled(enabled)
        ui.spinbox_fluid_nscalar_eq.setEnabled(enabled
                    and self.ui.checkbox_enable_fluid_scalar_eq.isChecked())

        # Solids Model selection tied to Solver
        # XXX What to do about solids that are already defined?
        self.setup_combobox_solids_model()

    def enable_energy_eq(self, state):
        # Additional callback on top of automatic keyword update,
        # since this has to change availabilty of several other GUI items
        self.ui.model_setup.checkbox_keyword_energy_eq.setChecked(state)
        ui = self.ui
        for item in (ui.combobox_fluid_specific_heat_model,
                     ui.combobox_fluid_conductivity_model,
                     # more ?
                     ):
            item.setEnabled(state)

        # TODO: these should probably be lineEdits not spinboxes
        spinbox = ui.spinbox_keyword_cp_g0 # cp_g0 == specific heat for fluid phase
        if state:
            spinbox.setEnabled(self.fluid_specific_heat_model == CONSTANT)
        else:
            spinbox.setEnabled(False)


    def set_subgrid_model(self, index):
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)

    def update_scalar_equations(self, prev_nscalar):
        # This is a little messy.  We may have reduced
        # nscalar, so we need to unset phase4scalar(i)
        # for any values of i > nscalar.
        nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
        if nscalar > 0:
            self.project.submit_change(None,{"nscalar": nscalar})
        else:
            self.unset_keyword("nscalar")

        for i in range(1,1+self.fluid_nscalar_eq):
            self.project.submit_change(None,{"phase4scalar":0},args=i)
        for i in range(1+self.fluid_nscalar_eq, 1+prev_nscalar):
            if self.project.get_value("phase4scalar", args=i)  == 0:
                self.unset_keyword("phase4scalar", i)


    ## Fluid phase
    def enable_fluid_species_eq(self, state):
        ui = self.ui
        for item in (ui.combobox_fluid_diffusion_model,
                     # more ?
                     ):
            item.setEnabled(state)
        spinbox = ui.spinbox_keyword_dif_g0 # dif_g0 == diffusion coeff model
        if state:
            spinbox.setEnabled(self.fluid_diffusion_model == CONSTANT)
        else:
            spinbox.setEnabled(False)

    def enable_fluid_scalar_eq(self, state):
        spinbox = self.ui.spinbox_fluid_nscalar_eq
        spinbox.setEnabled(state)
        if state:
            val = spinbox.value()
            self.set_fluid_nscalar_eq(val)
        else:
            # Don't call set_fluid_nscalar_eq(0) b/c that will clobber spinbox
            prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
            self.fluid_nscalar_eq = 0
            self.update_scalar_equations(prev_nscalar)

    def set_fluid_nscalar_eq(self, value):
        # This *sums into* nscalar - not a simple keyword
        prev_nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
        self.fluid_nscalar_eq = value
        spinbox = self.ui.spinbox_fluid_nscalar_eq
        if value != spinbox.value():
            spinbox.setValue(value)
            return
        self.update_scalar_equations(prev_nscalar)


    # molecular wt model only has 2 choices, and the key names don't
    # follow the same pattern, so create its setter specially
    def set_fluid_mol_weight_model(self, model):
        self.fluid_mol_weight_model = model
        combobox = self.ui.combobox_fluid_mol_weight_model
        prev_model = combobox.currentIndex()
        if model != prev_model:
            combobox.setCurrentIndex(model)
            return
        # Enable spinbox for constant mol_weight model
        spinbox = self.ui.spinbox_keyword_mw_avg
        spinbox.setEnabled(model==CONSTANT)
        if model == CONSTANT:
            value = spinbox.value() # Possibly re-enabled gui item
            if self.project.get_value("mw_avg") != value:
                self.set_keyword("mw_avg", value) # Restore keyword value
        else: # Mixture
            # TODO: validate, require mw for all component species
            self.unset_keyword("mw_avg")

    def set_fluid_phase_name(self, value):
        if value != self.ui.lineedit_fluid_phase_name.text():
            self.ui.lineedit_fluid_phase_name.setText(value)
            return
        value = None if value=='Fluid' else value # don't save default
        self.project.mfix_gui_comments['fluid_phase_name'] = value

    def disable_fluid_solver(self, state):
        enabled = not state # "disable"
        self.set_navigation_item_state("Fluid", enabled)
        ms = self.ui.model_setup
        checkbox = ms.checkbox_enable_turbulence
        checkbox.setEnabled(enabled)
        ms.combobox_turbulence_model.setEnabled(enabled and
                                                checkbox.isChecked())

    # Solids phase
    def setup_combobox_solids_model(self):
        """solids model combobox is tied to solver setting"""
        solver = self.project.solver
        if solver == SINGLE:
            # Note, if Single-Phase solver is enabled, this pane is disabled
            return
        return #tmp
        cb = self.ui.combobox_solids_model
        model = cb.model()
        #          TFM,  DEM,  PIC
        enabled = [False, False, False]
        enabled[0] = (solver==TFM or solver==HYBRID)
        enabled[1] = (solver==DEM or solver==HYBRID)
        enabled[2] = (solver==PIC)
        for (i, e) in enumerate(enabled):
            item = model.item(i, 0)
            flags = item.flags()
            if e:
                flags |= QtCore.Qt.ItemIsEnabled
            else:
                flags &= ~QtCore.Qt.ItemIsEnabled
            item.setFlags(flags)
        i = cb.currentIndex()
        if not enabled[i]:
            # Current selection no longer valid, so pick first valid choice
            # Don't leave a non-enabled item selected!
            j = enabled.index(True)
            cb.setCurrentIndex(j)

    def handle_combobox_solids_model(self, index):
        if self.current_solid_index is None:
            return # shouldn't get here
        models = ('TFM', 'DEM', 'PIC')
        self.update_keyword('solids_model', models[index], args=self.current_solid_index)


    # helper functions for __init__
    def __setup_other_widgets(self): # rename/refactor
        """setup widgets which are not tied to a simple keyword"""
        ui = self.ui
        model_setup = ui.model_setup
        combobox = model_setup.combobox_solver
        # activated: Only on user action, avoid recursive calls in set_solver
        combobox.activated.connect(self.set_solver)

        checkbox = model_setup.checkbox_disable_fluid_solver
        checkbox.stateChanged.connect(self.disable_fluid_solver)
        self.disable_fluid_solver(False)

        checkbox = model_setup.checkbox_keyword_energy_eq
        checkbox.stateChanged.connect(self.enable_energy_eq)

        checkbox = ui.checkbox_keyword_species_eq_args_0
        checkbox.stateChanged.connect(self.enable_fluid_species_eq)

        combobox = model_setup.combobox_subgrid_model
        combobox.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)

        #self.enable_energy_eq(False)

        # Fluid phase
        ui.lineedit_fluid_phase_name.textChanged.connect(
            self.set_fluid_phase_name)
        ui.checkbox_enable_fluid_scalar_eq.stateChanged.connect(
            self.enable_fluid_scalar_eq)
        ui.spinbox_fluid_nscalar_eq.valueChanged.connect(
            self.set_fluid_nscalar_eq)

        # Fluid phase models
        # Density
        for name in ('density', 'viscosity', 'mol_weight',
                     'specific_heat', 'mol_weight', 'conductivity', 'diffusion'):
            combobox = getattr(ui, 'combobox_fluid_%s_model' % name)
            setter = getattr(self,'set_fluid_%s_model' % name)
            combobox.currentIndexChanged.connect(setter)

        # Fluid species
        tb = ui.toolbutton_fluid_species_add
        tb.clicked.connect(self.fluid_species_add)
        tb = ui.toolbutton_fluid_species_copy # misnomer
        tb.clicked.connect(self.fluid_species_edit)
        tb.setEnabled(False)
        tb = ui.toolbutton_fluid_species_delete
        tb.setEnabled(False)
        tb.clicked.connect(self.fluid_species_delete)
        tw = ui.tablewidget_fluid_species
        tw.itemSelectionChanged.connect(self.handle_fluid_species_selection)

        # Solid phase
        tb = ui.solids.toolbutton_solids_add
        tb.clicked.connect(self.solids_add)
        tb = ui.solids.toolbutton_solids_delete
        tb.clicked.connect(self.solids_delete)
        tb.setEnabled(False)
        #tw = ui.solids.tablewidget_solids
        #tw.itemSelectionChanged.connect(self.handle_solids_table_selection)
        #cb = ui.solids.combobox_solids_model
        #cb.currentIndexChanged.connect(self.handle_combobox_solids_model)

        # connect solid tab btns
        for i, btn in enumerate([self.ui.solids.pushbutton_material,
                                 self.ui.solids.pushbutton_continuum,
                                 self.ui.solids.pushbutton_discrete,
                                 self.ui.solids.pushbutton_parcel]):
            btn.pressed.connect(
                make_callback(self.solids_change_tab, i, btn))

        # numerics
        ui.linear_eq_table = LinearEquationTable(ui.numerics)
        ui.numerics.gridlayout_leq.addWidget(ui.linear_eq_table)
        self.project.register_widget(ui.linear_eq_table,
                                     ['leq_method', 'leq_tol', 'leq_it',
                                      'leq_sweep', 'leq_pc', 'ur_fac'],
                                     args='*')


    def __setup_simple_keyword_widgets(self):
        """Look for and connect simple keyword widgets to the project manager.
        Keyword information from the namelist doc strings is added to each
        keyword widget. The widget must be named: *_keyword_<keyword> where
        <keyword> is the actual keyword. For example:
        lineedit_keyword_run_name"""

        def try_int(str):
            try:
                return int(str)
            except ValueError:
                return str

        # loop through all widgets looking for *_keyword_<keyword>
        for widget in widget_iter(self):
            name_list = str(widget.objectName()).lower().split('_')

            if 'keyword' in name_list:
                keyword_idx = name_list.index('keyword')
                args = None
                # Look for _args_ following <keyword>
                if 'args' in name_list:
                    args_idx = name_list.index('args')
                    args = map(try_int, name_list[args_idx+1:])
                    keyword = '_'.join(name_list[keyword_idx+1:args_idx])
                else:
                    keyword = '_'.join(name_list[keyword_idx+1:])

                # set the key attribute to the keyword
                widget.key = keyword
                widget.args = args

                # add info from keyword documentation
                if keyword in self.keyword_doc:
                    doc = self.keyword_doc[keyword]
                    widget.setdtype(doc['dtype'])
                    if 'required' in doc:
                        widget.setValInfo(
                            req=(doc['required'] == 'true'))
                    if 'validrange' in doc:
                        vr = doc['validrange']
                        if 'max' in vr:
                            widget.setValInfo(_max=vr['max'])
                        if 'min' in doc['validrange']:
                            widget.setValInfo(_min=vr['min'])

                    if 'initpython' in doc:
                        widget.default(doc['initpython'])

                    if isinstance(widget, QtWidgets.QComboBox) and widget.count() < 1:
                            widget.addItems(list(doc['valids'].keys()))

                # register the widget with the project manager
                self.project.register_widget(widget, keys=[keyword], args=args)

                # connect to unsaved method
                widget.value_updated.connect(self.unsaved)

    def __setup_vtk_widget(self):
        """initialize the vtk widget"""
        if not self.use_vtk:
            self.vtkwidget = None
            return

        from widgets.vtkwidget import VtkWidget

        self.vtkwidget = VtkWidget(self.project, parent=self)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkwidget)

        # register with project manager
        self.project.register_widget(self.vtkwidget,
                                     ['xmin', 'xlength', 'ymin', 'ylength',
                                      'zmin', 'zlength', 'imax', 'jmax',
                                      'kmax'])

        # add reference to other widgets
        self.ui.regions.vtkwidget = self.vtkwidget

    def __setup_workflow_widget(self):

        self.nodeChart = NodeWidget(showtoolbar=False)
        # Build default node library
        self.nodeChart.nodeLibrary.buildDefaultLibrary()
        self.ui.horizontalLayoutPyqtnode.addWidget(self.nodeChart)

    @classmethod
    def get_project_file(cls):
        """get the project filename, including full path"""
        last = cls.settings.value('project_file')
        return last if last else None

    @classmethod
    def set_project_file(cls, value):
        cls.settings.setValue('project_file', value)

    def get_project_dir(self):
        """get the current project directory"""
        project_file = self.get_project_file()
        return os.path.dirname(project_file) if project_file else None

    def mode_changed(self, mode):
        """change the Modeler, Workflow, Developer tab"""
        current_index = 0
        for i in range(self.ui.stackedwidget_mode.count()):
            widget = self.ui.stackedwidget_mode.widget(i)
            if mode == str(widget.objectName()):
                current_index = i
                break

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

        self.animate_stacked_widget(self.ui.stackedwidget_mode,
                                    self.ui.stackedwidget_mode.currentIndex(),
                                    current_index,
                                    'horizontal')

    # --- modeler pane navigation ---
    def change_pane(self, name):
        """change to the specified pane"""
        clist = self.ui.treewidget_model_navigation.findItems(
                    name,
                    Qt.MatchContains | Qt.MatchRecursive,
                    0)

        for item in clist:
            if str(item.text(0)).lower() == name.lower():
                break

        self.ui.treewidget_model_navigation.setCurrentItem(item)
        self.navigation_changed()

    def navigation_changed(self):
        """an item in the tree was selected, change panes"""
        current_selection = self.ui.treewidget_model_navigation.selectedItems()

        # Force any open popup to close
        # if dialog is modal we don't need this
        self.species_popup.done(0)

        if current_selection:
            text = str(current_selection[-1].text(0))
            text = '_'.join(text.lower().split(' '))
            current_index = 0
            for i in range(self.ui.stackedWidgetTaskPane.count()):
                widget = self.ui.stackedWidgetTaskPane.widget(i)
                #print(text, str(widget.objectName()))
                if text == str(widget.objectName()):
                    current_index = i
                    break
            self.animate_stacked_widget(
                self.ui.stackedWidgetTaskPane,
                self.ui.stackedWidgetTaskPane.currentIndex(),
                current_index)

    # --- animation methods ---
    def animate_stacked_widget(self, stackedwidget, from_, to,
                               direction='vertical', line=None, to_btn=None,
                               btn_layout=None):
        """animate changing of qstackedwidget"""

        # check to see if already animating
        if self.animating and self.stack_animation is not None:
            self.stack_animation.stop()

        from_widget = stackedwidget.widget(from_)
        to_widget = stackedwidget.widget(to)

        # get from geometry
        width = from_widget.frameGeometry().width()
        height = from_widget.frameGeometry().height()

        # offset
        # bottom to top
        if direction == 'vertical' and from_ < to:
            offsetx = 0
            offsety = height
        # top to bottom
        elif direction == 'vertical' and from_ > to:
            offsetx = 0
            offsety = -height
        elif direction == 'horizontal' and from_ < to:
            offsetx = width
            offsety = 0
        elif direction == 'horizontal' and from_ > to:
            offsetx = -width
            offsety = 0
        else:
            return

        # move to widget and show
        # set the geometry of the next widget
        to_widget.setGeometry(0 + offsetx, 0 + offsety, width, height)
        to_widget.show()
        to_widget.raise_()
        #to_widget.activateWindow() ? needed?

        # animate
        # from widget
        animnow = QtCore.QPropertyAnimation(from_widget, "pos".encode('utf-8'))
        animnow.setDuration(self.animation_speed)
        animnow.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnow.setStartValue(
            QtCore.QPoint(0,0))
        animnow.setEndValue(
            QtCore.QPoint(0 - offsetx,
                          0 - offsety))

        # to widget
        animnext = QtCore.QPropertyAnimation(to_widget, "pos".encode('utf-8'))
        animnext.setDuration(self.animation_speed)
        animnext.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnext.setStartValue(
            QtCore.QPoint(0 + offsetx,
                          0 + offsety))
        animnext.setEndValue(
            QtCore.QPoint(0,0))

        # line
        animline = None
        if line is not None and to_btn is not None:
            animline = QtCore.QPropertyAnimation(line, "pos".encode('utf-8'))
            animline.setDuration(self.animation_speed)
            animline.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
            animline.setStartValue(
                QtCore.QPoint(line.geometry().x(),
                              line.geometry().y()))
            animline.setEndValue(
                QtCore.QPoint(to_btn.geometry().x(),
                              line.geometry().y()))

        # animation group
        self.stack_animation = QtCore.QParallelAnimationGroup()
        self.stack_animation.addAnimation(animnow)
        self.stack_animation.addAnimation(animnext)
        if animline is not None:
            self.stack_animation.addAnimation(animline)
        self.stack_animation.finished.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, line)
            )
        self.stack_animation.stateChanged.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, line))

        self.animating = True
        self.stack_animation.start()

    def animate_stacked_widget_finished(self, widget, from_, to,
                                        btn_layout=None, line=None):
        """cleanup after animation"""
        if self.stack_animation.state() == QtCore.QAbstractAnimation.Stopped:
            widget.setCurrentIndex(to)
            from_widget = widget.widget(from_)
            from_widget.hide()
            from_widget.move(0, 0)

            if btn_layout is not None and line is not None:
                btn_layout.addItem(btn_layout.takeAt(
                    btn_layout.indexOf(line)), 1, to)

            self.animating = False

    # --- helper methods ---
    def message(self,
                title='Warning',
                icon='warning',
                text='This is a warning.',
                buttons=['ok'],
                default='ok',
                infoText=None,
                detailedtext=None,
                ):
        """Create a message box:
        title = 'title'
        icon = 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel',
            'discard'
        default = 'ok' the default selected button
        infotext = 'extended information text'
        detailedtext = 'Some details'"""

        # TODO: tie this in with logging & print_internal
        # TODO: suppress popups in 'test mode'

        msgBox = QtWidgets.QMessageBox(self)
        msgBox.setWindowTitle(title)

        # Icon
        if icon == 'warning':
            icon = QtWidgets.QMessageBox.Warning
        else:
            icon = QtWidgets.QMessageBox.Information

        msgBox.setIcon(icon)

        # Text
        msgBox.setText(text)

        if infoText:
            msgBox.setInformativeText(infoText)

        if detailedtext:
            msgBox.setDetailedText(detailedtext)

        # buttons
        qbuttonDict = {'ok':      QtWidgets.QMessageBox.Ok,
                       'yes':     QtWidgets.QMessageBox.Yes,
                       'no':      QtWidgets.QMessageBox.No,
                       'cancel':  QtWidgets.QMessageBox.Cancel,
                       'discard': QtWidgets.QMessageBox.Discard,
                       }
        for button in buttons:
            msgBox.addButton(qbuttonDict[button])

            if button == default:
                msgBox.setDefaultButton(qbuttonDict[button])

        ret = msgBox.exec_()

        for key, value in qbuttonDict.items():
            if value == ret:
                result = key
                break

        return result

    def print_internal(self, line, color=None, font=None): # basically same as handle_line
        qtextbrowser = self.ui.command_output
        if not line.endswith('\n'):
            line += '\n'
        lower = line.lower()
        msg = line.strip()
        # hack. TODO: real msg types, map to font/color
        strikeout = font and font.lower() == 'strikeout'
        if strikeout:
            msg = "unset " + msg
        if 'warning:' in lower or 'error:' in lower:
            log.warn(msg)
        else:
            log.info(msg)
        cursor = qtextbrowser.textCursor()
        cursor.movePosition(cursor.End)
        char_format = QtGui.QTextCharFormat()
        if color:
            char_format.setForeground(QtGui.QColor(color))
        if font:
            if strikeout: # hack
                char_format.setFontFamily("Monospace")
                char_format.setFontStrikeOut(True)
            else:
                char_format.setFontFamily(font)
        cursor.setCharFormat(char_format)
        cursor.insertText(line)
        scrollbar = qtextbrowser.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def handle_select_executable(self):
        """Enable/disable run options based on selected executable"""
        mfix_exe = self.ui.run.spinbox_mfix_executables.currentText()
        self.mfix_exe = mfix_exe
        if not mfix_exe:
            return
        config = self.monitor_thread.executables[mfix_exe]
        self.mfix_config = config
        self.smp_enabled = 'smp' in config
        self.dmp_enabled = 'dmp' in config

        self.pymfix_enabled = any(mfix_exe.lower().endswith(x)
                                  for x in ('pymfix', 'pymfix.exe'))

        self.ui.run.openmp_threads.setEnabled(self.smp_enabled)
        if not self.dmp_enabled:
            self.ui.run.spinbox_keyword_nodesi.setValue(1)
            self.ui.run.spinbox_keyword_nodesj.setValue(1)
            self.ui.run.spinbox_keyword_nodesk.setValue(1)
        self.ui.run.spinbox_keyword_nodesi.setEnabled(self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesj.setEnabled(self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesk.setEnabled(self.dmp_enabled)

    def remove_output_files(self, output_files=None, message_text=None):
        """ remove MFIX output files from current project directory

        :param output_files: List of patterns to be matched for file removal
        :type output_files: list
        :return: True for success, False for user cancel"""

        if not output_files:
            output_files = self.monitor_thread.get_outputs()
        if not message_text:
            message_text = 'Deleting output files %s' % '\n'.join(output_files)
        confirm = self.message(title='Warning',
                         icon='warning',
                         text=message_text,
                         buttons=['ok','cancel'],
                         default='cancel',
                         )
        if confirm != 'ok':
            return False

        for path in output_files:
            log.debug('Deleting file %s' % path)
            try:
                os.remove(path)
            except OSError as e:
                msg = 'Cannot delete %s: %s' % (path, e.strerror)
                self.print_internal("Warning: %s" % msg, color='red')
                log.warn(msg)
                self.message(title='Warning',
                                 icon='warning',
                                 text=msg,
                                 buttons=['ok'],
                                 default=['ok'])
                break
        return True

    def handle_run_stop(self):
        if self.run_thread.mfixproc is None:
            return self.run_mfix()
        else:
            return self.stop_mfix()

    def run_mfix(self):
        output_files = self.monitor_thread.get_outputs()
        if output_files:
            if self.monitor_thread.get_res():
                return self.resume_mfix()
            if not self.remove_output_files(output_files):
                log.info('output files exist and run was cancelled')
                return
        self.set_keyword('run_type', 'new')
        # FIXME only write it if updated
        self.project.writeDatFile(self.get_project_file()) #XXX
        self._start_mfix()

    def restart_mfix(self):
        """Restart MFIX. This will remove previous output and start a new run."""
        output_files = self.monitor_thread.get_outputs()
        if len(output_files) > 0:
            message = "Delete all output and resume files?"
            if not self.remove_output_files(output_files, message):
                log.debug('output or resume files exist and run was cancelled')
                return
        self.set_keyword('run_type', 'new')
        self.project.writeDatFile(self.get_project_file()) # XXX
        self._start_mfix()

    def pause_mfix(self):
        if not self.pymfix_enabled:
            return
        if requests:
            requests.put(self.pymfix_url + 'stop')
        self.update_run_options()

    def resume_mfix(self):
        """resume previously stopped mfix run"""
        if self.ui.run.use_spx_checkbox.isChecked():
            self.set_keyword('run_type', 'restart_1')
        else:
            # TODO: is it correct to remove all but *.RES ?
            spx_files = self.monitor_thread.get_outputs(['*.SP?', "*.pvd", "*.vtp"])
            if not self.remove_output_files(spx_files):
                # pass here reproduces QThread::wait error somewhat reliably
                log.debug('SP* files exist and run was cancelled')
                #pass
                return
            self.set_keyword('run_type', 'restart_2')
        self.project.writeDatFile(self.get_project_file()) # XXX
        self._start_mfix()

    def _start_mfix(self):
        """start a new local MFIX run, using pymfix, mpirun or mfix directly"""

        mfix_exe = self.mfix_exe

        if self.pymfix_enabled:
            # run pymfix.  python or python3, depending on sys.executable
            run_cmd = [sys.executable, mfix_exe]

        else:
            if self.dmp_enabled:
                np = (self.project.nodesi.value *
                      self.project.nodesj.value *
                      self.project.nodesk.value)

                run_cmd = ['mpirun', '-np', str(np), mfix_exe]

                # adjust environment for to-be called process
                # assume user knows what they are doing and don't override vars
                if not "OMP_NUM_THREADS" in environ:
                    os.environ["OMP_NUM_THREADS"] = str(np)
                log.info('Will start mpirun with OMP_NUM_THREADS=%d' % np)

            else:
                # no dmp support
                run_cmd = [mfix_exe]
                np = 1

        project_filename = os.path.basename(self.get_project_file())
        # Warning, not all versions of mfix support '-f' !
        run_cmd += ['-f', project_filename]

        msg = 'Running %s' % ' '.join(run_cmd)
        #log.info(msg) # print_internal logs
        self.print_internal(msg, color='blue')

        self.run_thread.start_command(
            cmd=run_cmd,
            cwd=self.get_project_dir(),
            env=os.environ)

        if self.pymfix_enabled:
            # get last pymfix url from config
            # compare last to ui.run.{something} to pick up user change
            # also set this elsewhere ...

            time.sleep(2) # this has to go away .. put pymfix api calls
            # into another thread, attach signals to update ui

            self.pymfix_url = 'http://127.0.0.1:5000/'
            if requests:
                requests.put(self.pymfix_url + 'start')

        self.update_run_options()

    def stop_mfix(self):
        """stop locally running instance of mfix"""
        if self.pymfix_enabled:
            if requests:
                requests.put(self.pymfix_url + 'exit')
        # check whether pymfix has stopped mfix then
        # do something here
        self.run_thread.stop_mfix()

        self.update_run_options()

    def update_residuals(self):
        self.ui.residuals.setText(str(self.update_residuals_thread.residuals))
        if self.update_residuals_thread.job_done:
            self.ui.mfix_browser.setHTML('')

    # --- open/save/new ---

    def export_project(self):
        # TODO: combine export_project and save_as
        if self.monitor_thread.get_outputs(["*.SP?"]):
            # TODO: copy SPx files - prompt user for confirmation
            pass
        output_files = self.monitor_thread.get_outputs(["*.RES", "*.STL"])
        project_file = self.get_save_filename()
        if not project_file:
            return
        project_dir = os.path.dirname(project_file)
        if not self.check_writable(project_dir):
            self.save_as()
            return
        for fn in output_files:
            #copy into new_project_directory
            shutil.copyfile(
                    fn, os.path.join(project_dir, os.path.basename(fn)))
        self.save_project(project_file)

    def save_project(self, filename=False):
        """save project, optionally as a new project.

        :param project_file: Filename of project (including path)
        :type project_file: str
        :return: None"""

        if filename:
            project_dir = os.path.dirname(filename)
            project_file = filename
        else:
            project_dir = self.get_project_dir()
            project_file = self.get_project_file()
        if self.use_vtk:
            self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))
        self.project.writeDatFile(project_file) # XXX
        #self.setWindowTitle('MFIX - %s' % project_file)

    def save_as(self):
        """Save As user dialog
        updates project.run_name with user-supplied data
        opens the new project"""
        # TODO: combine export_project and save_as
        project_file = self.get_save_filename()
        if not project_file:
            return
        project_dir = os.path.dirname(project_file)
        new_project_basename = os.path.basename(project_file)
        run_name = os.path.splitext(new_project_basename)[0]
        self.project.run_name.value = run_name

        if not self.check_writable(project_dir):
            self.save_as()
            return
        self.save_project(project_file)
        self.open_project(project_file)

    def get_save_filename(self, dialog_message=None):
        """wrapper for call to getSaveFileName for unit tests to override
        :return: filename
        :rtype: str"""
        if not dialog_message:
            dialog_message = 'Save Project As'
        filename = QtWidgets.QFileDialog.getSaveFileName(
                            self,
                            dialog_message,
                            os.path.join(
                                self.get_project_dir(),
                                self.project.run_name.value + ".mfx"),
                            "*.mfx")
        # User pressed "cancel"
        if not filename:
            return
        # qt4/qt5 compat hack
        #if type(filename) == tuple:
        if PYQT5:
            return filename[0]
        else:
            return filename

    def handle_save(self):
        return self.save_project()

    def handle_export(self):
        return self.export_project()

    def handle_save_as(self):
        return self.save_as()

    def unsaved(self):
        self.setWindowTitle('MFIX - %s *' % self.get_project_file())

    def check_writable(self, directory):
        """check whether directory is writable """
        try:
            import tempfile
            testfile = tempfile.TemporaryFile(dir=directory)
            testfile.close()
            return True

        except Exception as e:
            # maybe to debug, but not to user dialog
            #log.debug(e.message)
            self.message(
                title='Warning',
                icon='warning',
                text="The directory %s is not writable" % directory,
                buttons=['ok'],
                default='ok')

            return False

    def new_project(self, project_dir=None):
        if not project_dir:
            project_dir = str(
                QtWidgets.QFileDialog.getOpenFileName(
                    self, 'Create Project in Directory',
                    ""))
        if len(project_dir) < 1:
            return

        if os.path.isdir(project_path):
            project_dir = project_path
            project_file = os.path.join(project_dir, 'mfix.dat')
        else:
            project_dir = os.path.dirname(project_path)
            project_file = proj
        if not self.check_writable(project_dir):
            return

        shutil.copyfile('mfix.dat.template', project_path)

        self.open_project(project_file)

    def get_open_filename(self):
        """wrapper for call to getSaveFileName for unit tests to override"""
        project_dir = self.get_project_dir()
        return QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Project Directory', project_dir)

    def handle_open_action(self):
        """handler for toolbar Open button"""
        project_path = self.get_open_filename()
        # qt4/qt5 compat hack
        #if type(project_path) == tuple:
        if PYQT5:
            project_path = project_path[0]
        if not project_path:
            return # user pressed Cancel
        self.open_project(project_path)

    def open_project(self, project_path, auto_rename=True):
        """Open MFiX Project"""
        # see also project_manager.load_project_file

        # Make sure path is absolute
        if not os.path.isabs(project_path):
            project_path = os.path.abspath(project_path)
        # "path" may be a directory or a file
        if os.path.isdir(project_path):
            project_dir = project_path
            project_file = os.path.abspath(os.path.join(project_path, 'mfix.dat'))
        else:
            project_dir = os.path.dirname(project_path)
            project_file = project_path

        if not os.path.exists(project_file):
            self.message(title='Warning',
                         icon='warning',
                         text=('%s does not exist' % project_file),
                         buttons=['ok'],
                         default='ok')
            return

        # FIXME. move these into Project so they can get cleared there.
        # (It would be good to instantate a new Project instead of trying
        # to clear all data members)
        self.fluid_species.clear()

        self.ui.lineedit_fluid_phase_name.setText("Fluid") # default
        if self.saved_fluid_species:
            self.saved_fluid_species.clear()

        self.print_internal("Loading %s" % project_file, color='blue')
        try:
            self.project.load_project_file(project_file)
        except Exception as e:
            msg = 'Failed to load %s: %s: %s' % (project_file, e.__class__.__name__, e)
            self.print_internal("Warning: %s" % msg, color='red')
            self.message(title='Warning',
                         icon='warning',
                         text=msg,
                         buttons=['ok'],
                         default='ok')
            import traceback
            traceback.print_exception(*sys.exc_info())
            # Should we stick this in the output window?  no, for now.
            return

        if hasattr(self.project, 'run_name'):
            name = self.project.run_name.value
        else:
            name = 'new_file'
        for char in ('.', '"', "'", '/', '\\', ':'):
            name = name.replace(char, '_')
        runname_mfx = name + '.mfx'

        if auto_rename and not project_path.endswith(runname_mfx):
            ok_to_write = False
            save_msg = 'Saving %s as %s based on run name' % (project_path, runname_mfx)
            response = self.message(title='Info',
                                    icon='info',
                                    text=save_msg,
                                    buttons=['ok', 'cancel'],
                                    default='ok')
            if response == 'ok':
                ok_to_write = True
                renamed_project_file = os.path.join(project_dir, runname_mfx)
                if os.path.exists(renamed_project_file):
                    clobber_msg = '%s exists, replace?' % renamed_project_file
                    response = self.message(title='Warning',
                                    icon='warning',
                                    text=clobber_msg,
                                    buttons=['yes', 'no'],
                                    default='no')
                    if response == 'no':
                        ok_to_write = False

            if ok_to_write:
                project_file = renamed_project_file
                self.project.writeDatFile(project_file)
                self.print_internal(save_msg, color='blue')

        if not self.is_project_open(): # why?
            self.ui.stackedwidget_mode.setVisible(True)

        self.set_project_file(project_file)
        self.setWindowTitle('MFIX - %s' % project_file)

        # read the file (again)
        with open(project_file, 'r') as f:
            src = to_unicode_from_fs(f.read())
        if self.number_source_lines:
            src = '\n'.join('%4d:%s'%(lineno,line)
                            for (lineno,line) in enumerate(src.split('\n'),1))
        self.ui.mfix_dat_source.setPlainText(src)
        # self.mode_changed('developer')

        # Additional GUI setup based on loaded projects (not handled
        # by keyword updates)
        #    .... is there a way to verify that 'energy_eq' is boolean?
        #    should that get set from keyword doc?
        self.enable_energy_eq(bool(self.project.get_value('energy_eq')))

        # cgw - lots more model setup todo here.  Should we do this here or
        #  in ProjectManager.load_project_file (where we do guess/set_solver)

        # values that don't map to keywords, saved as #!MFIX-GUI params
        for (key, val) in self.project.mfix_gui_comments.items():
            if key == 'fluid_phase_name':
                self.set_fluid_phase_name(val)
            # Add more here


        # Fluid phase
        # fluid species table
        self.update_fluid_species_table()

        # fluid momentum and species eq. handled by _keyword_ widget

        # fluid scalar eq
        nscalar = self.project.get_value('nscalar', 0)
        self.fluid_nscalar_eq = sum(1 for i in range(1, nscalar+1)
            if self.project.get_value('phase4scalar', args=i) == 0)
        self.enable_fluid_scalar_eq(self.fluid_nscalar_eq > 0)

        # handle a bunch of items which are essentially the same
        for (setter, name) in ((self.set_fluid_density_model, 'ro'),
                               (self.set_fluid_viscosity_model, 'mu'),
                               (self.set_fluid_specific_heat_model, 'cp'),
                               (self.set_fluid_conductivity_model, 'k'),
                               (self.set_fluid_diffusion_model, 'dif')):
            name_g0 = name+'_g0'
            name_usr = 'usr_'+name+'g'
            val_g0 = self.project.get_value(name_g0)
            val_usr = self.project.get_value(name_usr)

            if (val_usr is not None and val_g0 is not None):
                self.print_internal('Warning: %s and %s are both set' % (name_g0, name_usr))
                # FIXME this is getting printed after error count ... should be included in # of errs
                # (another reason to move this to load_project_file)

            setter(CONSTANT if val_g0 is not None
                   else UDF if val_usr is not None
                   else 1)

        # molecular weight model is the odd one (only 2 settings)
        if self.project.get_value('mw_avg') is not None:
            self.set_fluid_mol_weight_model(CONSTANT)
        else:
            self.set_fluid_mol_weight_model(1)
        # requires molecular weights for all species components, should we valdate

        # Solids
        self.update_solids_table()

        # Look for geometry.stl and load automatically
        geometry = os.path.abspath(os.path.join(project_dir, 'geometry.stl'))
        if os.path.exists(geometry) and self.use_vtk:
            self.vtkwidget.add_stl(None, filename=geometry)

    # --- fluid species methods ---
    def fluid_species_revert(self):
        self.fluid_species = self.saved_fluid_species
        self.species_popup.defined_species = deepcopy(self.fluid_species)
        self.update_fluid_species_table()

    def fluid_species_save(self):
        self.fluid_species = deepcopy(self.species_popup.defined_species)
        self.update_fluid_species_table()

    def update_fluid_species_table(self):
        """Update table in fluids pane.  Also set nmax_g, species_g and species_alias_g keywords,
        which are not tied to a single widget"""

        hv = QtWidgets.QHeaderView
        table = self.ui.tablewidget_fluid_species
        if PYQT5:
            resize = table.horizontalHeader().setSectionResizeMode
        else:
            resize = table.horizontalHeader().setResizeMode
        for n in range(5):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)

        table.clearContents()
        if self.fluid_species is None:
            return
        nrows = len(self.fluid_species)
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem(str(val))
            set_item_noedit(item)
            return item
        old_nmax_g = self.project.get_value('nmax_g')
        nmax_g = len(self.fluid_species)
        if nmax_g > 0 and old_nmax_g is not None:
            self.update_keyword('nmax_g', nmax_g)
        for (row, (species,data)) in enumerate(self.fluid_species.items()):
            for (col, key) in enumerate(('alias', 'phase', 'mol_weight',
                                        'h_f', 'source')):
                alias = data.get('alias', species) # default to species if no alias
                data['alias'] = alias # for make_item
                table.setItem(row, col, make_item(data.get(key)))
                self.update_keyword('species_g', species, args=row+1)
                self.update_keyword('species_alias_g', alias, args=row+1)
                # We're avoiding mw_g in favor of the settings in THERMO DATA
                #self.update_keyword('mw_g', data['mol_weight'], args=row+1)#

        # Clear any keywords with indices above nmax_g
        if old_nmax_g is None:
            old_nmax_g = 0
        for i in range(nmax_g+1, old_nmax_g+1):
            self.unset_keyword('species_g', i)
            self.unset_keyword('species_alias_g', i)
            #self.unset_keyword('mw_g', i) # TBD

        self.project.update_thermo_data(self.fluid_species)

    def handle_fluid_species_selection(self):
        row = get_selected_row(self.ui.tablewidget_fluid_species)
        enabled = (row is not None)
        self.ui.toolbutton_fluid_species_delete.setEnabled(enabled)
        self.ui.toolbutton_fluid_species_copy.setEnabled(enabled)

    def fluid_species_add(self):
        sp = self.species_popup
        sp.phases='GL' # ? is this correct
        # how to avoid this if dialog open already?
        self.saved_fluid_species = copy.deepcopy(self.fluid_species) # So we can revert
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.update_defined_species()
        sp.setWindowTitle("Fluid Species")
        sp.show()
        sp.raise_()
        sp.activateWindow()

    def fluid_species_delete(self):
        # XXX FIXME this is potentially a big problem since
        # it results in species being renumbered, or a hole in
        # the sequence - either way is trouble.  Have to warn
        # user, if species is referenced elsewhere.
        table = self.ui.tablewidget_fluid_species
        row = get_selected_row(table)
        if row is None: # No selection
            return
        table.clearSelection()
        key = self.fluid_species.keys()[row]
        del self.fluid_species[key]
        self.update_fluid_species_table()
        # Sigh, we have to update the row in the popup too.
        # Should the popup just be modal, to avoid this?
        sp = self.species_popup
        sp.defined_species = self.fluid_species
        sp.update_defined_species()

    def fluid_species_edit(self):
        table = self.ui.tablewidget_fluid_species
        row = get_selected_row(table)
        sp = self.species_popup
        self.saved_fluid_species = copy.deepcopy(self.fluid_species) # So we can revert
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species
        sp.update_defined_species()
        if row is None:
            sp.tablewidget_defined_species.clearSelection()
        else:
            sp.tablewidget_defined_species.setCurrentCell(row, 0)
        sp.setWindowTitle("Fluid Species")
        sp.show()
        sp.raise_()
        sp.activateWindow()

    # --- solids phase methods ---
    def make_solids_name(self):
        n = 1
        while True:
            name = 'Solid %d' % n
            if name not in self.solids:
                break
            n += 1
        return name

    def solids_add(self):
        return
        tw = self.ui.tablewidget_solids
        nrows = tw.rowCount()
        tw.setRowCount(nrows + 1)
        name = self.make_solids_name()
        if self.project.solver == SINGLE: # Should not get here! this pane is disabled.
            return
        else:
            model = [None, 'TFM', 'DEM', 'PIC', 'TEM'][self.project.solver]
        diameter = 0.0
        density = 0.0
        self.solids[name] = {'model': model,
                             'diameter': diameter,
                             'density': density} # more?
        self.update_solids_table()
        tw.setCurrentCell(nrows, 0) # Select new item

    def handle_solids_table_selection(self):
        return
        tw = self.ui.tablewidget_solids
        row = get_selected_row(tw)
        enabled = (row is not None)
        self.ui.toolbutton_solids_delete.setEnabled(enabled)
        self.ui.toolbutton_solids_copy.setEnabled(enabled)
        name = None if row is None else tw.item(row,0).text()
        self.current_solid_index = (row+1) if row is not None else None
        self.update_solids_detail_pane(name)

    def update_solids_detail_pane(self, name):
        ui = self.ui
        sa = ui.scrollarea_solids_detail
        if name is None:
            sa.setEnabled(False)
            # Clear out all values?
        else:
            data = self.solids[name]
            sa.setEnabled(True)
            ui.lineedit_solids_name.setText(name)
            self.setup_combobox_solids_model()

    def update_solids_table(self):
        return
        hv = QtWidgets.QHeaderView
        table = self.ui.tablewidget_solids
        if PYQT5:
            resize = table.horizontalHeader().setSectionResizeMode
        else:
            resize = table.horizontalHeader().setResizeMode
        for n in range(4):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)

        table.clearContents()
        if self.solids is None:
            return
        nrows = len(self.solids)
        table.setRowCount(nrows)
        def make_item(val):
            item = QtWidgets.QTableWidgetItem(str(val))
            set_item_noedit(item)
            return item
        for (row,(k,v)) in enumerate(self.solids.items()):
            table.setItem(row, 0, make_item(k))
            for (col, key) in enumerate(('model', 'diameter', 'density'), 1):
                table.setItem(row, col, make_item(v[key]))

    def solids_delete(self):
        return
        tw = self.ui.tablewidget_solids
        row = get_selected_row(tw)
        if row is None: # No selection
            return
        name = tw.item(row, 0).text()
        del self.solids[name]
        tw.removeRow(row)
        tw.clearSelection()

    def solids_change_tab(self, tabnum, btn):
        """ switch mesh stacked widget based on selected """
        self.animate_stacked_widget(
            self.ui.solids.stackedwidget_solids,
            self.ui.solids.stackedwidget_solids.currentIndex(),
            tabnum,
            direction='horizontal',
            line=self.ui.solids.line_solids,
            to_btn=btn,
            btn_layout=self.ui.solids.gridlayout_solid_tab_btns,
            )



def Usage(name):
    print("""Usage: %s [directory|file] [-h, --help] [-l, --log=LEVEL] [-q, --quit]
    directory: open mfix.dat file in specified directory
    file: open mfix.dat or <RUN_NAME>.mfx project file
    -h, --help: display this help message
    -l, --log=LEVEL: set logging level (error,warning,info,debug)
    -n, --new:  open new project (do not autoload previous)
    -q, --quit: quit after opening file (for testing)"""  % name, file=sys.stderr)
    sys.exit(1)

def main(args):
    """Handle command line options and start the GUI"""
    name = args[0]
    try:
        opts, args = getopt.getopt(args[1:], "hqnl:", ["help", "quit", "new", "log="])
    except getopt.GetoptError as err:
        print(err)
        Usage(name)

    quit_after_loading = False
    project_file = None
    new_project = False
    log_level = 'WARN'

    for opt, arg in opts:
        if opt in ("-l", "--log"):
            log_level = arg
        elif opt in ("-h", "--help"):
            Usage(name)
            sys.exit(0)
        elif opt in ("-q", "--quit"):
            quit_after_loading = True
        elif opt in ("-n", "--new"):
            new_project = True
        else:
            Usage(name)

    logging.basicConfig(stream=sys.stdout,
                        filemode='w', level=getattr(logging, log_level.upper()),
                        format='%(name)s - %(levelname)s - %(message)s')

    if len(args) > 1:
        Usage(name)
    if args:
        project_file = args[0]

    qapp = QtWidgets.QApplication([])
    mfix = MfixGui(qapp, project_file=project_file)
    mfix.show()

    # --- print welcome message
    #mfix.print_internal("MFiX-GUI version %s" % mfix.get_version())

    if project_file is None and not new_project:
        # autoload last project
        project_file = mfix.get_project_file()

    if project_file:
        mfix.open_project(project_file, auto_rename=(not quit_after_loading))
    else:
        # disable all widgets except New and Open
        mfix.set_solver(SINGLE)
        mfix.ui.stackedwidget_mode.setVisible(False)
        mfix.ui.toolbutton_export.setEnabled(False)
        mfix.ui.toolbutton_run.setEnabled(False)
        mfix.ui.toolbutton_save.setEnabled(False)
        mfix.ui.toolbutton_save_as.setEnabled(False)
        # This gets set by guess_solver if we're loading a project, otherwise
        # we need to set the default.  (Do other defaults need to be set here?)

    # print number of keywords
    mfix.print_internal('Registered %d keywords' %
                        len(mfix.project.registered_keywords))

    # have to initialize vtk after the widget is visible!
    if mfix.use_vtk:
        mfix.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    if not quit_after_loading:
        qapp.exec_()

    qapp.deleteLater()
    sys.exit()

if __name__  == '__main__':
    main(sys.argv)
