""" MFIX GUI """


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals
import glob
import logging
import os
import shutil
import signal
import subprocess
import sys
import time
from collections import OrderedDict
import copy
import warnings

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

# import qt
from qtpy import QtCore, QtWidgets, QtGui
from qtpy.QtCore import (QObject, QThread, pyqtSignal, QUrl, QTimer, QSettings,
                         Qt)

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16
# TODO: cache ui file to a py file, use the ui file if newer, else py file
try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic


# Debugging hooks
def debug_trace():
  '''Set a tracepoint in the Python debugger that works with Qt'''
  from qtpy.QtCore import pyqtRemoveInputHook
  from pdb import set_trace
  pyqtRemoveInputHook()
  set_trace()


# local imports
from widgets.vtkwidget import VtkWidget
from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table)
from widgets.regions import RegionsWidget
from tools.mfixproject import Project, Keyword
from tools.general import (make_callback, get_icon, get_unique_string,
                           widget_iter, set_script_directory, CellColor)
from tools.namelistparser import buildKeywordDoc

# look for pyqtnode
try:
    from pyqtnode import NodeWidget
except ImportError:
    NodeWidget = None

logging.basicConfig(stream=sys.stdout,
                    filemode='w', level=logging.INFO,
                    format='%(name)s - %(levelname)s - %(message)s')

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
LOG = logging.getLogger(__name__)
LOG.debug(SCRIPT_DIRECTORY)
set_script_directory(SCRIPT_DIRECTORY)  # should this be in an __init__.py?

def get_mfix_home():
    " return the top level directory (since __file__ is mfix_home/gui/gui.py) "
    return os.path.dirname(
        os.path.dirname(os.path.realpath(__file__)))

def format_key_with_args(key, args=None):
    if args:
        return "%s(%s)" % (key, ','.join(str(a) for a in args))
    else:
        return str(key)

def plural(n, word):
    fmt = "%d %s" if n==1 else "%d %ss"
    return fmt % (n, word)

# Constants

# Solver types
# must match combobox_solver in model_setup.ui
SINGLE, TFM, DEM, PIC, HYBRID = range(5)
CONSTANT, AIR, UDF = 0, 1, 2

# --- Main Gui ---
class MfixGui(QtWidgets.QMainWindow):
    '''
    Main window class handling all gui interactions
    '''
    def __init__(self, app, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)

        # reference to qapp instance
        self.app = app

        # load ui file
        self.customWidgets = {'LineEdit':      LineEdit,
                              'CheckBox':      CheckBox,
                              'ComboBox':      ComboBox,
                              'DoubleSpinBox': DoubleSpinBox,
                              'SpinBox':       SpinBox,
                              'Table':         Table,
                              }

        self.ui = uic.loadUi(os.path.join('uifiles', 'gui.ui'), self)

        self.ui.general = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'general.ui'), self.ui.general)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.general)

        self.ui.geometry = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'geometry.ui'), self.ui.geometry)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.geometry)

        self.ui.mesh = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'mesh.ui'), self.ui.mesh)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.mesh)

        self.ui.regions = RegionsWidget()
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.regions)

        self.ui.model_setup = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'model_setup.ui'), self.ui.model_setup)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.model_setup)

        self.ui.numerics = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'numerics.ui'), self.ui.numerics)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.numerics)

        self.ui.output = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'output.ui'), self.ui.output)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.output)

        self.ui.monitors = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'monitors.ui'), self.ui.monitors)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.monitors)

        self.ui.run = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'run.ui'), self.ui.run)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.run)

        self.ui.interact = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'interact.ui'), self.ui.interact)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.interact)

        self.ui.post_processing = QtWidgets.QWidget()
        uic.loadUi(os.path.join('uifiles', 'post_processing.ui'), self.ui.post_processing)
        self.ui.stackedWidgetTaskPane.addWidget(self.ui.post_processing)

        # load settings
        self.settings = QSettings('MFIX', 'MFIX')

        # set title and icon
        self.setWindowTitle('MFIX')
        self.setWindowIcon(get_icon('mfix.png'))

        # build keyword documentation from namelist docstrings
        self.keyword_doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY,
                                                        os.pardir, 'model'))

        # create project manager
        self.project = ProjectManager(self, self.keyword_doc)

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
        self.solid_nscalar_eq = 0 # Infer these from phase4scalar?
        # Defaults
        self.solver = SINGLE
        self.fluid_density_model = CONSTANT
        self.fluid_viscosity_model = CONSTANT
        self.fluid_molecular_weight_model = CONSTANT
        self.fluid_specific_heat_model = CONSTANT
        self.fluid_conductivity_model = AIR
        self.fluid_diffusion_model = AIR

        # --- icons ---
        # loop through all widgets, because I am lazy
        for widget in widget_iter(self.ui):
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

        self.ui.toolbutton_compile.setIcon(get_icon('build.png'))
        self.ui.toolbutton_run.setIcon(get_icon('play.png'))
        self.ui.toolbutton_restart.setIcon(get_icon('restart.png'))
        self.ui.toolbutton_interact.setIcon(get_icon('flash.png'))

        self.ui.geometry.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        self.ui.geometry.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        self.ui.geometry.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        self.ui.geometry.toolbutton_geometry_intersect.setIcon(
            get_icon('intersect.png'))
        self.ui.geometry.toolbutton_geometry_difference.setIcon(
            get_icon('difference.png'))

        self.ui.toolButtonTFMSolidsDatabase.setIcon(get_icon('download.png'))

        # --- Connect Signals to Slots---
        # open/save/new project
        self.ui.toolbutton_open.clicked.connect(self.handle_open_action)
        self.ui.toolbutton_save.clicked.connect(self.save_project)

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.clicked.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        self.ui.treewidget_model_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # build/run/connect MFIX
        self.ui.run.run_mfix_button.clicked.connect(self.run_mfix)
        self.ui.run.connect_mfix_button.clicked.connect(self.connect_mfix)
        self.ui.run.clear_output_button.clicked.connect(self.clear_output)
        self.ui.run.mfix_executables.activated.connect(self.update_run)

        # --- Threads ---
        self.run_thread = MfixThread(self)
        self.clear_thread = MfixThread(self)
        self.monitor_executables_thread = MonitorExecutablesThread(self)

        def make_handler(qtextbrowser):
            " make a closure to read output from external process "

            def handle_line(line, color=None): # Combine with print_internal
                " closure to read output from external process "
                log = logging.getLogger(__name__)
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

        def update_run_executables(*args):
            "Updates the list of executables available"
            self.ui.run.mfix_executables.clear()
            output = self.monitor_executables_thread.get_output()
            for executable in output:
                self.ui.run.mfix_executables.addItem(executable)
            mfix_available = bool(output)
            self.ui.run.mfix_executables.setVisible(mfix_available)
            self.ui.run.mfix_executables_warning.setVisible(not mfix_available)

        self.run_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        self.clear_thread.line_printed.connect(
            make_handler(self.ui.command_output))
        update_run_executables()

        self.monitor_executables_thread.sig.connect(update_run_executables)
        self.monitor_executables_thread.start()

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

        # --- print welcome message
        self.print_internal("MFiX-GUI version %s" % self.get_version())

    def get_version(self):
        return "0.2x" # placeholder

    def set_navigation_item_state(self, item_name, state):
        on = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        off = Qt.ItemFlags(0)
        tree = self.ui.treewidget_model_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        items = tree.findItems(item_name, flags, 0)
        assert len(items) == 1, "Multiple menu items matching %s"%item_name
        items[0].setFlags(on if state else off)

    def set_solver(self, index):
        """ handler for "Solver" combobox in Model Setup """
        self.solver = index

        model_setup = self.ui.model_setup
        solver_name = model_setup.combobox_solver.currentText()
        self.print_internal("set solver to %d: %s" % (index, solver_name))

        item_names =  ("Solids", "Continuum Solids Model",
                       "Discrete Element Model", "Particle in Cell Model")

        item_states = {SINGLE: (False, False, False, False),
                       TFM: (True, True, False, False),
                       DEM: (True, False, True, False),
                       PIC: (True, False, False, True),
                       HYBRID: (True, True, True, False)}

        for item_name, item_state in zip(item_names, item_states[index]):
            self.set_navigation_item_state(item_name, item_state)

        # Options which require TFM, DEM, or PIC
        enabled = self.solver in (TFM, DEM, PIC)
        interphase = model_setup.groupbox_interphase
        interphase.setEnabled(enabled)

        # TFM only
        # use a groupbox here, instead of accessing combobox + label?
        enabled = (self.solver == TFM)
        model_setup.combobox_subgrid_model.setEnabled(enabled)
        model_setup.label_subgrid_model.setEnabled(enabled)
        model_setup.groupbox_subgrid_params.setEnabled(enabled and
                                                       self.subgrid_model > 0)

        self.ui.checkbox_enable_fluid_scalar_eq.setEnabled(enabled)
        self.ui.spinbox_fluid_nscalar_eq.setEnabled(enabled
            and self.ui.checkbox_enable_fluid_scalar_eq.isChecked())

    def set_keyword(self, key, value, args=None):
        self.project.submit_change(None, {key:value}, args)

    def unset_keyword(self, key, args=None):
        if self.project.removeKeyword(key, args, warn=False):
            self.print_internal("unset %s" %
                                format_key_with_args(key, args))

    def enable_energy_eq(self, state):
        # Additional callback on top of automatic keyword update,
        # since this has to change availabilty of a bunch of other GUI items
        ui = self.ui
        for item in (ui.combobox_fluid_specific_heat_model,
                     ui.combobox_fluid_conductivity_model,
                     # more ?
                     ):
            item.setEnabled(state)


        spinbox = ui.spinbox_keyword_cp_g0 # cp_g0 == specific heat for fluid phase
        if state:
            spinbox.setEnabled(self.fluid_specific_heat_model == CONSTANT)
        else:
            spinbox.setEnabled(False)

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
        self.ui.spinbox_fluid_nscalar_eq.setEnabled(state)
        if state:
            self.set_fluid_nscalar_eq(self.fluid_nscalar_eq)

    def set_fluid_nscalar_eq(self, value):
        # TODO:  load from mfix.dat (do we save this explicitly,
        # or infer from PHASE4SCALAR settings?
        self.fluid_nscalar_eq = value
        self.project.submit_change(None,{"nscalar":
                                         self.fluid_nscalar_eq + self.solid_nscalar_eq})
        for i in range(1,1+value):
            self.project.submit_change(None,{"phase4scalar":0},args=i)

    # note, the set_fluid_*_model methods have a lot of repeated code

    def set_fluid_density_model(self, value):
        self.fluid_density_model = value

        # Enable spinbox for constant density model
        spinbox = self.ui.spinbox_keyword_ro_g0
        spinbox.setEnabled(value==0)
        if value == CONSTANT:
            self.set_keyword("ro_g0", spinbox.value())
            self.unset_keyword("usr_rog")
        elif value == 1: # Ideal Gas Law
            self.unset_keyword("ro_g0")
            self.unset_keyword("usr_rog")
        elif value == UDF:
            self.unset_keyword("ro_g0")
            self.set_keyword("usr_rog", True)

    def set_fluid_viscosity_model(self, value):
        self.fluid_viscosity_model = value

        # Enable spinbox for constant viscosity model
        spinbox = self.ui.spinbox_keyword_mu_g0
        spinbox.setEnabled(value==0)
        if value == CONSTANT:
            self.set_keyword("mu_g0", spinbox.value())
            self.unset_keyword("usr_mug")
        elif value == 1: # Sutherland's Law
            self.unset_keyword("mu_g0")
            self.unset_keyword("usr_mug")
        elif value == UDF:
            self.unset_keyword("mu_g0")
            self.set_keyword("usr_mug", True)

    def set_fluid_molecular_weight_model(self, value):
        self.fluid_molecular_weight_model = value

        # Enable spinbox for constant molecular_weight model
        spinbox = self.ui.spinbox_keyword_mw_avg
        spinbox.setEnabled(value==0)
        if value == CONSTANT:
            self.set_keyword("mw_avg", spinbox.value())
        elif value == 1: # Mixture
            # TODO: require mw for all component species
            self.unset_keyword("mw_avg")

    def set_fluid_specific_heat_model(self, value):
        self.fluid_specific_heat_model = value

        # Enable spinbox for constant specific_heat model
        spinbox = self.ui.spinbox_keyword_cp_g0
        spinbox.setEnabled(value==0)
        if value == CONSTANT:
            self.set_keyword("cp_g0", spinbox.value())
            self.unset_keyword("usr_cpg")
        elif value == 1: # Mixture
            # TODO: require cp for all component species
            self.unset_keyword("cp_g0")
            self.unset_keyword("usr_mug")
        elif value == UDF:
            self.unset_keyword("cp_g0")
            self.set_keyword("usr_cpg", True)

    def set_fluid_conductivity_model(self, value):
        self.fluid_conductivity_model = value

        # Enable spinbox for constant (thermal) conductivity model
        spinbox = self.ui.spinbox_keyword_k_g0
        spinbox.setEnabled(value==0)
        if value == CONSTANT:
            self.set_keyword("k_g0", spinbox.value())
            self.unset_keyword("usr_kg")
        elif value == 1: # Temperature dep.
            # TODO: require cp for all component species
            self.unset_keyword("k_g0")
            self.unset_keyword("usr_kg")
        elif value == UDF:
            self.unset_keyword("k_g0")
            self.set_keyword("usr_kg", True)

    def set_fluid_diffusion_model(self, value):
        self.fluid_diffusion_model = value

        # Enable spinbox for constant diffusion model
        spinbox = self.ui.spinbox_keyword_dif_g0
        spinbox.setEnabled(value==CONSTANT)

        if value == CONSTANT:
            self.set_keyword("dif_g0", spinbox.value())
            self.unset_keyword("usr_difg")
        elif value == AIR: # Temperature dep.
            # TODO: require temperature field for full domain
            self.unset_keyword("dif_g0")
            self.unset_keyword("usr_difg")
        elif value == UDF:
            self.unset_keyword("dif_g0")
            self.set_keyword("usr_difg", True)

    def disable_fluid_solver(self, state):
        self.set_navigation_item_state("Fluid", not state)

    def set_subgrid_model(self, index):
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)

    def __setup_other_widgets(self): # rename/refactor
        """ setup widgets which are not tied to a simple keyword """

        model_setup = self.ui.model_setup
        combobox = model_setup.combobox_solver
        combobox.currentIndexChanged.connect(self.set_solver)
        self.set_solver(self.solver) # Default - Single Phase - (?)

        checkbox = model_setup.checkbox_disable_fluid_solver
        checkbox.stateChanged.connect(self.disable_fluid_solver)
        self.disable_fluid_solver(False)

        checkbox = model_setup.checkbox_keyword_energy_eq
        checkbox.stateChanged.connect(self.enable_energy_eq)

        checkbox = self.ui.checkbox_keyword_species_eq_args_0
        checkbox.stateChanged.connect(self.enable_fluid_species_eq)

        combobox = model_setup.combobox_subgrid_model
        combobox.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)

        self.enable_energy_eq(False)
        # something is wrong, we shouldn't have to explicitly uncheck this box
        self.model_setup.checkbox_keyword_energy_eq.setChecked(False)


        # Fluid phase
        self.ui.checkbox_enable_fluid_scalar_eq.stateChanged.connect(
            self.enable_fluid_scalar_eq)
        self.ui.spinbox_fluid_nscalar_eq.valueChanged.connect(
            self.set_fluid_nscalar_eq)

        # Fluid phase models
        # Density
        self.ui.combobox_fluid_density_model.currentIndexChanged.connect(
            self.set_fluid_density_model)
        self.set_fluid_density_model(self.fluid_density_model)
        # Viscosity
        self.ui.combobox_fluid_viscosity_model.currentIndexChanged.connect(
            self.set_fluid_viscosity_model)
        self.set_fluid_viscosity_model(self.fluid_viscosity_model)
        # Molecular Weight
        self.ui.combobox_fluid_molecular_weight_model.currentIndexChanged.connect(
            self.set_fluid_molecular_weight_model)
        self.set_fluid_molecular_weight_model(self.fluid_molecular_weight_model)
        # Specific Heat
        self.ui.combobox_fluid_specific_heat_model.currentIndexChanged.connect(
            self.set_fluid_specific_heat_model)
        self.set_fluid_specific_heat_model(self.fluid_specific_heat_model)
        # (Thermal) Conductivity
        self.ui.combobox_fluid_conductivity_model.currentIndexChanged.connect(
            self.set_fluid_conductivity_model)
        self.set_fluid_conductivity_model(self.fluid_conductivity_model)
        # Diffusion (Coefficient)
        self.ui.combobox_fluid_diffusion_model.currentIndexChanged.connect(
            self.set_fluid_diffusion_model)
        self.set_fluid_diffusion_model(self.fluid_diffusion_model)

        # Fluid species
        # Automate the connecting?
        toolbutton = self.ui.toolbutton_fluid_species_add
        toolbutton.clicked.connect(self.fluid_species_add)

        # lineedit = self.ui.lineedit_fluid_species_name
        # class FluidSpeciesNameValidator(QtGui.QValidator):
        #     def validate(self, text, pos):
        #         if text == "":
        #             return (QtGui.QValidator.Intermediate, text, pos)
        #         elif not isalpha(text[0]):
        #             return (QtGui.QValidator.Invalid, text, pos)
        #         # additional constraints, ie no spaces, etc? - cgw

        #         return (QtGui.QValidator.Acceptable, text, pos)

        #     # see http://stackoverflow.com/questions/26759623/why-is-the-return-of-qtgui-qvalidator-validate-so-inconsistent-robust-way-to
        # lineedit.setValidator(FluidSpeciesNameValidator())



    def __setup_simple_keyword_widgets(self):
        """
        Look for and connect simple keyword widgets to the project manager.
        Keyword information from the namelist doc strings is added to each
        keyword widget. The widget must be named: *_keyword_<keyword> where
        <keyword> is the actual keyword. For example:
        lineedit_keyword_run_name
        """

        def try_int(str):
            try:
                return int(str)
            except ValueError:
                return str

        # loop through all widgets looking for *_keyword_<keyword>
        for widget in widget_iter(self.ui):
            name_list = str(widget.objectName()).lower().split('_')

            if 'keyword' in name_list:
                keyword_idx = name_list.index('keyword')
                args = None
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
        " setup the vtk widget "

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

    def get_project_path(self):
        " get the current project directory "

        last = self.settings.value('project_path')
        if last:
            return last
        else:
            return None

    def mode_changed(self, mode):
        " change the Modeler, Workflow, Developer tab"

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
        """ change to the specified pane """

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
        """ an item in the tree was selected, change panes """
        current_selection = self.ui.treewidget_model_navigation.selectedItems()

        if current_selection:

            text = str(current_selection[-1].text(0))
            text = '_'.join(text.lower().split(' '))

            current_index = 0
            for i in range(self.ui.stackedWidgetTaskPane.count()):
                widget = self.ui.stackedWidgetTaskPane.widget(i)
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
        """ animate changing of qstackedwidget """

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

        # animate
        # from widget
        animnow = QtCore.QPropertyAnimation(from_widget, "pos")
        animnow.setDuration(self.animation_speed)
        animnow.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnow.setStartValue(
            QtCore.QPoint(0,
                          0))
        animnow.setEndValue(
            QtCore.QPoint(0 - offsetx,
                          0 - offsety))

        # to widget
        animnext = QtCore.QPropertyAnimation(to_widget, "pos")
        animnext.setDuration(self.animation_speed)
        animnext.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animnext.setStartValue(
            QtCore.QPoint(0 + offsetx,
                          0 + offsety))
        animnext.setEndValue(
            QtCore.QPoint(0,
                          0))

        # line
        animline = None
        if line is not None and to_btn is not None:
            animline = QtCore.QPropertyAnimation(line, "pos")
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
                          stackedwidget, from_, to, btn_layout, line)
            )

        self.animating = True
        self.stack_animation.start()

    def animate_stacked_widget_finished(self, widget, from_, to,
                                        btn_layout=None, line=None):
        """ cleanup after animation """
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
        '''
        Create a message box:
        title = 'title'
        icon = 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel',
            'discard'
        default = 'ok' the default selected button
        infotext = 'extended information text'
        detailedtext = 'Some details'
        '''
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
        LOG.info(str(line).strip())
        cursor = qtextbrowser.textCursor()
        cursor.movePosition(cursor.End)
        char_format = QtGui.QTextCharFormat()
        if color:
            char_format.setForeground(QtGui.QColor(color))
        if font:
            char_format.setFontFamily(font)
        cursor.setCharFormat(char_format)
        cursor.insertText(line)
        scrollbar = qtextbrowser.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def update_run(self):
        """Enable/disable run options based on selected executable """
        mfix_exe = self.ui.run.mfix_executables.currentText()
        if not mfix_exe:
            return
        config = self.monitor_executables_thread.output[mfix_exe]
        smp_enabled = 'smp' in config
        dmp_enabled = 'dmp' in config
        pymfix_enabled = 'pymfix' in mfix_exe
        self.ui.run.openmp_threads.setEnabled(smp_enabled)
        if not dmp_enabled:
            self.ui.run.spinbox_keyword_nodesi.setValue(1)
            self.ui.run.spinbox_keyword_nodesj.setValue(1)
            self.ui.run.spinbox_keyword_nodesk.setValue(1)
        self.ui.run.spinbox_keyword_nodesi.setEnabled(dmp_enabled)
        self.ui.run.spinbox_keyword_nodesj.setEnabled(dmp_enabled)
        self.ui.run.spinbox_keyword_nodesk.setEnabled(dmp_enabled)
        self.ui.run.mfix_host.setEnabled(pymfix_enabled)
        self.ui.run.mfix_port.setEnabled(pymfix_enabled)
        self.ui.run.connect_mfix_button.setEnabled(pymfix_enabled)

    def clear_output(self):
        self.run_thread.start_command(
            cmd='rm -v -f *.LOG *.OUT *.RES *.SP? *.pvd *vtp',
            cwd=self.get_project_path())

    def run_mfix(self):
        mfix_exe = self.ui.run.mfix_executables.currentText()
        config = self.monitor_executables_thread.output[mfix_exe]

        if mfix_exe.endswith('pymfix'):
            # run with python or python3, depending on sys.executable
            mfix_exe = '{} {}'.format(sys.executable, mfix_exe)

        if 'dmp' in config:
            nodesi = self.project.nodesi
            nodesj = self.project.nodesj
            nodesk = self.project.nodesk
            total = nodesi * nodesj * nodesk
            # FIXME: maybe we should save NODES* keywords to runname.mfx instead of passing them on command line?
            mfix_exe = 'mpirun -np {} {} NODESI={} NODESJ={} NODESK={}'.format(
                total, mfix_exe, nodesi, nodesj, nodesk)

        run_cmd = '{} -f {}'.format(mfix_exe, self.get_mfix_dat())
        project_dir = os.path.dirname(self.get_project_path())
        self.run_thread.start_command(cmd=run_cmd, cwd=project_dir)

    def update_residuals(self):
        self.ui.residuals.setText(str(self.update_residuals_thread.residuals))
        if self.update_residuals_thread.job_done:
            self.ui.mfix_browser.setHTML('')

    def connect_mfix(self):
        """ connect to running instance of mfix """
        url = "http://{}:{}".format(self.ui.run.mfix_host.text(),
                                    self.ui.run.mfix_port.text())
        log = logging.getLogger(__name__)
        log.debug("trying to connect to {}".format(url))
        qurl = QUrl(url)
        self.ui.interact.mfix_browser.load(qurl)

        self.update_residuals_thread = UpdateResidualsThread()
        self.update_residuals_thread.sig.connect(self.update_residuals)
        self.update_residuals_thread.start()

        self.change_pane('interact')

    # --- open/save/new ---
    def save_project(self):
        project_path = self.settings.value('project_path')
        project_dir = os.path.dirname(project_path)
        # export geometry
        self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))

        self.project.writeDatFile(project_path)
        with open(project_path, 'r') as mfx:
            src = mfx.read()
        self.ui.mfix_dat_source.setPlainText(src)
        self.setWindowTitle('MFIX - %s' % project_path)

    def unsaved(self):
        project_path = self.settings.value('project_path')
        self.setWindowTitle('MFIX - %s *' % project_path)

    def check_writable(self, directory):
        " check whether directory is writable "
        try:
            import tempfile
            testfile = tempfile.TemporaryFile(dir=directory)
            testfile.close()
            return True
        except IOError:
            writable = False

            self.message(title='Warning',
                         icon='warning',
                         text=('You do not have write access to this '
                               'directory.\n'),
                         buttons=['ok'],
                         default='ok',
                         )
            return False

    def new_project(self, project_path=None):
        if not project_path:
            project_path = str(
                QtWidgets.QFileDialog.getOpenFileName(
                    self, 'Create Project in Directory',
                    ""))
        if len(project_path) < 1:
            return

        if os.path.isdir(project_path):
            project_dir = project_path
            project_path = os.path.join(project_dir, 'mfix.dat')
        else:
            project_dir = os.path.dirname(project_path)

        if not self.check_writable(project_dir):
            return

        shutil.copyfile('mfix.dat.template',project_path)

        self.open_project(project_path)

    def handle_open_action(self):
        project_dir = os.path.dirname(self.get_project_path())
        project_path = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Project Directory', project_dir)

        if len(project_path) < 1:
            return

        project_path = project_path[0]

        if not os.path.exists(project_path):
            return

        self.open_project(project_path)

    def open_project(self, project_path):
        """
        Open MFiX Project
        """

        if os.path.isdir(project_path):
            project_dir = project_path
            project_path = os.path.abspath(os.path.join(project_path, 'mfix.dat'))
        else:
            project_dir = os.path.dirname(project_path)

        if not os.path.exists(project_path):
            return

        self.print_internal("Loading %s" % project_path, color='blue')
        self.project.load_mfix_dat(project_path)
        self.print_internal("Loaded %s" % project_path, color='blue')

        if hasattr(self.project, 'run_name'):
            name = self.project.run_name.value
        else:
            name = 'new_file'
        for char in ('.', '"', "'", '/', '\\', ':'):
            name = name.replace(char, '_')
        runname_mfx = name + '.mfx'

        if not project_path.endswith(runname_mfx):
            self.message(title='Warning',
                         icon='warning',
                         text=('Saving %s as %s based on run name\n' % (project_path, runname_mfx)),
                         buttons=['ok'],
                         default='ok',
            )
            project_path = os.path.join(project_dir, runname_mfx)
            self.project.writeDatFile(project_path)

        # project_dir = os.path.dirname(project_path)
        self.settings.setValue('project_path', project_path)
        self.setWindowTitle('MFIX - %s' % project_path)

        # read the file
        with open(project_path, 'r') as mfx:
            src = mfx.read()
        self.ui.mfix_dat_source.setPlainText(src)
        # self.mode_changed('developer')

        # Additional GUI setup based on loaded projects (not handled
        # by keyword updates)
        self.enable_energy_eq(self.project['energy_eq'])
        # cgw - lots more model setup todo here

    # --- fluid species methods ---
    def fluid_species_add(self):
        pass
        # table =  self.ui.tablewidget_fluid_species
        # n = table.rowCount()
        # table.setRowCount(n+1)
        # for i in range(table.columnCount()):
        #     item = QtWidgets.QTableWidgetItem("New" if i==0 else "")
        #     item.setFlags(item.flags() ^ Qt.ItemIsEditable)
        #     table.setItem(n, i, item)
        # table.setCurrentCell(n, 0)
        # lineedit =  self.ui.lineedit_fluid_species_name
        # lineedit.clear()
        # lineedit.insert("New")
        # lineedit.selectAll()
        # lineedit.setFocus()


# --- Threads ---
class MfixThread(QThread):
    line_printed = pyqtSignal(str, str)

    def __init__(self, parent):
        super(MfixThread, self).__init__(parent)
        self.cmd = None
        self.cwd = None

    def start_command(self, cmd, cwd):
        log = logging.getLogger(__name__)
        log.debug("Running in %s : %s", cwd, cmd)
        self.cmd = cmd
        self.cwd = cwd
        self.start()

    def run(self):
        if self.cmd:
            popen = subprocess.Popen(self.cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True, cwd=self.cwd)
            lines_iterator = iter(popen.stdout.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(str(line), None)
            lines_iterator = iter(popen.stderr.readline, b"")
            for line in lines_iterator:
                self.line_printed.emit(str(line), "red")


class MonitorExecutablesThread(QThread):

    sig = QtCore.pyqtSignal()

    def __init__(self, parent):
        QThread.__init__(self)
        self.parent = parent
        self.mfix_home = get_mfix_home()
        self.cache = {}
        self.output = self.get_output()

    def get_output(self):
        """ returns a dict mapping full [mfix|pymfix] paths
        to configuration options."""

        def mfix_print_flags(mfix_exe, cache=self.cache):
            """ Determine mfix configuration by running mfix --print-flags.  Cache results"""
            try: # Possible race, file may have been deleted/renamed since isfile check!
                stat = os.stat(mfix_exe)
            except OSError:
                return ''

            cached_stat, cached_flags = cache.get(mfix_exe, (None, None))
            if cached_stat and cached_stat == stat:
                return cached_flags

            popen = subprocess.Popen(mfix_exe + " --print-flags",
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True)
            out = popen.stdout.read()
            err = popen.stderr.read()
            flags = '' if err else out
            cache[mfix_exe] = (stat, flags)
            return flags

        config_options = {}

        # Check system PATH dirs first
        PATH = os.environ.get("PATH")
        if PATH:
            dirs =  PATH.split(os.pathsep)
        else:
            dirs = []

        # Look in subdirs of build dir
        build_dir = os.path.join(self.mfix_home,'build')
        if os.path.exists(build_dir):
            dirs.extend(os.path.join(build_dir, subdir)
                         for subdir in os.listdir(build_dir))
        # Check run_dir
        project_path = self.parent.get_project_path()
        if project_path:
            run_dir = os.path.dirname(project_path)
            if run_dir:
                dirs.append(run_dir)

        # Now look for mfix/pymfix in these dirs
        for dir in dirs:
            for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                exe = os.path.join(dir, name)
                if os.path.isfile(exe):
                    config_options[exe] = mfix_print_flags(exe)

        return config_options

    def run(self):
        while True:
            output = self.get_output()
            if output != self.output:
                self.output = output
                self.sig.emit()
            time.sleep(1)

class MonitorOutputFilesThread(QThread):

    sig = QtCore.pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            output = glob.glob(run_name+'*')
            if output:
                self.sig.emit('output')
            else:
                self.sig.emit('no_output')
            time.sleep(1)


class UpdateResidualsThread(QThread):

    sig = QtCore.pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            try:
                self.residuals = urlopen('http://localhost:5000/residuals').read()
            except Exception:
                log = logging.getLogger(__name__)
                log.debug("cannot retrieve residuals; pymfix process must have terminated.")
                self.job_done = True
                return
            time.sleep(1)
            self.sig.emit('update')


# --- Project Manager ---
class ProjectManager(Project):
    '''
    Manages the editing of the MFiX.dat file
    '''
    def __init__(self, parent=None, keyword_doc=None):
        Project.__init__(self, keyword_doc=keyword_doc)

        self.parent = parent
        self.keyword_and_args_to_widget = {}
        self.registered_keywords = set()
        self._widget_update_stack = [] # prevent circular updates

    def submit_change(self, widget, newValueDict, args=None,
                      forceUpdate=False): # forceUpdate unused (?)
        '''
        Submit a value change

        Examples:
        submitChange(lineEdit, {'run_name':'new run name'}, args)
        '''

        # Note, this may be a callback from Qt (in which case 'widget' is
        # the widget that the user activated), or from the initial mfix
        # loading (in which case widget == None)

        if isinstance(args, int):
            args = [args]

        for key, newValue in newValueDict.items():
            if isinstance(newValue, dict):
                if args:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=args+[ind],
                                     forceUpdate=forceUpdate)
                else:
                    for ind, value in newValue.items():
                        self._change(widget, key, value, args=[ind],
                                     forceUpdate=forceUpdate)
            else:
                self._change(widget, key, newValue, args=args,
                             forceUpdate=forceUpdate)

        self._cleanDeletedItems()


    def _change(self, widget, key, newValue, args=None, forceUpdate=False):
        # prevent circular updates, from this widget or any higher in stack
        if widget in self._widget_update_stack:
            return

        key = key.lower()
        updatedValue = None
        if isinstance(newValue, Keyword): # why needed?
            keyword = newValue
            newValue = keyword.value
        else:
            keyword = None
        if args is None:
            args = []

        updatedValue = self.updateKeyword(key, newValue, args)

        keytuple = tuple([key]+args)
        widgets_to_update = self.keyword_and_args_to_widget.get(keytuple)
        warn = False
        if widgets_to_update == None:
            warn = True
            widgets_to_update = []

        # Are we using the 'all' mechanism?
        #widgets_to_update.extend(
        #    self.keyword_and_args_to_widget.get("all", []))

        for w in widgets_to_update:
            # Prevent circular updates
            self._widget_update_stack.append(w)
            try:
                w.updateValue(key, updatedValue, args)
            except Exception, e:
                raise ValueError("%s = %s" % (format_key_with_args(key, args),
                                              updatedValue))
            finally:
                self._widget_update_stack.pop()

        self.parent.print_internal("%s = %s" % (format_key_with_args(key, args),
                                                updatedValue),
                                   font="Monospace", color='red' if warn else None)


    def load_mfix_dat(self, mfixDat):
        """Load an MFiX project file."""
        n_errs = 0
        errlist = []
        with warnings.catch_warnings(record=True) as ws:
            self.parsemfixdat(fname=mfixDat)
            # emit loaded keys
            # some of these changes may cause new keywords to be instantiated,
            # so iterate over a copy of the list, which may change
            kwlist = list(self.keywordItems())
            for keyword in kwlist:
                try:
                    self.submit_change(None, {keyword.key: keyword.value},
                                       args=keyword.args, forceUpdate=True)
                except ValueError, e:
                    errlist.append(e)

            # report any errors
            for w in errlist + ws:
                self.parent.print_internal("Warning: %s" % w.message, color='red')
            n_errs = len(errlist) + len(ws)
            if n_errs:
                self.parent.print_internal("Warning: %s" % plural(n_errs, "error"), color='red')
            else:
                self.parent.print_internal("0 errors", color='darkgreen')

    def register_widget(self, widget, keys=None, args=None):
        '''
        Register a widget with the project manager. The widget must have a
        value_updated signal to connect to.
        '''
        if args is None:
            args = []

        LOG.debug('ProjectManager: Registering {} with keys {}, args {}'.format(
            widget.objectName(),
            keys, args))

        # add widget to dictionary of widgets to update
        d = self.keyword_and_args_to_widget
        for key in keys:
            keytuple = tuple([key]+args)
            if keytuple not in d:
                d[keytuple] = []
            d[keytuple].append(widget)

        widget.value_updated.connect(self.submit_change)

        self.registered_keywords = self.registered_keywords.union(set(keys))

    def objectName(self):
        return 'Project Manager'


if __name__ == '__main__':
    qapp = QtWidgets.QApplication(sys.argv)
    mfix = MfixGui(qapp)
    mfix.show()
    # --- print welcome message
    mfix.print_internal("MFiX-GUI version %s" % mfix.get_version())

    if len(sys.argv) > 1:
        mfix.open_project(sys.argv[-1])
    else:
        # autoload last project
        project_path = mfix.get_project_path()
        if project_path:
            mfix.open_project(project_path)

    # print number of keywords
    mfix.print_internal('Registered %d keywords' %
                        len(mfix.project.registered_keywords), color='blue')

    # have to initialize vtk after the widget is visible!
    mfix.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()

    qapp.deleteLater()
    sys.exit()
