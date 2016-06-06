#!/usr/bin/env python
"""MFIX GUI"""

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import getopt
import logging
import os
import shutil
import signal
import sys
import time
import traceback
from collections import OrderedDict

# Initialize logger early
SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
log = logging.getLogger(__name__)
log.debug(SCRIPT_DIRECTORY)

# import qt
from qtpy import QtCore, QtWidgets, QtGui, PYQT4, PYQT5
from qtpy.QtCore import Qt, QFileSystemWatcher, QSettings, QUrl, QUrl, pyqtSignal

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16

PRECOMPILE_UI = False

if not PRECOMPILE_UI:
    try:
        from PyQt5 import uic
    except ImportError:
        from PyQt4 import uic


# local imports
from project_manager import ProjectManager
from job import Job
from monitor import Monitor

from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table, BaseWidget)
from widgets.regions import RegionsWidget
from widgets.linear_equation_table import LinearEquationTable
from widgets.species_popup import SpeciesPopup
from widgets.workflow import WorkflowWidget, PYQTNODE_AVAILABLE

from fluid_handler import FluidHandler
from solid_handler import SolidHandler

from tools.general import (make_callback, get_icon, get_mfix_home,
                           widget_iter, set_script_directory,
                           format_key_with_args, to_unicode_from_fs,
                           set_item_noedit, get_selected_row)

set_script_directory(SCRIPT_DIRECTORY)

from tools.namelistparser import buildKeywordDoc

from constants import *

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


class MfixGui(QtWidgets.QMainWindow, FluidHandler, SolidHandler):
    """Main window class handling all gui interactions"""

    settings = QSettings('MFIX', 'MFIX')

    stdout_signal = pyqtSignal(str)
    stderr_signal = pyqtSignal(str)
    update_run_options_signal = pyqtSignal(str)

    def __init__(self, app, parent=None, project_file=None):
        # load settings early so get_project_file returns the right thing.
        if project_file:
            self.set_project_file(project_file)

        QtWidgets.QMainWindow.__init__(self, parent)
        self.setWindowIcon(get_icon('mfix.png'))

        # reference to qapp instance (why?)
        #self.app = app

        # Initialize data members
        self.solver_name = None
        self.mfix_exe = None
        self.mfix_config = None
        self.smp_enabled = False
        self.dmp_enabled = False
        self.mfix_available = False
        self.pymfix_enabled = False
        self.open_succeeded = False
        self.unsaved_flag = False

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
            self.ui.panes = [] # Convenience
            def make_widget(cls):
                # Create an instance of a new class which is a subclass
                # of QWidget and the specified class
                class Widget(QtWidgets.QWidget, cls):
                    def __init__(self):
                        QtWidgets.QWidget.__init__(self, parent)
                        self.setupUi(self)
                return Widget()

            for cls in (Ui_general, Ui_geometry, Ui_mesh, RegionsWidget,
                        Ui_model_setup, Ui_fluid, Ui_solids,
                        Ui_numerics, Ui_output, Ui_vtk,
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
                self.ui.panes.append(widget)

        else:  # not precompiled
            uifiles = os.path.join(SCRIPT_DIRECTORY, 'uifiles')
            self.ui = uic.loadUi(os.path.join(uifiles, 'gui.ui'))
            self.ui.panes = [] # Convenience
            self.setCentralWidget(self.ui)
            assert self is not self.ui

            for name in ('general', 'geometry', 'mesh', 'regions',
                         'model_setup', 'fluid', 'solids', 'numerics',
                         'output', 'vtk','monitors', 'run'):
                if name == 'regions':  # not loaded from .ui file
                    widget = RegionsWidget()
                else:
                    widget = QtWidgets.QWidget()
                    try:
                        path = os.path.join(uifiles, name+'.ui')
                        uic.loadUi(path, widget)
                    except Exception:
                        # report which ui file it was, otherwise stack trace
                        # is too generic to be helpful.
                        print("Error loading", path)
                        raise

                # assign 'self.ui.general', etc
                setattr(self.ui, name, widget)
                self.ui.stackedWidgetTaskPane.addWidget(widget)
                self.ui.panes.append(widget)
        # end of ui loading

        #self.status_message("Loading")
        # build keyword documentation from namelist docstrings
        # TODO: pregenerate this?
        self.keyword_doc = buildKeywordDoc(os.path.join(SCRIPT_DIRECTORY,
                                                        os.pardir, 'model'))

        if False:
            keys = self.keyword_doc.keys()
            keys.sort()
            with open('/tmp/keys','w') as f:
                f.write('\n'.join(keys))


        self.species_popup = SpeciesPopup(QtWidgets.QDialog())
        #self.species_popup.setModal(True) # ?

        self.init_fluid_handler()
        self.init_solid_handler()

        # create project manager
        # NOTE.  it's a ProjectManager, not a Project.  But
        # ProjectManager is a subclass of Project.  Please
        # do not "fix" the code by renaming self.project to
        # self.project_manager
        self.project = ProjectManager(self, self.keyword_doc)

        # --- animation data
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper}
        self.animation_speed = 400
        self.animating = False
        self.stack_animation = None

        # --- icons ---
        # loop through all widgets & set icons for any ToolButton with add/delete/copy
        #  in the name
        for widget in widget_iter(self):
            if isinstance(widget, QtWidgets.QToolButton):
                name = str(widget.objectName())
                if 'add' in name:
                    widget.setIcon(get_icon('add.png'))
                elif 'delete' in name or 'remove' in name:
                    widget.setIcon(get_icon('remove.png'))
                elif 'copy' in name:
                    widget.setIcon(get_icon('copy.png'))

        # Toolbuttons at top of frame
        ui = self.ui
        for (button, icon_name, function) in (
                (ui.toolbutton_new, 'newfolder', self.new_project),
                (ui.toolbutton_open, 'openfolder', self.handle_open),
                (ui.toolbutton_save, 'save', self.handle_save),
                (ui.toolbutton_run_mfix, 'play', self.handle_run),
                (ui.toolbutton_pause_mfix, 'pause', self.handle_pause),
                (ui.toolbutton_stop_mfix, 'stop', self.handle_stop),
                (ui.toolbutton_reset_mfix, 'restart', self.remove_output_files)):
            button.setIcon(get_icon(icon_name+'.png'))
            button.clicked.connect(function)

        # Make sure lineedits lose focus so keywords update before save/run !!
        for button in (ui.toolbutton_run_mfix, ui.toolbutton_save, ui.toolbutton_more):
            button.setFocusPolicy(Qt.ClickFocus)

        # "More" submenu
        ui.toolbutton_more.setIcon(get_icon('more_vert_black_crop.png'))
        ui.menu_more = QtWidgets.QMenu()
        ui.toolbutton_more.setMenu(ui.menu_more)
        self.ui.action_save_as = self.ui.menu_more.addAction(
            get_icon('save.png'), 'Save As', self.handle_save_as)
        self.ui.action_export = self.ui.menu_more.addAction(
            get_icon('open_in_new.png'), 'Export', self.handle_export)
        self.ui.action_compile_tool = self.ui.menu_more.addAction(
            get_icon('build.png'), 'Compile', self.handle_compile)
        self.ui.menu_more.addSeparator()
        self.ui.action_about = self.ui.menu_more.addAction(
            get_icon('settings.png'), 'Settings', self.handle_settings)
        self.ui.action_about = self.ui.menu_more.addAction(
            get_icon('help.png'), 'Help', self.handle_help)
        self.ui.action_about = self.ui.menu_more.addAction(
            get_icon('infooutline.png'), 'About', self.handle_about)
        self.ui.menu_more.addSeparator()
        self.ui.action_close = self.ui.menu_more.addAction(
            get_icon('close.png'), 'Close', self.close)

        # Geometry toolbuttons
        geo = self.ui.geometry
        geo.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        geo.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        geo.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        geo.toolbutton_geometry_intersect.setIcon(get_icon('intersect.png'))
        geo.toolbutton_geometry_difference.setIcon(get_icon('difference.png'))


        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.clicked.connect(make_callback(self.mode_changed, mode))

        # navigation tree
        ui.treewidget_model_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        self.job = Job(parent=self)
        self.rundir_watcher = QFileSystemWatcher() # Move to monitor class
        self.rundir_watcher.directoryChanged.connect(self.slot_rundir_changed)

        self.exe_watcher = QFileSystemWatcher()
        self.exe_watcher.directoryChanged.connect(self.slot_exes_changed)

        self.monitor = Monitor(self)

        # buttons in 'run' pane
        run = ui.run
        run.button_run_mfix.clicked.connect(self.handle_run)
        run.button_pause_mfix.clicked.connect(self.handle_pause)
        run.button_pause_mfix.setVisible(self.pymfix_enabled)
        run.button_stop_mfix.clicked.connect(self.handle_stop)
        run.button_reset_mfix.clicked.connect(self.remove_output_files)
        run.combobox_mfix_exes.activated.connect(self.handle_select_exe)

        # Print welcome message.  Do this early so it appears before any
        # other messages that may occur during this __init__
        self.print_welcome()

        ## Run signals
        self.stdout_signal.connect(self.print_out)
        self.stderr_signal.connect(self.print_err)
        self.update_run_options_signal.connect(self.update_run_options)

        # --- setup widgets ---
        self.__setup_simple_keyword_widgets()
        self.__setup_other_widgets()  # refactor/rename - cgw

        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        self.__setup_workflow_widget()

        # --- default ---
        self.mode_changed('modeler')
        self.change_pane('general') #? start at the top?

        # some data fields that are not in Project
        self.fluid_species = OrderedDict()
        self.solids = OrderedDict()
        self.solids_current_phase = None

        # Update run options
        self.update_run_options()

        # Reset everything to default values
        self.reset() # Clear command_output too?
        self.last_line_blank = False
        # end of __init__ (hooray!)


    def set_no_project(self):
        """setup mode when no project is open"""
        self.open_succeeded = False
        self.set_solver(None)
        self.set_project_file(None)
        self.clear_unsaved_flag()
        self.update_window_title()
        self.enable_input(False)
        self.ui.toolbutton_new.setEnabled(True)
        self.ui.toolbutton_open.setEnabled(True)
        # This gets set by guess_solver if we're loading a project, otherwise
        # we need to set the default.  (Do other defaults need to be set here?)
        self.status_message("No project - open existing MFIX project or create a new one")

    def reset(self):
        """Reset all widgets to default values and set GUI to blank-slate"""
        #self.mfix_exe = None
        #self.mfix_config = None
        #self.smp_enabled = False
        #self.dmp_enabled = False
        #self.pymfix_enabled = False

        # ---- parameters which do not map neatly to keywords
        self.fluid_nscalar_eq = 0
        self.solid_nscalar_eq = 0 # Infer these from phase4scalar
        # Defaults

        self.solver_name = None

        self.project.reset() # Clears all keywords & collections

        #reset filesystem watchers: TODO: promote watchers to their own objects
        for w in (self.exe_watcher, self.rundir_watcher):
            for d in w.directories():
                w.removePath(d)
        self.slot_rundir_changed()
        # just the system dirs, no project dirs (TODO: reevaluate $PATH)
        PATH = os.environ.get("PATH")
        if PATH:
            dirs = set(PATH.split(os.pathsep))
        else:
            dirs = set()
        mfix_home = get_mfix_home()
        if mfix_home:
            dirs.add(mfix_home)
            dirs.add(os.path.join(mfix_home, 'bin'))
            dirs.add(os.path.join(mfix_home, 'build'))
        for d in dirs:
            # filter out empty strings and current directory from $PATH
            if d and d != os.path.curdir and os.path.isdir(d):
                if d not in self.exe_watcher.directories():
                    self.exe_watcher.addPath(d)
        self.slot_exes_changed()

        if self.mfix_exe: # Do we need to do this again here?
            cb = self.ui.run.combobox_mfix_exes
            if cb.findText(self.mfix_exe) == -1:
                cb.addItem(self.mfix_exe)
            cb.setCurrentText(self.mfix_exe)

        self.reset_fluids()
        self.reset_solids()

        # Set all custom widgets to default
        for w in widget_iter(self):
            if isinstance(w, BaseWidget):
                w.default()
            elif hasattr(w, 'default'):
                w.default()
            else:
                pass # What to do for rest of widgets?

        self.unsaved_flag = False
        #self.clear_unsaved_flag() - sets window title to MFIX - $project_file
        #self.set_project_file(None)  - do we want to do this?


    def confirm_close(self):
        # TODO : option to save
        msg = None
        if self.job.is_running():
            msg = "Stop running job?"
        elif self.job.is_paused():
            msg = "Stop paused job?"
        if msg:
            confirm = self.message(text=msg,
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm == 'no':
                return
            log.info("Stopping mfix at application exit")
            self.job.stop_mfix()

        if self.unsaved_flag:
            confirm = self.message(text="File not saved, really quit?",
                                   buttons=['yes', 'no'],
                                   default='no')
            return confirm == 'yes'
        else:
            return True

    def set_keyword(self, key, value, args=None):
        """convenience function to set keyword"""
        self.set_unsaved_flag()#?
        self.project.submit_change(None, {key:value}, args)

    def update_keyword(self, key, value, args=None):
        """like set_keyword but no action if value already set"""
        if self.project.get_value(key, args=args) == value:
            return
        self.set_keyword(key, value, args)

    def unset_keyword(self, key, args=None):
        """Undefine keyword.  Report to user, also catch and report any errors"""
        if isinstance(args, int):
            args = [args]
        elif args is None:
            args = []
        try:
            success = self.project.removeKeyword(key, args, warn=False)
            if success:
                self.set_unsaved_flag()
                self.print_internal("%s" % format_key_with_args(key, args),
                                    font='strikeout')
        except Exception as e:
            msg = 'Warning: Failed to unset %s: %s: %s' % (format_key_with_args(key, args),
                                                  e.__class__.__name__, e)
            self.print_internal(msg, color='red')
            traceback.print_exception(*sys.exc_info())

    def unimplemented(self):
        self.message(title='Unimplemented',
                     text='Feature not implemented')


    def new_project_X(self):
        # Make run name editable FIXME
        self.ui.general.lineedit_keyword_run_name.setEnabled(True)
        self.ui.stackedwidget_mode.setEnabled(True)
        self.reset()

    def status_message(self, message=''):
        self.ui.label_status.setText(message)

    def update_no_mfix_warning(self):
        ok = bool(self.mfix_exe or self.mfix_available)
        self.ui.run.label_mfix_exes_warning.setVisible(not ok)
        self.ui.run.combobox_mfix_exes.setVisible(ok)
        if not ok:
            self.print_internal("Warning: no MFIX executables available")

    def slot_rundir_changed(self):
        # Note: since log files get written to project dirs, this callback
        # is triggered frequently during a run.
        log.debug("rundir changed")
        # TODO figure out if we really need to do this update
        self.update_run_options()

    def slot_exes_changed(self):
        # The list of executables (maybe) changed.
        # Note: since log files get written to project dirs, this callback
        # is triggered frequently during a run
        running = self.job.is_running()
        res_file_exists = bool(self.monitor.get_res_files())

        exes = list(self.monitor.get_exes())
        mfix_exe = self.mfix_exe

        # Did the current exe go away?
        if mfix_exe and not os.path.exists(mfix_exe):
            self.print_internal("Warning: %s is gone" % mfix_exe)
            mfix_exe = self.mfix_exe = None

        # Make sure we don't loose the current selection
        if mfix_exe and mfix_exe not in exes:
            exes.insert(0, mfix_exe)
        self.mfix_available = bool(exes)
        self.update_no_mfix_warning()

        cb = self.ui.run.combobox_mfix_exes
        if not self.mfix_available:
            # Disable run/pause/stop buttons
            self.set_run_button(enabled=False)
            self.set_stop_button(enabled=False)
            self.set_reset_button(enabled=(res_file_exists and not running))
            # How did we get here if running? maybe somebody deleted the exe!
            return

        cb.clear()
        for exe in exes:
            cb.addItem(exe)

        if self.mfix_exe:
            cb.setCurrentText(self.mfix_exe)

    def set_run_button(self, text=None, enabled=None):
        if text is not None:
            self.ui.run.button_run_mfix.setText(text)
            self.ui.toolbutton_run_mfix.setToolTip('Resume previous MFIX run' if text=='Resume'
                                                   else text+' MFIX')
        if enabled is not None:
            for b in (self.ui.run.button_run_mfix, self.ui.toolbutton_run_mfix):
                b.setEnabled(enabled)

    def set_pause_button(self, text=None, enabled=None, visible=None):
        buttons = (self.ui.run.button_pause_mfix, self.ui.toolbutton_pause_mfix)

        if enabled is not None:
            for b in buttons:
                b.setEnabled(enabled)
        if visible is not None:
            for b in buttons:
                b.setVisible(visible)
        if text is not None:
            self.ui.run.button_pause_mfix.setText(text)
            self.ui.toolbutton_pause_mfix.setToolTip(text + ' MFIX')

    def set_stop_button(self, enabled):
        for b in (self.ui.run.button_stop_mfix, self.ui.toolbutton_stop_mfix):
            b.setEnabled(enabled)

    def set_reset_button(self, enabled):
        for b in (self.ui.run.button_reset_mfix, self.ui.toolbutton_reset_mfix):
            b.setEnabled(enabled)

    def enable_input(self, enabled):
        # Enable/disable all inputs (while job running, etc)
        # Stop/reset buttons are left enabled
        for pane in self.ui.panes:
            pane.setEnabled(enabled)


    # TODO:  separate this into different functions - this is called by
    # several different signals for different reasons
    # This function is called a lot, and it does too much work each time
    # 1) executables changed
    # 2) project directory changed
    # 3) process started
    # 4) process stopped
    def update_run_options(self, message=None):
        """Updates list of of mfix executables and sets run dialog options"""
        # This is the main state-transition handler

        if message is not None:
            # highlight for visibility, this is an important state chage
            self.print_internal(message, color='blue')

        # TODO: set this in __init__ or another early setup method
        # assemble list of available executables

        ui = self.ui
        project_file = os.path.basename(self.get_project_file() or '')

        project_open = bool(project_file and self.open_succeeded)
        paused = self.job.is_pymfix and self.job.is_paused()
        running = self.job.is_running() and not paused
        resumable = bool(self.monitor.get_res_files()) # overlaps with running & paused
        ready = project_open and not (running or paused or resumable)

        log.debug("UPDATE RUN OPTIONS", "running=", running, "paused=", paused,
                  "resumable=", resumable)

        self.update_window_title() # put run state in window titlebar

        self.enable_input(ready)
        self.ui.run.setEnabled(project_open)

        #handle buttons in order:  RESET RUN PAUSE STOP
        # Pause only available w/ pymfix
        if running:
            self.status_message("MFIX running, process %s" % self.job.mfix_pid)
            ui.run.combobox_mfix_exes.setEnabled(False)
            # also disable spinboxes for dt, tstop unless interactive
            self.set_reset_button(enabled=False)
            self.set_run_button(enabled=False)
            self.set_pause_button(enabled=self.job.is_pausable(), visible=self.pymfix_enabled)
            self.set_stop_button(enabled=True)
            self.change_pane('run')

        elif paused:
            self.status_message("MFIX paused, process %s" % self.job.mfix_pid)
            ui.run.combobox_mfix_exes.setEnabled(False)
            self.set_reset_button(enabled=False)
            self.set_pause_button(visible=True, enabled=False)
            self.set_run_button(text="Unpause", enabled=True)
            self.set_stop_button(enabled=False)
            self.change_pane('run')

        elif resumable:
            self.status_message("Previous MFIX run is resumable.  Reset job to edit model")
            ui.run.combobox_mfix_exes.setEnabled(False)
            self.set_reset_button(enabled=True)
            self.set_run_button(text='Resume', enabled=True)
            self.set_pause_button(enabled=False, visible=self.pymfix_enabled)
            self.set_stop_button(enabled=False)
            self.change_pane('run')

        else: # Not running, ready for input
            self.status_message("Ready" if project_open else "Loading %s"%project_file)
            ui.run.combobox_mfix_exes.setEnabled(True)
            self.set_reset_button(enabled=False)
            self.set_run_button(text="Run", enabled=self.mfix_available and project_open)
            self.set_pause_button(text="Pause", enabled=False, visible=self.pymfix_enabled)
            self.set_stop_button(enabled=False)

        ui.run.spinbox_openmp_threads.setEnabled(self.smp_enabled and ready)
        ui.run.spinbox_keyword_nodesi.setEnabled(self.dmp_enabled and ready)
        ui.run.spinbox_keyword_nodesj.setEnabled(self.dmp_enabled and ready)
        ui.run.spinbox_keyword_nodesk.setEnabled(self.dmp_enabled and ready)

        ui.run.use_spx_checkbox.setEnabled(resumable)
        ui.run.use_spx_checkbox.setChecked(resumable)


    def print_welcome(self):
        self.print_internal("Welcome to MFIX - https://mfix.netl.doe.gov",
                            color='blue')
        self.print_internal("MFIX-GUI version %s" % self.get_version(),
                            color='blue')

    def get_version(self):
        return "0.2x" # placeholder

    def closeEvent(self, event):
        if not self.confirm_close():
            event.ignore()
            return
        event.accept()

    def find_navigation_tree_item(self, item_name):
        tree = self.ui.treewidget_model_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        items = tree.findItems(item_name, flags, 0)
        assert len(items) == 1
        return items[0]


    # Top-level "Model Setup"
    def set_solver(self, solver):
        """handler for "Solver" combobox in Model Setup"""
        self.project.solver = solver
        if solver is None: #
            return

        ui = self.ui
        model_setup = ui.model_setup
        cb = model_setup.combobox_solver
        if cb.currentIndex != solver:
            cb.setCurrentIndex(solver)

        solver_name = {SINGLE:"MFIX Single-Phase",
                       TFM:"MFIX-TFM",
                       DEM:"MFIX-DEM",
                       PIC:"MFIX-PIC",
                       HYBRID:"MFIX-Hybrid"}.get(solver, "MFIX")

        self.print_internal("Solver: %s" % solver_name)
        self.solver_name = solver_name

        item_names =  ("Material", "TFM", "DEM", "PIC")

        item_states = {SINGLE: (False, False, False, False),
                       TFM: (True, True, False, False),
                       DEM: (True, False, True, False),
                       PIC: (True, False, False, True),
                       HYBRID: (True, True, True, False)}

        items = (self.find_navigation_tree_item("Solids"),
                 self.ui.solids.pushbutton_solids_tfm,
                 self.ui.solids.pushbutton_solids_dem,
                 self.ui.solids.pushbutton_solids_pic)

        for (item, state) in zip(items, item_states[solver]):
            item.setDisabled(not state)

        # Don't stay on a disabled tab!
        # Do we ever disable "Materials"?
        active_tab = self.ui.solids.stackedwidget_solids.currentIndex()
        if active_tab > 0 and not item_states[solver][active_tab]:
            if solver==SINGLE:
                i, p = 0, self.ui.solids.pushbutton_solids_materials  # This one should be open (?)
            elif solver in (TFM, HYBRID):
                i, p = 1, self.ui.solids.pushbutton_solids_tfm
            elif solver==DEM:
                i, p = 2, self.ui.solids.pushbutton_solids_dem
            elif solver==PIC:
                i, p = 3, self.ui.solids.pushbutton_solids_pic
            self.solids_change_tab(i, p)

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

        ui.fluid.checkbox_enable_fluid_scalar_eq.setEnabled(enabled)
        ui.fluid.spinbox_fluid_nscalar_eq.setEnabled(enabled
                    and self.ui.fluid.checkbox_enable_fluid_scalar_eq.isChecked())

        # Solids Model selection tied to Solver
        # FIXME XXX What to do about solids that are already defined?

        valid_models = (("DEM",) if solver==DEM
                        else ("TFM",) if solver==TFM
                        else ("PIC",) if solver==PIC
                        else ("TFM", "DEM"))

        for (i,(k,v)) in enumerate(self.solids.items(), 1):
            model = v.get('model')
            if model not in valid_models:
                model = valid_models[0]
                self.update_keyword('solids_model', model, args=[i])
                v['model'] = model
        self.update_solids_table()
        self.setup_combobox_solids_model()
        self.update_solids_detail_pane()
        self.update_window_title()

    def enable_energy_eq(self, state):
        # Additional callback on top of automatic keyword update,
        # since this has to change availabilty of several other GUI items
        self.ui.model_setup.checkbox_keyword_energy_eq.setChecked(state)
        ui = self.ui
        for item in (ui.fluid.combobox_fluid_specific_heat_model,
                     ui.fluid.combobox_fluid_conductivity_model,
                     # more ?
                     ):
            item.setEnabled(state)
        # c_pg0 == specific heat for fluid phase
        lineedit = ui.fluid.lineedit_keyword_c_pg0
        if state:
            lineedit.setEnabled(self.fluid_specific_heat_model == CONSTANT)
        else:
            lineedit.setEnabled(False)


    def set_subgrid_model(self, index):
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)

    def update_scalar_equations(self, prev_nscalar):
        """sets nscalar and phase4scalar(#) for all phases"""
        # Used by both fluid & solid han
        # This is a little messy.  We may have reduced
        # nscalar, so we need to unset phase4scalar(i)
        # for any values of i > nscalar.
        nscalar = self.fluid_nscalar_eq + self.solid_nscalar_eq
        if nscalar > 0:
            self.update_keyword("nscalar", nscalar)
        else:
            self.unset_keyword("nscalar")

        for i in range(1,1+self.fluid_nscalar_eq):
            self.update_keyword("phase4scalar", 0, args=i)
        i = 1+self.fluid_nscalar_eq
        for (phase, s) in enumerate(self.solids.values(), 1):
            n = s.get('nscalar_eq', 0)
            for j in range(n):
                self.update_keyword("phase4scalar", phase, args=i)
                i += 1
        while i <= prev_nscalar:
            self.unset_keyword("phase4scalar", i)
            i += 1


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

        combobox = model_setup.combobox_subgrid_model
        combobox.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)

        #self.enable_energy_eq(False)

        # Fluid phase - move to fluid_handlers.py
        checkbox = ui.fluid.checkbox_keyword_species_eq_args_0
        checkbox.stateChanged.connect(self.enable_fluid_species_eq)

        ui.fluid.lineedit_fluid_phase_name.editingFinished.connect(
            self.handle_fluid_phase_name)
        ui.fluid.checkbox_enable_fluid_scalar_eq.stateChanged.connect(
            self.enable_fluid_scalar_eq)
        ui.fluid.spinbox_fluid_nscalar_eq.valueChanged.connect(
            self.set_fluid_nscalar_eq)

        # Fluid phase models
        # Density
        for name in ('density', 'viscosity', 'specific_heat', 'mol_weight',
                     'conductivity', 'diffusion'):
            combobox = getattr(ui.fluid, 'combobox_fluid_%s_model' % name)
            setter = getattr(self,'set_fluid_%s_model' % name)
            combobox.currentIndexChanged.connect(setter)

        # Fluid species
        f = ui.fluid
        tb = f.toolbutton_fluid_species_add
        tb.clicked.connect(self.fluid_species_add)
        tb = f.toolbutton_fluid_species_copy # misnomer
        tb.clicked.connect(self.fluid_species_edit)
        tb.setEnabled(False)
        tb = f.toolbutton_fluid_species_delete
        tb.setEnabled(False)
        tb.clicked.connect(self.fluid_species_delete)
        tw = f.tablewidget_fluid_species
        tw.itemSelectionChanged.connect(self.handle_fluid_species_selection)

        # Solid phase
        s = ui.solids
        tb = s.toolbutton_solids_add
        tb.clicked.connect(self.solids_add)
        tb = s.toolbutton_solids_delete
        tb.clicked.connect(self.solids_delete)
        tb.setEnabled(False)
        tw = s.tablewidget_solids
        # Hack - force summary table to update  on kw updates
        class TableWidgetProxy:
            def objectName(self):
                return "proxy"
            def updateValue(*args):
                self.update_solids_table()

        self.project.register_widget(TableWidgetProxy(),
                                     ['solids_model', 'd_p0', 'ro_s0'], args='*')
        tw.itemSelectionChanged.connect(self.handle_solids_table_selection)
        cb = s.combobox_solids_model
        cb.currentIndexChanged.connect(self.handle_combobox_solids_model)
        s.lineedit_solids_phase_name.editingFinished.connect(self.handle_solids_phase_name)
        s.checkbox_enable_scalar_eq.stateChanged.connect(self.enable_solid_scalar_eq)
        s.spinbox_nscalar_eq.valueChanged.connect(self.set_solid_nscalar_eq)

        # connect solid tab btns
        for i, btn in enumerate((s.pushbutton_solids_materials,
                                 s.pushbutton_solids_tfm,
                                 s.pushbutton_solids_dem,
                                 s.pushbutton_solids_pic)):
            btn.pressed.connect(
                make_callback(self.solids_change_tab, i, btn))

        self.solids_change_tab(0, s.pushbutton_solids_materials) # ?

        # numerics
        ui.linear_eq_table = LinearEquationTable(ui.numerics)
        ui.numerics.gridlayout_leq.addWidget(ui.linear_eq_table)
        self.project.register_widget(ui.linear_eq_table,
                                     ['discretize', 'leq_method', 'leq_tol',
                                      'leq_it', 'leq_sweep', 'leq_pc',
                                      'ur_fac'],
                                     args='*')


    def __setup_simple_keyword_widgets(self):
        """Look for and connect simple keyword widgets to the project manager.
        Keyword information from the namelist doc strings is added to each
        keyword widget. The widget must be named: *_keyword_<keyword> where
        <keyword> is the actual keyword.
        Args are also supported via widgetname_keyword_KEY_args_ARGS"""

        def try_int(str):
            try:
                return int(str)
            except ValueError:
                return str

        # loop through all widgets looking for *_keyword_<keyword>
        for widget in widget_iter(self):
            name_list = str(widget.objectName()).split('_')

            if 'keyword' in name_list:
                key_idx = name_list.index('keyword')
                args = None
                # Look for _args_ following <keyword>
                if 'args' in name_list:
                    args_idx = name_list.index('args')
                    args = map(try_int, name_list[args_idx+1:])
                    key = '_'.join(name_list[key_idx+1:args_idx])
                else:
                    key = '_'.join(name_list[key_idx+1:])

                # sometimes multiple widgets point to the same key ...
                # name them widget_keyword_KW_1, _2, etc
                if key not in self.keyword_doc:
                    name_list = key.split('_')
                    if name_list[-1].isdigit():
                        base_key = '_'.join(name_list[:-1])
                        if base_key not in self.keyword_doc:
                            log.error("UNKNOWN KEYWORD %s: not registering %s (also tried %s)" %
                                      (key, widget.objectName(), base_key))
                            continue
                        key = base_key

                # set the key attribute to the keyword
                widget.key = key
                widget.args = args

                # add info from keyword documentation
                if key in self.keyword_doc:
                    doc = self.keyword_doc[key]
                    widget.setdtype(doc['dtype'])
                    req = doc.get('required')
                    if req is not None:
                        widget.setValInfo(req)
                    vr = doc.get('validrange')
                    if vr is not None:
                        if 'max' in vr:
                            widget.setValInfo(_max=vr['max'])
                        if 'min' in doc['validrange']:
                            widget.setValInfo(_min=vr['min'])

                    default = doc.get('initpython') # "Initial Python Value"
                    if default is not None:
                        widget.default(default)

                    if isinstance(widget, QtWidgets.QComboBox) and widget.count() < 1:
                            widget.addItems(list(doc['valids'].keys()))
                else:
                    log.error("UNKNOWN KEYWORD %s: not registering %s" % (key, widget.objectName()))
                    continue

                # register the widget with the project manager
                self.project.register_widget(widget, keys=[key], args=args)

    def __setup_vtk_widget(self):
        """initialize the vtk widget"""

        if 'MFIX_NO_VTK' in os.environ:
            log.info("MFIX_NO_VTK set, creating fake VTK")
            # Create a dummy object, so we don't have to test for 'if use_vtk' all over
            class FakeVtk:
                def noop(self, *args, **kwargs):
                    return None
                def __getattr__(self, key):
                    return self if key=='vtkiren' else self.noop
            self.vtkwidget = FakeVtk()
            self.ui.regions.vtkwidget = self.vtkwidget
            return

        try:
            from widgets.vtkwidget import VtkWidget
        except ImportError:
            print("vtk not available.  Set MFIX_NO_VTK in environment to skip vtk for testing")
            sys.exit(-1)
        self.vtkwidget = VtkWidget(self.project, parent=self)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkwidget)

        # register with project manager
        self.project.register_widget(self.vtkwidget,
                                     ['xmin', 'xlength', 'ymin', 'ylength',
                                      'zmin', 'zlength', 'imax', 'jmax',
                                      'kmax', 'no_k'])

        # add reference to other widgets
        self.ui.regions.vtkwidget = self.vtkwidget

    def __setup_workflow_widget(self):
        """set up the workflow widgets if pyqtnode is available"""
        if PYQTNODE_AVAILABLE:
            self.ui.workflow_widget = WorkflowWidget(self.project, self)
            self.ui.verticallayout_workflow_mode.addWidget(
                self.ui.workflow_widget)
        else:
            self.ui.pushButtonWorkflow.setEnabled(False)
            self.ui.pushButtonWorkflow.setToolTip(
                "Workflow disabled, can't import pyqtnode")

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
                    Qt.MatchFixedString | Qt.MatchRecursive, 0)
        assert len(clist) == 1
        item = clist[0]
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
        try:
            if self.stack_animation.state() == QtCore.QAbstractAnimation.Stopped:
                widget.setCurrentIndex(to)
                from_widget = widget.widget(from_)
                from_widget.hide()
                from_widget.move(0, 0)

                if btn_layout is not None and line is not None:
                    btn_layout.addItem(btn_layout.takeAt(
                        btn_layout.indexOf(line)), 1, to)
        except AttributeError: # Happends during shutdown. TODO: remove this hack
            pass
        finally:
            self.animating = False

    # --- helper methods ---
    def message(self,
                title='Warning',
                icon='warning',
                text='This is a warning.',
                buttons=['ok'],
                default='ok',
                infoText=None,
                detailedText=None,
                ):
        """Create and display a modal message box:
        title = 'title'
        icon = 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel',
            'discard'
        default = 'ok' the default selected button
        infoText = 'extended information text'
        detailedText = 'Some details'

        Returns the pressed button.  Also prints & logs message"""
        # TODO: disable this when running in QTest
        self.print_internal(title + ": " + text)
        if infoText:
            self.print_internal(infoText)
        if detailedText:
            self.print_internal(detailedText)

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

        if detailedText:
            msgBox.setDetailedText(detailedText)

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

    def print_out(self, text):
        for line in text.split('\n'):
            self.print_internal(line, color=None, font='Courier')

    def print_err(self, text):
        for line in text.split('\n'):
            self.print_internal(line, color='red', font='Courier') # Bold fond?

    def print_internal(self, line, color=None, font=None):
        qtextbrowser = self.ui.command_output
        if not line.strip():
            if self.last_line_blank:
                return
            self.last_line_blank=True
        else:
            self.last_line_blank=False

        if not line.endswith('\n'):
            line += '\n'
        lower = line.lower()
        msg = line.strip()
        # hack. TODO: real msg types, map to font/color
        strikeout = font and font.lower() == 'strikeout'
        if strikeout:
            msg = "unset " + msg
        if 'error:' in lower:
            log.error(msg)
        if 'warning:' in lower:
            log.warn(msg)
        else:
            log.info(msg)
        cursor = qtextbrowser.textCursor()
        cursor.movePosition(cursor.End)
        char_format = QtGui.QTextCharFormat()
        # HACK is this going too far?  we should define message types, not infer from messages
        if any(x in lower for x in ('error', 'warn', 'fail')) and not 'error%' in lower:
            color='red'
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


    def handle_select_exe(self):
        """Enable/disable run options based on selected executable"""
        mfix_exe = self.ui.run.combobox_mfix_exes.currentText()
        if mfix_exe == self.mfix_exe:
            return

        self.mfix_exe = mfix_exe

        if not mfix_exe:
            self.update_run_options()
            return

        self.settings.setValue('mfix_exe', mfix_exe)
        config = self.monitor.exes.get(mfix_exe)
        self.mfix_config = config
        self.smp_enabled = 'smp' in config if config else False
        self.dmp_enabled = 'dmp' in config if config else False

        self.pymfix_enabled = any(mfix_exe.lower().endswith(x)
                                  for x in ('pymfix', 'pymfix.exe'))

        self.ui.run.spinbox_openmp_threads.setEnabled(self.smp_enabled)
        if not self.dmp_enabled:
            self.ui.run.spinbox_keyword_nodesi.setValue(1)
            self.ui.run.spinbox_keyword_nodesj.setValue(1)
            self.ui.run.spinbox_keyword_nodesk.setValue(1)
        self.ui.run.spinbox_keyword_nodesi.setEnabled(self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesj.setEnabled(self.dmp_enabled)
        self.ui.run.spinbox_keyword_nodesk.setEnabled(self.dmp_enabled)
        self.update_run_options()

    def remove_output_files(self, output_files=None, message_text=None):
        """ remove MFIX output files from current project directory

        :param output_files: List of patterns to be matched for file removal
        :type output_files: list
        :return: True for success, False for user cancel"""

        if not output_files:
            output_files = self.monitor.get_outputs()
        if not message_text:
            message_text = 'Deleting output files:\n %s' % '\n'.join(output_files)

        confirm = self.message(title="Info",
                               icon="info",
                               text=message_text,
                               buttons=['ok','cancel'],
                               default='cancel')

        if confirm != 'ok':
            return False

        for path in output_files:
            log.debug('Deleting file %s' % path)
            try:
                os.remove(path)
            except OSError as err:
                msg = 'Cannot delete %s: %s' % (path, err.strerror)
                self.print_internal("Error: %s" % msg, color='red')
                self.message(text=msg,
                             buttons=['ok'],
                             default=['ok'])
                break
        return True

    # Don't make these depend on current state, since (esp for pymfix)
    # the state variables are cached and potentially outdated
    def handle_run(self):
        name = 'Run'
        try:
            if not self.job.is_running():
                self.run_mfix()
            else:
                name='unpause'
                self.job.unpause()
        except Exception as e:
            self.print_internal("%s: error %s" % (name, e))

    def handle_pause(self):
        try:
            self.job.pause()
        except Exception as e:
            self.print_internal("Pause: error %s" % e)

    def handle_stop(self):
        try:
            self.job.stop_mfix()
        except Exception as e:
            self.print_internal("Stop: error %s" % e)

    def run_mfix(self):
        output_files = self.monitor.get_outputs()
        res_files = self.monitor.get_res_files()
        if output_files:
            if res_files:
                return self.resume_mfix() # ugh
            if not self.remove_output_files(output_files):
                log.info('output files exist and run was cancelled')
                return

        if self.unsaved_flag:
            response = self.message(title="Save?",
                                    icon="question",
                                    text="Save current project?",
                                    buttons=['yes', 'no'])
            if response=='yes': #FIXME need to catch/report errors, writeDatFile is too low-level
                self.project.writeDatFile(self.get_project_file())
            else:
                return
        self.update_keyword('run_type', 'new')
        self.project.writeDatFile(self.get_project_file())
        self.clear_unsaved_flag()
        self.update_source_view()
        self._start_mfix()

    def restart_mfix(self):
        """Restart MFIX. This will remove previous output and start a new run."""
        output_files = self.monitor.get_outputs()
        if len(output_files) > 0:
            message = "Delete all output and resume files?"
            if not self.remove_output_files(output_files, message):
                log.debug('output or resume files exist and run was cancelled')
                return
        self.update_keyword('run_type', 'new')
        #FIXME need to catch/report errors, writeDatFile is too low-level
        self.project.writeDatFile(self.get_project_file()) # XXX
        self.clear_unsaved_flag()
        self.update_source_view()
        self._start_mfix()

    def resume_mfix(self):
        """resume previously stopped mfix run"""
        if self.ui.run.use_spx_checkbox.isChecked():
            self.update_keyword('run_type', 'restart_1')
        else:
            # TODO: is it correct to remove all but *.RES ?
            spx_files = self.monitor.get_outputs(['*.SP?', "*.pvd", "*.vtp"])
            if not self.remove_output_files(spx_files):
                log.debug('SP* files exist and run was cancelled')
                return
            self.update_keyword('run_type', 'restart_2')
        #FIXME need to catch/report errors, writeDatFile is too low-level
        self.project.writeDatFile(self.get_project_file()) # XXX
        self.clear_unsaved_flag()
        self.update_source_view()
        self._start_mfix()

    def _start_mfix(self):
        """start a new local MFIX run, using pymfix, mpirun or mfix directly"""

        if not self.mfix_exe:
            self.print_internal("ERROR: MFIX not available")
            return

        mfix_exe = self.mfix_exe

        if self.dmp_enabled:
            mpiranks = (self.project.nodesi.value *
                        self.project.nodesj.value *
                        self.project.nodesk.value)

            run_cmd = ['mpirun', '-np', str(mpiranks), mfix_exe]
        else:
            # no dmp support
            run_cmd = [mfix_exe]

        if self.smp_enabled:
            if not "OMP_NUM_THREADS" in os.environ:
                os.environ["OMP_NUM_THREADS"] = self.spinbox_openmp_threads.value()
            log.info('SMP enabled with OMP_NUM_THREADS=%d', self.spinbox_openmp_threads.value())

        project_filename = os.path.basename(self.get_project_file())
        # Warning, not all versions of mfix support '-f' !
        run_cmd += ['-f', project_filename]

        msg = 'Starting %s' % ' '.join(run_cmd)
        #log.info(msg) # print_internal logs
        self.print_internal(msg, color='blue')

        self.job.start_command(
            is_pymfix=self.pymfix_enabled,
            cmd=run_cmd,
            cwd=self.get_project_dir(),
            env=os.environ)


    def export_project(self):
        """Copy project files to new directory, but do not switch to new project"""
        # Note, this does not do a save first.  Just copies existing files.
        project_file = self.get_project_file()
        if not project_file:
            self.message(text="Nothing to export",
                         buttons=['ok'])
            return

        export_file = self.get_save_filename()
        if not export_file: # User bailed out
            return
        export_dir = os.path.dirname(export_file)
        if not self.check_writable(export_dir):
            return

        files_to_copy = [project_file]

        sp_files = self.monitor.get_outputs(["*.SP?"])
        if (sp_files):
            yes = self.message(text="Copy .SP files?\n%s" % '\n'.join(output_files),
                               buttons=['yes', 'no'],
                               default='yes')
            if yes=='yes':
                files_to_copy.extend(sp_files)

        files_to_copy.extend(self.monitor.get_outputs(["*.RES", "*.STL"]))
        #copy project files into new_project_directory
        for f in files_to_copy:
            try:
                shutil.copyfile(f, os.path.join(export_dir, os.path.basename(f)))
            except Exception as e:
                self.message(text="Error copying file:\n%s" % e,
                             buttons=['ok'])


    def save_project(self, filename=None):
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

        # save geometry
        self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))

        # save regions
        self.project.mfix_gui_comments['regions_dict'] = self.ui.regions.regions_to_str()

        project_base = os.path.basename(project_file)
        run_name = os.path.splitext(project_base)[0]
        self.update_keyword('run_name', run_name)
        self.project.writeDatFile(project_file)
        self.clear_unsaved_flag()

    def save_as(self):
        """Prompt user for new filename, save project to that file and make
        it the active project"""
        new_file = self.get_save_filename()
        if not new_file:
            return
        new_dir = os.path.dirname(new_file)
        if not self.check_writable(new_dir):
            return
        # Force run name to file name.  Is this a good idea?
        basename = os.path.basename(new_file)
        run_name = os.path.splitext(basename)[0]
        self.set_project_file(new_file)
        self.update_keyword('run_name', run_name)
        self.save_project()

    def get_save_filename(self, dialog_message=None):
        """wrapper for call to getSaveFileName, override in unit tests"""

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
        project_file = self.get_project_file()
        try:
            self.save_project()
        except Exception as e:
            msg = 'Failed to save %s: %s: %s' % (project_file, e.__class__.__name__, e)
            self.print_internal("Error: %s" % msg, color='red')
            self.message(title='Error',
                         icon='error',
                         text=msg,
                         buttons=['ok'],
                         default='ok')
            traceback.print_exception(*sys.exc_info())
            return

    def handle_export(self):
        self.export_project()

    def handle_save_as(self):
        self.save_as()

    def handle_settings(self):
        """handle user settings"""
        # TODO: implement
        self.unimplemented()

    def handle_compile(self):
        """compiling tool"""
        # TODO: implement
        self.unimplemented()

    def handle_help(self):
        """show help popup"""
        # TODO: implement
        self.unimplemented()

    def handle_about(self):
        """show about popup"""
        # TODO: implement
        self.unimplemented()

    def update_window_title(self):
        title = self.solver_name or 'MFIX'
        project_file = self.get_project_file()
        if project_file:
            title += " - " + os.path.basename(project_file)
        if self.unsaved_flag:
            title += '*'

        if self.job.is_running() and not self.job.is_paused():
            title += ', RUNNING'
            if self.job.mfix_pid is not None:
                title += ', process %s'% self.job.mfix_pid
        elif self.job.is_paused():
            title += ', PAUSED'
        elif self.monitor.get_res_files():
            title += ', STOPPED, resumable'
        else:
            pass
            #title += ', EDITING'
        self.setWindowTitle(title)


    def set_save_button(self, enabled):
        # enable/disable the save button and save_as menu item
        self.ui.toolbutton_save.setEnabled(enabled)
        self.ui.action_save_as.setEnabled(enabled)


    def set_unsaved_flag(self):
        if not self.unsaved_flag:
            log.info("Project is unsaved")
        self.unsaved_flag = True
        self.update_window_title()
        self.set_save_button(enabled=True)


    def clear_unsaved_flag(self):
        if self.unsaved_flag:
            log.info("Project is saved")
        self.unsaved_flag = False
        self.update_window_title()
        self.set_save_button(enabled=False)

        # This should only happen when we're not already running
        # or paused
        runnable = bool(self.get_project_file()) and self.mfix_available
        self.set_run_button(enabled=runnable)

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
                QtWidgets.QFileDialog.getExistingDirectory(
                    self, 'Create Project in Directory',
                    ""))
        if not project_dir:
            return

        # TODO: allow user to set run name
        project_file = os.path.join(project_dir, 'mfix.dat')
        if not self.check_writable(project_dir):
            return
        # Start with a nice template - note, there's too much set in this file.
        # FIXME, this can clobber files
        template = os.path.join(get_mfix_home(), 'gui', 'mfix.dat.template')
        shutil.copyfile(template, project_file)
        self.open_project(project_file)

    def get_open_filename(self):
        """wrapper for call to getOpenFileName, override in for unit tests"""
        project_dir = self.get_project_dir()
        return QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Project Directory', project_dir)

    def handle_open(self):
        """handler for toolbar Open button"""
        if self.unsaved_flag:
            confirm = self.message(text="Project not saved\nData will be lost!\nProceed?",
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm != 'yes':
                return
            self.clear_unsaved_flag()

        project_path = self.get_open_filename()
        # qt4/qt5 compat hack
        #if type(project_path) == tuple:
        if PYQT5:
            project_path = project_path[0]
        if not project_path:
            return # user pressed Cancel
        self.open_project(project_path)


    def update_source_view(self, number_lines=True):
        project_file = self.get_project_file()
        if not project_file:
            src = ''
        else:
            try:
                with open(project_file, 'rb') as f:
                    src = to_unicode_from_fs(f.read())
            except Exception as e:
                log.error("error opening %s: %s" % (project_file, e))
                src = ''

        if number_lines:
            lines = src.split('\n')
            # Avoid extra blank lines at end
            while lines and lines[-1] == '':
                lines.pop(-1)
            src = '\n'.join('%4d:%s'%(lineno,line)
                            for (lineno,line) in enumerate(lines,1))

        self.ui.mfix_dat_source.setPlainText(src)

    def force_default_settings(self):
        self.update_keyword('chk_batchq_end', True)

    def open_project(self, project_path, auto_rename=True):
        """Open MFiX Project"""

        self.open_succeeded = False  # set to true on success

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
            self.message(title='Error',
                         icon='error',
                         text=('%s does not exist' % project_file),
                         buttons=['ok'],
                         default='ok')
            self.set_no_project()
            return

        self.reset() # resets gui, keywords, file system watchers, etc

        self.print_internal("Loading %s" % project_file, color='blue')
        try:
            self.project.load_project_file(project_file)
        except Exception as e:
            msg = 'Failed to load %s: %s: %s' % (project_file, e.__class__.__name__, e)
            self.print_internal("Error: %s" % msg, color='red')
            self.message(title='Error',
                         icon='error',
                         text=msg,
                         buttons=['ok'],
                         default='ok')
            traceback.print_exception(*sys.exc_info())
            # Should we stick this in the output window?  no, for now.

            self.set_no_project()
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
                try:
                    self.force_default_settings()
                    self.project.writeDatFile(project_file) #XX
                    self.print_internal(save_msg, color='blue')
                    self.clear_unsaved_flag()
                except Exception as ex:
                    msg = 'Failed to save %s: %s: %s' % (project_file, ex.__class__.__name__, ex)
                    self.print_internal("Error: %s" % msg, color='red')
                    self.message(title='Error',
                                 icon='error',
                                text=msg,
                                 buttons=['ok'],
                                 default='ok')
                    traceback.print_exception(*sys.exc_info())
                    return

        self.set_project_file(project_file)
        self.clear_unsaved_flag()
        self.update_source_view()

        # set up rundir watcher
        self.rundir_watcher.addPath(project_dir)
        self.slot_rundir_changed()

        # set up exe watcher (it got cleared in 'reset')
        dirs = (project_dir, ) # any others?
        for d in dirs:
            if d not in self.exe_watcher.directories():
                self.exe_watcher.addPath(d)
        self.slot_exes_changed()

        # Additional GUI setup based on loaded projects (not handled
        # by keyword updates)
        #    .... is there a way to verify that 'energy_eq' is boolean?
        #    should that get set from keyword doc?
        self.enable_energy_eq(bool(self.project.get_value('energy_eq')))

        # cgw - lots more model setup todo here.  Should we do this here or
        #  in ProjectManager.load_project_file (where we do guess/set_solver)
        # make sure exceptions are handled & reported
        # values that don't map to keywords, saved as #!MFIX-GUI params
        solids_phase_names = {}
        for (key, val) in self.project.mfix_gui_comments.items():
            if key == 'fluid_phase_name':
                self.set_fluid_phase_name(val)
            if key.startswith('solids_phase_name('):
                n = int(key.split('(')[1][:-1])
                solids_phase_names[n] = val
            if key == 'regions_dict':
                self.ui.regions.regions_from_str(val)
            # Add more here

        # hack, copy ordered dict to modify keys w/o losing order
        if solids_phase_names:
            s = OrderedDict()
            for (i, (k, v)) in enumerate(self.solids.items(), 1):
                s[solids_phase_names.get(i, k)] = v
            self.solids = s

        #### Fluid phase
        # fluid species table
        self.update_fluid_species_table()
        # fluid momentum and species eq. handled by _keyword_ widget
        # fluid scalar eq
        nscalar = self.project.get_value('nscalar', default=0)

        self.fluid_nscalar_eq = sum(1 for i in range(1, nscalar+1)
                                    if self.project.get_value('phase4scalar', args=i) == 0)
        self.solid_nscalar_eq = sum(1 for i in range(1, nscalar+1)
                                    if self.project.get_value('phase4scalar', args=i) != 0)

        self.enable_fluid_scalar_eq(self.fluid_nscalar_eq > 0)
        # solid scalar eq checkbox will be handled in update_solids_detail_pane

        #Move to fluid_handler
        # handle a bunch of items which are essentially the same
        for (setter, name) in ((self.set_fluid_density_model, 'ro'),
                               (self.set_fluid_viscosity_model, 'mu'),
                               (self.set_fluid_specific_heat_model, 'c_p'), # inconsistent
                               (self.set_fluid_conductivity_model, 'k'),
                               (self.set_fluid_diffusion_model, 'dif')):
            name_g0 = 'c_pg0' if name=='c_p' else name+'_g0'
            name_usr = 'usr_cpg' if name=='c_p' else 'usr_'+name+'g'
            val_g0 = self.project.get_value(name_g0)
            val_usr = self.project.get_value(name_usr)

            if (val_usr is not None and val_g0 is not None):
                self.print_internal('Warning: %s and %s are both set' % (name_g0, name_usr))
                # FIXME this is getting printed after error count ... should be included in # of errs
                # (another reason to move this to load_project_file)

            # XXX FIXME conflicts with default fluid models
            setter(CONSTANT if val_g0 is not None
                   else UDF if val_usr is not None
                   else 1)

        # molecular weight model is the odd one (only 2 settings)
        if self.project.get_value('mw_avg') is not None:
            self.set_fluid_mol_weight_model(CONSTANT)
        else:
            self.set_fluid_mol_weight_model(1)
        # requires molecular weights for all species components, when should we validate?

        ### Solids
        self.update_solids_table()
        self.update_solids_detail_pane()

        ### Geometry
        # Look for geometry.stl and load automatically
        geometry_file = os.path.abspath(os.path.join(project_dir, 'geometry.stl'))
        if os.path.exists(geometry_file):
            self.vtkwidget.add_stl(None, filename=geometry_file)

        # Look for regions in IC, BC, PS, etc.
        self.ui.regions.extract_regions(self.project)

        self.force_default_settings()

        if self.unsaved_flag: # Settings changed after loading
            self.save_project()

        self.open_succeeded = True
        self.update_run_options()

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
    if new_project and project_file: # Should we allow this - create new proj. by name?
        Usage(name)
    if project_file and not os.path.exists(project_file):
        print("%s: no such file or directory" % project_file)
        Usage(name)

    qapp = QtWidgets.QApplication([])
    mfix = MfixGui(qapp, project_file=project_file)
    mfix.show()

    # --- print welcome message
    #mfix.print_internal("MFiX-GUI version %s" % mfix.get_version())

    saved_exe = mfix.settings.value('mfix_exe') #
    cb =  mfix.ui.run.combobox_mfix_exes
    if saved_exe is not None and os.path.exists(saved_exe):
        if cb.findText(saved_exe) == -1:
            cb.addItem(saved_exe)
        cb.setCurrentText(saved_exe)
    #mfix.mfix_exe = saved_exe
    mfix.handle_select_exe()

    mfix.update_no_mfix_warning()

    if project_file is None and not new_project:
        # autoload last project
        project_file = mfix.get_project_file()

    if project_file:
        mfix.open_project(project_file, auto_rename=(not quit_after_loading))
    else:
        mfix.set_no_project()

    # print number of keywords
    mfix.print_internal('Registered %d keywords' %
                        len(mfix.project.registered_keywords))

    # have to initialize vtk after the widget is visible!
    mfix.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    if not quit_after_loading:
        qapp.exec_()

    qapp.deleteLater()
    sys.exit()

if __name__  == '__main__':
    main(sys.argv)
