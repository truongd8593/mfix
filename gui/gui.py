#!/usr/bin/env python
from __future__ import print_function, absolute_import, unicode_literals, division

"""MFIX GUI"""

import getopt
import logging
import multiprocessing
import os
import re
import shutil
import signal
import socket
import sys
import traceback
from collections import OrderedDict

# Initialize logger early
SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), ))
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))
log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)
log.debug(SCRIPT_DIRECTORY)

# import qt
from qtpy import QtCore, QtWidgets, QtGui, PYQT5
from qtpy.QtCore import Qt, QFileSystemWatcher, QSettings, pyqtSignal
UserRole = QtCore.Qt.UserRole

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16

PRECOMPILE_UI = False

if not PRECOMPILE_UI:
    from qtpy import uic


# local imports
from project_manager import ProjectManager
from project import Equation
from job import JobManager, get_dict_from_pidfile
from monitor import Monitor

from widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table, BaseWidget)
from widgets.regions import RegionsWidget
from widgets.linear_equation_table import LinearEquationTable
from widgets.species_popup import SpeciesPopup
from widgets.regions_popup import RegionsPopup
from widgets.run_popup import RunPopup
from widgets.workflow import WorkflowWidget, PYQTNODE_AVAILABLE
from widgets.parameter_dialog import ParameterDialog

from fluid_handler import FluidHandler
from solids_handler import SolidsHandler
from ics import ICS
from bcs import BCS

from interpreter import Interpreter

from tools.general import (make_callback, get_icon, get_mfix_home,
                           widget_iter, set_script_directory,
                           format_key_with_args, to_unicode_from_fs)


set_script_directory(SCRIPT_DIRECTORY)

from tools.namelistparser import buildKeywordDoc

from constants import *

if PRECOMPILE_UI:
    try:
        from uifiles.fluid import Ui_fluid
        from uifiles.geometry import Ui_geometry
        from uifiles.gui import Ui_MainWindow
        from uifiles.mesh import Ui_mesh
        from uifiles.model import Ui_model
        from uifiles.monitors import Ui_monitors
        from uifiles.numerics import Ui_numerics
        from uifiles.output import Ui_output
        from uifiles.post_processing import Ui_post_processing
        from uifiles.run import Ui_run
        from uifiles.solids import Ui_solids
        from uifiles.vtk import Ui_vtk
    except ImportError:
        print("You must compile ui files! (run 'make'")
        sys.exit(1)

# --- Main Gui ---

class MfixGui(QtWidgets.QMainWindow,
              FluidHandler,
              SolidsHandler,
              ICS,
              BCS,
              Interpreter):
    """Main window class for MFIX-GUI"""
    settings = QSettings('MFIX', 'MFIX')

    stdout_signal = pyqtSignal(str)
    stderr_signal = pyqtSignal(str)
    signal_update_runbuttons = pyqtSignal(str)

    # Allow LineEdit widgets to report out-of-bounds values.
    def popup_value_error(self, exc):
        self.message(title='Error', text=str(exc))

    def error(self, msg, popup=False):
        """Convenience function to show the user a warning & log it"""
        self.print_internal('Error: %s' % msg)
        # No popup

    def warn(self, msg, popup=False):
        """Convenience function to show the user a warning & log it"""
        if not popup:
            self.print_internal("Warning: %s" % msg)
            # print_internal will call log.warn if message starts with "Warning"
        else:
            self.message("Warning: %s" % msg)
            # Will also print-internal and log
    warning = warn

    def __init__(self, app, parent=None, project_file=None):
        # load settings early so get_project_file returns the right thing.
        if project_file:
            self.set_project_file(project_file)

        self.app = app
        LineEdit.value_error = self.popup_value_error

        QtWidgets.QMainWindow.__init__(self, parent)
        self.setWindowIcon(get_icon('mfix.png'))

        self.message_box = None # for tests to access


        # Initialize data members - make sure these values match 'reset'!
        self.solver_name = None
        self.fluid_solver_disabled = False
        self.mfix_exe = None
        self.mfix_exe_flags = {}
        self.commandline_option_exe = None
        self.mfix_available = False
        self.open_succeeded = False
        self.unsaved_flag = False
        self.run_dialog = None

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

            for cls in (Ui_model_setup,
                        Ui_geometry,
                        Ui_mesh,
                        RegionsWidget,
                        Ui_Fluid,
                        Ui_Solids,
                        Ui_initial_conditions,
                        # Ui_boundary_conditions,
                        # Ui_point_sources,
                        # Ui_internal_surfaces,
                        # Ui_chemistry,
                        Ui_numerics,
                        Ui_output,
                        Ui_vtk,
                        Ui_monitors,
                        Ui_post_processing,
                        Ui_run):
                if cls == RegionsWidget: # not loaded from ui file
                    widget = RegionsWidget(parent=self)
                    name = 'regions'
                else:
                    widget = make_widget(cls)
                    name = cls.__name__.split('_',1)[1] # part after "Ui_"
                # assign 'self.ui.model_setup', etc
                setattr(self.ui, name, widget)
                self.ui.stackedWidgetTaskPane.addWidget(widget)
                self.ui.panes.append(widget)

        else:  # not precompiled
            uifiles = os.path.join(SCRIPT_DIRECTORY, 'uifiles')
            self.ui = uic.loadUi(os.path.join(uifiles, 'gui.ui'))
            self.ui.panes = [] # Convenience
            self.setCentralWidget(self.ui)
            assert self is not self.ui

            for name in ('model_setup',
                         'geometry',
                         'mesh',
                         'regions',
                         'fluid',
                         'solids',
                         'initial_conditions',
                         #'boundary_conditions',
                         #'point_sources',
                         #'internal_surfaces',
                         #'chemistry',
                         'numerics',
                         'output',
                         'vtk',
                         'monitors',
                         'run',
                         #'post-processing'
                         ):
                if name == 'regions':  # not loaded from .ui file
                    widget = RegionsWidget(parent=self)
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

        # A few more ranges etc not mentioned in namelist doc
        self.add_extra_keyword_doc()

        if False:
            keys = self.keyword_doc.keys()
            keys.sort()
            with open('/tmp/keys','w') as f:
                f.write('\n'.join(keys))


        # Setup the navigation tree widget
        tw = self.ui.treewidget_navigation
        self.max_label_len = tw.fontMetrics().width('Boundary Conditions') + 10

        self.nav_labels = [("Model Setup", "Model"),
                           ("Post-processing", "Post"),
                           ("Boundary Conditions", "BCs"),
                           ("Initial Conditions", "ICs"),
                           ("Point Sources", "Points"), ### Only under "Model Setup", not "Monitors"
                           ("Internal Surfaces", "Surfaces")]

        # Set tooltips for nav tree & set a data property on
        # Monitors / Points to distinguish it
        root = tw.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            item.setToolTip(0, item.text(0))
            for j in range(item.childCount()):
                subitem = item.child(j)
                subitem.setToolTip(0, subitem.text(0))
                if item.text(0)=='Monitors' and subitem.text(0) == 'Points':
                    subitem.setData(UserRole, 0, True) # Mark this item

        # Intercept the resize event
        tw.resizeEvent = (lambda old_method:
                          (lambda event:
                           (self._on_resized(event),
                            old_method(event))[-1]))(tw.resizeEvent)

        # Disable items that are not yet implemented
        for name in ('Boundary Conditions',
                     'Point Sources',
                     'Internal Surfaces',
                     'Chemistry',
                     'Monitors',
                     'Points',
                     'Planes',
                     'Volumes',
                     'Post-processing',
                     'Export',
                     'Plugins'):
            self.find_navigation_tree_item(name).setDisabled(True)


        # Initialize popup dialogs
        self.species_popup = SpeciesPopup(QtWidgets.QDialog())
        #self.species_popup.setModal(True) # ?
        self.regions_popup = RegionsPopup(QtWidgets.QDialog())

        # Create project manager
        # NOTE.  it's a ProjectManager, not a Project.  But
        # ProjectManager is a subclass of Project.  Please
        # do not "fix" the code by renaming self.project to
        # self.project_manager
        self.project = ProjectManager(self, self.keyword_doc)

        # Extra setup for fluid & solids panes.  Needs to happen
        # after we create ProjectManager, because widgets get registered
        self.init_fluid_handler()
        self.init_solids_handler()
        self.init_ics()
        self.init_bcs()

        # In-process REPL (for development, should we enable this for users?)
        self.init_interpreter()

        # --- animation data
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper,
                               'interpreter': self.ui.pushButtonInterpreter}
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
                (ui.toolbutton_settings, 'settings', self.toggle_nav_menu),
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
        self.ui.action_compile_tool = self.ui.menu_more.addAction(
            get_icon('functions.png'), 'Parameters', self.handle_parameters)
        self.ui.menu_more.addSeparator()
        #self.ui.action_about = self.ui.menu_more.addAction(
        #    get_icon('settings.png'), 'Settings', self.handle_settings)
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
        ui.treewidget_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # Make tree fully open & non-closable
        # We expect "rootIsDecorated" has been set False in the .ui file
        tree = ui.treewidget_navigation
        tree.expandAll()
        tree.setExpandsOnDoubleClick(False)
        tree.setMaximumWidth(tree.fontMetrics().width('Boundary Conditions') + 10)
        tree.setMinimumWidth(tree.fontMetrics().width('Chemistry') + 10)

        # Job manager / monitor
        self.job_manager = JobManager(self)
        self.job_manager.sig_change_run_state.connect(self.slot_update_runbuttons)
        self.job_manager.sig_update_job_status.connect(self.slot_update_residuals)
        self.rundir_watcher = QFileSystemWatcher() # Move to monitor class
        self.rundir_watcher.directoryChanged.connect(self.slot_rundir_changed)

        self.monitor = Monitor(self)

        # buttons in 'run' pane
        run = ui.run
        run.button_run_mfix.clicked.connect(self.handle_run)
        run.button_pause_mfix.clicked.connect(self.handle_pause)
        run.button_pause_mfix.setVisible(True)
        run.button_stop_mfix.clicked.connect(self.handle_stop)
        run.button_reset_mfix.clicked.connect(self.remove_output_files)
        run.checkbox_pymfix_output.stateChanged.connect(self.handle_set_pymfix_output)

        # Print welcome message.  Do this early so it appears before any
        # other messages that may occur during this __init__
        self.print_welcome()

        ## Run signals
        self.stdout_signal.connect(self.handle_stdout)
        self.stderr_signal.connect(self.handle_stderr)
        self.signal_update_runbuttons.connect(self.slot_update_runbuttons)

        # --- setup widgets ---
        self.__setup_simple_keyword_widgets()
        self.__setup_other_widgets()  # refactor/rename - cgw

        # --- vtk setup ---
        self.__setup_vtk_widget()

        # --- workflow setup ---
        self.__setup_workflow_widget()

        # --- parameter dialog
        self.parameter_dialog = ParameterDialog(self)

        # --- default ---
        self.mode_changed('modeler')
        self.change_pane('model setup')

        # Update run options
        log.info('init update_runbuttons')
        self.signal_update_runbuttons.emit('')

        # Reset everything to default values
        self.reset() # Clear command_output too?
        # end of __init__ (hooray!)

    def add_extra_keyword_doc(self):
        # Add a little extra knowledge ...
        # These are all fractions, must be btwn 0 and 1, not documented as such
        for key in ('des_em', 'eps_f_min'):
            self.keyword_doc[key]['validrange'] = {'min':0.0, 'max':1.0}
        self.keyword_doc['particles']['validrange'] = {'min':0.0}

        # IC_EP_G should be btwn 0 and 1.
        self.keyword_doc['ic_ep_g']['validrange'] = {'min':0.0, 'max':1.0}
        self.keyword_doc['ic_ep_s']['validrange'] = {'min':0.0, 'max':1.0}

        # All temperatures > 0 ?

    def set_no_project(self):
        """setup mode when no project is open"""
        self.open_succeeded = False
        self.set_solver(None)
        self.set_project_file(None)
        self.clear_unsaved_flag() # sets save button
        self.set_save_as_action(enabled=False)
        self.update_window_title()
        self.enable_input(False)
        self.ui.toolbutton_new.setEnabled(True)
        self.ui.toolbutton_open.setEnabled(True)
        # This gets set by guess_solver if we're loading a project, otherwise
        # we need to set the default.  (Do other defaults need to be set here?)
        self.status_message("No project - open existing MFIX project or create a new one")
        self.change_pane("model setup")

    def reset(self):
        """Reset all widgets to default values and set GUI to blank-slate"""

        self.change_pane("model setup") # Default pane

        # ---- parameters which do not map neatly to keywords
        self.fluid_nscalar_eq = 0
        self.solids_nscalar_eq = 0 # Infer these from phase4scalar

        # Defaults - see __init__
        self.solver_name = None
        self.fluid_solver_disabled = False  # TODO: infer at load time

        self.project.reset() # Clears all keywords & collections

        self.slot_rundir_changed()

        self.reset_fluids()
        self.reset_solids()
        self.reset_ics()
        self.reset_bcs()

        # Set all custom widgets to default
        for w in widget_iter(self):
            if isinstance(w, BaseWidget):
                w.default()
            elif hasattr(w, 'default'):
                w.default()
            else:
                pass # What to do for rest of widgets?

        # reset parameters
        base_parameters = OrderedDict([(key, 0.0) for key in SPECIAL_PARAMETERS])
        PARAMETER_DICT.clear()
        PARAMETER_DICT.update(base_parameters)

        self.unsaved_flag = False

        self.update_nav_tree()
        #self.clear_unsaved_flag() - sets window title to MFIX - $project_file
        #self.set_project_file(None)  - do we want to do this?


    def confirm_close(self):
        """before closing, ask user whether to end job and save project"""
        if self.job_manager.job:
            confirm = self.message(text="Stop running job?",
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm == 'yes':
                log.info("Stopping mfix at application exit")
                self.job_manager.job.stop_mfix()

        if self.unsaved_flag:
            confirm = self.message(text="Save project before quitting?",
                                   buttons=['yes', 'no', 'cancel'],
                                   default='Cancel')
            if confirm == 'yes':
                self.save_project()
            return confirm != 'cancel'
        else:
            return True

    def set_keyword(self, key, value, args=None):
        """convenience function to set keyword"""
        self.set_unsaved_flag()#?
        self.project.submit_change(None, {key:value}, args)

    def update_keyword(self, key, value, args=None):
        """like set_keyword but no action if value already set"""
        if not isinstance(value, Equation):  # Always update equations
            if self.project.get_value(key, args=args) == value:
                return
        self.set_keyword(key, value, args=args)

    def unset_keyword(self, key, args=None):
        """Undefine keyword.  Report to user, also catch and report any errors"""
        #  Note - keyword is still registered!  This method does not deregister
        #   keywords with project manager
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

    def _on_resized(self, ev):
        ui = self.ui
        w = ev.size().width()
        if w < self.max_label_len:
            self.short_labels()
        else:
            self.long_labels()

        g = ui.treewidget_navigation.geometry()

    def short_labels(self):
        tree = self.ui.treewidget_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        for (long, short) in self.nav_labels:
            items = tree.findItems(long, flags, 0)
            for item in items:
                if item.data(UserRole, 0): # Avoid toggling
                    continue
                item.setText(0, short)

    def long_labels(self):
        tree = self.ui.treewidget_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        for (long, short) in self.nav_labels:
            items = tree.findItems(short, flags, 0)
            for item in items:
                if item.data(UserRole,0): # Avoid toggling
                    continue
                item.setText(0, long)


    def unimplemented(self):
        self.message(title='Unimplemented',
                     text='Feature not implemented')

    def update_nav_tree(self):
        n_regions = len(self.ui.regions.get_region_dict()) # causes an extra copy
        enabled = (n_regions > 0) and (self.fluid_nscalar_eq > 0
                                       or self.solids_nscalar_eq > 0
                                       or len(self.solids) > 0
                                       or not(self.fluid_solver_disabled) )
        for name in ("Initial Conditions",): # "Boundary Conditions"):
            item = self.find_navigation_tree_item(name)
            item.setDisabled(not enabled)

    def toggle_nav_menu(self):
        nav_menu = self.ui.treewidget_navigation
        nav_menu.setVisible(not nav_menu.isVisible())

    def status_message(self, message=''):
        self.ui.label_status.setText(message)

    def slot_rundir_changed(self):
        # Note: since log files get written to project dirs, this callback
        # is triggered frequently during a run.
        runname = self.get_runname()
        runname_mfx, runname_pid = runname + '.mfx', runname + '.pid'
        if self.get_project_dir() and not self.job_manager.job:
            log.debug("SLOT_RUNDIR_CHANGED was called")
            full_runname_pid = os.path.join(self.get_project_dir(), runname_pid)
            self.job_manager.try_to_connect(full_runname_pid)

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

        # Selected rows look bad with input disabled
        if not enabled:
            self.ui.fluid.tablewidget_fluid_species.clearSelection()
            self.ui.solids.tablewidget_solids.clearSelection()
            self.ui.solids.tablewidget_solids_species.clearSelection()
        else:
            self.update_fluid_species_table() # reenable sel?
            self.update_solids_table() # ?
            #self.update_solids_detail_pane()

    def slot_update_residuals(self):
        """Get job status from JobManager and update residuals pane"""
        #if not self.job_manager:
        #    return
        if self.job_manager and self.job_manager.job:
            log.debug('update_residuals')
            self.ui.residuals.setText(self.job_manager.job.pretty_status)
        else:
            log.debug('no Job object (update_residuals)')

    # TODO:  separate this into different functions - this is called by
    # several different signals for different reasons
    # This function is called a lot, and it does too much work each time
    # 1) executables changed
    # 2) project directory changed
    # 3) process started
    # 4) process stopped
    def slot_update_runbuttons(self, message=None):
        """Updates list of of mfix executables and sets run dialog options"""
        # This is the main state-transition handler

        if message is not None:
            # highlight for visibility, this is an important state chage
            self.print_internal(message, color='blue')

        # TODO: set this in __init__ or another early setup method
        # assemble list of available executables

        ui = self.ui
        project_file = os.path.basename(self.get_project_file() or '')

        log.debug('accessing job_manager object: %s (slot_update_runbuttons)', self.job_manager)
        if self.job_manager.job:
            log.debug('accessing job_manager.job object: %s (slot_update_runbuttons)', self.job_manager.job)
        project_open = bool(project_file and self.open_succeeded)
        pending = self.job_manager.is_job_pending()
        # why both paused and unpaused states?
        paused = self.job_manager.job and self.job_manager.job.is_paused()
        unpaused = self.job_manager.job and not paused
        resumable = bool(self.monitor.get_res_files()) and not self.job_manager.job
        editable = project_open and not (pending or paused or unpaused or resumable)

        log.debug("UPDATE RUN OPTIONS: pending=%s paused=%s resumable=%s",
                   pending, paused, resumable)

        self.update_window_title() # put run state in window titlebar

        self.enable_input(editable)
        self.ui.run.setEnabled(project_open)

        #handle buttons in order:  RESET RUN PAUSE STOP
        pause_visible = bool(self.job_manager.job)
        if pending:
            self.status_message("MFIX starting up, process %s" % self.job_manager.job.mfix_pid)
            # also disable spinboxes for dt, tstop unless interactive
            self.set_reset_button(enabled=False)
            self.set_run_button(enabled=False)
            self.set_pause_button(text="Pause", enabled=False, visible=pause_visible)
            self.set_stop_button(enabled=True)
            self.change_pane('run')

        elif unpaused:
            self.status_message("MFIX running, process %s" % self.job_manager.job.mfix_pid)
            # also disable spinboxes for dt, tstop unless interactive
            self.set_reset_button(enabled=False)
            self.set_run_button(enabled=False)
            self.set_pause_button(text="Pause", enabled=True, visible=pause_visible)
            self.set_stop_button(enabled=True)
            self.change_pane('run')

        elif paused:
            self.status_message("MFIX paused, process %s" % self.job_manager.job.mfix_pid)
            self.set_reset_button(enabled=False)
            self.set_run_button(text="Unpause", enabled=True)
            self.set_pause_button(text="Pause", enabled=False, visible=pause_visible)
            self.set_stop_button(enabled=True)
            self.change_pane('run')

        elif resumable:
            self.status_message("Previous MFIX run is resumable.  Reset job to edit model")
            self.set_reset_button(enabled=True)
            self.set_run_button(text='Resume', enabled=True)
            self.set_pause_button(text="Pause", enabled=False, visible=pause_visible)
            self.set_stop_button(enabled=False)
            self.change_pane('run')

        else: # Not running, ready for input
            self.status_message("Ready" if project_open else "Loading %s"%project_file)
            self.set_reset_button(enabled=False)
            self.set_run_button(text="Run", enabled=project_open)
            self.set_pause_button(text="Pause", enabled=False, visible=pause_visible)
            self.set_stop_button(enabled=False)

        ui.run.use_spx_checkbox.setEnabled(resumable)
        ui.run.use_spx_checkbox.setChecked(resumable)
        ui.run.checkbox_pymfix_output.setEnabled(True)

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
        tree = self.ui.treewidget_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        items = tree.findItems(item_name, flags, 0)
        if len(items) == 1:
            return items[0]
        else:
            for (long,short) in self.nav_labels:
                if (item_name == long):
                    items = tree.findItems(short, flags, 0)
                    if len(items) == 1:
                        return items[0]


    # Top-level "Model"
    def set_solver(self, solver):
        """handler for "Solver" combobox in Model pane"""
        self.project.solver = solver
        if solver is None: #
            return

        ui = self.ui
        m = ui.model_setup
        cb = m.combobox_solver
        if cb.currentIndex != solver:
            cb.setCurrentIndex(solver)

        solver_name = {SINGLE:"MFIX Single-Phase",
                       TFM:"MFIX-TFM",
                       DEM:"MFIX-DEM",
                       PIC:"MFIX-PIC",
                       HYBRID:"MFIX-Hybrid"}.get(solver, "MFIX")

        self.print_internal("Solver: %s" % solver_name)
        self.solver_name = solver_name

        enabled = (solver != SINGLE)
        self.find_navigation_tree_item("Solids").setDisabled(not enabled)
        if enabled:
            self.solids_update_tabs()

        # Options which require TFM, DEM, or PIC
        enabled = solver in (TFM, DEM, PIC)
        interphase = m.groupbox_interphase
        interphase.setEnabled(enabled)

        # TFM only
        # use a groupbox here, instead of accessing combobox + label?
        enabled = (solver == TFM)
        m.combobox_subgrid_model.setEnabled(enabled)
        m.label_subgrid_model.setEnabled(enabled)
        m.groupbox_subgrid_params.setEnabled(enabled and
                                                       self.subgrid_model > 0)

        enabled = (self.fluid_nscalar_eq > 0)
        ui.fluid.checkbox_enable_scalar_eq.setChecked(self.fluid_nscalar_eq>0)
        ui.fluid.spinbox_nscalar_eq.setEnabled(enabled)

        # Equiv for solids is done in update_solids_detail_pane

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
        self.update_solids_table() # some of these settings are dependent on solver
        self.setup_combobox_solids_model()
        #self.update_solids_detail_pane()
        self.update_window_title()

    def disable_fluid_solver(self, disabled):
        self.fluid_solver_disabled = disabled
        m = self.ui.model_setup
        enabled = not disabled
        item = self.find_navigation_tree_item("Fluid")
        item.setDisabled(disabled)
        if disabled:
            self.enable_turbulence(False)
        m.checkbox_enable_turbulence.setEnabled(enabled)
        m.combobox_turbulence_model.setEnabled(enabled and
                        m.checkbox_enable_turbulence.isChecked())

        self.ics_update_enabled()
        self.update_nav_tree()
        # TODO update nscalar

    def enable_energy_eq(self, enabled):
        # Additional callback on top of automatic keyword update,
        # since this has to change availabilty of several other GUI items

        # TODO : move this to fluids_handler, like the way we do for solids

        self.ui.model_setup.checkbox_keyword_energy_eq.setChecked(enabled)

        # It might not be necessary to do all this - will the fluid or
        # solid panes get updated before we display them?
        ui = self.ui
        f = ui.fluid
        for item in (f.label_fluid_specific_heat_model,
                     f.combobox_fluid_specific_heat_model,
                     f.label_fluid_conductivity_model,
                     f.combobox_fluid_conductivity_model,
                     # more ?
                     ):
            item.setEnabled(enabled)

        # c_pg0 == specific heat for fluid phase
        lineedit = f.lineedit_keyword_c_pg0
        label = f.label_c_pg0_units
        for item in (lineedit, label):
            item.setEnabled(enabled and (self.fluid_specific_heat_model == CONSTANT))

        # k_g0 == thermal conductivity fluid phase
        lineedit = f.lineedit_keyword_k_g0
        label = f.label_k_g0_units
        for item in (lineedit, label):
            item.setEnabled(enabled and (self.fluid_conductivity_model == CONSTANT))


    def set_subgrid_model(self, index):
        self.subgrid_model = index
        groupbox_subgrid_params = self.ui.model_setup.groupbox_subgrid_params
        groupbox_subgrid_params.setEnabled(index > 0)

    def enable_turbulence(self, enabled):
        m = self.ui.model_setup
        if enabled != m.checkbox_enable_turbulence.isChecked():
            m.checkbox_enable_turbulence.setChecked(enabled)
        for item in (m.label_turbulence_model,
                     m.combobox_turbulence_model,
                     m.label_mu_gmax,
                     m.lineedit_keyword_mu_gmax,
                     m.label_mu_gmax_units):
            item.setEnabled(enabled)
        if not enabled:
            self.unset_keyword('turbulence_model')
            self.unset_keyword('mu_gmax')
        else:
            self.set_turbulence_model(m.combobox_turbulence_model.currentIndex())

    def set_turbulence_model(self, val):
        m = self.ui.model_setup
        cb = m.combobox_turbulence_model
        if cb.currentIndex() != val:
            cb.setCurrentIndex(val)
            return
        self.update_keyword('turbulence_model',
                            ['MIXING_LENGTH', 'K_EPSILON'][val])
        val = m.lineedit_keyword_mu_gmax.value
        if val=='' or val is None:
            val = '1.e+03'
        self.update_keyword('mu_gmax', val)


    def update_scalar_equations(self, prev_nscalar):
        """sets nscalar and phase4scalar(#) for all phases"""
        # Used by both fluid & solid han
        # We may have reduced nscalar, so we need to unset phase4scalar(i)
        # for any values of i > nscalar.
        nscalar = self.fluid_nscalar_eq + self.solids_nscalar_eq
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
            # TODO implement a way to unset keys with wildcard
            for ic in range(1, 1+len(self.ics)):
                self.unset_keyword('ic_scalar', args=[ic, i])

            i += 1

        # ICs enabled/disabled depends on nscalar
        self.ics_update_enabled()
        self.update_nav_tree()

    # helper functions for __init__
    def __setup_other_widgets(self): # rename/refactor
        """setup widgets which are not tied to a simple keyword"""
        ui = self.ui
        model = ui.model_setup
        combobox = model.combobox_solver
        # activated: Only on user action, avoid recursive calls in set_solver
        combobox.activated.connect(self.set_solver)

        checkbox = model.checkbox_disable_fluid_solver
        checkbox.clicked.connect(self.disable_fluid_solver)
        self.disable_fluid_solver(self.fluid_solver_disabled)

        checkbox = model.checkbox_keyword_energy_eq
        checkbox.clicked.connect(self.enable_energy_eq)

        checkbox = model.checkbox_enable_turbulence
        checkbox.clicked.connect(self.enable_turbulence)

        combobox = model.combobox_turbulence_model
        combobox.currentIndexChanged.connect(self.set_turbulence_model)

        combobox = model.combobox_subgrid_model
        combobox.currentIndexChanged.connect(self.set_subgrid_model)
        self.set_subgrid_model(0)

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
                    args = [try_int(name) for name in name_list[args_idx+1:]]
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
                    try:
                        widget.setdtype(doc['dtype'])
                    except:
                        print(widget, widget.objectName())
                        raise
                    vr = doc.get('validrange', {})
                    widget.setValInfo(min=vr.get("min"), max=vr.get("max"),
                                      required=doc.get("required"))

                    default = doc.get('initpython') # "Initial Python Value"
                    if default is not None:
                        widget.default(default)

                    description = doc.get('description')
                    if description is not None:
                        widget.help_text = description
                        widget.setToolTip('%s<br/>' % doc['description'])

                    if isinstance(widget, QtWidgets.QLineEdit) and widget.dtype in [int, float]:
                        widget.allow_parameters = True

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

        from widgets.vtkwidget import VtkWidget
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
        self.capture_output(mode=='interpreter')
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
        """set current pane to the one matching 'name'.  Must be the long
        (non-abbreviated) navigation label.  Case-insensitive"""

        items = self.ui.treewidget_navigation.findItems(
                    name,
                    Qt.MatchFixedString | Qt.MatchRecursive, 0)
        if not items: # Nav menu may be in abbreviated mode.  Might be better
            # to identify navigation items by something other than text, since
            # that can change (long/short) and is possibly non-unique (eg "points")
            for (long, short) in self.nav_labels:
                if name.lower() == long.lower():
                    items = self.ui.treewidget_navigation.findItems(
                        short,
                        Qt.MatchFixedString | Qt.MatchRecursive, 0)
                    if items:
                        break
        assert len(items) == 1
        self.ui.treewidget_navigation.setCurrentItem(items[0])
        self.navigation_changed()

    def navigation_changed(self):
        """an item in the tree was selected, change panes"""

        current_selection = self.ui.treewidget_navigation.selectedItems()
        if not current_selection:
            return
        text = str(current_selection[-1].text(0))
        # Translate from short to long name
        for (long, short) in self.nav_labels:
            if text==short:
                text=long
                break
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

        self.setup_current_tab()

    def setup_current_tab(self):
        # Force any open popup to close
        # (if dialog is modal we don't need this)
        self.species_popup.done(0)
        self.regions_popup.done(0)
        current_selection = self.ui.treewidget_navigation.selectedItems()
        if not current_selection:
            return
        text = str(current_selection[-1].text(0))
        text = '_'.join(text.lower().split(' '))
        if text == 'solids': # Special helper for setting up subpanes,
            # since params may have changed
            self.setup_solids_tab(self.solids_current_tab)
        # TODO: improve handling of long/short names
        elif text in ('initial_conditions', 'ics'):
            self.setup_ics()
        elif text in ('boundary_conditions', 'bcs'):
            self.setup_bcs()


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
            offsetx = 0
            offsety = 0

        self.stack_animation = QtCore.QParallelAnimationGroup()

        if (from_ != to):
            # move to widget and show
            # set the geometry of the next widget
            to_widget.setGeometry(0 + offsetx, 0 + offsety, width, height)
            to_widget.show()
            to_widget.raise_()
            #to_widget.activateWindow() ? needed?


            # animate
            # from widget
            self.animation_setup(from_widget, 0, 0, -offsetx, -offsety)

            # to widget
            self.animation_setup(to_widget, offsetx, offsety, 0, 0)

        # line
        line_to = None
        if line is not None and to_btn is not None:
            self.animation_setup(line, line.geometry().x(), line.geometry().y(), to_btn.geometry().x(), line.geometry().y())
            if btn_layout is not None:
                for i in range(0, btn_layout.columnCount()):
                    if btn_layout.itemAtPosition(0,i) == to_btn:
                        line_to = i
                        break

        # animation group FIXME why 2 callbacks?
        self.stack_animation.finished.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, to_btn, line, line_to))

        self.stack_animation.stateChanged.connect(
            make_callback(self.animate_stacked_widget_finished,
                          stackedwidget, from_, to, btn_layout, to_btn, line, line_to))

        self.animating = True
        self.stack_animation.start()

    def animation_setup(self, target, x_start, y_start, x_end, y_end):
        """setup an animation widget"""
        animation = QtCore.QPropertyAnimation(target, "pos".encode('utf-8'))
        animation.setDuration(self.animation_speed)
        animation.setEasingCurve(QtCore.QEasingCurve.InOutQuint)
        animation.setStartValue(QtCore.QPoint(x_start, y_start))
        animation.setEndValue(QtCore.QPoint(x_end,y_end))
        self.stack_animation.addAnimation(animation)

    def animate_stacked_widget_finished(self, widget, from_, to,
                                        btn_layout=None, to_btn=None, line=None, line_to=None):
        """cleanup after animation"""
        try:
            if self.stack_animation.state() == QtCore.QAbstractAnimation.Stopped:
                widget.setCurrentIndex(to)
                if (from_ != to):
                    from_widget = widget.widget(from_)
                    from_widget.hide()
                    from_widget.move(0, 0) #why?
                if btn_layout is not None and line is not None:
                    btn_layout.addItem(btn_layout.takeAt(
                        btn_layout.indexOf(line)), 1, line_to or to)
        except AttributeError: # Happens during shutdown. TODO: remove this hack
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
        icon = 'question' 'warning' or 'info'
        text = 'test to show'
        buttons = ['ok',...] where value is 'ok', 'yes', 'no', 'cancel',
            'discard'
        default = 'ok' the default selected button
        infoText = 'extended information text'
        detailedText = 'Some details'

        Returns the pressed button.  Also prints & logs message"""

        self.print_internal(title + ": " + text)
        if infoText:
            self.print_internal(infoText)
        if detailedText:
            self.print_internal(detailedText)

        message_box = QtWidgets.QMessageBox(self)
        self.message_box = message_box # Make it accessible to tests
        message_box.setWindowTitle(title)

        # Icon
        if icon == 'warning':
            icon = QtWidgets.QMessageBox.Warning
        elif icon == 'question':
            icon = QtWidgets.QMessageBox.Question
        else:
            icon = QtWidgets.QMessageBox.Information

        message_box.setIcon(icon)

        # Text
        message_box.setText(text)

        if infoText:
            message_box.setInformativeText(infoText)

        if detailedText:
            message_box.setDetailedText(detailedText)

        # TODO: support more standard buttons from #
        #   http://doc.qt.io/qt-4.8/qmessagebox.html#StandardButton-enum
        # or drop this translation dict
        # buttons
        qbuttonDict = {'ok':      QtWidgets.QMessageBox.Ok,
                       'yes':     QtWidgets.QMessageBox.Yes,
                       'no':      QtWidgets.QMessageBox.No,
                       'cancel':  QtWidgets.QMessageBox.Cancel,
                       'discard': QtWidgets.QMessageBox.Discard,
                       }
        for button in buttons:
            message_box.addButton(qbuttonDict[button])

            if button == default:
                message_box.setDefaultButton(qbuttonDict[button])
        ret = message_box.exec_()

        for key, value in qbuttonDict.items():
            if value == ret:
                return key


    def scan_errors(self, lines):
        ### "Error 1000: A keyword pair on line 129"
        ### "Error 2000: Unable to process line 185"
        # TODO: capture more of the error text and produce a fuller message
        # in the popup
        lineno = bad_line = err = None
        re_err_1000 = re.compile("Error 1000: A keyword pair on line (\d+)")
        re_err_2000 = re.compile("Error 2000: Unable to process line (\d+)")
        for (i, line) in enumerate(lines):
            for (re_err, err_type) in ((re_err_1000, 'deprecated'),
                                       (re_err_2000, 'invalid')):
                match = re_err.search(line)
            if match:
                try:
                    lineno = int(match.group(1))
                except ValueError:
                    return
                try:
                    bad_line = self.datfile_lines[lineno-1]
                except IndexError:
                    return
                break
        # TODO:  colorize errors in source view (red)
        if bad_line:
            try:
                key = bad_line.split("=", 1)[0].strip()
                self.report_keyword_error(key, bad_line, err_type)
            except:
                pass # Don't introduce additional errors in error handler


    def report_keyword_error(self, key, line, err_type='deprecated'):
        """Give the user a chance to omit or edit deprecated keys"""
        #  This is a first implementation.  There's a lot more
        #  we could do here - suggest fixes, link to documentation,
        #  attempt to retain order in dat_file_list, etc.
        key = key.lower()
        line = line.strip()

        message_box = QtWidgets.QMessageBox(self)
        self.message_box = message_box
        message_box.setWindowTitle("%s keyword" % err_type.title())
        message_box.setIcon(QtWidgets.QMessageBox.Warning)
        text="'%s' is %s" % (key, err_type)
        message_box.setText(text)
        buttons = ['Drop Key', 'Edit', 'Ignore']
        for b in buttons:
            role = QtWidgets.QMessageBox.AcceptRole # Seems to be OK to use for all buttons
            message_box.addButton(QtWidgets.QPushButton(b),  role)

        def drop(key):
            if '(' in key: # This is a little crude, use parser instead?
                k, a = key.split('(', 1)
                a = a.split(')')[0]
                a = [int(x) for x in a.split(',')]
            else:
                a = None
            self.unset_keyword(k,a)

        resp = buttons[message_box.exec_()].lower()
        if not resp or resp=='ignore': # User bailed out
            return

        elif resp == 'drop key':
            drop(key)

        elif resp == 'edit':
            q = QtWidgets.QInputDialog(self)#,
            text, ok = q.getText(self, 'Edit keyword', "Deprecated keyword: %s\n\nEnter replacement text"%line, text=line)
            if not ok:
                return
            text = text.strip()
            if not text:
                drop(key) #User replaced it with blank
            for (new_key, new_args, new_value) in self.project.parseKeywordLine(text):
                if new_key:
                    drop(key) # Only drop the old one once we have a valid replacemnent
                    self.update_keyword(new_key, new_value, args=new_args) # validate key?
                else:
                    self.print_internal("Error:  cannot parse %s" % text)



    def can_skip(self, line,
                 boilerplate=set(['Program Terminated.',
                                  'Fatal error reported on one or more processes. The .LOG file',
                                  'may contain additional information about the failure.',
                                  'ERROR STOP 1',
                                  'Please see the user documentation and update the mfix.dat file.',
                                  ])):

        # These are routine messages that we are not going to trouble the user with
        # Note, the set is only constructed once (at load time)
        # Also skip lines containing only '*'
        stripped = line.strip()
        return all(c=='*' for c in stripped) or stripped in boilerplate


    def handle_stdout(self, text):
        """collect stderr from mfix/pymfix process, format and display
        to user.  Note that some errors are printed to stdout"""
        color = 'red' if "Error" in text else None
        lines = text.split('\n')
        # Scanning for errors may trigger popups, etc, so get the output to
        # the user first.
        for line in lines:
            if self.can_skip(line):
                continue

            self.print_internal(line, color=color, font='Courier')
        self.scan_errors(lines)

    def handle_stderr(self, text):
        """collect stderr from mfix/pymfix process, format and display
        to user."""
        # Scanning for errors may trigger popups, etc, so get the output to
        # the user first.

        lines = text.split('\n')
        for line in lines:
            if self.can_skip(line):
                continue
            self.print_internal(line, color='red', font='Courier') # Bold font?
        self.scan_errors(lines)

    def print_internal(self, line, color=None, font=None):
        qtextbrowser = self.ui.command_output
        stripped = line.strip()
        if not stripped:
            # Let's just skip blank lines completely
            return
        if not line.endswith('\n'):
            line += '\n'
        lower = line.lower()
        logmsg = stripped
        # hack. TODO: real msg types, map to font/color
        strikeout = font and font.lower() == 'strikeout'
        if strikeout:
            logmsg = "unset " + logmsg
        if lower.startswith("error:"):
            log.error(logmsg[6:])
        elif lower.startswith("warning:"):
            log.warning(logmsg[8:])
        elif lower.startswith("info:"):
            log.info(logmsg[5:])
            color='blue'
        else:
            log.info(logmsg)
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
        scrollbar = qtextbrowser.verticalScrollBar()
        scrolled_to_end = (scrollbar.value() == scrollbar.maximum())
        cursor.insertText(line)
        if scrolled_to_end:
            scrollbar.setValue(scrollbar.maximum())


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
                msg = 'Error: cannot delete %s: %s' % (path, err.strerror)
                self.message(title='File error',
                             icon='error',
                             text=msg,
                             buttons=['ok'],
                             default=['ok'])
                break
        self.signal_update_runbuttons.emit('')
        return True

    # Don't make these depend on current state, since (esp for pymfix)
    # the state variables are cached and potentially outdated
    def handle_run(self):
        # name?
        name = 'Run'
        try:
            if not self.job_manager.job:
                # open the run dialog for job options
                self.open_run_dialog()
                return
            else:
                name='unpause' #?
                self.job_manager.job.unpause()
        except Exception:
            log.exception('problem in handle_run')
            self.print_internal('problem in handle_run')
        self.signal_update_runbuttons.emit('')

    def handle_set_pymfix_output(self):
        try:
            self.job_manager.job.set_pymfix_output(
              self.ui.run.checkbox_pymfix_output.isChecked())
        except Exception:
            log.exception('problem in handle_set_pymfix_output')
            self.print_internal('problem in handle_set_pymfix_output')

    def handle_pause(self):
        try:
            self.job_manager.job.pause()
        except Exception:
            log.exception('problem in handle_pause')
            self.print_internal('problem in handle_pause')

    def handle_stop(self):
        try:
            self.job_manager.stop_mfix()
        except Exception:
            log.exception('problem in handle_stop')
            self.print_internal('problem in handle_stop')

    def check_save(self):
        if self.unsaved_flag:
            response = self.message(title="Save?",
                                    icon="question",
                                    text="Save current project?",
                                    buttons=['yes', 'no'])
            if response != 'yes':
                return False
        #FIXME need to catch/report errors, writeDatFile is too low-level
        self.project.writeDatFile(self.get_project_file())
        self.clear_unsaved_flag()
        self.update_source_view()
        return True

    def open_run_dialog(self):
        """Open run popup dialog"""
        if not self.check_save():
            return
        self.run_dialog = RunPopup(self.commandline_option_exe, self)
        self.run_dialog.set_run_mfix_exe.connect(self.handle_exe_changed)
        self.run_dialog.label_cores_detected.setText("Running with %d cores" % multiprocessing.cpu_count())
        self.run_dialog.setModal(True)
        self.run_dialog.popup()

    def handle_exe_changed(self):
        """callback from run dialog when combobox is changed"""
        self.mfix_exe = self.run_dialog.mfix_exe
        self.settings.setValue('mfix_exe', self.mfix_exe)
        log.debug('exe changed signal recieved: %s' % self.mfix_exe)
        self.signal_update_runbuttons.emit('')

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
            resp = self.message(text="Copy .SP files?\n%s" % '\n'.join(sp_files),
                               icon='question',
                               buttons=['yes', 'no'],
                               default='yes')
            if resp=='yes':
                files_to_copy.extend(sp_files)

        files_to_copy.extend(self.monitor.get_outputs(["*.RES", "*.STL"]))
        #copy project files into new_project_directory
        for f in files_to_copy:
            try:
                shutil.copyfile(f, os.path.join(export_dir, os.path.basename(f)))
            except Exception as e:
                self.message(title='File error',
                             text="Error copying file:\n%s" % e,
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
        self.project.mfix_gui_comments['geometry'] = self.vtkwidget.geometry_to_str()

        # save regions
        self.project.mfix_gui_comments['regions_dict'] = self.ui.regions.regions_to_str()

        # save IC regions
        if self.ics:
            self.project.mfix_gui_comments['ic_regions'] = self.ics_to_str()
        else:
            self.project.mfix_gui_comments.pop('ic_regions', None)

        project_base = os.path.basename(project_file)
        run_name = os.path.splitext(project_base)[0]
        self.update_keyword('run_name', run_name)
        self.print_internal("Info: Saving %s" % project_file)
        self.project.writeDatFile(project_file)

        # save workflow
        if PYQTNODE_AVAILABLE:
            self.ui.workflow_widget.save(
                os.path.abspath(os.path.join(project_dir, 'workflow.nc')))

        self.clear_unsaved_flag()
        self.update_source_view()

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
        if PYQT5:
            filename = filename[0]
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

#    def handle_settings(self):
#        """handle user settings"""
#        # TODO: implement
#        self.unimplemented()

    def handle_parameters(self):
        """add/change parameters"""
        changed_params = self.parameter_dialog.get_parameters()
        self.set_unsaved_flag()
        self.update_parameters(changed_params)

    def update_parameters(self, changed_params):
        """update the changed parameters"""
        self.ui.regions.update_parameters(changed_params)
        self.vtkwidget.update_parameters(changed_params)
        self.project.update_parameters(changed_params)

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

        if self.job_manager.job:
            if not self.job_manager.job.is_paused():
                title += ', RUNNING'
                if self.job_manager.job.mfix_pid is not None:
                    title += ', process %s'% self.job_manager.job.mfix_pid
            elif self.job_manager.job.is_paused():
                title += ', PAUSED'
            elif self.monitor.get_res_files():
                title += ', STOPPED, resumable'
        else:
            pass
            #title += ', EDITING'
        self.setWindowTitle(title)

    def set_save_button(self, enabled):
        self.ui.toolbutton_save.setEnabled(enabled)

    def set_save_as_action(self, enabled):
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

    def new_project(self):
        # Why not use get_save_filename here?
        project_file = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Create Project File',
            "", "MFX files (*.mfx)")

        if PYQT5:
            project_file = project_file[0]
        if not project_file:
            return

        if not project_file.endswith(".mfx"):
            project_file = project_file + ".mfx"

        project_filename = os.path.basename(project_file)
        project_dir = os.path.dirname(project_file)
        run_name = os.path.splitext(project_filename)[0]
        if not self.check_writable(project_dir):
            self.message(text='Unable to write to %s' % project_dir,
                         buttons=['ok'],
                         default='ok')
            return
        # Start with a nice template - note, there's too much set in this file.
        # FIXME, this can clobber files
        template = os.path.join(get_mfix_home(), 'gui', 'mfix.dat.template')
        shutil.copyfile(template, project_file)
        self.open_project(project_file, auto_rename=False)
        self.project.run_name.updateValue(run_name)
        self.save_project()

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

        self.datfile_lines = src.split('\n')

        if number_lines:
            lines = src.split('\n')
            # Avoid extra blank lines at end
            while lines and lines[-1] == '':
                lines.pop(-1)
            src = '\n'.join('%4d:%s'%(lineno,line)
                            for (lineno,line) in enumerate(lines,1))

        self.ui.mfix_dat_source.setPlainText(src)

    def force_default_settings(self):
        # Should these just be in the template file? TODO
        self.update_keyword('chk_batchq_end', True)

    def get_runname(self):
        name = self.project.get_value('run_name', default='new_file')
        for char in ('.', '"', "'", '/', '\\', ':'):
            name = name.replace(char, '_')
        return name

    def open_project(self, project_path, auto_rename=True):
        """Open MFiX Project"""

        self.open_succeeded = False  # set to true on success
        self.vtkwidget.defer_render = True # defer rendering vtk until load finished

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

        # --- read the mfix.dat or *.mfx file
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

        runname = self.get_runname()
        runname_mfx, runname_pid = runname + '.mfx', runname + '.pid'

        if os.path.exists(runname_pid):
            # previously started job may be running, try to reconnect
            log.debug('attempting to connect to running job %s' % runname_pid)
            self.job_manager.try_to_connect(runname_pid)

        if auto_rename and not project_path.endswith(runname_mfx):
            ok_to_write = False
            save_msg = 'Renaming %s to %s based on run name' % (project_file, runname_mfx)
            response = self.message(title='Info',
                                    icon='question',
                                    text=save_msg,
                                    buttons=['ok', 'cancel'],
                                    default='ok')
            if response == 'ok':
                ok_to_write = True
                renamed_project_file = os.path.join(project_dir, runname_mfx)
                if os.path.exists(renamed_project_file):
                    clobber_msg = '%s exists, replace?' % renamed_project_file
                    response = self.message(title='Warning',
                                    icon='question',
                                    text=clobber_msg,
                                    buttons=['yes', 'no'],
                                    default='no')
                    if response == 'no':
                        ok_to_write = False
                        log.info("Rename cancelled at user request")

            if ok_to_write:
                project_file = renamed_project_file
                try:
                    self.force_default_settings()
                    self.print_internal("Info: Saving %s" % project_file)
                    self.project.writeDatFile(project_file) #XX
                    #self.print_internal(save_msg, color='blue')
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
            else:
                self.print_internal("Rename canceled at user request")

        self.set_project_file(project_file)
        self.clear_unsaved_flag()
        self.set_save_as_action(enabled=True)

        self.setup_current_tab() # update vals in any open tabs
        self.update_source_view()

        # set up rundir watcher
        self.rundir_watcher.addPath(project_dir)
        self.slot_rundir_changed()

        ### Geometry
        # Look for geometry.stl and load automatically
        geometry_file = os.path.abspath(os.path.join(project_dir, 'geometry.stl'))
        if os.path.exists(geometry_file) and 'geometry' not in self.project.mfix_gui_comments:
            self.vtkwidget.add_stl(None, filename=geometry_file)
            msg = '%s will be overwritten when the project is saved' % os.path.basename(geometry_file)
            self.message(title='Warning', text=msg)

        # TODO: load more geometry

        # Additional GUI setup based on loaded projects (not handled
        # by keyword updates)
        #    .... is there a way to verify that 'energy_eq' is boolean?
        #    should that get set from keyword doc?

        # Note that energy_eq is TRUE by default according to MFIX but
        # FALSE by default according to SRS document.  We have to respect
        # MFIX.  The GUI defaults are enforced in the template file used
        # to create new projects.
        self.enable_energy_eq(bool(self.project.get_value('energy_eq', default=True)))

        # cgw - lots more model setup todo here.  Should we do this here or
        #  in ProjectManager.load_project_file (where we do guess/set_solver)
        # make sure exceptions are handled & reported
        # values that don't map to keywords, saved as #!MFIX-GUI params

        # Model parameters
        # TODO: set 'disable fluid solver' checkbox if appropriate

        turbulence_model = self.project.get_value('turbulence_model')
        if turbulence_model is None:
            self.enable_turbulence(False)
        else:
            turbulence_models = ['MIXING_LENGTH', 'K_EPSILON']
            if turbulence_model not in turbulence_models:
                self.message("Error: Invalid turbulence model %s" % turbulence_model)
                self.unset_keyword('turbulence_model')
            else:
                # set model first to avoid extra keyword settings
                self.set_turbulence_model(turbulence_models.index(turbulence_model))
                self.enable_turbulence(True)


        #  Non-keyword params stored as !#MFIX_GUI comments
        # TODO: add exception handling here
        solids_phase_names = {}
        for (key, val) in self.project.mfix_gui_comments.items():
            if key == 'fluid_phase_name':
                self.set_fluid_phase_name(val)
            if key.startswith('solids_phase_name('):
                n = int(key.split('(')[1][:-1])
                solids_phase_names[n] = val
            if key == 'regions_dict':
                self.ui.regions.regions_from_str(val)
            if key == 'geometry':
                self.vtkwidget.geometry_from_str(val)
            if key == 'ic_regions':
                self.ics_regions_from_str(val)
            # Add more here

        # Copy ordered dict to modify keys w/o losing order
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
        self.solids_nscalar_eq = sum(1 for i in range(1, nscalar+1)
                                    if self.project.get_value('phase4scalar', args=i) != 0)


        # solid scalar eq checkbox will be handled in update_solids_detail_pane
        # need to initialize per-solid 'nscalar_eq' and 'saved_nscalar_eq' so
        # that update_scalar_equations will do the right thing
        for (i, (k,v)) in enumerate(self.solids.items(), 1):
            v['nscalar_eq'] = sum(1 for j in range(1, nscalar+1)
                                  if self.project.get_value('phase4scalar', args=j) == i)

            v['saved_nscalar_eq'] = v['nscalar_eq']

        # setValue will trigger update_scalar_equations (ok)
        self.ui.fluid.spinbox_nscalar_eq.setValue(self.fluid_nscalar_eq)
        self.enable_fluid_scalar_eq(self.fluid_nscalar_eq > 0)


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
        self.solids_update_tabs()
        #self.update_solids_detail_pane()

        ### Regions
        # Look for regions in IC, BC, PS, etc.
        self.ui.regions.extract_regions(self.project)
        self.ics_extract_regions()

        # Take care of updates we deferred during extract_region
        self.ui.regions.tablewidget_regions.fit_to_contents()

        ### Workflow
        workflow_file = os.path.abspath(os.path.join(project_dir, 'workflow.nc'))
        if PYQTNODE_AVAILABLE and os.path.exists(workflow_file):
            self.ui.workflow_widget.clear()
            self.ui.workflow_widget.load(workflow_file)

        # FIXME: is this a good idea?  it means we can't open a file read-only
        #self.force_default_settings()
        #if self.unsaved_flag: # Settings changed after loading
        #    self.save_project()

        self.vtkwidget.reset_view()
        self.vtkwidget.render(defer_render=False)
        self.open_succeeded = True
        self.signal_update_runbuttons.emit('')

def Usage(name):
    print("""Usage: %s [directory|file] [-h, --help] [-l, --log=LEVEL] [-q, --quit]
    directory: open mfix.dat file in specified directory
    file: open mfix.dat or <RUN_NAME>.mfx project file
    -h, --help: display this help message
    -e, --exe:  specify MFIX executable (full path)
    -l, --log=LEVEL: set logging level (error,warning,info,debug)
    -n, --new:  open new project (do not autoload previous)
    -q, --quit: quit after opening file (for testing)"""  % name, file=sys.stderr)
    sys.exit(1)


def main(args):
    global gui
    args = sys.argv
    name = args[0]
    try:
        opts, args = getopt.getopt(args[1:], "hqnl:e:", ["help", "quit", "new", "log=", "exe="])
    except getopt.GetoptError as err:
        print(err)
        Usage(name)

    quit_after_loading = False
    project_file = None
    new_project = False
    log_level = 'WARN'
    mfix_exe_option = None

    for opt, arg in opts:
        if opt in ("-l", "--log"):
            log_level = arg
        elif opt in ("-h", "--help"):
            Usage(name)
        elif opt in ("-q", "--quit"):
            quit_after_loading = True
        elif opt in ("-n", "--new"):
            new_project = True
        elif opt in ("-e", "--exe"):
            mfix_exe_option = arg
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
    gui = MfixGui(qapp, project_file=project_file)
    gui.show()

    if mfix_exe_option:
        #print('exe option passed: %s' % mfix_exe_option)
        gui.commandline_option_exe = mfix_exe_option

    if project_file is None and not new_project:
        # autoload last project
        project_file = gui.get_project_file()

    if project_file and os.path.isfile(project_file):
        gui.open_project(project_file, auto_rename=(not quit_after_loading))
    else:
        gui.set_no_project()

    # print number of keywords
    gui.print_internal('Registered %d keywords' %
                        len(gui.project.registered_keywords))

    # have to initialize vtk after the widget is visible!
    gui.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    # This makes it too easy to skip the exit confirmation dialog.  Should we allow it?
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    if not quit_after_loading:
        qapp.exec_()

    qapp.deleteLater()
    sys.exit()

if __name__  == '__main__':
    main(sys.argv)
