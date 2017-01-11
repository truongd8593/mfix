#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, unicode_literals, division

"""MFIX GUI"""
__version__ = [2017, 1, 0, 'a']
__version_str__ = '.'.join([str(i) for i in __version__])

import argparse
import glob
import logging
import os
import re
import shutil
import signal
import sys
import time
from collections import OrderedDict
import json
import traceback
import copy

# Initialize logger early
log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

from mfixgui.tools.general import SCRIPT_DIRECTORY
sys.path.append(os.path.join(SCRIPT_DIRECTORY, 'pyqtnode'))

# import qt
from qtpy import QtCore, QtWidgets, QtGui, PYQT5
from qtpy.QtCore import Qt, QFileSystemWatcher, QSettings, Signal
UserRole = QtCore.Qt.UserRole

# TODO: add pyside? There is an issue to add this to qtpy:
# https://github.com/spyder-ide/qtpy/issues/16

PRECOMPILE_UI = False

if not PRECOMPILE_UI:
    from qtpy import uic

# local imports
from mfixgui.project_manager import ProjectManager
from mfixgui.job import JobManager
from mfixgui.monitor import Monitor

from mfixgui.widgets.base import (LineEdit, CheckBox, ComboBox, SpinBox, DoubleSpinBox,
                          Table, BaseWidget)
from mfixgui.widgets.regions import RegionsWidget

from mfixgui.widgets.species_popup import SpeciesPopup
from mfixgui.widgets.regions_popup import RegionsPopup
from mfixgui.widgets.run_popup import RunPopup
from mfixgui.widgets.parameter_dialog import ParameterDialog
from mfixgui.widgets.output_selection_popup import OutputSelectionPopup

from mfixgui.model_setup import ModelSetup
from mfixgui.fluid_handler import FluidHandler
from mfixgui.solids_handler import SolidsHandler
from mfixgui.ics import ICS
from mfixgui.bcs import BCS
from mfixgui.pss import PSS
from mfixgui.iss import ISS
from mfixgui.mesh import Mesh
from mfixgui.main_menu import MainMenu
from mfixgui.chemistry import Chemistry
from mfixgui.numerics import Numerics
from mfixgui.output import Output
from mfixgui.graphic_tabs import GraphicTabs

from mfixgui.interpreter import Interpreter

from mfixgui.tools.general import (get_icon, widget_iter,
                           is_text_string, is_unicode, get_image_path,
                           format_key_with_args, to_unicode_from_fs,
                           get_username, convert_string_to_python)

from mfixgui.tools.namelistparser import getKeywordDoc
from mfixgui.tools.keyword_args import keyword_args

from mfixgui.constants import *

if PRECOMPILE_UI:
    try:
        from uifiles.model_setup import Ui_model_setup
        from uifiles.mesh import Ui_mesh
        from uifiles.fluid import Ui_fluid
        from uifiles.solids import Ui_solids
        from uifiles.geometry import Ui_geometry
        from uifiles.gui import Ui_MainWindow
        from uifiles.initial_conditions import Ui_initial_conditions
        from uifiles.boundary_conditions import Ui_boundary_conditions
        from uifiles.point_sources import Ui_point_sources
        from uifiles.internal_surfaces import Ui_internal_surfaces
        from uifiles.chemistry import Ui_chemistry
        from uifiles.numerics import Ui_numerics
        from uifiles.output import Ui_output
        from uifiles.monitors import Ui_monitors
        from uifiles.post_processing import Ui_post_processing
        from uifiles.run import Ui_run

    except ImportError as e:
        print(e)
        print("You must compile ui files! (run 'make')")
        sys.exit(1)

SETTINGS = QSettings('MFIX', 'MFIX')
ANIMATION_SPEED = int(SETTINGS.value('animation_speed', 400))

# --- Main Gui ---

class MfixGui(QtWidgets.QMainWindow,
              ModelSetup,
              Mesh,
              FluidHandler,
              SolidsHandler,
              ICS, BCS, PSS, ISS,
              Chemistry,
              Numerics,
              Output,
              GraphicTabs,
              Interpreter,
              MainMenu):
    # Main window class for MFIX-GUI

    settings = SETTINGS

    stdout_signal = Signal(str)
    stderr_signal = Signal(str)
    signal_update_runbuttons = Signal(str)

    # Allow LineEdit widgets to report out-of-bounds values.
    def popup_value_error(self, exc):
        self.message(title='Error', text=str(exc))

    def error(self, msg, popup=False):
        # Show the user a warning & log it - use this instead of log.error
        if not popup:
            self.print_internal('Error: %s' % msg)
        else:
            self.message(title='Error',
                         icon='error',
                         text=msg,
                         buttons=['ok'],
                         default='ok')

    def warn(self, msg, popup=False):
        # Show the user a warning & log it - use instead of log.warn
        if not popup:
            self.print_internal("Warning: %s" % msg)
            # print_internal will call log.warn if message starts with "Warning"
        else:
            self.message(text=msg)
            # Will also print-internal and log
    warning = warn


    def __init__(self, app, parent=None, project_file=None, loadworkflow=True,
                 loadvtk=True):
        self.app = app
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setObjectName('mfixgui')
        if project_file is not None:
            self.set_project_file(project_file)

        LineEdit.report_value_error = self.popup_value_error

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

        # Hack -remove these once better 'restore value' framework exists
        self.saved_ro_g0 = None

        # load UI

        self.customWidgets = {'LineEdit':      LineEdit,
                              'CheckBox':      CheckBox,
                              'ComboBox':      ComboBox,
                              'DoubleSpinBox': DoubleSpinBox,
                              'SpinBox':       SpinBox,
                              'Table':         Table,
                              }

        if PRECOMPILE_UI: # Slightly faster startup (?)
            self.ui = Ui_MainWindow()
            self.ui.setupUi(self)
            self.ui.panes = []
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
                        Ui_fluid,
                        Ui_solids,
                        Ui_initial_conditions,
                        Ui_boundary_conditions,
                        Ui_point_sources,
                        Ui_internal_surfaces,
                        Ui_chemistry,
                        Ui_numerics,
                        Ui_output,
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

        else:  # not precompiled - for developers
            uifiles = os.path.join(SCRIPT_DIRECTORY, 'uifiles')
            self.ui = uic.loadUi(os.path.join(uifiles, 'gui.ui'))
            self.ui.panes = []
            self.setCentralWidget(self.ui)
            assert self is not self.ui

            for name in ('model_setup',
                         'geometry',
                         'mesh',
                         'regions',
                         'fluid',
                         'solids',
                         'initial_conditions',
                         'boundary_conditions',
                         'point_sources',
                         'internal_surfaces',
                         'chemistry',
                         'numerics',
                         'output',
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
                    except Exception as e:
                        # report which ui file it was, otherwise stack trace
                        # is too generic to be helpful.
                        print("Error loading", path)
                        raise

                # assign 'self.ui.general', etc
                setattr(self.ui, name, widget)
                self.ui.stackedWidgetTaskPane.addWidget(widget)
                self.ui.panes.append(widget)
        # end of ui loading

        # build keyword documentation from namelist docstrings
        self.keyword_doc = getKeywordDoc()

        # Fixes
        doc = self.keyword_doc['turbulence_model']
        doc['valids'] = OrderedDict((('NONE', {'note': 'No turbulence model'}),
                                     (TURBULENCE_MODELS[1], doc['valids'][TURBULENCE_MODELS[1]]),
                                     (TURBULENCE_MODELS[2], doc['valids'][TURBULENCE_MODELS[2]])))

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
                           ("Point Sources", "PSs"),
                           ("Internal Surfaces", "ISs")]

        # Set tooltips for nav tree & set a data property on
        # Monitors / Points to distinguish it
        root = tw.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            item.setToolTip(0, item.text(0))
            for j in range(item.childCount()):
                subitem = item.child(j)
                subitem.setToolTip(0, subitem.text(0))
                #if item.text(0)=='Monitors' and subitem.text(0) == 'Points':
                #    subitem.setData(UserRole, 0, True) # Mark this item

        # Intercept the resize event
        tw.resizeEvent = (lambda old_method:
                          (lambda event:
                           (self._on_resized(event),
                            old_method(event))[-1]))(tw.resizeEvent)

        # Disable items that are not yet implemented
        for name in ('Monitors',
                     'Post-processing'):

            self.find_navigation_tree_item(name).setDisabled(True)

        # Initialize popup dialogs
        self.species_popup = SpeciesPopup(QtWidgets.QDialog(), self)
        #self.species_popup.setModal(True) # ?
        self.regions_popup = RegionsPopup(QtWidgets.QDialog(), self)

        # Create project manager
        # NOTE.  it's a ProjectManager, not a Project.  But
        # ProjectManager is a subclass of Project.  Please
        # do not "fix" the code by renaming self.project to
        # self.project_manager
        self.project = ProjectManager(self, self.keyword_doc)

        # Extra setup for fluid & solids panes.  Needs to happen
        # after we create ProjectManager, because widgets get registered
        self.init_model_setup()
        self.init_fluid_handler()
        self.init_solids_handler()
        self.init_mesh()
        self.init_ics()
        self.init_bcs()
        self.init_pss()
        self.init_iss()
        self.init_chemistry()
        self.init_numerics()
        self.init_output()
        self.init_graphic_tabs()

        # In-process REPL (for development, should we enable this for users?)
        self.init_interpreter()

        # delete popup
        self.output_selection_popup = OutputSelectionPopup(self)

        # --- animation data
        self.modebuttondict = {'modeler':   self.ui.pushButtonModeler,
                               'workflow':  self.ui.pushButtonWorkflow,
                               'developer': self.ui.pushButtonDeveloper,
                               'interpreter': self.ui.pushButtonInterpreter}
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
                (ui.toolbutton_file, 'menu', self.handle_main_menu),
                (ui.toolbutton_save, 'save', self.handle_save),
                (ui.toolbutton_run_mfix, 'play', self.handle_run),
                (ui.toolbutton_pause_mfix, 'pause', self.handle_pause),
                (ui.toolbutton_stop_mfix, 'stop', self.handle_stop),
                (ui.toolbutton_reset_mfix, 'delete', self.remove_output_files),
                (ui.toolbutton_compile, 'build', self.handle_compile),
                (ui.toolbutton_parameters, 'functions', self.handle_parameters),
                ):
            button.setIcon(get_icon(icon_name+'.png'))
            button.clicked.connect(function)

        # Make sure lineedits lose focus so keywords update before save/run !!
        for button in (ui.toolbutton_run_mfix, ui.toolbutton_save):
            button.setFocusPolicy(Qt.ClickFocus)

        # Geometry toolbuttons
        geo = self.ui.geometry
        geo.toolbutton_add_geometry.setIcon(get_icon('geometry.png'))
        geo.toolbutton_add_filter.setIcon(get_icon('filter.png'))
        geo.toolbutton_wizard.setIcon(get_icon('wand.png'))
        geo.toolbutton_geometry_union.setIcon(get_icon('union.png'))
        geo.toolbutton_geometry_intersect.setIcon(get_icon('intersect.png'))
        geo.toolbutton_geometry_difference.setIcon(get_icon('difference.png'))

        # mode (modeler, workflow, developer)
        for mode, btn in self.modebuttondict.items():
            btn.pressed.connect(lambda mode=mode: self.change_mode(mode))
            btn.released.connect(lambda btn=btn:  btn.setChecked(True)) # Force button to remain in 'checked' state

        # navigation tree
        ui.toolbutton_collapse_navigation.clicked.connect(self.toggle_nav_menu)
        icon = self.icon_collapse_right = QtGui.QIcon()
        icon.addFile(get_image_path('right_light.png'), QtCore.QSize(), icon.Normal, icon.Off)
        icon.addFile(get_image_path('right.png'), QtCore.QSize(), icon.Active, icon.Off)
        icon = self.icon_collapse_left = QtGui.QIcon()
        icon.addFile(get_image_path('left_light.png'), QtCore.QSize(), icon.Normal, icon.Off)
        icon.addFile(get_image_path('left.png'), QtCore.QSize(), icon.Active, icon.Off)
        ui.toolbutton_collapse_navigation.setIcon(icon)

        ui.treewidget_navigation.itemSelectionChanged.connect(
            self.navigation_changed)

        # Make tree fully open & non-closable
        # We expect "rootIsDecorated" has been set False in the .ui file
        tree = ui.treewidget_navigation
        tree.expandAll()
        tree.setExpandsOnDoubleClick(False)
        tree.setMaximumWidth(tree.fontMetrics().width('Boundary Conditions') + 10)
        tree.setMinimumWidth(tree.fontMetrics().width('Chemistry') + 10)

        # Make splitters non-collapsing
        for widget in widget_iter(self):
            if isinstance(widget, QtWidgets.QSplitter):
                widget.setChildrenCollapsible(False)

        # Job manager / monitor
        self.job_manager = JobManager(self)
        self.job_manager.sig_change_run_state.connect(self.slot_update_runbuttons)
        self.job_manager.sig_update_job_status.connect(self.slot_update_residuals)
        self.rundir_watcher = QFileSystemWatcher() # Move to monitor class
        self.rundir_watcher.directoryChanged.connect(self.slot_rundir_changed)

        self.monitor = Monitor(self)

        # setup Run pane #TODO move to 'setup_run'
        ui.run.checkbox_pymfix_output.stateChanged.connect(self.handle_set_pymfix_output)

        # Print welcome message.  Do this early so it appears before any
        # other messages that may occur during this __init__
        self.print_welcome()

        ## Run signals
        self.stdout_signal.connect(self.handle_stdout)
        self.stderr_signal.connect(self.handle_stderr)
        self.signal_update_runbuttons.connect(self.slot_update_runbuttons)

        # --- Register widgets ---
        self.register_keyword_widgets()

        # --- Create main menu ---
        self.init_main_menu()

        # --- vtk setup ---
        self.init_vtk_widget(loadvtk)

        # --- workflow setup ---
        self.init_workflow_widget(loadworkflow)

        # --- parameter dialog
        self.parameter_dialog = ParameterDialog(self)

        # --- default --- # do we need this?  note this gets reset in main anyhow
        self.change_mode('modeler')
        self.change_pane('model setup')

        # Update run options
        log.info('init update_runbuttons')
        self.signal_update_runbuttons.emit('')

        # Reset everything to default values
        # This is done in 'load_project'.  so why do it now?
        #self.reset() # Clear command_output too?

    def add_extra_keyword_doc(self):
        # Add a little extra knowledge ...
        # TODO , move all of this to *.f
        # These are all fractions, must be btwn 0 and 1, not documented as such

        for key in ('des_em',
                    'eps_f_min',
                    'bc_xw_g',
                    'phip0'):
            self.keyword_doc[key]['validrange'] = {'min':0.0, 'max':1.0}

        # MAX_INLET_VEL_FAC
        # DEFAULT value of 1.0
        # Error check: Value greater than or equal to 1.0
        self.keyword_doc['max_inlet_vel_fac']['validrange'] = {'min': 1.0}
        # Error check: value bounded between 0 and 1
        self.keyword_doc['ur_kth_sml']['validrange'] = {'min':0.0, 'max':1.0}

        # of particles must be non-negative
        self.keyword_doc['particles']['validrange'] = {'min':0.0}

        # Times must be nonnegative
        self.keyword_doc['res_dt']['validrange'] = {'min':0.0}

        self.keyword_doc['res_backups']['validrange'] = {'min': 0}

        self.keyword_doc['nrr']['validrange'] = {'min': 0}

        # Remove mention of 'cylindrical' since we don't support it
        self.keyword_doc['no_k']['description'] = 'Flag to disable the third dimension (i.e., 2D simulation).'

        # Remove this docstring completely - it refers to cylindrical coordinates (annluar region)
        del self.keyword_doc['xmin']['description']

        # All temperatures > 0 ?


    def set_no_project(self):
        """setup mode when no project is open"""
        self.open_succeeded = False
        self.set_solver(None)
        self.set_project_file(None)
        self.clear_unsaved_flag() # sets save button
        self.update_window_title()
        self.enable_input(False)
        # This gets set by guess_solver if we're loading a project, otherwise
        # we need to set the default.  (Do other defaults need to be set here?)
        self.status_message("No project - open existing MFIX project or create a new one")
        self.change_pane("model setup")

    def reset(self):
        """Reset all widgets to default values and set GUI to blank-slate"""

        self.change_pane("model setup") # Default pane

        # Defaults - see __init__
        self.solver_name = None
        self.fluid_solver_disabled = False  # TODO: infer at load time
        self.saved_ro_g0 = None # Hack

        self.project.reset() # Clears all keywords & collections

        self.slot_rundir_changed()

        self.reset_model_setup()
        self.reset_fluid()
        self.reset_solids()
        self.ui.regions.reset_regions()
        self.reset_ics()
        self.reset_bcs()
        self.reset_iss()
        self.reset_pss()
        self.reset_chemistry()
        self.reset_output()

        self.saved_ro_g0 = None # hack

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

        self.update_nav_tree()
        self.clear_unsaved_flag()
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
        """set_keyword(key, value, args)
             sets key(args) = value and updates GUI state.
             Always use this function,
             do not assign to keyword object directly! """
        self.set_unsaved_flag()
        self.project.submit_change(None, {key:value}, args)


    def update_keyword(self, key, value, args=None):
        """like set_keyword but no action if value already set"""

        expected_args = keyword_args.get(key)
        if expected_args is not None:
            if isinstance(args, int):
                args = [args]
            elif args is None:
                args = []
            if len(args) != len(expected_args):
                self.error("keyword %s: argument mismatch, expected %s, got %s" %
                           (key, str(expected_args), str(args)))
                return

        if value is None or value=='':
            self.unset_keyword(key, args)
            return

        #check whether it's actually changed
        v = self.project.get_value(key, args=args)

        #we might have updated 3 with 3.0 or with @(2+1), or @(1+2) with @(1+1+1)
        def typematch(v1, v2):
            if is_text_string(v1) or is_unicode(v1):
                return is_text_string(v2) or is_unicode(v2)
            return (type(v1) == type(v2))

        if typematch(v, value) and str(v)==str(value):
                return
        self.set_keyword(key, value, args=args)


    def unset_keyword(self, key, args=None):
        """Undefine keyword.  Report to user, also catch and report any errors"""
        #  Note - keyword is still registered!  This method does not deregister
        #   keywords with project manager
        if isinstance(args, int):
            args = [args]
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
                #if item.data(UserRole, 0): # Avoid toggling
                #    continue
                item.setText(0, short)


    def long_labels(self):
        tree = self.ui.treewidget_navigation
        flags =  Qt.MatchFixedString | Qt.MatchRecursive
        for (long, short) in self.nav_labels:
            items = tree.findItems(short, flags, 0)
            for item in items:
                #if item.data(UserRole,0): # Avoid toggling
                #    continue
                item.setText(0, long)


    def unimplemented(self):
        self.message(title='Unimplemented',
                     text='Feature not implemented')


    def update_nav_tree(self):
        self.ics_update_enabled()
        self.bcs_update_enabled()
        self.pss_update_enabled()
        self.iss_update_enabled()
        self.chemistry_update_enabled()


    def update_region(self, name, data):
        for update in (self.ics_update_region,
                       self.bcs_update_region,
                       self.iss_update_region,
                       self.pss_update_region,
                       self.output_update_region):
            update(name, data)


    def check_region_in_use(self, name):
        return any(check(name) for check in (self.ics_check_region_in_use,
                                             self.bcs_check_region_in_use,
                                             self.pss_check_region_in_use,
                                             self.iss_check_region_in_use,
                                             self.output_check_region_in_use))
                                             # any more places region can be used?


    def change_region_name(self, name, new_name):
        self.bcs_change_region_name(name, new_name)
        self.ics_change_region_name(name, new_name)
        self.iss_change_region_name(name, new_name)
        self.pss_change_region_name(name, new_name)
        self.output_change_region_name(name, new_name)
        # any more places region can be used?


    def toggle_nav_menu(self):
        nav_menu = self.ui.treewidget_navigation
        if self.mode != 'modeler':
            self.change_mode('modeler')
            self.change_pane('model setup')
            nav_menu.setVisible(True)
        else:
            nav_menu.setVisible(not nav_menu.isVisible())

        if nav_menu.isVisible():
            self.ui.toolbutton_collapse_navigation.setIcon(self.icon_collapse_left)
        else:
            self.ui.toolbutton_collapse_navigation.setIcon(self.icon_collapse_right)


    def status_message(self, message=''):
        self.ui.label_status.setText(message)
        if message != 'Ready': # Don't clutter the console with unimportant msgs
            self.print_internal(message, color='blue')


    def slot_rundir_changed(self):
        # Note: since log files get written to project dirs, this callback
        # is triggered frequently during a run.
        log.debug("entering slot_rundir_changed")
        runname_pid = self.get_pid_name()
        log.debug('job_manager.job: %s' % self.job_manager.job)
        if self.get_project_dir() and not self.job_manager.job:
            log.debug("slot_rundir_changed was called, pid {}".format(runname_pid))
            full_runname_pid = os.path.join(self.get_project_dir(), runname_pid)
            self.job_manager.try_to_connect(full_runname_pid)


    def set_run_button(self, text=None, enabled=None):
        if text is not None:
            self.ui.toolbutton_run_mfix.setToolTip('Resume previous MFIX run' if text=='Resume'
                                                   else text+' MFIX')
        if enabled is not None:
            b = self.ui.toolbutton_run_mfix
            b.setEnabled(enabled)


    def set_pause_button(self, text=None, enabled=None):
        b = self.ui.toolbutton_pause_mfix
        if enabled is not None:
            b.setEnabled(enabled)
        if text is not None:
            b.setToolTip(text + ' MFIX')


    def set_stop_button(self, enabled):
        b = self.ui.toolbutton_stop_mfix
        b.setEnabled(enabled)
        # tooltip?


    def set_reset_button(self, enabled):
        b =  self.ui.toolbutton_reset_mfix
        b.setEnabled(enabled)


    def enable_input(self, editable=True, partial=False):
        # Enable/disable all inputs (while job running, etc)
        # enable some inouts if partial==True (paused, resumable)
        # Stop/reset buttons are left enabled
        for pane in self.ui.panes:
            name = pane.objectName()
            pane.setEnabled(editable or partial and name in [
                'boundary_conditions', 'point_sources', 'internal_surfaces',
                'numerics', 'output', 'monitors', 'run'])


    def slot_update_residuals(self):
        """Get job status from JobManager and update residuals pane"""
        if self.job_manager and self.job_manager.job:
            log.debug('update_residuals')
            self.ui.residuals.setText(self.job_manager.job.pretty_status)

            # update plot data
            if not self.job_manager.job.is_paused():
                self.update_plots(self.job_manager.job.status)
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
            # highlight for visibility, this is an important state change
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

        log.debug("UPDATE RUN OPTIONS: pending=%s paused=%s resumable=%s",
                   pending, paused, resumable)

        self.update_window_title() # put run state in window titlebar

        self.enable_input(editable=project_open and not (pending or unpaused or paused or resumable),
                          partial=project_open and (paused or resumable))

        #handle buttons in order:  RESET RUN PAUSE STOP

        if pending:
            self.status_message("MFIX starting up, process %s" % self.job_manager.job.mfix_pid)
            # also disable spinboxes for dt, tstop unless interactive
            self.set_reset_button(enabled=False)
            self.set_run_button(enabled=False)
            self.set_pause_button(text="Pause", enabled=False)
            self.set_stop_button(enabled=True)

        elif unpaused:
            self.status_message("MFIX running, process %s" % self.job_manager.job.mfix_pid)
            # also disable spinboxes for dt, tstop unless interactive
            self.set_reset_button(enabled=False)
            self.set_run_button(enabled=False)
            self.set_pause_button(text="Pause", enabled=True)
            self.set_stop_button(enabled=True)

        elif paused:
            self.status_message("MFIX paused, process %s" % self.job_manager.job.mfix_pid)
            self.set_reset_button(enabled=False)
            self.set_run_button(text="Unpause", enabled=True)
            self.set_pause_button(text="Pause", enabled=False)
            self.set_stop_button(enabled=True)

        elif resumable:
            self.status_message("Previous MFIX run is resumable.  Reset job to edit model")
            self.set_reset_button(enabled=True)
            self.set_run_button(text='Resume', enabled=True)
            self.set_pause_button(text="Pause", enabled=False)
            self.set_stop_button(enabled=False)

        else: # Not running, ready for input
            self.status_message("Ready")
            self.set_reset_button(enabled=False)
            self.set_run_button(text="Run", enabled=project_open)
            self.set_pause_button(text="Pause", enabled=False)
            self.set_stop_button(enabled=False)

        ui.run.checkbox_pymfix_output.setEnabled(bool(paused or unpaused))


    def print_welcome(self):
        self.print_internal("Welcome to MFIX - https://mfix.netl.doe.gov",
                            color='blue')
        self.print_internal("MFIX-GUI version %s" % self.get_version(),
                            color='blue')

    def get_version(self):
        return __version_str__

    def resizeEvent(self, event):
        '''over-ride of qt resize event'''
        if self.main_menu.isVisible():
            self.main_menu.setGeometry(0, 0, self.width(), self.height())

        QtWidgets.QMainWindow.resizeEvent(self, event)

    def closeEvent(self, event):
        '''over-ride of qt close event'''
        if not self.confirm_close():
            event.ignore()
            return
        # save the state
        self.settings.setValue('geometry', self.saveGeometry())
        self.settings.setValue('splitter_left_right', self.ui.splitter_left_right.sizes())
        self.settings.setValue('splitter_graphics_cmd_output', self.ui.splitter_graphics_cmd_output.sizes())
        self.settings.setValue('mode', self.mode)

        name = self.curr_nav_label
        for (long, short) in self.nav_labels:
            if name==short:
                name=long
                break
        self.settings.setValue('navigation', name)

        # clean up
        self.main_menu.close()

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


    # move to 'scalar_handler.py'
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
            # Use "keyword_args" to delete all keys with deleted scalar index
            for IC in range(1, 1+len(self.ics)):
                self.unset_keyword('ic_scalar', args=[IC, i])
            for BC in range(1, 1+len(self.bcs)):
                self.unset_keyword('bc_scalar', args=[BC, i])
                # there are more bc_*_scalar keywords!
            i += 1

        # ICs enabled/disabled depends on nscalar
        self.update_nav_tree()



    def register_keyword_widgets(self):
        """Look for and connect keyword widgets to the project manager.
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

            if name_list[0] == 'label': # Note no 'keyword'
                if 'args' in name_list:
                    args_idx = name_list.index('args')
                    args = [try_int(name) for name in name_list[args_idx+1:]]
                    key = '_'.join(name_list[1:args_idx])
                else:
                    if name_list[-1].isdigit(): # strip suffix
                        name_list = name_list[:-1]
                    args = None
                    key = '_'.join(name_list[1:])
                if args:
                    self.add_tooltip(widget, key, value=args[0])
                else:
                    self.add_tooltip(widget, key)

            elif 'keyword' in name_list:
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

                    self.add_tooltip(widget, key)

                    # NB not all widgets get set up this way
                    if isinstance(widget, QtWidgets.QLineEdit) and widget.dtype in (int, float):
                        widget.allow_parameters = True
                    elif isinstance(widget, QtWidgets.QComboBox) and widget.count() < 1:
                            widget.addItems(list(doc['valids'].keys()))
                    elif isinstance(widget, QtWidgets.QCheckBox):
                        # TODO: infer which value is true/false?
                        values = sorted([convert_string_to_python(v) for v in doc['valids'].keys()])
                        if len(values) == 2:
                            widget.false_value = values[0]
                            widget.true_value = values[1]

                else:
                    log.error("UNKNOWN KEYWORD %s: not registering %s" % (key, widget.objectName()))
                    continue

                # register the widget with the project manager
                self.project.register_widget(widget, keys=[key], args=args)

    def init_vtk_widget(self, load_vtk=True):
        #initialize the vtk widget
        disable_vtk = False
        if not 'MFIX_NO_VTK' in os.environ and load_vtk: # Avoid importing vtkwidget if MFIX_NO_VTK set
            from mfixgui.widgets.vtkwidget import VTK_AVAILABLE
            disable_vtk = not VTK_AVAILABLE
        else: # env var set
            disable_vtk = True

        if disable_vtk:
            log.info("MFIX_NO_VTK set or vtk not importable, creating fake VTK")
            # Create a dummy object, so we don't have to test for 'if use_vtk' all over
            class FakeVtk:
                def noop(self, *args, **kwargs):
                    return None
                def __getattr__(self, key):
                    return self if key=='vtkiren' else self.noop
            self.vtkwidget = FakeVtk()
            self.ui.regions.vtkwidget = self.vtkwidget
            label = QtWidgets.QLabel('Could not import vtk, please check your installation.')
            label.setAlignment(QtCore.Qt.AlignCenter);
            self.ui.horizontalLayoutModelGraphics.addWidget(label)
            return

        from mfixgui.widgets.vtkwidget import VtkWidget
        self.vtkwidget = VtkWidget(self.project, parent=self)
        self.ui.horizontalLayoutModelGraphics.addWidget(self.vtkwidget)

        # register with project manager
        self.project.register_widget(
            self.vtkwidget, ['xmin', 'xlength', 'ymin', 'ylength', 'zmin',
                             'zlength', 'imax', 'jmax', 'kmax', 'no_k',
                             ])

        # add reference to other widgets
        self.ui.regions.vtkwidget = self.vtkwidget

    def init_workflow_widget(self, load_workflow=True):
        # initialize the workflow widgets if pyqtnode is available
        if load_workflow:
            from mfixgui.widgets.workflow import WorkflowWidget, PYQTNODE_AVAILABLE
            if PYQTNODE_AVAILABLE:
                self.ui.workflow_widget = WorkflowWidget(self.project, self)
                self.ui.workflow_widget.PYQTNODE_AVAILABLE = True
                self.ui.verticallayout_workflow_mode.addWidget(
                    self.ui.workflow_widget)
            else:
                self.make_fake_workflow()
        else:
            self.make_fake_workflow()

    def make_fake_workflow(self):
        class FakeWorkflow:
            PYQTNODE_AVAILABLE = False
            def noop(self, *args, **kwargs):
                return None
            def __getattr__(self, key):
                return self if key=='nodeChart' else self.noop
        self.ui.workflow_widget = FakeWorkflow()
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

    def change_mode(self, mode):
        """change the Modeler, Workflow, Developer tab"""
        self.mode = mode
        to_index = None
        if mode == 'interpreter':
            self.capture_output(True)
            self.setup_interpreter()
        else:
            self.capture_output(False)
        for i in range(self.ui.stackedwidget_mode.count()):
            widget = self.ui.stackedwidget_mode.widget(i)
            if mode == str(widget.objectName()):
                to_index = i
                break
        if to_index is None:
            self.error("Invalid mode %s" % mode)
            return

        for key, btn in self.modebuttondict.items():
            btn.setChecked(mode == key)

        workflow = (mode=='workflow')
        if self.ui.workflow_widget.PYQTNODE_AVAILABLE:
            nc = self.ui.workflow_widget.nodeChart
            for btn in [nc.runToolButton, nc.autorunToolButton, nc.stopToolButton, nc.stepToolButton]:
                btn.setVisible(workflow)

        ui = self.ui
        for btn in (ui.toolbutton_reset_mfix, ui.toolbutton_run_mfix,
                    ui.toolbutton_pause_mfix, ui.toolbutton_stop_mfix):
            btn.setVisible(not workflow)

        self.animate_stacked_widget(self.ui.stackedwidget_mode,
                                    self.ui.stackedwidget_mode.currentIndex(),
                                    to_index,
                                    'horizontal')

        if mode == 'modeler': # open navigation menu whenever we go back to Modeler mode
            ui.treewidget_navigation.setVisible(True)
            self.ui.toolbutton_collapse_navigation.setIcon(self.icon_collapse_left)


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
        #self.navigation_changed() # not needed, since setCurrentItem triggers callback

    def navigation_changed(self):
        """an item in the tree was selected, change panes"""
        current_selection = self.ui.treewidget_navigation.selectedItems()
        if not current_selection:
            return None
        self.curr_nav_label = name = str(current_selection[0].text(0))
        # Translate from short to long name
        for (long, short) in self.nav_labels:
            if name==short:
                name=long
                break
        name = '_'.join(name.lower().split(' '))

        matches = [i
                   for i in range(self.ui.stackedWidgetTaskPane.count())
                   if self.ui.stackedWidgetTaskPane.widget(i).objectName() == name]

        assert len(matches) == 1

        to_index = matches[0]

        self.animate_stacked_widget(
            self.ui.stackedWidgetTaskPane,
            self.ui.stackedWidgetTaskPane.currentIndex(),
            to_index)

        self.setup_current_pane()


    def setup_current_pane(self):
        # Force any open popup to close
        # (if dialog is modal we don't need this)
        self.species_popup.done(0)
        self.regions_popup.done(0)
        current_selection = self.ui.treewidget_navigation.selectedItems()
        if not current_selection:
            return
        text = str(current_selection[-1].text(0))
        text = '_'.join(text.lower().split(' '))
        # Make sure panes are properly initialized as we change to them.  This way
        # we do not have to worry so much about inter-pane state dependencies
        if text == 'model_setup':
            self.setup_model_setup()
        if text == 'fluid' : #
            self.setup_fluid()
        elif text == 'solids':
            self.setup_solids()
        elif text in ('initial_conditions', 'ics'):
            self.setup_ics()
        elif text in ('boundary_conditions', 'bcs'):
            self.setup_bcs()
        elif text in ('point_sources', 'pss'):
            self.setup_pss()
        elif text in ('internal_surfaces', 'iss'):
            self.setup_iss()
        elif text == 'chemistry':
            self.setup_chemistry()
        elif text == 'numerics':
            self.setup_numerics()
        elif text == 'output':
            self.setup_output()

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
            return #?
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

        # animation group
        # call back for when the animation is finsished/canceled.
        self.stack_animation.stateChanged.connect(lambda: self.animate_stacked_widget_finished(
            stackedwidget, from_, to, btn_layout, to_btn, line, line_to))

        self.animating = True
        self.stack_animation.start()

    def animation_setup(self, target, x_start, y_start, x_end, y_end):
        """setup an animation widget"""
        animation = QtCore.QPropertyAnimation(target, "pos".encode('utf-8'))
        animation.setDuration(ANIMATION_SPEED)
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
                text='',
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
            self.unset_keyword(k, a) #TODO FIXME does not weed out malformed keywords

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
            logmsg = logmsg[5:].strip()
            line = line[5:].strip() # Supress 'Info:' in console window (?)
            log.info(logmsg)
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


    def remove_output_files(self, output_files=None, message_text=None, force_remove=False):
        """ remove MFIX output files from current project directory

        output_files: List of patterns to be matched for file removal
        return: True for success, False for user cancel"""

        if not output_files:
            output_files = self.monitor.get_outputs()
        if not message_text:
            message_text = 'Deleting output files:\n %s' % '\n'.join(output_files)

        ok_to_delete = self.output_selection_popup.exec_(output_files, force_remove)

        if not ok_to_delete:
            return False

        for path in self.output_selection_popup.get_output_files():
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
        try:
            # no job, create a new one
            if not self.job_manager.job:
                # open the run dialog for job options
                self.open_run_dialog()
                return
            # paused
            else:
                if self.unsaved_flag and self.job_manager.job.is_paused():
                    self.handle_reinit()
                self.job_manager.job.unpause()
        except Exception as e:
            self.error('handle_run: %s' % str(e))
        self.signal_update_runbuttons.emit('')

    def handle_set_pymfix_output(self):
        try:
            self.job_manager.job.set_pymfix_output(
              self.ui.run.checkbox_pymfix_output.isChecked())
        except Exception as e:
            self.error('handle_set_pymfix_output: %s' % e)


    def handle_pause(self):
        try:
            self.job_manager.job.pause()
        except Exception as e:
            self.error('handle_pause: %s' % e)

    def handle_reinit(self):
        """Save current project with copy <run_name>.reinit.NN, then
        reinitialize through job manager"""
        if ( self.unsaved_flag
             and self.job_manager.job
             and self.job_manager.job.is_paused()
             and not self.job_manager.is_job_pending()):

            project_dir = self.get_project_dir()
            project_file = self.get_project_file()
            run_name = self.get_runname()
            # save reinit.X version
            try:
                rfiles = glob.glob(os.path.join(
                           project_dir, run_name+'.reinit.[0-9]*'))
                rfiles.sort()
                filenum = rfiles[-1].split('.')
                filenum = int(filenum[-1]) + 1
            except IndexError:
                filenum = 1
            newname = "%s.reinit.%s" % (run_name, filenum)
            reinit_file = os.path.join(project_dir, newname)
            self.project.writeDatFile(reinit_file)
            self.print_internal('reinitialize file saved %s' % reinit_file)

            try:
                self.job_manager.job.reinit(''.join(self.project.to_string()).encode('utf-8'))
                log.debug('gui returned from job reinit call')
                #TODO: save project file upon reinit success
            except Exception as e:
                self.error('handle_reinit: %s' % e)
        else:
            log.debug('reinitialize called in invalid state')
        log.debug('gui leaving handle_reinit')

    def handle_stop(self):
        try:
            self.job_manager.stop_mfix()
        except Exception as e:
            self.error('handle_stop: %s' % e)

    def check_save(self, text="Save current project?"):
        if self.unsaved_flag:
            response = self.message(title="Save?",
                                    icon="question",
                                    text=text,
                                    buttons=['yes', 'no'])
            if response != 'yes':
                return False
        self.save_project()
        return True

    def open_run_dialog(self):
        """Open run popup dialog"""
        if not self.check_save():
            return
        self.run_dialog = RunPopup(self.commandline_option_exe, self)
        self.run_dialog.set_run_mfix_exe.connect(self.handle_exe_changed)
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
        self.check_save(text='Save current project before exporting?')
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

        # new project name
        new_project_file = os.path.splitext(os.path.basename(export_file))[0]

        # copy the project file
        export_prj = copy.deepcopy(self.project)
        export_prj.updateKeyword('run_name', new_project_file)
        self.print_internal("Info: Exporting %s" % new_project_file)
        export_prj.writeDatFile(os.path.join(export_dir, new_project_file+'.mfx'))

        # collect files to copy
        files_to_copy = []
        # output files *.RES, *.SP?, *.vtu, etc...
        output_files = self.monitor.get_outputs()
        ok_to_copy = self.output_selection_popup.exec_(output_files, heading='Copy the following files?', title='Copy files?')
        if not ok_to_copy: # user bailed out
            return
        files_to_copy.extend(self.output_selection_popup.get_output_files())

        # stl files
        files_to_copy.extend(self.monitor.get_outputs(["*.[Ss][Tt][Ll]"]))

        #copy project files into new_project_directory
        for f in files_to_copy:
            try:
                shutil.copyfile(f, os.path.join(export_dir, os.path.basename(f)))
            except Exception as e:
                self.message(title='File error',
                             text="Error copying file:\n%s" % e,
                             buttons=['ok'])
        self.handle_main_menu_hide()

    def save_project(self, filename=None):
        """save project, optionally as a new project.

        project_file: string, full filename of project (including path)
        return: None"""

        if filename:
            project_dir = os.path.dirname(filename)
            project_file = filename
        else:
            project_dir = self.get_project_dir()
            project_file = self.get_project_file()

        if project_dir is None or project_file is None:
            return

        # save version
        v = self.project.mfix_gui_comments.get('project_version', 0)
        self.project.mfix_gui_comments['project_version'] = str(int(v) + 1)
        self.project.mfix_gui_comments['gui_version'] = self.get_version()

        self.project.mfix_gui_comments['project_notes'] = json.dumps(self.ui.main_menu_project_notes.toPlainText())

        # save users
        users = self.project.mfix_gui_comments.get('modified_by', [])
        if not isinstance(users, list):
            users = users.split('|')
        users.append(get_username())
        users = set(users)
        self.project.mfix_gui_comments['modified_by'] = '|'.join(users)

        # save geometry
        self.vtkwidget.export_stl(os.path.join(project_dir, 'geometry.stl'))
        self.project.mfix_gui_comments['geometry'] = self.vtkwidget.geometry_to_str()
        self.project.mfix_gui_comments['visual_props'] = self.vtkwidget.visual_props_to_str()

        # save regions
        self.project.mfix_gui_comments['regions_dict'] = self.ui.regions.regions_to_str()

        for (data, key) in ((self.bcs_to_str(), 'bc_regions'),
                            (self.ics_to_str(), 'ic_regions'),
                            (self.iss_to_str(), 'is_regions'),
                            (self.pss_to_str(), 'ps_regions'),
                            (self.vtks_to_str(), 'vtk_regions'),
                            (self.chemistry_to_str(), 'chemistry')):
            if data:
                self.project.mfix_gui_comments[key] = data
            else:
                self.project.mfix_gui_comments.pop(key, None)


        project_base = os.path.basename(project_file)
        run_name = os.path.splitext(project_base)[0]
        self.update_keyword('run_name', run_name)
        self.print_internal("Info: Saving %s\n" % project_file)
        self.project.writeDatFile(project_file)

        # save workflow
        if self.ui.workflow_widget.PYQTNODE_AVAILABLE:
            self.ui.workflow_widget.save(
                os.path.abspath(os.path.join(project_dir, 'workflow.nc')))

        self.clear_unsaved_flag()
        self.update_source_view()
        self.save_recent_projects()

    def save_as(self):
        """Prompt user for new filename, save project to that file and make
        it the active project"""
        new_file = self.get_save_filename()
        if not new_file:
            return
        new_dir = os.path.dirname(new_file)
        if not self.check_writable(new_dir):
            return

        old_dir = self.get_project_dir()

        # Force run name to file name.  Is this a good idea?
        basename = os.path.basename(new_file)
        run_name = os.path.splitext(basename)[0]
        self.set_project_file(new_file)
        self.update_keyword('run_name', run_name)
        self.save_project()

        # change file watcher
        self.rundir_watcher.removePath(old_dir)
        self.rundir_watcher.addPath(new_dir)
        self.slot_rundir_changed()
        self.signal_update_runbuttons.emit('')
        self.handle_main_menu_hide()

    def get_save_filename(self, message=None):
        """wrapper for call to getSaveFileName, override in unit tests"""
        if message is None:
            message = 'Save Project As'
        default = os.path.join(self.get_project_dir(), self.get_runname()+".mfx",)
        filename = QtWidgets.QFileDialog.getSaveFileName(self, message, default, "*.mfx")
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
        self.handle_main_menu_hide()

    def handle_export(self):
        self.export_project()

    def handle_save_as(self):
        self.save_as()

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
            # add entire path to title, abbreviate user dir
            title += " - " + project_file.replace(os.path.expanduser('~'), '~')
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
        self.setWindowTitle(title)

    def set_save_button(self, enabled):
        self.ui.toolbutton_save.setEnabled(enabled)
        # set the one in the main menu too!
        self.disable_main_menu_items([] if enabled else ['save'])

    def set_unsaved_flag(self):
        if not self.unsaved_flag:
            # For debugging problems where flag gets set during load
            #import traceback
            #traceback.print_stack()
            log.info("Project is unsaved")
        self.unsaved_flag = True
        self.update_window_title()
        self.set_save_button(enabled=True)

    def clear_unsaved_flag(self):
        if self.unsaved_flag:
            log.info("Project is saved")
        self.unsaved_flag = False
        self.set_save_button(enabled=False)
        # reinit support
        self.signal_update_runbuttons.emit('')

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

    def open_new_from_template(self, mfx_file):

        """ This function takes the path mfx_file, saves it to a new
        location specified by the user, and opens the new location"""

        project_file = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Create Project File',
            "", "MFX files (*.mfx)")

        if PYQT5:
            project_file = project_file[0]
        if not project_file:
            return

        if not project_file.endswith(".mfx"):
            project_file = project_file + ".mfx"

        if os.path.exists(project_file):
            if not self.confirm_clobber(project_file):
                return

        basename = os.path.basename(project_file)
        project_dir = os.path.dirname(project_file)
        run_name = os.path.splitext(basename)[0]
        if not self.check_writable(project_dir):
            self.message(text='Unable to write to %s' % project_dir,
                         buttons=['ok'],
                         default='ok')
            return
        # Start from template
        template = mfx_file
        creator = get_username()
        creation_time = time.strftime('%Y-%m-%d %H:%M')
        try:
            with open(template) as infile:
                with open(project_file,'w') as outfile:
                    for line in infile:
                        if line.startswith('# Generic'):
                            # Put something better in the comment field, otherwise
                            # "generic template" lingers there forever
                            line = '# Created by %s on %s\n' % (creator, creation_time)
                        elif '${run_name}' in line:
                            line = line.replace('${run_name}', run_name)
                        outfile.write(line)
        except Exception as e:
            self.message(text="Error %s creating new project" % e,
                         buttons=['ok'],
                         default=['ok'])
            self.set_no_project()
            return

        self.open_project(project_file)

        # add some info
        self.project.mfix_gui_comments['author'] = creator
        self.project.mfix_gui_comments['created_date'] = creation_time

        self.save_project()


    def get_open_filename(self):
        """wrapper for call to getOpenFileName, override in for unit tests"""
        project_dir = self.get_project_dir()
        project_path = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open Project Directory', project_dir,
            'MFIX Project (*.mfx *.dat);; All Files (*)')
        if PYQT5:
            # qt4/qt5 compat hack
            project_path = project_path[0]
        return project_path

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

    def get_runname(self, default='new_file'):
        name = str(self.project.get_value('run_name', default=default))
        for char in ('.', '"', "'", '/', '\\', ':'):
            name = name.replace(char, '_')
        return name

    def get_pid_name(self, check_exists=False):
        '''return the pid name if it exists, else an empty string if
        check_exists==True'''
        # look for all caps!
        run_name = self.get_runname().upper() + '.pid'
        if not os.path.exists(run_name) and check_exists:
            run_name = ''
        return run_name

    def save_recent_projects(self):
        '''Save recent projects'''

        rec_prjs = self.settings.value('recent_projects')
        if rec_prjs is None:
            rec_prjs = []
        else:
            rec_prjs = rec_prjs.split('|')

        new_rec_prjs = list(set([self.get_project_file()] + rec_prjs))[:MAX_RECENT_PROJECTS]
        self.settings.setValue('recent_projects', '|'.join(new_rec_prjs))


    def open_project(self, project_path, interactive=True):
        """Open MFiX Project"""
        # Too much going on in this method, maybe split some of this out

        if self.main_menu.isVisible():
            self.handle_main_menu_hide()

        #self.change_pane('model setup')

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
        # (may not be needed on initial load)

        basename, pathname = os.path.split(project_file)

        self.print_internal("Loading %s from %s" % (pathname, basename), color='blue')
        try:
            self.project.load_project_file(project_file)
        except Exception as e:
            msg = 'Failed to load %s: %s: %s' % (os.path.basename(project_file),
                                                 e.__class__.__name__, e)
            self.error(msg, popup=interactive) # don't popup in -q mode

            traceback.print_exception(*sys.exc_info())
            # Should we stick this in the output window?  no, for now.

            self.set_no_project()
            return

        default_runname = os.path.splitext(os.path.basename(project_file))[0]
        runname = self.get_runname(default=default_runname)
        runname_mfx = runname + '.mfx'
        runname_pid = self.get_pid_name(True)

        if runname_pid:
            # previously started job may be running, try to reconnect
            log.debug('attempting to connect to running job %s' % runname_pid)
            self.job_manager.try_to_connect(runname_pid)

        self.setup_current_pane() # update vals in any open tabs
        self.update_source_view()

        # set up rundir watcher
        self.rundir_watcher.addPath(project_dir)
        self.slot_rundir_changed()

        ### Geometry
        # Look for geometry.stl and load automatically
        geometry_file = os.path.abspath(os.path.join(project_dir, 'geometry.stl'))
        if os.path.exists(geometry_file) and 'geometry' not in self.project.mfix_gui_comments:
            self.vtkwidget.add_stl(None, filename=geometry_file)
        else:
            # order needs to be visual_props -> geometry -> regions (loaded below)
            # load props first
            props = self.project.mfix_gui_comments.get('visual_props')
            if props:
                self.vtkwidget.visual_props_from_str(props)
            geo = self.project.mfix_gui_comments.get('geometry')
            if geo:
                self.vtkwidget.geometry_from_str(geo)

        #  Non-keyword params stored as !#MFIX-GUI comments
        solids_phase_names = {}
        try:
            for (key, val) in self.project.mfix_gui_comments.items():
                if key == 'fluid_phase_name':
                    self.set_fluid_phase_name(val)
                elif key.startswith('solids_phase_name('):
                    n = int(key.split('(')[1][:-1])
                    solids_phase_names[n] = val
                elif key == 'regions_dict':
                    self.ui.regions.regions_from_str(val)
                elif key == 'ic_regions':
                    self.ics_regions_from_str(val)
                elif key == 'bc_regions':
                    self.bcs_regions_from_str(val)
                elif key == 'is_regions':
                    self.iss_regions_from_str(val)
                elif key == 'ps_regions':
                    self.pss_regions_from_str(val)
                elif key == 'vtk_regions':
                    self.vtk_regions_from_str(val)
                elif key == 'chemistry':
                    self.chemistry_from_str(val)
                elif key == 'geometry':
                    pass # handled in 'geometry' section above
                elif key == 'visual_props':
                    pass # handled in 'geometry' section above
                elif key == 'parameters':
                    pass # handled in project
                # Add more here
                #else:  # Too many warnings!
                #    self.warn("Unknown mfix-gui setting '%s'" % key)

        except Exception as e:
            self.error("%s: %s" % (key, e))

        # Copy ordered dict to modify keys w/o losing order
        if solids_phase_names:
            s = OrderedDict()
            for (i, (k, v)) in enumerate(self.solids.items(), 1):
                s[solids_phase_names.get(i, k)] = v
            self.solids = s

        #### Fluid phase
        # fluid species table # move this stuff to setup_fluid!
        self.update_fluid_species_table()

        # fluid momentum and species eq. handled by _keyword_ widget
        # fluid scalar eq
        nscalar = self.project.get_value('nscalar', default=0)
        self.fluid_nscalar_eq = sum(1 for i in range(1, nscalar+1)
                                    if self.project.get_value('phase4scalar', args=i) == 0)
        self.solids_nscalar_eq = sum(1 for i in range(1, nscalar+1)
                                    if self.project.get_value('phase4scalar', args=i) != 0)


        # solid scalar eq checkbox will be handled in update_solids_detail_pane
        # need to initialize per-solid 'nscalar_eq' so
        # that update_scalar_equations will do the right thing
        for (i, (k,v)) in enumerate(self.solids.items(), 1):
            v['nscalar_eq'] = sum(1 for j in range(1, nscalar+1)
                                  if self.project.get_value('phase4scalar', args=j) == i)


        # setValue will trigger update_scalar_equations (ok)
        self.ui.fluid.spinbox_nscalar_eq.setValue(self.fluid_nscalar_eq)
        self.enable_fluid_scalar_eq(self.fluid_nscalar_eq > 0)

        #Move to fluid_handler!
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
                # this is getting printed after error count ... should be included in # of errs

            # XXX FIXME conflicts with default fluid models (?)
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
        # Needed?  will this get done when we switch to solids tab?
        self.update_solids_table()
        self.solids_update_tabs()
        self.update_solids_detail_pane()

        ### Regions
        # Look for regions in IC, BC, PS, etc.
        self.ui.regions.extract_regions(self.project, project_dir)
        # Take care of updates we deferred during extract_region
        # FIXME Do this when switching to regions pane
        self.ui.regions.tablewidget_regions.fit_to_contents()

        # background mesh
        self.init_background_mesh()

        # "Extract" pulls out info from non-GUI project files, which don't have MFIX_GUI section
        # Initial conditions
        self.ics_extract_regions()

        # Boundary conditions
        self.bcs_extract_regions()

        # Point sources
        self.pss_extract_regions()

        # Internal surfaces
        self.iss_extract_regions()

        # VTK output regions
        self.vtk_extract_regions()

        # Chemistry
        self.chemistry_extract_info()

        ### Workflow
        if self.ui.workflow_widget.PYQTNODE_AVAILABLE:
            workflow_file = os.path.abspath(os.path.join(project_dir, 'workflow.nc'))
            if os.path.exists(workflow_file):
                self.ui.workflow_widget.clear()
                self.ui.workflow_widget.load(workflow_file)
            self.ui.workflow_widget.look_for_projects(basename)


        if interactive and not project_path.endswith(runname_mfx):
            ok_to_write =  self.confirm_rename(project_file, runname_mfx)
            if ok_to_write:
                renamed_project_file = os.path.join(project_dir, runname_mfx)
                if os.path.exists(renamed_project_file):
                    ok_to_write = self.confirm_clobber(renamed_project_file)
            if not ok_to_write:
                self.print_internal("Rename canceled at user request")
                return

            project_file = renamed_project_file
            try:
                self.print_internal("Info: Saving %s" % project_file)
                self.project.writeDatFile(project_file) #XX
                #self.print_internal(save_msg, color='blue')
                self.clear_unsaved_flag()
            except Exception as e:
                msg = 'Failed to save %s: %s: %s' % (os.path.basename(project_file),
                                                     e.__class__.__name__, e)
                self.print_internal("Error: %s" % msg, color='red')
                self.message(title='Error',
                             icon='error',
                             text=msg,
                             buttons=['ok'],
                             default='ok')
                traceback.print_exception(*sys.exc_info())
                return

        self.set_project_file(project_file)

        self.vtkwidget.reset_view()
        self.vtkwidget.render(defer_render=False)
        self.open_succeeded = True
        self.signal_update_runbuttons.emit('')
        self.update_nav_tree()

        #if self.unsaved_flag: #
        # Settings changed during loading
        #    self.save_project()- let the user do this

        self.save_recent_projects()


    def add_tooltip(self, widget, key, description=None, value=None):
        if description is None:
            doc = self.keyword_doc.get(key)
            if not doc:
                return
            description = doc.get('description')
            if description is None:
                return
            if value is not None and 'valids' in doc:
                vkey = '.FALSE.' if value is False else '.TRUE.' if value is True else str(value)
                if vkey in doc['valids']:
                    description = doc['valids'][vkey].get('note', str(value))

        #if 'cylindrical' in description.lower():
        #    print('CYLINDRICAL found! key=%s tooltip=%s' % (key, description))
        #    return

        # Clean it up a little
        description = description.strip()

        # '>' and '<' will mess up HTML
        description = description.replace('<', '&lt;')
        description = description.replace('>', '&gt;')

        # non-breaking hyphen
        description = description.replace('-', '&#8209;')

        ### epsilon
        description = description.replace('epsilon', 'ε')

        # Don't split diff. eq's over multiple lines
        pat = re.compile('boundary condition: *([^,]*, )')
        match = pat.search(description)
        if match:
            text = match.group(1)
            description = description.replace(text, '<br/>%s<br/>'%text[:-2].replace(' ', '&nbsp;'))

        # Default
        #pat = re.compile(r'\[[^]]+\]')
        #while True:
        #    match = pat.search(description)
        #    if not match:
        #        break
        #    text = match.group(0)
        #    description = description.replace(text, '<i>Default: %s</i>'%text[1:-1])

        # Bullets
        description = description.replace(' o ', '<br/>&bull;')

        args = widget.args if hasattr(widget, 'args') else None
        if args is None:
            nargs = 0
        elif isinstance(args, int):
            nargs = 1
        else:
            nargs = len(args)

        # This is a weird place to be doing this check
        if getattr(widget, 'key', None) == key and len(keyword_args.get(key, [])) != nargs:
            if not key.startswith('des_'):
                self.error("keyword args mismatch: key=%s: expected %s, got %s" %
                           (key, keyword_args.get(key), args))

        if isinstance(args, int):
            key = '%s(%s)' % (key, args)
        elif args:
            key = '%s(%s)' % (key, ','.join(map(
                lambda x: str(x[0] if isinstance(x, (tuple, list)) else str(x)), args)))
        if value is not None:
            key = '%s=%s' % (key, value)
        if key is None:
            msg = '<b></b>%s</br>' % description
        else:
            msg = '<b>%s</b>: %s</br>' % (key, description)
        widget.setToolTip(msg)
        widget.help_text = msg # Can we get more info here, so help_text is not just repeating tooltip?


    # Following functions are overrideable for test runner
    def confirm_rename(self, project_file, runname_mfx):
        rename_msg = 'Renaming %s to %s based on run name' % (project_file, runname_mfx)
        response = self.message(title='Info',
                                icon='question',
                                text=rename_msg,
                                buttons=['ok', 'cancel'],
                                default='ok')
        return response == 'ok'

    def confirm_clobber(self, renamed_project_file):
        clobber_msg = '%s exists, replace?' % renamed_project_file
        response = self.message(title='Warning',
                        icon='question',
                        text=clobber_msg,
                        buttons=['yes', 'no'],
                        default='no')
        return response == 'yes'

    def confirm_delete_files(self, message_text):
        response = self.message(title="Info",
                               icon="info",
                               text=message_text,
                               buttons=['ok','cancel'],
                               default='cancel')
        return response == 'ok'


def main():
    global gui

    # build the arg parser
    av_styles = [s.lower() for s in QtWidgets.QStyleFactory.keys()]
    parser = argparse.ArgumentParser(description='MFIX GUI')
    parser.add_argument('project', action='store', nargs='?', default=None,
                        help='open mfix.dat or <RUN_NAME>.mfx project file or search a specified directory for project files')
    parser.add_argument('-e', '--exe',  metavar='EXE', action='store', default=None,
                        help='specify MFIX executable (full path)')
    parser.add_argument('-l', '--log', metavar='LOG', action='store', default='WARN',
                        choices=['error', 'warning', 'info', 'debug'],
                        help='set logging level (error, warning, info, debug)')
    parser.add_argument('-s', '--style', metavar='STYLE', action='store', default=None,
                        choices=av_styles,
                        help='specify app style (windowsvista, fusion, cleanlooks,...)')
    parser.add_argument('-n', '--noload', action='store_true',
                        help='do not autoload previous project')
    parser.add_argument('-w', '--noworkflow', action='store_false',
                        help='do not load the workflow environment')
    parser.add_argument('-k', '--novtk', action='store_false',
                         help='do not load vtk')
    parser.add_argument('-q', '--quit', action='store_true',
                        help='quit after opening file (for testing)')
    parser.add_argument('-d', '--default_geo', action='store_true',
                        help="Use default geometry, don't restore previous state.")
    parser.add_argument('-v', '--version', action='version', version=__version_str__)

    # parse the args
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                        filemode='w', level=getattr(logging, args.log.upper()),
                        format='%(name)s - %(levelname)s - %(message)s')

    project_file = args.project
    if project_file and os.path.isdir(project_file):
        mfx_files = glob.glob(os.path.join(project_file, '*.mfx'))
        if mfx_files:
            project_file = mfx_files[0]
        else:
            dat_files = glob.glob(os.path.join(project_file, 'mfix.dat'))
            if dat_files:
                project_file = dat_files[0]
            else:
                print("Can't find *.mfx or mfix.dat in directory: %s" % project_file)
                parser.print_help()
                sys.exit()

    elif project_file and not os.path.isfile(project_file):
        print("%s: is not a file " % project_file)
        parser.print_help()
        sys.exit()

    # create the QApplication
    qapp = QtWidgets.QApplication([]) # TODO pass args to qt
    # set style
    app_style_arg = args.style
    app_style_settings = SETTINGS.value('app_style')
    if app_style_arg is not None:
        qapp.setStyle(app_style_arg.lower())
        SETTINGS.setValue('app_style', app_style_arg)
    elif app_style_settings is not None:
        qapp.setStyle(app_style_settings.lower())

    # create the gui
    gui = MfixGui(qapp, project_file=project_file, loadworkflow=args.noworkflow,
                  loadvtk=args.novtk)

    # set previous geometry
    geo = SETTINGS.value('geometry')
    if geo is not None and not args.default_geo:
        gui.restoreGeometry(geo)
        left_right = SETTINGS.value('splitter_left_right')
        if left_right is not None:
            gui.ui.splitter_left_right.setSizes([int(num) for num in left_right])
        cmd_output = SETTINGS.value('splitter_graphics_cmd_output')
        if cmd_output is not None:
            gui.ui.splitter_graphics_cmd_output.setSizes([int(num) for num in cmd_output])
    gui.show()

    if args.exe:
        #print('exe option passed: %s' % mfix_exe_option)
        gui.commandline_option_exe = args.exe

    last_prj = None
    if project_file is None and not args.noload:
        # autoload last project
        last_prj = project_file = gui.get_project_file()

    if project_file and not args.noload:
        gui.open_project(project_file, interactive=(not args.quit))
        # change mode and navigation if loaded last project
        if last_prj == project_file:
            m = SETTINGS.value('mode')
            if m is not None:
                gui.change_mode(m)
            n = SETTINGS.value('navigation')
            if n is not None:
                gui.change_pane(n)
    else:
        gui.set_no_project()
        gui.handle_main_menu()

    # have to initialize vtk after the widget is visible!
    gui.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    # This makes it too easy to skip the exit confirmation dialog.  Should we allow it?
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    if not args.quit:
        qapp.exec_()

    qapp.deleteLater()
    sys.exit()

if __name__  == '__main__':
    main()
