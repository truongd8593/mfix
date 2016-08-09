# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
This module contains the work flow widget.
'''
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import pyqtSignal
from collections import OrderedDict
import copy
import os
import shutil
import glob

try:
    from pyqtnode import NodeWidget, Node, tools
    PYQTNODE_AVAILABLE = True
except ImportError:
    NodeWidget = None
    PYQTNODE_AVAILABLE = False
    Node = object

# local imports
from widgets.base import Table
from widgets.run_popup import RunPopup
from tools.general import get_icon
from constants import PARAMETER_DICT
from job import JobManager
from project import Project


# --- Custom MFIX GUI Nodes ---
class TestNode(Node):
    name = 'test'

    def __init__(self):
        self.terminalOpts = OrderedDict([
            ('used parameters', {'widget': 'pushbutton',
                      'in': False,
                      'out': False,
                      'showlabel': False,
                      'dtype': bool,
                      }),
            ('export', {'widget': 'pushbutton',
                        'in': False,
                        'out': False,
                        'showlabel': False,
                        'dtype': bool,
                        }),
            ('run', {'widget': 'pushbutton',
                        'in': False,
                        'out': False,
                        'showlabel': False,
                        'dtype': bool,
                        }),
             ])

        Node.__init__(self)

        self.terminals['used parameters'].valueChanged.connect(self.used_parameters)
        self.terminals['export'].valueChanged.connect(self.export)
        self.terminals['run'].valueChanged.connect(self.run_project)

        self.terminals['run'].widget.setEnabled(False)
        self.test_num = 1

    def used_parameters(self):
        print(self.parent.workflow_widget.used_parameters)

    def export(self):

        test_name = 'test'+str(self.test_num)
        self.test_num+=1
        curr_proj_dir = self.parent.mfixgui.get_project_dir()
        self.exp_path = os.path.join(curr_proj_dir, test_name)
        if not os.path.exists(self.exp_path):
            os.mkdir(self.exp_path)

        self.proj_file = self.parent.workflow_widget.export_project(
            self.exp_path, # export path
            {'x': 1, 'y': 2, 'z': 1}, # parameters
            {'drag_type': 'WEN_YU', 'BC_V_g,1': 5.0}, # keywords
            )

        self.terminals['run'].widget.setEnabled(True)

    def run_project(self):
        self.parent.workflow_widget.run_project(self.proj_file)


# --- Mock Parent for job submission ---
class Mock(object):
    def noop(self, *args, **kwargs):
        return False
    def __getattr__(self, name):
        return self.noop


class MockParent(QtWidgets.QWidget):
    stderr_signal = pyqtSignal(object)
    stdout_signal = pyqtSignal(object)
    signal_update_runbuttons = pyqtSignal(object)
    project = None
    project_dir = None
    project_file = None
    settings = None
    monitor = Mock()
    ui = Mock()
    ui.tabWidgetGraphics = Mock()

    def __init__(self, parent=None, mfix_gui = None):
        QtWidgets.QWidget.__init__(self, parent)
        self.mfix_gui = mfix_gui

    def noop(self, *args, **kwargs):
        return None

    def __getattr__(self, name):
        return self.noop

    def get_project_dir(self):
        return self.project_dir

    def get_project_file(self):
        return self.project_file

    def print_internal(self, text, color=None):
        if self.mfix_gui is not None:
            self.mfix_gui.print_internal(text, color)

    def message(self, title='Warning', icon='warning',
                text='This is a warning.', buttons=['ok'], default='ok',
                infoText=None, detailedText=None, ):
        if self.mfix_gui is not None:
            self.mfix_gui.message(title, icon, text, buttons, default,
                                  infoText, detailedText)

# --- custom run_pop ---
class WorkflowRunPopup(RunPopup):
    def __init__(self, mfix_exe, parent):
        RunPopup.__init__(self, mfix_exe, parent)

    # over-rides
    def handle_run(self):
        self.finish_with_dialog()

    def handle_submit(self):
        self.finish_with_dialog()

# --- Workflow Widget ---
class WorkflowWidget(QtWidgets.QWidget):
    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.project = project
        self.job_dict = {}
        self.update_timer = QtCore.QTimer()
        self.update_timer.timeout.connect(self.update_job_status)
        self.run_cmd = None
        self.file_watcher = QtCore.QFileSystemWatcher()
        self.file_watcher.directoryChanged.connect(self.create_job_manager)
        self.mock_parents = []

        # --- initalize the node widget ---
        self.nodeChart = NodeWidget(showtoolbar=True)
        if hasattr(self.nodeChart, 'needsSavedEvent'):
            self.nodeChart.needsSavedEvent.connect(self.set_save_btn)
        self.nodeChart.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Preferred)
        self.nodeChart.setGeometry(0, 0, 100, 1000)

        # modify toolbar
        i = 2
        for tool, icon, callback in [('import', 'import.png', self.handle_import),
                                     ('export', 'open_in_new.png', self.handle_import)]:
            btn = QtWidgets.QToolButton()
            btn.setIcon(get_icon(icon))
            btn.pressed.connect(callback)
            btn.setAutoRaise(True)
            btn.setToolTip(tool)
            self.nodeChart.toolbarLayout.insertWidget(i, btn)
            i += 1

        # add an attribute for the project manager
        self.nodeChart.project = project

        # add an attribute for the mfixgui
        self.mfixgui = parent
        self.nodeChart.workflow_widget = self
        self.nodeChart.mfixgui = parent

        # Build defualt node library
        self.nodeChart.nodeLibrary.buildDefaultLibrary()

        # Add custom Nodes
        for node in [TestNode]:
            self.nodeChart.nodeLibrary.addNode(node, ['MFIX', ])

        # --- initialize job status table ---
        self.job_frame = QtWidgets.QWidget()
        self.job_layout = QtWidgets.QVBoxLayout(self.job_frame)
        self.job_layout.setContentsMargins(0, 0, 0, 0)
        self.job_layout.setSpacing(0)
        self.job_frame.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Preferred)
        self.job_frame.setLayout(self.job_layout)

        self.job_toolbar = QtWidgets.QWidget()
        self.job_toolbar.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                       QtWidgets.QSizePolicy.Fixed)
        self.job_toolbar_layout = QtWidgets.QHBoxLayout(self.job_toolbar)
        self.job_toolbar_layout.setContentsMargins(0, 0, 0, 0)
        self.job_toolbar_layout.setSpacing(0)
        self.job_toolbar.setLayout(self.job_toolbar_layout)
        self.job_layout.addWidget(self.job_toolbar)

        self.tool_btn_dict = OrderedDict()
        for tool, icon, callback in [
                ('play', 'play.png', self.handle_play),
                ('stop', 'stop.png', self.handle_stop),
                ('pause', 'pause.png', self.handle_pause),
                ('restart', 'restart.png', self.handle_restart),
                ('auto restart', 'autorenew.png', self.handle_renew),
                ('remove from queue', 'removefromqueue.png', self.handle_remove_from_queue),
                ('submit to queue', 'addtoqueue.png', self.handle_add_to_queue),
                ('open', 'folder.png', self.handle_open)]:
            btn = QtWidgets.QToolButton()
            btn.setIcon(get_icon(icon))
            btn.pressed.connect(callback)
            btn.setAutoRaise(True)
            btn.setEnabled(False)
            btn.setToolTip(tool)
            self.tool_btn_dict[tool] = btn
            self.job_toolbar_layout.addWidget(btn)
        self.job_toolbar_layout.addStretch()

        self.job_status_table = Table(
            dtype=OrderedDict,
            columns=['status', 'progress', 'dt', 'time remaining', 'path'],
            column_delegate={1: {'widget': 'progressbar'}},
            multi_selection=True
            )
        self.job_status_table.set_value(OrderedDict())
        self.job_status_table.show_vertical_header(True)
        self.job_status_table.auto_update_rows(True)
        self.job_status_table.default_value = OrderedDict()
        self.job_status_table.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                            QtWidgets.QSizePolicy.Preferred)
        self.job_status_table.new_selection.connect(self.update_btns)
        self.job_layout.addWidget(self.job_status_table)

        # splitter
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.splitter.addWidget(self.nodeChart)
        self.splitter.addWidget(self.job_frame)

        # main layout
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        self.setLayout(self.layout)
        self.layout.addWidget(self.splitter)

    def set_save_btn(self):
        self.mfixgui.set_unsaved_flag()

    @property
    def parameter_dict(self):
        return PARAMETER_DICT

    @property
    def used_parameters(self):
        return self.mfixgui.project.parameter_key_map.keys()

    def export_project(self, path=None, param_dict={}, keyword_dict={}):
        """
        export a mfix project

        :path: directory to export the project to
        :param_dict: dictionary of parameters and values to use {'x':1.3}
        :keyword_dict: dictionary of keywords and values {'BC_V_g,1': 5.0}
        """
        # copy parameters
        param_copy = copy.deepcopy(PARAMETER_DICT)

        # copy project
        proj = copy.deepcopy(self.mfixgui.project)

        # copy files
        proj_dir = self.mfixgui.get_project_dir()
        files_to_copy = glob.glob(os.path.join(proj_dir, '*.stl'))
        for f in files_to_copy:
            shutil.copyfile(f, os.path.join(path, os.path.basename(f)))

        # change parameters
        PARAMETER_DICT.update(param_dict)

        # change keywords
        for key_args, value in keyword_dict.items():
            key_args = key_args.split(',')
            key = key_args[0]
            if len(key_args) > 1:
                args = [int(arg) for arg in key_args[1:]]
            else:
                args = []
            proj.updateKeyword(key, value, args=args)

        copied_proj = os.path.join(path, proj.run_name.value+'.mfx')
        self.mfixgui.print_internal("Exporting to: %s" % copied_proj,
                                    color='green')
        proj.writeDatFile(copied_proj)

        # reset parameters
        PARAMETER_DICT.update(param_copy)
        for key in set(PARAMETER_DICT.keys())-set(param_copy.keys()):
            PARAMETER_DICT.pop(key)

        return copied_proj

    def run_popup(self, mfx_file):

        proj_dir = os.path.dirname(mfx_file)

        parent = MockParent()
        parent.mfix_gui = self.mfixgui
        parent.project_dir = proj_dir
        parent.project_file = mfx_file
        parent.project = self.project
        parent.settings = self.mfixgui.settings
        run_dialog = WorkflowRunPopup(None, parent)
        run_dialog.exec_()

        self.run_cmd = run_dialog.get_run_cmd()

    def run_project(self, mfx_file):
        """
        Run the mfix project

        :mfx_file: path to mfx project file
        """

        proj_dir = os.path.dirname(mfx_file)
        dir_base = os.path.basename(proj_dir)

        if self.run_cmd is None:
            self.run_popup(mfx_file)

        data = self.job_status_table.value
        data[dir_base] = {'status':'waiting for pid', 'progress':0, 'path':proj_dir, 'dt':'None', 'time remaining':'None'}
        self.job_status_table.set_value(data)
        self.mfixgui.print_internal("Starting: %s" % str(proj_dir), color='green')

        if not os.path.exists(mfx_file):
            self.mfixgui.print_internal("Error: No project file: %s" % proj_dir)
            return

        parent = MockParent()
        parent.mfix_gui = self.mfixgui
        parent.project_dir = proj_dir
        parent.project_file = mfx_file
        parent.project = Project(mfx_file)
        parent.settings = self.mfixgui.settings
        parent.stderr_signal.connect(self.job_error)
        self.mock_parents.append(parent)
        run_dialog = WorkflowRunPopup(None, parent)

        # Queue
        if False:
            msg = 'Submitting to Queue: %s' % ' '.join(self.run_cmd)
            run_dialog.submit_command(self.run_cmd)
        # local
        else:
            msg = 'Run Command: %s' % ' '.join(self.run_cmd)
            run_dialog.start_command(self.run_cmd, proj_dir, os.environ)

        self.file_watcher.addPath(proj_dir)
        self.mfixgui.print_internal(msg, color='green')

        if not self.update_timer.isActive():
            self.update_timer.start(1000)

    def create_job_manager(self, proj_dir):
        pid_files = glob.glob(os.path.join(proj_dir, '*.pid'))
        if pid_files:
            if len(pid_files) > 1:
                self.mfixgui.print_internal('more than one pid file', color='red')

            dir_base = os.path.basename(proj_dir)
            mfx_files = glob.glob(os.path.join(proj_dir, '*.mfx'))

            parent = MockParent()
            parent.mfix_gui = self.mfixgui
            parent.project_dir = proj_dir
            parent.project_file = mfx_files[0]
            parent.project = Project(mfx_files[0])
            parent.settings = self.mfixgui.settings

            full_runname_pid = os.path.join(proj_dir, pid_files[0])
            job = JobManager(parent)
            job.try_to_connect(full_runname_pid)

            self.job_dict[dir_base] = job
            self.file_watcher.removePath(proj_dir)

    def job_error(self, error):
        self.mfixgui.print_internal(error, color='red')

    def save(self, fname):
        """save a node chart file at the given path"""
        self.nodeChart.save(path=fname)

    def load(self, fname):
        """load node chart file"""
        self.nodeChart.open(path=fname)

    def clear(self):
        """clear all nodes"""
        self.nodeChart.deleteAllNodes(confirm=False)

    # --- job update ---
    def update_job_status(self):
        """update the current job status"""

#        {   u'dt': 0.0001234567901234568,
#    u'nit': 18,
#    u'paused': True,
#    u'profiling': [   [   u'time_step_init',
#                          0.0012989044189453125],
#                      [   u'do_iteration',
#                          18,
#                          0.01675891876220703],
#                      [   u'time_step_end',
#                          1.2874603271484375e-05],
#                      [   u'des_time_init',
#                          u'unknown'],
#                      [   u'des_time_steps',
#                          0,
#                          u'unknown'],
#                      [   u'des_time_end',
#                          u'unknown']],
#    u'residuals': [   [u'', u'0.0'],
#                      [   u'P0  ',
#                          u'0.000817730624263'],
#                      [   u'P1  ',
#                          u'9.72339499832e-06'],
#                      [   u'U0  ',
#                          u'3.51669989309e-08'],
#                      [   u'V0  ',
#                          u'1.85474211159e-06'],
#                      [   u'U1  ',
#                          u'5.18117609376e-08'],
#                      [   u'V1  ',
#                          u'1.85322488506e-06'],
#                      [u'    ', u'0.0']],
#    u'time': 0.0012666666666666666,
#    u'tstop': 2.0,
#    u'walltime_elapsed': 3.575956106185913,
#    u'walltime_remaining': u'5642.67052735'}

        data = self.job_status_table.value

        for job_name in data.keys():
            if job_name in self.job_dict:

                job = self.job_dict[job_name].job

                if data[job_name]['status'] != 'stopped':
                    if job.is_paused():
                        status = 'paused'
                    else:
                        status = 'running'

                    if job.status:
                        data[job_name]['status'] = status
                        data[job_name]['progress'] = job.status['time']/job.status['tstop']*100
                        data[job_name]['dt'] = '{0:.2e}'.format(job.status['dt'])
                        try:
                            data[job_name]['time remaining'] = '{0:.0f}'.format(float(job.status['walltime_remaining']))
                        except:
                            pass

        self.job_status_table.set_value(data)

    # --- job managment ---
    def get_selected_jobs(self):
        """get the currently selected jobs"""
        projs = list(self.job_status_table.value.keys())
        return [self.job_dict[projs[i]] for i in self.job_status_table.current_rows()]

    def get_selected_projects(self):
        """get the currently selected project names"""
        projs = list(self.job_status_table.value.keys())
        return [projs[i] for i in self.job_status_table.current_rows()]

    def update_btns(self):
        """enable/diable btns"""

        enable_list = [False]*len(self.tool_btn_dict)

        projs = self.get_selected_projects()

        if projs:
            enable_list[:3] = [True]*4
            enable_list[-1] = True
#        elif len(projs) == 1:
#

        for enable, btn in zip(enable_list, self.tool_btn_dict.values()):
            btn.setEnabled(enable)

    def handle_play(self):
        """play the selected job"""
        jobs = self.get_selected_jobs()
        for job in jobs:
            job.job.unpause()

    def handle_stop(self):
        """stop the selected job"""
        projs = self.get_selected_projects()
        data = self.job_status_table.value
        for proj in projs:
            job = self.job_dict[proj]
            job.stop_mfix()
            data[proj]['status'] = 'stopped'
        self.job_status_table.set_value(data)

    def handle_pause(self):
        """pause the selected job"""
        jobs = self.get_selected_jobs()
        for job in jobs:
            job.job.pause()

    def handle_restart(self):
        """restart the selected job"""
        jobs = self.get_selected_jobs()
        for job in jobs:
            job.stop_mfix()

    def handle_remove_from_queue(self):
        """remove job from queue"""
        projs = self.get_selected_projects()
        print('remove')

    def handle_add_to_queue(self):
        """add job to queue"""
        projs = self.get_selected_projects()
        print('add')

    def handle_renew(self):
        """auto restart the selected job"""
        projs = self.get_selected_projects()
        print('auto restart')

    def handle_open(self):
        """open the selected job"""
        projs = self.get_selected_projects()
        print('open')

    def handle_import(self):
        """immport a nc file"""
        print('import')

    def handle_export(self):
        """export a nc file"""
        print('export')
