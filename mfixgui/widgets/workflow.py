# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
This module contains the work flow widget.
'''
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtCore, QtWidgets, PYQT5
from qtpy.QtCore import Signal
from collections import OrderedDict
import copy
import os
import shutil
import glob
import subprocess
import json
import sys
import re

try:
    from pyqtnode import NodeWidget, Node, tools
    PYQTNODE_AVAILABLE = True
except ImportError:
    NodeWidget = None
    PYQTNODE_AVAILABLE = False
    Node = object

# local imports
from mfixgui.widgets.base import Table
from mfixgui.widgets.run_popup import RunPopup
from mfixgui.tools.general import get_icon, SCRIPT_DIRECTORY, is_vnc, replace_with_dict
from mfixgui.constants import PARAMETER_DICT
from mfixgui.job import JobManager
from mfixgui.project import Project


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


class ExportNode(Node):
    name = 'Export'

    def __init__(self):
        self.terminalOpts = OrderedDict([
            ('path', {'widget': 'browse',
                      'in': True,
                      'out': True,
                      'showlabel': True,
                      'dtype': str,
                      }),
            ('parameters', {'widget': None,
                        'in': True,
                        'out': False,
                        'showlabel': True,
                        'dtype': dict,
                        'default':{}
                        }),
             ])

        Node.__init__(self)

    def process(self):

        exp_path = self.terminals['path'].value
        if not os.path.exists(exp_path):
            os.mkdir(exp_path)

        self.proj_file = self.parent.workflow_widget.export_project(
            exp_path, # export path
            self.terminals['parameters'].value, # parameters
            )


# --- Mock Parent for job submission ---
class Mock(object):
    def noop(self, *args, **kwargs):
        return False
    def __getattr__(self, name):
        return self.noop


class MockParent(QtWidgets.QWidget):
    stderr_signal = Signal(object)
    stdout_signal = Signal(object)
    signal_update_runbuttons = Signal(object)
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
        self.submit_queue = False

    # over-rides
    def handle_run(self):
        self.submit_queue = False
        self.finish_with_dialog()

    def handle_submit(self):
        self.submit_queue = True
        self.finish_with_dialog()

class FakeJob(object):
    job=None

# --- Workflow Widget ---
class WorkflowWidget(QtWidgets.QWidget):
    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.project = project
        self.job_dict = {}
        self.update_timer = QtCore.QTimer()
        self.update_timer.timeout.connect(self.update_job_status)
        self.run_cmd = None
        self.submit_cmd = None
        self.queue = False
        self.file_timer = QtCore.QTimer()
        self.file_timer.timeout.connect(self.look_for_pid)
        self.watch_dir_paths = []
        self.mock_parents = []
        self.running_in_vnc = is_vnc()

        # --- initalize the node widget ---
        nc = self.nodeChart = NodeWidget(showtoolbar=False)
        if hasattr(self.nodeChart, 'needsSavedEvent'):
            self.nodeChart.needsSavedEvent.connect(self.set_save_btn)
        self.nodeChart.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Preferred)
        self.nodeChart.setGeometry(0, 0, 100, 1000)

        # add btns to gui toolbar
        ly = parent.ui.horizontallayout_menu_bar
        i = ly.count() - 5

        for btn in [nc.runToolButton, nc.autorunToolButton, nc.stopToolButton, nc.stepToolButton]:
            btn.setIconSize(QtCore.QSize(24, 24))
            ly.insertWidget(i, btn)
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
        for node in [TestNode, ExportNode]:
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
                ('delete', 'delete.png', self.handle_delete),
                ('restart', 'restart.png', self.handle_restart),
                ('auto restart', 'autorenew.png', self.handle_renew),
                ('remove from queue', 'removefromqueue.png', self.handle_remove_from_queue),
                ('submit to queue', 'addtoqueue.png', self.handle_add_to_queue),
                ('open', 'folder.png', self.handle_open),
                ('settings', 'settings.png', self.handle_settings)]:
            btn = QtWidgets.QToolButton()
            btn.setIcon(get_icon(icon))
            btn.pressed.connect(callback)
            btn.setAutoRaise(True)
            btn.setEnabled(tool == 'settings')
            btn.setToolTip(tool)
            self.tool_btn_dict[tool] = btn
            self.job_toolbar_layout.addWidget(btn)
        self.job_toolbar_layout.addStretch()

        self.job_status_table = Table(
            dtype=OrderedDict,
            columns=['status', 'job id', 'progress', 'dt', 'time remaining', 'path'],
            column_delegate={2: {'widget': 'progressbar'}},
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

    def look_for_projects(self, path):
        data = self.job_status_table.value
        for root, dirs, files in os.walk(path):
            for name in files:
                if root != path and name.endswith('.mfx'):
                    dir_base = os.path.basename(root)
                    data[dir_base] = {'status':'waiting for pid', 'progress':0,
                                      'path':root, 'dt':'None', 'time remaining':'None',
                                      'job id': None}
                    run_data_path = os.path.join(root, '.workflow_run_cmd.json')
                    if os.path.exists(run_data_path):
                        with open(run_data_path) as f:
                            d = json.load(f)
                        data[dir_base].update(d)
                    self.create_job_manager(root)
                    self.watch_dir_paths.append(root)

        self.job_status_table.set_value(data)
        if not self.update_timer.isActive():
            self.update_timer.start(1000)
        if not self.file_timer.isActive():
            self.file_timer.start(1000)

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

        run_name = proj.get_value('run_name')
        if run_name is None:
            self.mfixgui.error('The project does not have a run_name, is it a valid proejct?')
            return
        copied_proj = os.path.join(path, run_name+'.mfx')
        self.mfixgui.print_internal("Exporting to: %s" % copied_proj,
                                    color='green')
        proj.writeDatFile(copied_proj)

        # reset parameters
        PARAMETER_DICT.update(param_copy)
        for key in set(PARAMETER_DICT.keys())-set(param_copy.keys()):
            PARAMETER_DICT.pop(key)

        return copied_proj

    def run_popup(self):

        parent = MockParent()
        parent.mfix_gui = self.mfixgui
        parent.project_dir = self.mfixgui.get_project_dir()
        parent.project_file = self.mfixgui.get_project_file()
        parent.project = self.project
        parent.settings = self.mfixgui.settings
        run_dialog = WorkflowRunPopup(None, parent)
        run_dialog.exec_()

        self.queue = run_dialog.submit_queue
        if self.queue:
            self.submit_cmd = run_dialog.get_submit_command()
        else:
            self.submit_cmd = None
        self.run_cmd = run_dialog.get_run_command()

    def run_project(self, mfx_file):
        """
        Run the mfix project

        :mfx_file: path to mfx project file
        """
        proj_dir = os.path.dirname(mfx_file)
        dir_base = os.path.basename(proj_dir)

        data = self.job_status_table.value
        if dir_base in data:
            self.mfixgui.print_internal("Error: Project already submited: %s" % proj_dir)
            return

        if self.run_cmd is None:
            self.run_popup()

        data[dir_base] = {'status':'waiting for pid', 'progress':0,
                          'path':proj_dir, 'dt':'None', 'time remaining':'None',
                          'job id':None}
        self.job_status_table.set_value(data)

        if not os.path.exists(mfx_file):
            self.mfixgui.print_internal("Error: No project file: %s" % proj_dir)
            return

        # run it
        job_id = self._run(proj_dir, mfx_file)

        # save some info on the run
        c = data[dir_base]['cmd'] = copy.deepcopy(self.run_cmd)
        s = data[dir_base]['submit'] = copy.deepcopy(self.submit_cmd)
        data[dir_base]['queue'] = self.queue
        data[dir_base]['file'] = mfx_file
        self.job_status_table.set_value(data)

        save_data = {
            'cmd': c,
            'submit': s,
            'queue': self.queue,
            'file': mfx_file,
            'job id': job_id,
        }

        with open(os.path.join(proj_dir, '.workflow_run_cmd.json'), 'w') as f:
            json.dump(save_data, f)

        if not self.update_timer.isActive():
            self.update_timer.start(1000)
        if not self.file_timer.isActive():
            self.file_timer.start(1000)

    def _run(self, proj_dir, mfx_file):
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
        job_id = None
        if self.queue:
            msg = 'Submitting to Queue'
            job_id = self.submit_to_queue(proj_dir)
        # local
        else:
            msg = 'Running: %s' % ' '.join(self.run_cmd)
            run_dialog.start_command(self.run_cmd, proj_dir, os.environ)

        self.watch_dir_paths.append(proj_dir)
        self.mfixgui.print_internal(msg, color='green')
        return job_id

    def submit_to_queue(self, project_dir):
        script_name = '.qsubmit_script'

        script, sub_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict = self.submit_cmd

        # write the script
        with open(os.path.join(project_dir, script_name), 'w') as f:
            f.write(script)

        replace_dict['SCRIPT'] = script_name

        submit_cmd = replace_with_dict(sub_cmd, replace_dict)

        # submit the job
        self.mfixgui.print_internal("Job submit CMD: {}".format(submit_cmd),
                                   color='blue')

        proc = subprocess.Popen(submit_cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=project_dir)
        out, err = proc.communicate()
        if job_id_regex is not None and out:
            job_id = re.findall(job_id_regex, str(out))
        else:
            job_id = []
        if job_id:
            job_id = job_id[0]
            self.mfixgui.print_internal("Job successfully submitted with job id: {}".format(job_id),
                                       color='blue')
        else:
            self.mfixgui.error('Could not determine job id')
            job_id = None
        if err:
            self.mfixgui.error('Error with submission:\n{}'.format(err))

        dir_base = os.path.basename(project_dir)
        data = self.job_status_table.value
        data[dir_base]['job id'] = job_id
        data[dir_base]['status'] = 'submitted to queue'

        return job_id

    def look_for_pid(self):

        for d in copy.deepcopy(self.watch_dir_paths):
            pid_files = glob.glob(os.path.join(d, '*.pid'))
            if pid_files:
                self.create_job_manager(d)
        if len(self.watch_dir_paths) == 0:
            self.file_timer.stop()

    def create_job_manager(self, proj_dir):
        pid_files = glob.glob(os.path.join(proj_dir, '*.pid'))
        dir_base = os.path.basename(proj_dir)
        if pid_files:
            if len(pid_files) > 1:
                self.mfixgui.print_internal('more than one pid file', color='red')

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
            try:
                self.watch_dir_paths.remove(proj_dir)
            except ValueError:
                # already removed?
                pass
        else:
            self.job_dict[dir_base] = FakeJob()

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

        data = self.job_status_table.value

        for job_name in data.keys():
            if job_name in self.job_dict:

                job = self.job_dict[job_name].job

                if job is None:
                    continue

                if job.is_paused():
                    status = 'paused'
                else:
                    status = 'running'

                if job.status:

                    p = data[job_name]['progress'] = job.status['time']/job.status['tstop']*100
                    if p >= 99:
                        status = 'finished'
                    data[job_name]['dt'] = '{0:.2e}'.format(job.status['dt'])
                    try:
                        data[job_name]['time remaining'] = '{0:.0f}'.format(float(job.status['walltime_remaining']))
                    except:
                        pass
                    data[job_name]['status'] = status

        self.job_status_table.set_value(data)

    # --- job managment ---
    def get_selected_jobs(self):
        """get the currently selected jobs"""
        projs = list(self.job_status_table.value.keys())
        return [self.job_dict[projs[i]] for i in self.job_status_table.current_rows() if self.job_dict[projs[i]].job is not None]

    def get_selected_projects(self):
        """get the currently selected project names"""
        projs = list(self.job_status_table.value.keys())
        return [projs[i] for i in self.job_status_table.current_rows()]

    def update_btns(self):
        """enable/diable btns"""

        n_btns = len(self.tool_btn_dict)
        enable_list = [False]*(n_btns-1)

        projs = self.get_selected_projects()

        if projs:
            enable_list[:4] = [True]*5
        if len(projs) == 1:
            enable_list[n_btns-2] = True

        for enable, btn in zip(enable_list, list(self.tool_btn_dict.values())[:-1]):
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
            if not isinstance(job, FakeJob):
                job.stop_mfix()
            data[proj]['status'] = 'stopped'
        self.job_status_table.set_value(data)

    def handle_pause(self):
        """pause the selected job"""
        jobs = self.get_selected_jobs()
        for job in jobs:
            job.job.pause()

    def handle_delete(self):
        """delete the selected job"""
        projs = self.get_selected_projects()

        btn = self.mfixgui.message(
            text='The selectect directories will be delete.\nContinue?',
            buttons=['yes', 'no'],
            default='no',)

        if btn != 'yes':
            return

        data = self.job_status_table.value
        for proj in projs:

            # make sure the job is stopped
            job = self.job_dict[proj]
            if not isinstance(job, FakeJob):
                job.stop_mfix()
                data[proj]['status'] = 'stopped'

            # remove the dir
            path = data[proj]['path']
            shutil.rmtree(path)

            data.pop(proj)
            self.job_dict.pop(proj)

        self.job_status_table.set_value(data)

    def handle_restart(self):
        """restart the selected job"""
        projs = self.get_selected_projects()
        data = self.job_status_table.value
        for proj in projs:
            p = data[proj]

            # make sure the job is stopped
            job = self.job_dict[proj]
            if not isinstance(job, FakeJob):
                job.stop_mfix()
                p['status'] = 'stopped'

            # look for *.SP? files
            spx_files = glob.glob(os.path.join(p['path'], '*SP?'))

            cmd = copy.deepcopy(p['cmd'])
            if spx_files:
                cmd.append('run_type=restart_1')
            else:
                cmd.append('run_type=restart_2')

            self._run(p['path'], p['file'], cmd, p['queue'])
            p['status'] = 'waiting for pid'

        self.job_status_table.set_value(data)


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
        if not projs:
            return

        data = self.job_status_table.value
        path = data[projs[0]]['path']

        #TODO: this is hard coded to > python gui.py project
        gui_path = os.path.join(SCRIPT_DIRECTORY, 'gui.py')
        cmd = []
        if self.running_in_vnc:
            cmd.append('vglrun')
        cmd += ['python', gui_path, path]
        subprocess.Popen(cmd)

    def handle_settings(self):
        """open the run settings dialog"""
        self.run_popup()

    def handle_import(self):
        """immport a nc file"""
        new_nc = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open a node chart",
            self.mfixgui.get_project_dir(),
            'Node Chart (*.nc);; All Files (*)')
        if not new_nc:
            return
        if PYQT5:
            new_nc = new_nc[0]

        self.nodeChart.open(path=new_nc)
        self.handle_main_menu_hide()

    def handle_export(self):
        """export a nc file"""

        new_nc = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save the node chart",
            self.mfixgui.get_project_dir(),
            'Node Chart (*.nc)')
        if not new_nc:
            return
        if PYQT5:
            new_nc = new_nc[0]

        self.nodeChart.save(path=new_nc)
        self.handle_main_menu_hide()
