# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
This module contains the work flow widget.
'''
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtCore, QtWidgets
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
    Node = None

# local imports
from widgets.base import Table
from tools.general import make_callback, get_icon
from constants import PARAMETER_DICT

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
             ])

        Node.__init__(self)

        self.terminals['used parameters'].valueChanged.connect(self.used_parameters)
        self.terminals['export'].valueChanged.connect(self.export)

    def used_parameters(self):
        print(self.parent.workflow_widget.used_parameters)

    def export(self):

        curr_proj_dir = self.parent.mfixgui.get_project_dir()
        exp_path = os.path.join(curr_proj_dir, 'test')
        if not os.path.exists(exp_path):
            os.mkdir(exp_path)

        self.parent.workflow_widget.export_project(
            exp_path, # export path
            {'x': 1, 'y': 2, 'z': 1}, # parameters
            {'drag_type': 'WEN_YU', 'BC_V_g,1': 5.0}, # keywords
            )


# --- Workflow Widget ---
class WorkflowWidget(QtWidgets.QWidget):
    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.project = project

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

        self.tool_btn_dict = {}
        for tool, icon, callback in [('stop', 'stop.png', self.handle_stop_job),
                                     ('open', 'folder.png', self.handle_open_job)]:
            btn = QtWidgets.QToolButton()
            btn.setIcon(get_icon(icon))
            btn.pressed.connect(callback)
            btn.setAutoRaise(True)
            btn.setToolTip(tool)
            self.tool_btn_dict[tool] = btn
            self.job_toolbar_layout.addWidget(btn)
        self.job_toolbar_layout.addStretch()

        self.job_status_table = Table(dtype=OrderedDict,
                                      columns=['job', 'status', 'progress'])
        self.job_status_table.set_value(OrderedDict())
        self.job_status_table.show_vertical_header(True)
        self.job_status_table.auto_update_rows(True)
        self.job_status_table.default_value = OrderedDict()
        self.job_status_table.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                            QtWidgets.QSizePolicy.Preferred)
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
        self.mfixgui.print_internal("Exporting to: %s" % copied_proj, color='blue')
        proj.writeDatFile(copied_proj)

        # reset parameters
        PARAMETER_DICT.update(param_copy)

    def save(self, fname):
        """save a node chart file at the given path"""
        self.nodeChart.save(path=fname)

    def load(self, fname):
        """load node chart file"""
        self.nodeChart.open(path=fname)

    def clear(self):
        """clear all nodes"""
        self.nodeChart.deleteAllNodes(confirm=False)

    def handle_stop_job(self):
        """stop the selected job"""
        print('stop')

    def handle_open_job(self):
        """open the selected job"""
        print('open')

    def handle_import(self):
        """immport a nc file"""
        print('import')

    def handle_export(self):
        """export a nc file"""
        print('export')
