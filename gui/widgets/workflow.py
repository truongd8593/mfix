# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
This module contains the work flow widget.
'''
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtCore, QtWidgets
from collections import OrderedDict

try:
    from pyqtnode import NodeWidget
    PYQTNODE_AVAILABLE = True
except ImportError:
    NodeWidget = None
    PYQTNODE_AVAILABLE = False

# local imports
from widgets.base import Table
from tools.general import make_callback, get_icon


# --- Workflow Widget ---
class WorkflowWidget(QtWidgets.QWidget):
    def __init__(self, project, parent=None):
        QtWidgets.QWidget.__init__(self, parent)

        self.project = project

        # --- initalize the node widget ---
        self.nodeChart = NodeWidget(showtoolbar=True)

        # add an attribute for the project manager
        self.nodeChart.project = project

        # add an attribute for the mfixgui
        self.nodeChart.mfixgui = parent

        # Build defualt node library
        self.nodeChart.nodeLibrary.buildDefaultLibrary()

        # Add custom Nodes
        for node in []:
            self.nodeChart.nodeLibrary.addNode(node, ['MFIX', ])

        # --- initialize job status table ---
        self.job_frame = QtWidgets.QWidget()
        self.job_layout = QtWidgets.QVBoxLayout(self.job_frame)
        self.job_layout.setContentsMargins(0, 0, 0, 0)
        self.job_layout.setSpacing(0)
        self.job_frame.setLayout(self.job_layout)

        self.job_toolbar = QtWidgets.QWidget()
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
