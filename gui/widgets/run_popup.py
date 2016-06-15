#!/usr/bin/env python

import os
import sys
import signal

from qtpy import QtWidgets, QtCore

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

class RunPopup(QtWidgets.QDialog):

    run = QtCore.Signal()
    cancel = QtCore.Signal()

    def __init__(self, app, parent=None):
        super(RunPopup, self).__init__(parent)

        self.app = app
        thisdir = os.path.abspath(os.path.dirname(__file__))
        datadir = thisdir
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = self.ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(lambda: self.cancel.emit())

    def handle_abort(self):
        pass

    def handle_run(self):
        os.environ['OMP_NUM_THREADS'] = str(self.ui.spinbox_threads.value())
        # set NODES* here
        self.run.emit()

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()
