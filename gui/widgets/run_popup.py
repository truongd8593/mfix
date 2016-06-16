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

    def __init__(self, project, title, parent=None):
        super(RunPopup, self).__init__(parent)

        self.project = project
        thisdir = os.path.abspath(os.path.dirname(__file__))
        datadir = thisdir
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = self.ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)

        self.setWindowTitle(title)

        if "OMP_NUM_THREADS" in os.environ:
            self.ui.spinbox_threads.setValue(int(os.environ['OMP_NUM_THREADS']))
        self.ui.spinbox_keyword_nodesi.setValue(project.get_value('nodesi', 1))
        self.ui.spinbox_keyword_nodesj.setValue(project.get_value('nodesj', 1))
        self.ui.spinbox_keyword_nodesk.setValue(project.get_value('nodesk', 1))

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(lambda: self.cancel.emit())

    def handle_abort(self):
        pass

    def handle_run(self):
        os.environ['OMP_NUM_THREADS'] = str(self.ui.spinbox_threads.value())
        self.project.updateKeyword('nodesi', self.ui.spinbox_keyword_nodesi.value())
        self.project.updateKeyword('nodesj', self.ui.spinbox_keyword_nodesj.value())
        self.project.updateKeyword('nodesk', self.ui.spinbox_keyword_nodesk.value())
        self.run.emit()

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()

if __name__ == '__main__':
    args = sys.argv
    qapp = QtWidgets.QApplication(args)
    dialog = QtWidgets.QDialog()
    species_popup = SpeciesPopup(dialog, phases='GL')
    species_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    qapp.exec_()
    qapp.deleteLater()

    sys.exit()
