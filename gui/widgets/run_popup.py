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

    def __init__(self, project, settings, title, parent=None):

        super(RunPopup, self).__init__(parent)

        self.project = project
        self.settings = settings
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

        self.ui.button_browse_exe.clicked.connect(self.handle_browse_exe)

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(lambda: self.cancel.emit())

    def handle_browse_exe(self):
        """Handle file open dialog for user specified exe"""
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable")
        self.update_combobox_mfix_exes(new_exe)

    def set_saved_exe_list(self, exe_list):
        self.settings.setValue('saved_mfix_exes', ','.join(exe_list))

    def get_saved_exe_list(self):
        saved_exes = self.settings.value('saved_mfix_exes')
        if saved_exes is not None:
            saved_exes = saved_exes.split(',')
        return saved_exes if saved_exes else []

    def update_combobox_mfix_exes(self, new_exe):
        saved_exes = self.get_saved_exe_list()
        saved_exes.append(new_exe)
        self.set_saved_exe_list(saved_exes)
        self.ui.combobox_mfix_exes.addItem(new_exe)
        #TODO: change the global mfix_exe in gui

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
