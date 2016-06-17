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
    mfix_exe_changed = QtCore.Signal()
    SAVED_EXE_LIMIT = 5

    def __init__(self, project, settings, title, parent=None):

        super(RunPopup, self).__init__(parent)

        self.project = project
        self.settings = settings
        self.mfix_exe = ''
        if parent is not None:
            self.parent = parent
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = self.ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)

        self.setWindowTitle(title)
        self.populate_combobox_mfix_exes()

        if "OMP_NUM_THREADS" in os.environ:
            self.ui.spinbox_threads.setValue(int(os.environ['OMP_NUM_THREADS']))
        try:
            self.ui.spinbox_keyword_nodesi.setValue(project.get_value('nodesi'))
            self.ui.spinbox_keyword_nodesj.setValue(project.get_value('nodesj'))
            self.ui.spinbox_keyword_nodesk.setValue(project.get_value('nodesk'))
            self.NODES_SET = True
        except:
            self.NODES_SET = False

        self.ui.button_browse_exe.clicked.connect(self.handle_browse_exe)

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(self.handle_abort)

    def populate_combobox_mfix_exes(self):
        exe_list = self.get_saved_exe_list()
        try:
            project_exe = self.project.get_value('mfix_exe')
            exe_list.append(project_exe)
        except:
            pass
        self.ui.combobox_mfix_exes.clear()
        self.ui.combobox_mfix_exes.addItems(exe_list)

    def set_saved_exe_list(self, exe_list):
        new_list = []
        for exe in exe_list:
            if type(exe) is not 'str':
                continue
            if len(exe) == 0:
                continue
            new_list.append(exe)
        self.settings.setValue('saved_mfix_exes', ','.join(new_list))

    def get_saved_exe_list(self):
        saved_exes = self.settings.value('saved_mfix_exes')
        if saved_exes is not None:
            saved_exes = saved_exes.split(',')
        print(saved_exes)
        return saved_exes if saved_exes else []

    def saved_exe_append(self, new_exe):
        saved_exes = get_saved_exe_list()
        if not len(saved_exes) <= SAVED_EXE_LIMIT:
            saved_exes = saved_exes[:4]
        if new_exe not in saved_exes:
            saved_exes.append(new_exe)
        set_saved_exe_list(saved_exes)

    def update_combobox_mfix_exes(self, new_exe):
        saved_exes = self.get_saved_exe_list()
        saved_exes.append(new_exe)
        self.set_saved_exe_list(saved_exes)
        self.ui.combobox_mfix_exes.addItem(new_exe)

    def handle_browse_exe(self):
        """Handle file open dialog for user specified exe"""
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable")
        self.mfix_exe = new_exe
        self.update_combobox_mfix_exes(new_exe)
        # update mfix_exe in main loop
        #self.parent.mfix_exe = new_exe
        # connected to gui.py:handle_mfix_exe_changed()
        self.mfix_exe_changed.emit()

    def handle_abort(self):
        self.cancel.emit()

    def handle_run(self):
        os.environ['OMP_NUM_THREADS'] = str(self.ui.spinbox_threads.value())
        if self.NODES_SET:
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
    app = QtWidgets.QApplication(args)
    run_popup = Run_Popup(QtWidgets.QDialog())
    run_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    app.exec_()
    app.deleteLater()

    sys.exit()
