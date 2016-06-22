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

    def __init__(self, project, settings, title, parent=None):

        super(RunPopup, self).__init__(parent)

        self.SAVED_EXE_LIMIT = 5
        self.project = project
        self.settings = settings
        self.mfix_exe = ''
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
        self.ui.combobox_mfix_exes.currentIndexChanged.connect(self.handle_exe_change)

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(self.handle_abort)

    # event handlers
    def handle_abort(self):
        self.cancel.emit()

    def handle_run(self):
        thread_count = str(self.ui.spinbox_threads.value())
        os.environ['OMP_NUM_THREADS'] = thread_count
        self.project.mfix_gui_comments['OMP_NUM_THREADS'] = thread_count
        if self.NODES_SET:
            self.project.updateKeyword('nodesi', self.ui.spinbox_keyword_nodesi.value())
            self.project.updateKeyword('nodesj', self.ui.spinbox_keyword_nodesj.value())
            self.project.updateKeyword('nodesk', self.ui.spinbox_keyword_nodesk.value())
        self.mfix_exe_changed.emit()
        self.run.emit()

    def handle_exe_change(self):
        """ emit signals when exe combobox changes """
        # get exe config features (smp/dmp)
        # update run options (enable/disable NODES*)
        self.mfix_exe = new_exe = self.ui.combobox_mfix_exes.currentText()

        # save new exe list, set mfix_exe in settings, set mfix_exe in project
        self.saved_exe_append(new_exe)
        self.set_config_mfix_exe(new_exe)
        self.set_project_exe(new_exe)
        self.mfix_exe_changed.emit()

    def handle_browse_exe(self):
        """Handle file open dialog for user specified exe"""
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable", directory=self.parent.get_project_dir())
        self.mfix_exe = new_exe
        self.update_combobox_mfix_exes(new_exe)
        # save new exe list, set mfix_exe in settings, set mfix_exe in project
        self.saved_exe_append(new_exe)
        self.set_config_mfix_exe(new_exe)
        self.set_project_exe(new_exe)
        # connected to gui.py:handle_mfix_exe_changed()
        print('changed exe: %s' % new_exe)
        self.mfix_exe_changed.emit()

    # utils
    def set_saved_exe_list(self, exe_list):
        new_list = self.verify_exe_list(exe_list)
        print('set saved exes: %s' % new_list)
        self.settings.setValue('saved_mfix_exes', ','.join(new_list))

    def set_project_exe(self, new_exe):
        self.project.mfix_gui_comments['mfix_exe'] = new_exe

    def set_config_mfix_exe(self, new_exe):
        self.settings.setValue('mfix_exe', new_exe)

    def verify_exe_list(self, exe_list):
        """ check that all items in list are found in filesystem """
        new_list = []
        for exe in exe_list:
            if os.path.exists(exe):
                new_list.append(exe)
        return new_list

    def get_exe_list(self):
        """ assemble list of executables from:
        - config item 'saved_mfix_exes' (ordered by last used)
        - project file 'mfix_exe'
        ? command line
        ? 
        """
        exe_list = self.settings.value('saved_mfix_exes')
        if exe_list is not None:
            exe_list = self.verify_exe_list(exe_list.split(','))
        else:
            exe_list = []
        project_exe = self.project.get_value('mfix_exe')
        if project_exe:
            exe_list.append(project_exe)
        print('get exes: %s' % exe_list)
        exe_list = self.verify_exe_list(exe_list)
        return exe_list if exe_list else []

    def saved_exe_append(self, new_exe):
        print('append new: %s' % new_exe)
        saved_exes = self.get_exe_list()
        if len(saved_exes) >= self.SAVED_EXE_LIMIT:
            saved_exes = saved_exes[:3]
        if new_exe not in saved_exes:
            saved_exes.append(new_exe)
        self.set_saved_exe_list(saved_exes)

    # UI update functions
    def populate_combobox_mfix_exes(self):
        exe_list = self.get_exe_list()
        print('current saved: %s' % exe_list)
        self.ui.combobox_mfix_exes.clear()
        self.ui.combobox_mfix_exes.addItems(exe_list)
        print('parent.mfix_exe: %s' % self.parent.mfix_exe)

    def set_combobox_active_item(self, selected):
        pass

    def update_combobox_mfix_exes(self, new_exe):
        self.ui.combobox_mfix_exes.addItem(new_exe)

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
