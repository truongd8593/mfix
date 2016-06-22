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
    set_run_mfix_exe = QtCore.Signal()

    def __init__(self, project, settings, title, parent=None):

        super(RunPopup, self).__init__(parent)

        self.saved_exe_limit = 5
        self.mfix_exe = ''

        self.project = project
        self.settings = settings
        self.parent = parent
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = self.ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)

        self.setWindowTitle(title)

        # set initial UI element values
        self.populate_combobox_mfix_exes()

        # set OMP_NUM_THREADS
        project_threads = self.project.mfix_gui_comments.get('OMP_NUM_THREADS', None)
        env_threads = os.environ.get('OMP_NUM_THREADS', None)
        if project_threads:
            self.ui.spinbox_threads.setValue(int(project_threads))
        elif env_threads:
            self.ui.spinbox_threads.setValue(int(env_threads))
        else:
            self.ui.spinbox_threads.setValue(1)

        # FIXME: enable/disable based on smp/dmp
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


    # UI update functions
    def populate_combobox_mfix_exes(self):
        """ add items from get_exe_list() to combobox """
        # FIXME: handle get_exe_list being empty (ie during first run)
        #  - if empty, search for exe
        #  - if no exe found, disable combobox and present popup warning
        exe_list = self.get_exe_list()
        self.ui.combobox_mfix_exes.addItems(exe_list)
        self.set_combobox_active_item(exe_list[-1])

    def set_combobox_active_item(self, selected):
        self.ui.combobox_mfix_exes.setCurrentText(selected)
        self.set_run_mfix_exe.emit()

    def update_combobox_mfix_exes(self, new_exe):
        self.ui.combobox_mfix_exes.addItem(new_exe)
        self.set_combobox_active_item(new_exe)

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()


    # event handlers
    def handle_abort(self):
        self.cancel.emit()

    def handle_run(self):
        """ persist run options in project file, then emit run signal """
        thread_count = str(self.ui.spinbox_threads.value())
        os.environ['OMP_NUM_THREADS'] = thread_count
        self.project.mfix_gui_comments['OMP_NUM_THREADS'] = thread_count
        if self.NODES_SET:
            self.project.updateKeyword('nodesi', self.ui.spinbox_keyword_nodesi.value())
            self.project.updateKeyword('nodesj', self.ui.spinbox_keyword_nodesj.value())
            self.project.updateKeyword('nodesk', self.ui.spinbox_keyword_nodesk.value())
        self.set_run_mfix_exe.emit()
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
        self.set_run_mfix_exe.emit()

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable", directory=self.parent.get_project_dir())
        self.mfix_exe = new_exe
        self.update_combobox_mfix_exes(new_exe)
        # save new exe list, set mfix_exe in settings, set mfix_exe in project
        self.saved_exe_append(new_exe)
        self.set_config_mfix_exe(new_exe)
        self.set_project_exe(new_exe)
        print('changed exe: %s' % new_exe)
        self.set_run_mfix_exe.emit()



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
        if exe_list:
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
        if len(saved_exes) >= self.saved_exe_limit:
            saved_exes = saved_exes[:3]
        if new_exe not in saved_exes:
            saved_exes.append(new_exe)
        else:
            saved_exes.pop(saved_exes.index(new_exe))
            saved_exes.append(new_exe)
        self.set_saved_exe_list(saved_exes)



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
