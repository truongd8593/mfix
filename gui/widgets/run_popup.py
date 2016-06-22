#!/usr/bin/env python

import os
import sys
import signal
from glob import glob

from qtpy import QtWidgets, QtCore

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

class RunPopup(QtWidgets.QDialog):

    run = QtCore.Signal()
    cancel = QtCore.Signal()
    set_run_mfix_exe = QtCore.Signal()

    def __init__(self, title, parent):

        super(RunPopup, self).__init__(parent)
        print(parent.__repr__())

        self.recent_exe_limit = 5
        self.mfix_exe = ''
        self.mfix_exe_flags = {}
        self.parent = parent
        self.project = parent.project
        self.settings = parent.settings
        self.project_dir = parent.get_project_dir()
        self.gui_comments = self.project.mfix_gui_comments

        # load ui
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)
        self.setWindowTitle(title)

        # set initial UI element values
        self.populate_combobox_mfix_exe()

        # set OMP_NUM_THREADS
        project_threads = self.gui_comments.get('OMP_NUM_THREADS', None)
        env_threads = os.environ.get('OMP_NUM_THREADS', None)
        if project_threads:
            ui.spinbox_threads.setValue(int(project_threads))
        elif env_threads:
            ui.spinbox_threads.setValue(int(env_threads))
        else:
            ui.spinbox_threads.setValue(1)

        # FIXME: enable/disable based on smp/dmp
        try:
            ui.spinbox_keyword_nodesi.setValue(project.get_value('nodesi'))
            ui.spinbox_keyword_nodesj.setValue(project.get_value('nodesj'))
            ui.spinbox_keyword_nodesk.setValue(project.get_value('nodesk'))
            self.NODES_SET = True
        except:
            self.NODES_SET = False

        ui.button_browse_exe.clicked.connect(self.handle_browse_exe)
        ui.combobox_mfix_exe.currentIndexChanged.connect(self.handle_exe_change)

        buttons = self.ui.buttonBox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(self.handle_abort)



    # UI update functions

    def populate_combobox_mfix_exe(self):
        """ add items from generate_exe_list() to combobox """
        # FIXME: handle generate_exe_list being empty (ie during first run)
        #  - if empty, search for exe
        #  - if no exe found, disable combobox and present popup warning
        exe_list = self.generate_exe_list()
        print(exe_list)
        self.ui.combobox_mfix_exe.addItems(exe_list)
        self.ui.combobox_mfix_exe.setCurrentText(exe_list[-1])
        self.set_run_mfix_exe.emit()

    def update_combobox_mfix_exe(self, new_exe):
        self.ui.combobox_mfix_exe.addItem(new_exe)
        self.ui.combobox_mfix_exe.setCurrentText(new_exe)
        self.set_run_mfix_exe.emit()

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
        self.gui_comments['OMP_NUM_THREADS'] = thread_count
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
        self.mfix_exe = new_exe = self.ui.combobox_mfix_exe.currentText()
        self.persist_selected_exe(new_exe)

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable", directory=self.project_dir)
        self.mfix_exe = new_exe
        self.update_combobox_mfix_exe(new_exe)
        self.persist_selected_exe(new_exe)



    # utils

    def persist_selected_exe(self, new_exe):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        # save new exe list, set mfix_exe in settings, set mfix_exe in project
        self.append_recent_exe(new_exe)
        self.set_config_mfix_exe(new_exe)
        self.set_project_exe(new_exe)
        self.set_run_mfix_exe.emit()

    def set_saved_exe_list(self, exe_list):
        new_list = self.verify_exe_list(exe_list)
        self.settings.setValue('recent_executables', ','.join(new_list))

    def set_project_exe(self, new_exe):
        self.gui_comments['mfix_exe'] = new_exe

    def set_config_mfix_exe(self, new_exe):
        self.settings.setValue('mfix_exe', new_exe)

    def generate_exe_list(self, new_exe=None):
        """ assemble list of executables from:
        - argument 'new_exe'
        - config item 'recent_executables'
        - project file 'mfix_exe'
        - project dir
        - default install location
        ? command line
        """
        exe_list = []

        def append(exe, exe_list):
            print('append %s' % exe)
            new_list = exe_list
            """ verify exe exists and is executable, append to list """
            #if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
            #    return exe_list
            if exe in new_list:
                new_list.pop(new_list.index(exe))
                new_list.append(exe)
            if len(exe_list) >= self.recent_exe_limit:
                new_list = new_list[-4:]
            return new_list

        # default install location
        # TODO: where will the default binaries be installed?

        # project dir executable
        for name in ['mfix', 'mfix.exe', 'pymfix', 'pymfix.exe']:
            for exe in glob(os.path.join(self.project_dir, name)):
                exe_list = append(os.path.abspath(exe), exe_list)

        # recently used executables
        recent_list = self.settings.value('recent_executables')
        if recent_list:
            for exe in recent_list.split(','):
                exe_list = append(exe, exe_list)

        # project executable
        project_exe = self.project.get_value('mfix_exe')
        if project_exe:
            exe_list = append(project_exe, exe_list)

        # append new exe if provided
        if new_exe:
            exe_list = append(new_exe, exe_list)

        return exe_list

    def get_exe_features(self, exe):
        """ run mfix -f to get executable features (like MP support) """
        pass

    def append_recent_exe(self, new_exe):
        """ Get exe list, append new_exe. Ensure list items are unique. """
        exe_list = self.generate_exe_list(new_exe)
        self.set_saved_exe_list(exe_list)


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
