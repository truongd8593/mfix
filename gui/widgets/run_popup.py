#!/usr/bin/env python

import os
import sys
import signal
from subprocess import Popen, PIPE
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

        self.recent_exe_limit = 5
        self.mfix_exe = False
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

        # create initial executable list
        self.generate_exe_list()
        self.set_run_mfix_exe.emit()

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

        try:
            ui.spinbox_keyword_nodesi.setValue(project.get_value('nodesi'))
            ui.spinbox_keyword_nodesj.setValue(project.get_value('nodesj'))
            ui.spinbox_keyword_nodesk.setValue(project.get_value('nodesk'))
            self.NODES_SET = True
        except:
            self.NODES_SET = False

        ui.button_browse_exe.clicked.connect(self.handle_browse_exe)
        ui.combobox_mfix_exe.currentIndexChanged.connect(self.handle_exe_change)

        buttons = self.ui.buttonbox.buttons()
        self.ok_button = buttons[0]
        self.ok_button.clicked.connect(self.handle_run)
        buttons[1].clicked.connect(self.handle_abort)



    # UI update functions

    def populate_combobox_mfix_exe(self):
        """ Add items from self.exe_list to combobox, select the last item """
        if not self.exe_list:
            self.ui.combobox_mfix_exe.setEnabled(False)
            return
        self.ui.combobox_mfix_exe.clear()
        self.ui.combobox_mfix_exe.addItems(self.exe_list)
        self.ui.combobox_mfix_exe.setCurrentIndex(0)
        self.update_run_options()
        # TODO: handle anything gui.py is doing with mfix_exe
        self.set_run_mfix_exe.emit()

    def update_run_options(self):
        """ Enable or disable run options based on self.mfix_exe features,
        local or remote settings """

        mfix_exe_available = self.mfix_exe if not self.mfix_exe else True
        self.ui.buttonbox.buttons()[0].setEnabled(mfix_exe_available)
        self.ui.combobox_mfix_exe.setEnabled(mfix_exe_available)

        mode = 'local' # TODO: base element status on local vs queue
        queue_enabled = mode == 'queue'
        self.ui.groupbox_queue_options.setEnabled(queue_enabled)

        cfg = self.mfix_exe_flags.get(self.mfix_exe, None) 
        dmp = 'dmp' in cfg if cfg else False
        smp = 'smp' in cfg if cfg else False
        #dmp = smp = True
        self.ui.spinbox_keyword_nodesi.setEnabled(dmp)
        self.ui.spinbox_keyword_nodesj.setEnabled(dmp)
        self.ui.spinbox_keyword_nodesk.setEnabled(dmp)
        self.ui.spinbox_threads.setEnabled(smp)

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
        # TODO: handle anything gui.py is doing with mfix_exe
        self.set_run_mfix_exe.emit()
        self.run.emit()

    def handle_exe_change(self):
        """ emit signals when exe combobox changes """
        # get exe config features (smp/dmp)
        # update run options (enable/disable NODES*)
        self.mfix_exe = new_exe = self.ui.combobox_mfix_exe.currentText()
        self.append_to_exe_list(new_exe)
        self.persist_selected_exe(new_exe)
        # TODO: handle anything gui.py is doing with mfix_exe
        self.set_run_mfix_exe.emit()

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable", directory=self.project_dir)
        self.mfix_exe = new_exe
        self.append_to_exe_list(new_exe)
        self.persist_selected_exe(new_exe)
        self.populate_combobox_mfix_exe()
        # TODO: handle anything gui.py is doing with mfix_exe
        self.set_run_mfix_exe.emit()



    # utils

    def persist_selected_exe(self, new_exe):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        # save new exe list, set mfix_exe in settings, set mfix_exe in project
        self.settings.setValue('mfix_exe', new_exe)
        self.gui_comments['mfix_exe'] = new_exe
        self.settings.setValue('recent_executables', ','.join(self.exe_list))

    def append_to_exe_list(self, exe):
        """ Verify exe exists, is executable, and appears only once in list.
        Truncate to 5 items """
        exe_list = self.exe_list
        exe_list.reverse()
        if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
            return exe_list
        if exe in exe_list:
            exe_list.pop(exe_list.index(exe))

        exe_list.append(exe)
        self.update_exe_flags(exe)

        if len(exe_list) >= self.recent_exe_limit:
            drop_exe_list = exe_list[5:]
            exe_list = exe_list[:4]
            # cull dropped exes from self.mfix_exe_flags
            for exe in drop_exe_list:
                self.mfix_exe_flags.pop(exe)

        exe_list.reverse()
        self.mfix_exe = exe_list[0]
        self.exe_list = exe_list

    def generate_exe_list(self):
        """ assemble list of executables from:
        ? command line
        - config item 'recent_executables'
        - project file 'mfix_exe'
        - project dir
        - default install location
        """
        self.exe_list = []

        # default install location(s)
        # TODO: where will the default binaries be installed?
        #for location in default_install_dirs:
        #    for name in ['mfix', 'mfix.exe', 'pymfix', 'pymfix.exe']:
        #        for exe in glob(os.path.join(self.project_dir, name)):
        #            exe = os.path.abspath(exe)
        #            self.append_to_exe_list(exe)

        # project dir executable
        for name in ['mfix', 'mfix.exe', 'pymfix', 'pymfix.exe']:
            for exe in glob(os.path.join(self.project_dir, name)):
                exe = os.path.abspath(exe)
                self.append_to_exe_list(exe)
                self.update_exe_flags(exe)

        # recently used executables
        recent_list = self.settings.value('recent_executables')
        if recent_list:
            for recent_exe in recent_list.split(','):
                self.append_to_exe_list(recent_exe)
                self.update_exe_flags(recent_exe)

        # project executable
        project_exe = self.project.get_value('mfix_exe')
        if project_exe:
            self.append_to_exe_list(project_exe)
            self.update_exe_flags(project_exe)

        # no exe found
        if not self.exe_list:
            self.parent.message(
                icon='warning',
                text='MFIX executable not found. Please browse for an executable.',
                buttons=['ok','cancel'],
                default='ok')


    def update_exe_flags(self, mfix_exe):
        """ run mfix to get executable features (like dmp/smp support) """
        cache = self.mfix_exe_flags
        try: # Possible race, file may have been deleted/renamed since isfile check!
            stat = os.stat(mfix_exe)
        except OSError as err:
            log.exception("could not run %s --print-flags", mfix_exe)

        cached_stat, cached_flags = cache.get(mfix_exe, (None, None))
        if cached_stat and cached_stat == stat:
            return cached_flags

        exe_dir = os.path.dirname(mfix_exe)
        popen = Popen(mfix_exe + " --print-flags",
                      cwd=exe_dir,
                      stdout=PIPE,
                      stderr=PIPE,
                      shell=True) # needed?
        (out, err) = popen.communicate()
        flags = '' if err else out.strip()
        cache[mfix_exe] = (stat, flags)
        return flags




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
