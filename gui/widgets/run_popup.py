#!/usr/bin/env python

import os
import sys
import signal
import logging
from collections import OrderedDict
from subprocess import Popen, PIPE
from glob import glob

from qtpy import QtWidgets, QtCore, PYQT4, PYQT5

from tools.general import get_mfix_home

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

RECENT_EXE_LIMIT = 5
MFIX_EXE_NAMES = ['mfix', 'mfix.exe', 'pymfix', 'pymfix.exe']

class RunPopup(QtWidgets.QDialog):

    run = QtCore.Signal()
    cancel = QtCore.Signal()
    set_run_mfix_exe = QtCore.Signal()

    def __init__(self, title, mfix_exe, parent):

        super(RunPopup, self).__init__(parent)

        self.commandline_option_exe = mfix_exe if mfix_exe else None
        self.mfix_exe = None
        self.mfix_exe_list = []
        self.mfix_exe_flags = {}
        self.title = title
        self.parent = parent
        self.project = parent.project
        self.settings = parent.settings
        self.project_dir = parent.get_project_dir()
        self.gui_comments = self.project.mfix_gui_comments

        # load ui
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)
        self.initialize_ui()

        ui.button_browse_exe.clicked.connect(self.handle_browse_exe)
        ui.combobox_mfix_exe.currentIndexChanged.connect(self.handle_exe_change)

        buttons = self.ui.buttonbox.buttons()
        buttons[0].clicked.connect(self.handle_run)
        buttons[1].clicked.connect(self.handle_abort)

    # UI update functions

    def initialize_ui(self):

        ui = self.ui

        self.setWindowTitle(self.title)

        # create initial executable list
        self.generate_exe_list()
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

    def populate_combobox_mfix_exe(self):
        """ Add items from self.exe_list to combobox, select the last item """
        if not self.mfix_exe_list:
            self.ui.combobox_mfix_exe.setEnabled(False)
            return
        self.ui.combobox_mfix_exe.clear()
        self.ui.combobox_mfix_exe.addItems(self.mfix_exe_list)
        self.ui.combobox_mfix_exe.setCurrentIndex(0)
        self.update_dialog_options()
        self.set_run_mfix_exe.emit()

    def update_dialog_options(self):
        """ Enable or disable options based on self.mfix_exe features,
        local or remote settings """


        self.update_no_mfix_warning()

        mode = 'local' # TODO: base element status on local vs queue
        queue_enabled = mode == 'queue'
        self.ui.groupbox_queue_options.setEnabled(queue_enabled)

        cfg = self.mfix_exe_flags.get(self.mfix_exe, None)
        dmp = 'dmp' in cfg['flags'] if cfg else False
        smp = 'smp' in cfg['flags'] if cfg else False
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
        self.set_run_mfix_exe.emit()
        self.run.emit()

    def handle_exe_change(self):
        """ emit signals when exe combobox changes """
        self.mfix_exe = new_exe = self.ui.combobox_mfix_exe.currentText()
        self.prepend_to_exe_list(new_exe)
        self.persist_selected_exe(new_exe)
        self.update_dialog_options()
        self.set_run_mfix_exe.emit()

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Executable", directory=self.project_dir, options=QtWidgets.QFileDialog.DontResolveSymlinks)
        if not new_exe:
            return
        if PYQT5:
            new_exe = new_exe[0]
        self.mfix_exe = new_exe
        if not self.update_exe_flags(new_exe):
            self.parent.message(
                icon='warning',
                title='Warning',
                text='The selected file is not executable or is not an MFIX binary')
            return
        self.prepend_to_exe_list(new_exe)
        self.persist_selected_exe(new_exe)
        self.populate_combobox_mfix_exe()
        self.update_dialog_options()
        self.set_run_mfix_exe.emit()



    # utils

    def persist_selected_exe(self, new_exe):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        self.settings.setValue('mfix_exe', new_exe)
        self.gui_comments['mfix_exe'] = new_exe
        exe_list = self.mfix_exe_list
        #exe_list.reverse()
        self.settings.setValue('recent_executables', str(os.pathsep).join(exe_list))

    def prepend_to_exe_list(self, exe):
        """ Verify exe exists, is executable, and appears only once in list.
        Truncate to 5 items """
        exe_list = self.mfix_exe_list
        exe_list.reverse()

        if not self.update_exe_flags(exe):
            return exe_list

        if exe in exe_list:
            exe_list.pop(exe_list.index(exe))

        exe_list.append(exe)
        exe_list.reverse()
        self.mfix_exe_list = exe_list

    def update_no_mfix_warning(self):
        ok = bool(self.mfix_exe)
        self.ui.label_mfix_exe_warning.setVisible(not ok)
        self.ui.combobox_mfix_exe.setEnabled(ok)
        self.ui.buttonbox.buttons()[0].setEnabled(ok)
        if not ok:
            self.parent.print_internal("Warning: no MFIX executables available")

    def generate_exe_list(self):
        """ assemble list of executables from:
        - command line
        - project file 'mfix_exe'
        - project dir
        - config item 'recent_executables'
        - default install location
        """

        def default_install_location():
            # TODO? default install location(s)
            # ... where will the default binaries be installed?
            #for location in default_install_dirs:
            #    for name in ['mfix', 'mfix.exe', 'pymfix', 'pymfix.exe']:
            #        for exe in glob(os.path.join(self.project_dir, name)):
            #            exe = os.path.abspath(exe)
            #            self.prepend_to_exe_list(exe)
            pass

        def recently_used_executables():
            recent_list = self.settings.value('recent_executables')
            if recent_list:
                # limit recently used exes to RECENT_EXE_LIMIT
                recent_list = recent_list.split(os.pathsep)[:RECENT_EXE_LIMIT]
                for recent_exe in recent_list:
                    yield recent_exe

        def project_directory_executables():
            for name in MFIX_EXE_NAMES:
                for exe in glob(os.path.join(self.project_dir, name)):
                    yield os.path.abspath(exe)

        def project_file_executable():
            project_exe = self.project.get_value('mfix_exe')
            if project_exe:
                yield project_exe

        def os_path():
            PATH = os.environ.get("PATH")
            if PATH:
                # using OrderedDict to preserve PATH order
                dirs = OrderedDict.fromkeys(PATH.split(os.pathsep))
            else:
                dirs = OrderedDict()
            for d in dirs.keys():
                # filter out empty strings and current directory from $PATH
                if d and d != os.path.curdir and os.path.isdir(d):
                    for name in MFIX_EXE_NAMES:
                        for exe in glob(os.path.join(d, name)):
                            yield exe

        def mfix_build_directories():
            mfix_home = get_mfix_home()
            build_dir_name = 'build-aux'
            if mfix_home:
                dir_list = [
                    os.path.join(mfix_home, 'bin'),
                    os.path.join(mfix_home, 'build')]
                for d in dir_list:
                    if not os.path.isdir(d):
                        continue
                    for child in os.listdir(d):
                        build_dir = os.path.join(d, child, build_dir_name)
                        if os.path.isdir(build_dir):
                            for name in MFIX_EXE_NAMES:
                                for exe in glob(os.path.join(build_dir, name)):
                                    yield exe

        def command_line_option():
            if self.commandline_option_exe:
                yield self.commandline_option_exe

        exe_list_order = [
            recently_used_executables,
            project_directory_executables,
            project_file_executable,
            os_path,
            mfix_build_directories,
            command_line_option]

        # look for executables in the order listed in exe_list_order
        for exe_spec in exe_list_order:
            for exe in exe_spec():
                self.prepend_to_exe_list(exe)

        # no exe found
        if not self.mfix_exe_list:
            self.parent.message(
                icon='warning',
                text='MFIX executable not found. Please browse for an executable.',
                buttons=['ok','cancel'],
                default='ok')
            return
        self.mfix_exe = self.mfix_exe_list[0]

    def exe_exists(self, exe):
        """ Verify exe exists and is executable """
        return (os.path.isfile(exe) and os.access(exe, os.X_OK))

    def update_exe_flags(self, mfix_exe):
        """ run mfix to get executable features (like dmp/smp support) """
        if not self.exe_exists(mfix_exe):
            return False
        cache = self.mfix_exe_flags
        try: # Possible race, file may have been deleted/renamed since isfile check!
            stat = os.stat(mfix_exe)
        except OSError as err:
            log.exception("could not run %s --print-flags", mfix_exe)
            return False

        if any(mfix_exe.lower().endswith(x)
               for x in ('pymfix', 'pymfix.exe')):
            cache[mfix_exe] = {'stat': stat, 'flags': 'dmp smp'}

        cached = cache.get(mfix_exe, None)
        if cached and cached['stat'] == stat:
            return True

        exe_dir = os.path.dirname(mfix_exe)
        popen = Popen(mfix_exe + " --print-flags",
                      cwd=exe_dir,
                      stdout=PIPE,
                      stderr=PIPE,
                      shell=True) # needed?
        (out, err) = popen.communicate()
        flags = '' if err else str(out.strip())
        cache[mfix_exe] = {'stat': stat, 'flags': flags}
        return True




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
