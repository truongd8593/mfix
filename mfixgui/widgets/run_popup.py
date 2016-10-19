#!/usr/bin/env python

import logging
import os
import signal
import sys
import tempfile
import time
import StringIO

from collections import OrderedDict
from subprocess import Popen, PIPE
from glob import glob

from qtpy import PYQT5
from qtpy.QtCore import Signal, QProcess, QProcessEnvironment, QTimer
from qtpy.QtWidgets import (QDialog, QApplication, QFileDialog,
                            QDialogButtonBox, QLabel, QComboBox, QSpinBox,
                            QDoubleSpinBox, QCheckBox)

from mfixgui.tools.general import get_mfix_home, clear_layout, extract_config, replace_with_dict
from mfixgui.widgets.base import BASE_WIDGETS

try:
    from PyQt5 import uic
except ImportError:
    from PyQt4 import uic

try: #2.7
    import ConfigParser as configparser
except: # 3
    import configparser


log = logging.getLogger('mfix-gui' if __name__=='__main__' else __name__)

RECENT_EXE_LIMIT = 5
MFIX_EXE_NAMES = ['mfix', 'mfix.exe']


class RunPopup(QDialog):

    signal_run = Signal()
    signal_submit = Signal()
    signal_cancel = Signal()
    set_run_mfix_exe = Signal()

    def __init__(self, mfix_exe, parent):

        super(RunPopup, self).__init__(parent)

        self.commandline_option_exe = mfix_exe if mfix_exe else None
        self.mfix_available = False
        self.mfix_exe = None
        self.mfix_exe_cache = {}
        self.mfix_exe_list = []
        self.cmdline = None
        self.parent = parent
        self.project = parent.project
        self.settings = parent.settings
        self.project_dir = parent.get_project_dir()
        self.gui_comments = self.project.mfix_gui_comments

        # load ui
        thisdir = os.path.abspath(os.path.dirname(__file__))
        uidir = os.path.join(os.path.dirname(thisdir), 'uifiles')
        self.ui = ui = uic.loadUi(os.path.join(uidir, 'run_popup.ui'), self)

        ui.button_browse_exe.clicked.connect(self.handle_browse_exe)
        ui.button_browse_exe_2.clicked.connect(self.handle_browse_exe)
        ui.combobox_mfix_exe.currentIndexChanged.connect(self.handle_exe_change)

        if bool(self.parent.monitor.get_res_files()):
            self.title = 'Resume'
        else:
            self.title = 'Run'

        self.signal_run.connect(self.handle_run)
        self.signal_submit.connect(self.handle_submit)
        self.signal_cancel.connect(self.close)

        self.ui.button_local_run.clicked.connect(self.handle_run)
        self.ui.button_queue_submit.clicked.connect(self.handle_submit)
        self.ui.button_local_cancel.clicked.connect(self.handle_abort)
        self.ui.button_queue_cancel.clicked.connect(self.handle_abort)
        self.ui.pushbutton_browse_template.clicked.connect(self.handle_browse_template)

        self.initialize_ui()
        self.init_templates()

    # UI update functions

    def initialize_ui(self):

        ui = self.ui
        self.setWindowTitle(self.title)
        # create initial executable list
        self.generate_exe_list()

        # set OMP_NUM_THREADS
        project_threads = self.gui_comments.get('OMP_NUM_THREADS', None)
        env_threads = os.environ.get('OMP_NUM_THREADS', None)
        if project_threads:
            ui.spinbox_threads.setValue(int(project_threads))
        elif env_threads:
            ui.spinbox_threads.setValue(int(env_threads))
        else:
            ui.spinbox_threads.setValue(1)

        # issues/149
        #ui.spinbox_nodesi.setValue(self.project.get_value('nodesi', default=1))
        #ui.spinbox_nodesj.setValue(self.project.get_value('nodesj', default=1))
        #ui.spinbox_nodesk.setValue(self.project.get_value('nodesk', default=1))

        if self.mfix_exe_list:
            self.mfix_available = True
            self.mfix_exe = self.mfix_exe_list[0]
            self.populate_combobox_mfix_exe()
        else:
            self.mfix_available = False
            self.mfix_exe = None
            self.parent.message(
                icon='warning',
                text='MFIX not found. Please browse for an executable.',
                buttons=['ok','cancel'],
                default='ok')

        self.update_dialog_options()

    def init_templates(self):

        search_p = os.path.join(get_mfix_home(), 'mfixgui', 'queue_templates')

        self.templates = {}
        for root, dirs, files in os.walk(search_p):
            for f in files:
                p = os.path.join(root, f)
                self.add_queue_template(p)

        self.ui.combobox_template.addItems(list(self.templates.keys()))
        self.ui.combobox_template.currentIndexChanged.connect(self.update_queue_widgets)
        self.update_queue_widgets()

    def set_current_template(self, name):
        '''set the template file combobox'''
        cb = self.ui.combobox_template
        for itm in range(cb.count()):
            if str(name).lower() == str(cb.itemText(itm)).lower():
                cb.setCurrentIndex(itm)
                break

    def add_queue_template(self, path, select=False):
        config, script = extract_config(path)
        c = configparser.ConfigParser()
        c.readfp(StringIO.StringIO(config))

        d = OrderedDict([(s, dict(c.items(s))) for s in c.sections()])
        d['path'] = path
        d['script'] = script

        name = os.path.basename(path)
        if 'options' in d:
            name = d['options'].get('name', name)
        self.templates[name] = d

        if select:
            self.set_current_template(name)


    def update_queue_widgets(self):

        l = self.ui.groupbox_queue_options_gridlayout
        clear_layout(l)
        tp =  self.ui.combobox_template.currentText()

        wids_data = self.templates[tp]

        # add the widgets
        for i, wid in enumerate(list(wids_data.keys())):
            wd = wids_data[wid]
            if not isinstance(wd, dict) or wid == 'options': continue

            label = QLabel(wd.get('label', wid))
            l.addWidget(label, i, 0)
            widget = BASE_WIDGETS.get(wd.get('widget', 'lineedit'), BASE_WIDGETS['lineedit'])()
            if isinstance(widget, QComboBox):
                widget.addItems(wd.get('items','').split(';'))
            widget.updateValue('', wd.get('value'))
            widget.help_text = wd.get('help', 'No help avaliable.')
            l.addWidget(widget, i, 1)
            wd['widget_obj'] = widget

    def populate_combobox_mfix_exe(self):
        """ Add items from self.mfix_exe_list to combobox,
        select the first item """
        self.ui.combobox_mfix_exe.clear()
        self.ui.combobox_mfix_exe_2.clear()
        self.ui.combobox_mfix_exe.addItems(self.mfix_exe_list)
        self.ui.combobox_mfix_exe_2.addItems(self.mfix_exe_list)

    def update_dialog_options(self):
        """ Enable or disable options based on self.mfix_exe features,
        local or remote settings """

        self.ui.combobox_mfix_exe.setEnabled(self.mfix_available)
        self.ui.button_local_run.setEnabled(self.mfix_available)
        self.ui.button_queue_submit.setEnabled(self.mfix_available)

        self.update_no_mfix_warning()

        self.ui.groupbox_run_options.setEnabled(self.mfix_available)
        cfg = self.get_exe_flags(self.mfix_exe)
        dmp = 'dmp' in cfg['flags'] if cfg else False # why not use dmp_enabled
        smp = 'smp' in cfg['flags'] if cfg else False
        dmp = self.dmp_enabled()
        smp = self.smp_enabled()
        self.ui.spinbox_nodesi.setEnabled(dmp)
        self.ui.spinbox_nodesj.setEnabled(dmp)
        self.ui.spinbox_nodesk.setEnabled(dmp and not self.parent.project.get_value('no_k'))
        self.ui.spinbox_threads.setEnabled(smp)

    def update_no_mfix_warning(self):
        ok = bool(self.mfix_exe)
        self.ui.label_mfix_exe_warning.setVisible(not ok)
        self.ui.combobox_mfix_exe.setEnabled(ok)
        self.ui.button_local_run.setEnabled(ok)
        self.ui.button_queue_submit.setEnabled(ok)
        if not ok:
            self.parent.print_internal("Warning: no MFIX executables available")

    def popup(self):
        self.show()
        self.raise_()
        self.activateWindow()


    # event handlers

    def handle_abort(self):
        self.signal_cancel.emit()

    def finish_with_dialog(self):
        """ save run options in project file, then emit run signal """
        thread_count = str(self.ui.spinbox_threads.value())
        os.environ['OMP_NUM_THREADS'] = thread_count
        log.info('SMP enabled with OMP_NUM_THREADS=%s' % \
                 os.environ["OMP_NUM_THREADS"])
        self.gui_comments['OMP_NUM_THREADS'] = thread_count

        #self.project.updateKeyword('nodesi',
        #                            self.ui.spinbox_keyword_nodesi.value())
        #self.project.updateKeyword('nodesj',
        #                            self.ui.spinbox_keyword_nodesj.value())
        #self.project.updateKeyword('nodesk',
        #                            self.ui.spinbox_keyword_nodesk.value())
        self.save_selected_exe(self.mfix_exe)
        self.set_run_mfix_exe.emit() # possible duplication, but needed
                                     # in case signal has not yet been fired

        if self.title == 'Run':
            self.parent.update_keyword('run_type', 'new')
            output_files = self.parent.monitor.get_outputs()
            if output_files:
                if not self.parent.remove_output_files(output_files):
                    log.info('output files exist and run was canceled')
                    return False
        elif self.ui.use_spx_checkbox.isChecked():
            self.parent.update_keyword('run_type', 'restart_1')
        else:
            # TODO: is it correct to remove all but *.RES ?
            spx_files = self.parent.monitor.get_outputs(['*.SP?', "*.pvd", "*.vtp"])
            if not self.parent.remove_output_files(spx_files):
                log.debug('SP* files exist and run was canceled')
                return False
            self.parent.update_keyword('run_type', 'restart_2')

        self.parent.save_project()
        self.parent.update_source_view()
        self.close()
        self.parent.signal_update_runbuttons.emit('')
        return True

    def handle_run(self):
        if not self.finish_with_dialog():
            return

        run_cmd = self.get_run_command()
        msg = 'Starting %s' % ' '.join(run_cmd)
        self.parent.print_internal(msg, color='blue')

        self.start_command(
            cmd=run_cmd,
            cwd=self.parent.get_project_dir(),
            env=os.environ)

    def handle_submit(self):
        if not self.finish_with_dialog():
            return

        run_cmd = self.get_run_command()
        msg = 'Submitting %s' % ' '.join(run_cmd)
        self.parent.print_internal(msg, color='blue')

        self.submit_command(cmd=run_cmd)

    def handle_resume(self):
        """resume previously stopped mfix run"""
        #FIXME need to catch/report errors, writeDatFile is too low-level
        self.project.writeDatFile(self.get_project_file()) # XXX
        self._start_mfix(False)

    def handle_exe_change(self):
        """ emit signals when exe combobox changes """
        self.mfix_exe = new_exe = self.ui.combobox_mfix_exe.currentText()
        log.debug('selected new exe %s' % new_exe)
        self.update_dialog_options()
        self.set_run_mfix_exe.emit()

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe = QFileDialog.getOpenFileName(
            self, "Select Executable",
            directory=self.project_dir,
            options=QFileDialog.DontResolveSymlinks)
        if not new_exe:
            return
        if PYQT5:
            new_exe = new_exe[0]
        if not self.prepend_to_exe_list(new_exe):
            self.parent.message(
                icon='warning',
                title='Warning',
                text='The selected file is not an executable MFIX binary')
            return
        self.mfix_available = True
        self.mfix_exe = new_exe
        self.prepend_to_exe_list(new_exe)
        self.populate_combobox_mfix_exe()
        log.debug('selected new exe %s' % new_exe)
        self.set_run_mfix_exe.emit()

    def handle_browse_template(self):
        """ Handle file open dialog for user specified exe """
        new_temp = QFileDialog.getOpenFileName(
            self, "Select a Template",
            directory=self.project_dir,)
        if not new_temp:
            return
        if PYQT5:
            new_temp = new_temp[0]

        self.add_queue_template(new_temp, select=True)

    # utils

    def save_selected_exe(self, new_exe):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        self.settings.setValue('mfix_exe', new_exe)
        self.gui_comments['mfix_exe'] = new_exe
        recent_list = self.settings.value('recent_executables')
        if recent_list:
            recent_list = recent_list.split(os.pathsep)[:RECENT_EXE_LIMIT]
        else:
            recent_list = []
        if new_exe in recent_list:
            recent_list.pop(recent_list.index(new_exe))
        recent_list.insert(0, new_exe)
        self.settings.setValue(
                        'recent_executables',
                        str(os.pathsep).join(recent_list))

    def prepend_to_exe_list(self, exe):
        """ Verify exe exists, is executable, and appears only once in list."""
        if not (self.exe_exists(exe) and self.get_exe_flags(exe)):
            return False
        if exe in self.mfix_exe_list:
            self.mfix_exe_list.pop(self.mfix_exe_list.index(exe))
        self.mfix_exe_list.insert(0, exe)
        return True

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
            #    for name in ['mfix', 'mfix.exe']:
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
            bin_dir = os.path.join(mfix_home, 'bin')
            dir_list = set([mfix_home])
            if os.path.isdir(bin_dir):
                dir_list.add(bin_dir)
            for d in dir_list:
                for name in MFIX_EXE_NAMES:
                    for exe in glob(os.path.join(d, name)):
                        yield exe

        def get_saved_exe():
            last_exe = self.settings.value('mfix_exe')
            if last_exe:
                yield last_exe

        def command_line_option():
            if self.commandline_option_exe:
                yield self.commandline_option_exe

        exe_list_order = [
            recently_used_executables,
            project_directory_executables,
            os_path,
            mfix_build_directories,
            get_saved_exe,
            project_file_executable,
            command_line_option]

        # look for executables in the order listed in exe_list_order
        for exe_spec in exe_list_order:
            for exe in exe_spec():
                self.prepend_to_exe_list(exe)

    def exe_exists(self, exe):
        """ Verify exe exists and is executable """
        return (os.path.isfile(exe) and os.access(exe, os.X_OK))

    def get_exe_flags(self, mfix_exe):
        """ get and cache (and update) executable features """
        if mfix_exe is None:
            return None
        try:
            stat = os.stat(mfix_exe)
        except OSError as e:
            log.debug(str(e))
            return None

        # stat will have changed if the exe has been modified since last check
        if (stat, mfix_exe) in self.mfix_exe_cache:
            _, flags = self.mfix_exe_cache[(stat, mfix_exe)]
            return flags
        try:
            log.debug('Feature testing MFIX %s' % mfix_exe)
            exe_dir = os.path.dirname(mfix_exe)
            popen = Popen(mfix_exe + " --print-flags",
                        cwd=exe_dir, stdout=PIPE, stderr=PIPE, shell=True)
            (out, err) = popen.communicate()
            if err:
                log.error('MFIX %s' % str(err))
        except:
            log.error("could not run %s --print-flags", mfix_exe)
            return None

        flags = str(out.strip())
        mfix_exe_flags = {'flags': flags}
        self.mfix_exe_cache[(stat, mfix_exe)] = stat, mfix_exe_flags
        return mfix_exe_flags

    def dmp_enabled(self):
        config = self.get_exe_flags(self.mfix_exe)
        flags = config['flags'] if config else ''
        return 'dmp' in flags

    def smp_enabled(self):
        config = self.get_exe_flags(self.mfix_exe)
        flags = config['flags'] if config else ''
        return 'smp' in flags

    def get_run_command(self):

        nodesi = self.project.get_value('nodesi', 1)
        nodesj = self.project.get_value('nodesj', 1)
        nodesk = self.project.get_value('nodesk', 1)

        if self.dmp_enabled():
            dmp = ['mpirun', '-np', str(nodesi * nodesj * nodesk)]
        else:
            dmp = []

        if self.smp_enabled():
            smp = ['env', 'OMP_NUM_THREADS=%s' % str(self.ui.spinbox_threads.value())]
        else:
            smp = []

        run_cmd = smp + dmp + [self.mfix_exe]

        project_filename = os.path.basename(self.parent.get_project_file())
        # Warning, not all versions of mfix support '-f' !
        run_cmd += ['-f', project_filename]

        # Add key=value flags at end
        if self.dmp_enabled:
            run_cmd += ['nodesi=%s'%nodesi,
                        'nodesj=%s'%nodesj]
            if not self.parent.project.get_value('no_k'):
                run_cmd += ['nodesk=%s'%nodesk]
        return run_cmd

    def submit_command(self, cmd):

        tp =  self.ui.combobox_template.currentText()
        template = self.templates[tp]

        # collect widget values
        replace_dict = {
            'PROJECT_NAME': self.parent.project.get_value('run_name', default=''),
            'COMMAND': ' '.join(cmd),
        }
        for name, wid in template.items():
            if not isinstance(wid, dict): continue

            if 'widget_obj' in wid:
                wid_obj = wid['widget_obj']
                if isinstance(wid_obj, (QSpinBox, QDoubleSpinBox)):
                    v = wid_obj.value()
                elif isinstance(wid_obj, QCheckBox):
                    v = wid_obj.value
                    if v:
                        v = wid.get('true', '')
                    else:
                        v = wid.get('false', '')
                else:
                    v = wid_obj.value
                replace_dict[name] = v

        # replace twice ti make sure that any references added the first time
        # get replaced
        script = replace_with_dict(template['script'], replace_dict)
        script = replace_with_dict(script, replace_dict)

        print(script)

        sub_cmd = template['options'].get('submit', False)
        if not sub_cmd:
            self.parent.error('The template file at: {}\ndoes not have a submit_cmd defined'.format(tempfile['path']))
            return
        delete_cmd = template['options'].get('delete', False)
        status_cmd = template['options'].get('status', False)
        job_id_regex = template['options'].get('job_id_regex', None)


        self.parent.job_manager.submit_command(script, sub_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict)

    def start_command(self, cmd, cwd, env):
        """Start MFIX in QProcess"""

        mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        if os.path.exists(mfix_stop_file):
            try:
                os.remove(mfix_stop_file)
            except OSError:
                log.error("Cannot remove", mfix_stop_file)
                return

        self.mfixproc = QProcess()
        if not self.mfixproc:
            log.warn("QProcess creation failed")
            return
        self.mfixproc.setWorkingDirectory(cwd)
        process_env = QProcessEnvironment()
        for key, val in env.items():
            process_env.insert(key, val)
        self.mfixproc.setProcessEnvironment(process_env)

        def slot_start():
            # Keep a copy because it gets reset
            msg = "MFIX process %d is running" % self.mfixproc.pid()
            self.parent.signal_update_runbuttons.emit(msg)
            log.debug("Full MFIX startup parameters: %s", self.cmdline)

        def slot_read_out():
            out_str = bytes(self.mfixproc.readAllStandardOutput()).decode('utf-8')
            self.parent.stdout_signal.emit(out_str)

        def slot_read_err():
            err_str = bytes(self.mfixproc.readAllStandardError()).decode('utf-8')
            self.parent.stderr_signal.emit(err_str)

        def slot_finish(status):
            msg = "MFIX process has stopped"
            self.parent.signal_update_runbuttons.emit(msg)

        def slot_error(error):
            if error == QProcess.FailedToStart:
                msg = "Process failed to start "+self.cmdline
            elif error == QProcess.Crashed:
                msg = "Process exit "+self.cmdline
            elif error == QProcess.Timedout:
                msg = "Process timeout "+self.cmdline
            elif error in (QProcess.WriteError, QProcess.ReadError):
                msg = "Process communication error "+self.cmdline
            else:
                msg = "Unknown error "+self.cmdline
            log.warn(msg)
            # make the message print in red
            self.parent.stderr_signal.emit(msg)

        self.cmdline = ' '.join(cmd) # cmd is a list
        self.mfixproc.started.connect(slot_start)
        self.mfixproc.readyReadStandardOutput.connect(slot_read_out)
        self.mfixproc.readyReadStandardError.connect(slot_read_err)
        self.mfixproc.finished.connect(slot_finish)
        self.mfixproc.error.connect(slot_error)
        self.mfixproc.start(cmd[0], cmd[1:])

if __name__ == '__main__':

    args = sys.argv
    app = QApplication(args)
    run_popup = Run_Popup(QDialog())
    run_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    app.exec_()
    app.deleteLater()

    sys.exit()
