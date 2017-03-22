""" Dialog in the GUI for starting an MFIX solver job """

import logging
import os
import signal
import sys
import tempfile

try: #2.7
    from StringIO import StringIO
except ImportError: # 3
    from io import StringIO

import errno
import multiprocessing
import json

from collections import OrderedDict
from subprocess import Popen, PIPE
from glob import glob

try: #2.7
    import ConfigParser as configparser
except: # 3
    import configparser

from qtpy import PYQT5, uic
from qtpy.QtCore import (
    QProcess,
    QProcessEnvironment,
    Signal,
)
from qtpy.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QLabel,
    QSpinBox,
)

from mfixgui.tools.general import (
    clear_layout,
    extract_config,
    replace_with_dict
    )
from mfixgui.tools.util import (
    get_mfix_home,
)

from mfixgui.widgets.base import BASE_WIDGETS

log = logging.getLogger('mfix-gui' if __name__ == '__main__' else __name__)

RECENT_EXE_LIMIT = 5
MFIXSOLVER_GLOB_NAMES = ['mfixsolver', 'mfixsolver.exe']

class RunPopup(QDialog):

    signal_run = Signal()
    signal_submit = Signal()
    signal_cancel = Signal()

    def __init__(self, solver, parent):
        super(RunPopup, self).__init__(parent)
        self.commandline_option_exe = solver if solver else None
        self.mfix_available = False
        self.mfix_exe_cache = {}
        self.solver_list = []
        self.template_values = {}
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
        ui.combobox_solver_local.currentIndexChanged.connect(self.handle_exe_change)
        ui.combobox_solver_queue.currentIndexChanged.connect(self.handle_exe_change)

        if bool(self.parent.monitor.get_res_files()):
            self.title = 'Resume'
        else:
            self.title = 'Run'
        self.signal_run.connect(self.handle_run)
        self.signal_submit.connect(self.handle_submit)
        self.signal_cancel.connect(self.close)

        ui.button_local_run.clicked.connect(self.handle_run)
        ui.button_queue_submit.clicked.connect(self.handle_submit)
        ui.button_local_cancel.clicked.connect(self.handle_abort)
        ui.button_queue_cancel.clicked.connect(self.handle_abort)
        ui.pushbutton_browse_template.clicked.connect(self.handle_browse_template)

        ui.label_cores_detected.setText("%d cores available locally." % multiprocessing.cpu_count())

        self.initialize_ui()
        self.init_templates()

    @property
    def solver(self):
        """The currently selected solver, depends on current tab (queue/local)"""
        idx = self.ui.tabWidget.currentIndex()
        # local
        if idx == 1:
            solver = self.ui.combobox_solver_local.currentText()
        # queue
        else:
            solver = self.ui.combobox_solver_queue.currentText()

        if not solver:
            solver = None

        return solver

    # UI update functions
    def initialize_ui(self):

        ui = self.ui
        self.setWindowTitle(self.title)
        # create initial executable list
        self.solver_list = self.get_solver_list()

        # set OMP_NUM_THREADS
        project_threads = self.gui_comments.get('OMP_NUM_THREADS', None)
        env_threads = os.environ.get('OMP_NUM_THREADS', None)
        if project_threads:
            ui.spinbox_threads.setValue(int(project_threads))
        elif env_threads:
            ui.spinbox_threads.setValue(int(env_threads))
        else:
            ui.spinbox_threads.setValue(1)

        # local/queue
        self.ui.tabWidget.setCurrentIndex(int(self.gui_comments.get('run_location', 1)))

        if self.solver_list:
            self.mfix_available = True
            self.populate_combobox_solver()
        else:
            self.mfix_available = False
            self.parent.message(
                icon='warning',
                text='MFiX not found. Please browse for an executable.',
                buttons=['ok', 'cancel'],
                default='ok')

        self.update_dialog_options()

    def init_templates(self):

        # look for templates in MFIX_HOME/queue_templates
        search_p = os.path.join(get_mfix_home(), 'queue_templates')
        self.templates = {}
        for root, dirs, files in os.walk(search_p):
            for f in files:
                p = os.path.join(root, f)
                self.add_queue_template(p)

        # look for recent templates
        temp_paths = self.settings.value('queue_templates')
        if temp_paths:
            for temp_path in temp_paths.split('|'):
                if os.path.exists(temp_path):
                    self.add_queue_template(temp_path)

        self.ui.combobox_template.currentIndexChanged.connect(self.update_queue_widgets)

        temp = self.gui_comments.get('queue_template')
        if temp:
            self.template_values = json.loads(temp)
            t_name = self.template_values.get('template')
            if t_name:
                self.set_current_template(t_name)

        self.update_queue_widgets()

    def save_template(self):
        '''Save the current template data'''
        self.collect_template_values()
        template_txt = self.ui.combobox_template.currentText()
        self.template_values['template'] = template_txt
        self.gui_comments['queue_template'] = json.dumps(self.template_values)

    def collect_template_values(self):
        template_txt = self.ui.combobox_template.currentText()
        template = self.templates.get(template_txt, {})
        replace_dict = {}
        for name, wid in template.items():
            if not isinstance(wid, dict):
                continue

            if 'widget_obj' in wid:
                wid_obj = wid['widget_obj']
                if isinstance(wid_obj, (QSpinBox, QDoubleSpinBox)):
                    self.template_values[name] = v = wid_obj.value()
                elif isinstance(wid_obj, QCheckBox):
                    self.template_values[name] = v = wid_obj.value
                    if v:
                        v = wid.get('true', '')
                    else:
                        v = wid.get('false', '')
                else:
                    self.template_values[name] = v = wid_obj.value
                replace_dict[name] = v
        return replace_dict

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
        c.readfp(StringIO(config))

        d = OrderedDict([(s, dict(c.items(s))) for s in c.sections()])
        d['path'] = path
        d['script'] = script

        name = os.path.basename(path)
        if 'options' in d:
            name = d['options'].get('name', name)

        self.templates[name] = d

        self.ui.combobox_template.clear()
        self.ui.combobox_template.addItems(list(self.templates.keys()))

        if select:
            self.set_current_template(name)

    def update_queue_widgets(self):

        l = self.ui.groupbox_queue_options_gridlayout
        clear_layout(l)
        tp = self.ui.combobox_template.currentText()

        wids_data = self.templates.get(tp, None)
        if wids_data is None:
            return

        # add the widgets
        for i, wid in enumerate(list(wids_data.keys())):
            wd = wids_data[wid]
            if not isinstance(wd, dict) or wid == 'options':
                continue

            label = QLabel(wd.get('label', wid))
            l.addWidget(label, i, 0)
            widget = BASE_WIDGETS.get(wd.get('widget', 'lineedit'), BASE_WIDGETS['lineedit'])()
            items = [it.strip() for it in wd.get('items', '').split('|')]
            v = self.template_values.get(wid)
            if not v or self.template_values.get('template') != tp:
                v = wd.get('value')
            if isinstance(widget, QComboBox) and items:
                widget.addItems(items)
                if v not in items:
                    v = items[0]
            widget.updateValue('', v)
            widget.help_text = wd.get('help', 'No help avaliable.')
            l.addWidget(widget, i, 1)
            wd['widget_obj'] = widget

    def populate_combobox_solver(self):
        """ Add items from self.solver_list to combobox,
        select the first item """
        ui = self.ui
        for combo in [ui.combobox_solver_local, ui.combobox_solver_queue]:
            combo.clear()
            combo.addItems(self.solver_list)

    def update_dialog_options(self):
        """ Enable or disable options based on self.solver features,
        local or remote settings """
        ui = self.ui

        # Enable/disable widgets
        enable = self.mfix_available and self.solver is not None
        ui.combobox_solver_local.setEnabled(enable)
        ui.combobox_solver_queue.setEnabled(enable)
        ui.button_local_run.setEnabled(enable)
        ui.button_queue_submit.setEnabled(enable)
        ui.label_mfix_exe_warning.setVisible(not enable)
        ui.groupbox_run_options.setEnabled(enable)

        dmp = self.dmp_enabled()
        smp = self.smp_enabled()
        ui.spinbox_nodesi.setEnabled(dmp)
        ui.spinbox_nodesj.setEnabled(dmp)
        ui.spinbox_nodesk.setEnabled(dmp and not self.parent.project.get_value('no_k'))
        ui.spinbox_threads.setEnabled(smp)

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
        log.info('SMP enabled with OMP_NUM_THREADS=%s', os.environ["OMP_NUM_THREADS"])
        self.gui_comments['OMP_NUM_THREADS'] = thread_count
        self.gui_comments['run_location'] = self.ui.tabWidget.currentIndex()

        self.save_selected_exe()
        self.save_template()

        if self.title == 'Run':
            self.parent.update_keyword('run_type', 'new')
            output_files = self.parent.monitor.get_outputs()
            if output_files:
                if not self.parent.remove_output_files(output_files, force_remove=True):
                    log.info('output files exist and run was canceled')
                    return False
        elif self.ui.use_spx_checkbox.isChecked():
            self.parent.update_keyword('run_type', 'restart_1')
        else:
            # TODO: is it correct to remove all but *.RES ?
            spx_files = self.parent.monitor.get_outputs(['*.SP?', "*.pvd", "*.vtp"])
            if not self.parent.remove_output_files(spx_files, force_remove=True):
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

        msg = 'Submitting to queue'
        self.parent.print_internal(msg, color='blue')

        (script,
         sub_cmd,
         delete_cmd,
         status_cmd,
         job_id_regex,
         replace_dict) = self.get_submit_command()

        self.submit_command(script, sub_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict)

    def handle_exe_change(self):
        """emit signals when exe combobox changes"""
        log.debug('selected new solver %s', self.solver)
        self.update_dialog_options()

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

        self.save_selected_exe(new_exe) # non-exe get saved
        self.mfix_available = True
        self.solver_list = self.get_solver_list()
        self.populate_combobox_solver()
        log.debug('selected new exe %s', new_exe)

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

        # add it to the recent settings
        temp_paths = self.settings.value('queue_templates')
        good_paths = [os.path.abspath(new_temp)]
        if temp_paths:
            for temp_path in temp_paths.split('|'):
                if os.path.exists(temp_path):
                    good_paths.append(temp_path)
        self.settings.setValue('queue_templates',
                               '|'.join(list(set(good_paths))[:RECENT_EXE_LIMIT]))

    # utils
    def save_selected_exe(self, new_solver=None):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        if new_solver is None:
            new_solver = self.solver
        if new_solver is None:
            self.parent.error('No solver selected')
        self.settings.setValue('mfix_exe', new_solver)
        self.gui_comments['mfix_exe'] = new_solver
        recent_list = self.settings.value('recent_executables', None)
        if recent_list is not None:
            recent_list = recent_list.split(os.pathsep)[:RECENT_EXE_LIMIT]
        else:
            recent_list = []
        if new_solver in recent_list:
            recent_list.pop(recent_list.index(new_solver))
        recent_list.insert(0, new_solver)
        self.settings.setValue(
            'recent_executables',
            str(os.pathsep).join(recent_list))

    def get_solver_list(self):
        """ assemble list of executables from:
        - command line
        - project file 'mfix_exe'
        - project dir
        - config item 'recent_executables'
        - default install location
        """

        def recently_used_executables():
            recent_list = self.settings.value('recent_executables')
            if recent_list:
                # limit recently used exes to RECENT_EXE_LIMIT
                recent_lim = recent_list.split(os.pathsep)[:RECENT_EXE_LIMIT]
                recent_list = [exe for exe in recent_lim if os.path.exists(exe)]
                for recent_exe in recent_list:
                    yield recent_exe

        def project_directory_executables():
            for name in MFIXSOLVER_GLOB_NAMES:
                for exe in glob(os.path.join(self.project_dir, name)):
                    yield os.path.abspath(exe)

        def project_file_executable():
            project_exe = self.project.get_value('mfix_exe')
            if project_exe:
                yield project_exe

        def python_path():
            for d in sys.path:
                # filter out empty strings and current directory from $PATH
                if d and d != os.path.curdir and os.path.isdir(d):
                    for name in MFIXSOLVER_GLOB_NAMES:
                        for exe in glob(os.path.join(d, name)):
                            yield exe

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
                    for name in MFIXSOLVER_GLOB_NAMES:
                        for exe in glob(os.path.join(d, name)):
                            yield exe

        def mfix_build_directories():
            mfix_home = get_mfix_home()
            bin_dir = os.path.join(mfix_home, 'bin')
            dir_list = set([mfix_home])
            if os.path.isdir(bin_dir):
                dir_list.add(bin_dir)
            for d in dir_list:
                for name in MFIXSOLVER_GLOB_NAMES:
                    for exe in glob(os.path.join(d, name)):
                        yield exe

        def get_saved_exe():
            last_exe = self.settings.value('mfix_exe')
            if last_exe and os.path.exists(last_exe):
                yield last_exe

        def command_line_option():
            if self.commandline_option_exe and os.path.exists(self.commandline_option_exe):
                yield self.commandline_option_exe

        exe_list_order = [
            recently_used_executables,
            project_directory_executables,
            python_path,
            os_path,
            mfix_build_directories,
            get_saved_exe,
            project_file_executable,
            command_line_option]

        od = OrderedDict()
        # look for executables in the order listed in exe_list_order
        for exe_spec in exe_list_order:
            for exe in exe_spec():
                exe_exists = os.path.isfile(exe) and os.access(exe, os.X_OK)
                if exe_exists and self.get_exe_flags(exe):
                    od[exe] = True
                else:
                    self.parent.error('{} is not an executable'.format(exe))

        return list(od.keys())

    def get_exe_flags(self, solver):
        """ get and cache (and update) executable features """

        # let non-exe solvers through
        if solver and os.path.splitext(solver)[1] in ['.so', '.pyd']:
            return {'flags': 'python'}

        if solver is None:
            return None
        try:
            stat = os.stat(solver)
        except OSError as e:
            log.debug(str(e))
            return None

        # stat will have changed if the exe has been modified since last check
        if (stat, solver) in self.mfix_exe_cache:
            _, flags = self.mfix_exe_cache[(stat, solver)]
            return flags
        try:
            log.debug('Feature testing MFiX %s', solver)
            exe_dir = os.path.dirname(solver)
            popen = Popen(solver + " --print-flags",
                          cwd=exe_dir, stdout=PIPE, stderr=PIPE, shell=True)
            (out, err) = popen.communicate()
            if err:
                log.error('MFiX %s', str(err))
        except:
            log.error("could not run %s --print-flags", solver)
            return None

        flags = str(out.strip())
        mfix_exe_flags = {'flags': flags}
        self.mfix_exe_cache[(stat, solver)] = stat, mfix_exe_flags
        return mfix_exe_flags

    def dmp_enabled(self):
        config = self.get_exe_flags(self.solver)
        flags = config['flags'] if config else ''
        return 'dmp' in flags

    def smp_enabled(self):
        config = self.get_exe_flags(self.solver)
        flags = config['flags'] if config else ''
        return 'smp' in flags

    def get_run_command(self):

        nodesi = self.ui.spinbox_nodesi.value()
        nodesj = self.ui.spinbox_nodesj.value()
        nodesk = self.ui.spinbox_nodesk.value()

        if self.dmp_enabled():
            dmp = ['mpirun', '-np', str(nodesi * nodesj * nodesk)]
        else:
            dmp = []

        if self.smp_enabled():
            smp = ['env', 'OMP_NUM_THREADS=%s' % str(self.ui.spinbox_threads.value())]
        else:
            smp = []

        # FIXME: code-cleanup; mfix_exe really should be renamed to mfixsolver everywhere
        run_cmd = smp + dmp + [self.solver]

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

    def get_submit_command(self):

        cmd = self.get_run_command()

        template_txt = self.ui.combobox_template.currentText()
        template = self.templates[template_txt]

        # collect widget values
        replace_dict = self.collect_template_values()
        replace_dict.update({
            'PROJECT_NAME': self.parent.project.get_value('run_name', default=''),
            'COMMAND': ' '.join(cmd),
            'MFIX_HOME': get_mfix_home(),
        })

        # replace twice to make sure that any references added the first time
        # get replaced
        script = replace_with_dict(template['script'], replace_dict)
        script = replace_with_dict(script, replace_dict)

        sub_cmd = template['options'].get('submit', False)
        delete_cmd = template['options'].get('delete', False) # XXX
        status_cmd = template['options'].get('status', False)
        job_id_regex = template['options'].get('job_id_regex', None)

        ## FIXME, return something nicer than this 6-tuple
        return script, sub_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict

    def submit_command(self, script, sub_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict):

        self.remove_mfix_stop()

        if not sub_cmd:
            self.parent.error(('The template file at: {}\n'
                               'does not have a submit_cmd defined').format(tempfile['path']))
            return

        self.parent.job_manager.submit_command(script,
                                               sub_cmd,
                                               delete_cmd,
                                               status_cmd,
                                               job_id_regex,
                                               replace_dict)

    def remove_mfix_stop(self):
        mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        if os.path.exists(mfix_stop_file):
            try:
                os.remove(mfix_stop_file)
            except OSError:
                log.error("Cannot remove %s", mfix_stop_file)
                return

    def start_command(self, cmd, cwd, env):
        """Start MFIX in QProcess"""

        self.remove_mfix_stop()

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
            msg = "MFiX process %d is running" % self.mfixproc.processId()
            self.parent.signal_update_runbuttons.emit(msg)

        def slot_read_out():
            out_str = bytes(self.mfixproc.readAllStandardOutput()).decode('utf-8')
            self.parent.stdout_signal.emit(out_str)

        def slot_read_err():
            err_str = bytes(self.mfixproc.readAllStandardError()).decode('utf-8')
            self.parent.stderr_signal.emit(err_str)

        def slot_finish(status):

            if self.parent.job_manager.pidfile:
                try:
                    os.unlink(self.parent.job_manager.pidfile)
                    self.parent.job_manager.pidfile = None
                except OSError as e:
                    if e.errno != errno.ENOENT:
                        raise
                finally:
                    # Turn off timers!
                    if self.parent.job_manager.job:
                        self.parent.job_manager.job.cleanup_and_exit()
                    self.parent.job_manager.job = None
                    msg = "MFiX process has stopped"
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
