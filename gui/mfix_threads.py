"""classes to manage external MFIX process and monitor files"""

import os
import glob
import json
import logging
import time
import subprocess

try:
    # For Python 3.0 and later
    from urllib.request import urlopen, Request
    from urllib.parse import urlencode
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib import urlencode
    from urllib2 import urlopen, Request

from tools.general import get_mfix_home

from qtpy.QtCore import QProcess, QTimer

class MfixJobManager(object):
    """class for monitoring MFIX jobs"""

    def __init__(self, parent):
        self.parent = parent
        self.cmd = None
        self.cwd = None
        self.env = None
        self.is_pymfix = False
        self.mfixproc = None
        self.mfix_pid = None # Keep track of pid separately
        self.status = None
        self.timer = None

    def is_running(self):
        """indicate whether an MFIX job is running"""
        ret = (self.mfixproc is not None
               and self.mfixproc.state() == QProcess.Running)
        return ret

    def is_paused(self):
        """indicate whether pymfix is paused"""
        if not self.is_pymfix:
            return False
        return getattr(self.status, 'paused', False)

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        if self.is_pymfix:
            self.terminate_pymfix()
        else:
            #  GUI is nonresponsive while this runs.  Maybe it needs to be in another thread after all
            mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
            try:
                open(mfix_stop_file, "ab").close()
            except OSError:
                pass

        def force_kill():
            """attempt to kill job; return True if further force_kill() calls may be necessary"""
            if not self.is_running():
                return
            if not 'ok' in self.parent.message(text="MFIX is not responding. Force kill?",
                                               buttons=['ok', 'cancel'],
                                               default='cancel'):
                return

            mfixproc = self.mfixproc
            if not mfixproc:
                return
            if mfixproc.state() != QProcess.Running:
                return

            try:
                mfixproc.terminate()
            except OSError as err:
                log = logging.getLogger(__name__)
                log.error("Error terminating process: %s", err)
            # retry force kill if necessary
            QTimer.singleShot(100, force_kill)

        # initial force kill attempt
        QTimer.singleShot(500, force_kill)


    def start_command(self, is_pymfix, cmd, cwd, env):
        """Start MFIX in QProcess"""

        self.is_pymfix, self.cmd, self.cwd, self.env = is_pymfix, cmd, cwd, env
        self.mfixproc = QProcess()
        self.mfixproc.setWorkingDirectory(self.cwd)
        mfix_stop = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        try:
            os.remove(mfix_stop)
        except OSError:
            pass

        def slot_start():
            self.mfix_pid = self.mfixproc.pid() # Keep a copy because it gets reset
            msg = "MFIX process %d is running" % self.mfix_pid
            self.parent.update_run_options_signal.emit(msg)
            log = logging.getLogger(__name__)
            log.debug("Full MFIX startup parameters: %s", ' '.join(self.cmd))
            log.debug("starting mfix output monitor threads")
            if is_pymfix:
                self.timer = QTimer()
                self.timer.setInterval(1000)
                self.timer.timeout.connect(self.update_status)
                self.timer.start()
                self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)

        self.mfixproc.started.connect(slot_start)

        def slot_read_out():
            out_str = bytes(self.mfixproc.readAllStandardOutput()).decode('utf-8')
            self.parent.stdout_signal.emit(out_str)
        self.mfixproc.readyReadStandardOutput.connect(slot_read_out)

        def slot_read_err():
            err_str = bytes(self.mfixproc.readAllStandardError()).decode('utf-8')
            self.parent.stderr_signal.emit(err_str)
        self.mfixproc.readyReadStandardError.connect(slot_read_err)

        def slot_finish(status):
            msg = "MFIX process %s has stopped"%self.mfix_pid # by now mfixproc.pid() is 0
            self.mfix_pid = 0
            #self.parent.stdout_signal.emit("MFIX (pid %s) has stopped" % self.mfixproc.pid())
            self.mfixproc = None
            if is_pymfix:
                self.timer.stop()
                self.timer = None
            self.parent.update_run_options_signal.emit(msg)
        self.mfixproc.finished.connect(slot_finish)

        def slot_error(error):
            if error == QProcess.FailedToStart:
                msg = "Process failed to start"
            elif error == QProcess.Crashed:
                msg = "Process exit"
            elif error == QProcess.Timedout:
                msg = "Process timeout"
            elif error in (QProcess.WriteError, QProcess.ReadError):
                msg = "Process communication error"
            else:
                msg = "Unknown error"
            log = logging.getLogger(__name__)
            log.warn(msg)
            self.mfixproc = None
            self.parent.stderr_signal.emit(msg) # make the message print in red
            self.parent.update_run_options_signal.emit('')


        self.mfixproc.error.connect(slot_error)

        self.mfixproc.start(self.cmd[0], self.cmd[1:])

    def terminate_pymfix(self):
        """update the status of  the pymfix monitor"""
        url = 'http://localhost:5000/exit'
        values = {'timeout' : '1',}

        data = urlencode(values)
        req = Request(url, data)
        response = urlopen(req)
        response_str = response.read()
        log = logging.getLogger(__name__)
        log.debug("response_str is %s", response_str)
        # self.status = json.loads(status_str)
        # log.debug("status is %s", self.status)

    def update_status(self):
        """update the status of  the pymfix monitor"""
        status_str = urlopen('http://localhost:5000/status').read()
        log = logging.getLogger(__name__)
        log.debug("status_str is %s", status_str)
        self.status = json.loads(status_str)
        log.debug("status is %s", self.status)
        # FIXME: print JSON for now, plot it later
        import pprint
        status_str = pprint.PrettyPrinter(indent=4, width=40).pformat(self.status)
        self.parent.ui.residuals.setText(status_str)

class Monitor(object):
    """class for monitoring available MFIX executables and output files"""

    def __init__(self, parent):
        self.parent = parent
        self.cache = {}
        self.executables = None
        self.outputs = None
        self.res_exists = False
        self._update_executables()

    def get_executables(self):
        """returns a dict mapping full [mfix|pymfix] paths
        to configuration options."""

        return self.executables

    def _update_executables(self):
        """update self.executables"""

        def mfix_print_flags(mfix_exe, cache=self.cache):
            """Determine mfix configuration by running mfix --print-flags.  Cache results"""
            try: # Possible race, file may have been deleted/renamed since isfile check!
                stat = os.stat(mfix_exe)
            except OSError as err:
                log = logging.getLogger(__name__)
                log.exception("could not run %s --print-flags", mfix_exe)
                return ''

            cached_stat, cached_flags = cache.get(mfix_exe, (None, None))
            if cached_stat and cached_stat == stat:
                return cached_flags

            popen = subprocess.Popen(mfix_exe + " --print-flags",
                                     cwd = self.parent.get_project_dir(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True)
            (out, err) = popen.communicate()
            flags = '' if err else out
            cache[mfix_exe] = (stat, flags)
            return flags

        config_options = {}

        for d in self.parent.exe_watcher.directories():
            for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                exe = os.path.abspath(os.path.join(d, name))
                if os.path.isfile(exe):
                    log = logging.getLogger(__name__)
                    log.debug("found %s executable in %s", name, d)
                    config_options[exe] = str(mfix_print_flags(exe))

        self.executables = config_options

    def get_res(self):
        if not self.parent.get_project_dir():
            return
        pattern = os.path.join(self.parent.get_project_dir(), '*.RES')
        return glob.glob(pattern)

    def get_outputs(self, patterns=[]):
        project_dir = self.parent.get_project_dir()
        if project_dir is None:
            return
        if len(patterns) == 0:
            patterns = [
                '*.LOG', '*.OUT', '*.RES', '*.SP?', 'MFIX.STOP',
                '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        outputs = []
        for pat in patterns:
            outputs.extend(glob.glob(os.path.join(project_dir, pat)))
        return outputs
