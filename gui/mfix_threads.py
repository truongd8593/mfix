"""classes to manage external MFIX process and monitor files"""

import glob
import json
import logging
import os
import subprocess

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

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
        self.status = None
        self.timer = None

    def is_running(self):
        """indicate whether an MFIX job is running"""
        return self.mfixproc is not None

    def is_paused(self):
        """indicate whether pymfix is paused"""
        if not self.is_pymfix:
            return False
        return getattr(self.status, 'paused', False)

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""

        mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        try:
            open(mfix_stop_file, "ab").close()
        except OSError:
            pass

        def force_kill():
            confirm = self.parent.message(text="MFIX is not responding. Force kill?",
                                   buttons=['ok','cancel'],
                                   default='cancel')

            if confirm != 'ok':
                return

            mfixproc = self.mfixproc
            if not mfixproc:
                return
            try:
                mfixproc.terminate()
            except OSError as err:
                log = logging.getLogger(__name__)
                log.error("Error terminating process: %s", err)
            self.parent.stdout_signal.emit("Terminating MFIX (pid %s)" % mfixproc.pid())

            # python >= 3.3 has subprocess.wait(timeout), which would be good to loop wait
            # os.waitpid has a nohang option, but it's not available on Windows
            mfixproc.waitForFinished() # FIXME timeout
            self.mfixproc = None
        QTimer.singleShot(100, force_kill)

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
            msg = "MFIX (pid %d) is running" % self.mfix_pid
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
            self.parent.stdout_signal.emit(err_str)
        self.mfixproc.readyReadStandardError.connect(slot_read_err)

        def slot_finish(status):
            msg = "MFIX (pid %s) has stopped"%self.mfix_pid # by now mfixproc.pid() is 0
            self.mfix_pid = 0
            #self.parent.stdout_signal.emit("MFIX (pid %s) has stopped" % self.mfixproc.pid())
            self.mfixproc = None
            if is_pymfix:
                self.timer.stop()
                self.timer = None
            self.parent.update_run_options_signal.emit(msg)
        self.mfixproc.finished.connect(slot_finish)

        self.mfixproc.start(self.cmd[0], self.cmd[1:])

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

        # Check system PATH dirs first
        PATH = os.environ.get("PATH")
        if PATH:
            dirs = set(PATH.split(os.pathsep))
        else:
            dirs = set()

        # Look in subdirs of build dir
        build_dir = os.path.join(get_mfix_home(), 'build')
        if os.path.exists(build_dir):
            for subdir in os.listdir(build_dir):
                dirs.add(os.path.join(build_dir, subdir, 'build-aux'))

        # Check run_dir
        project_dir = self.parent.get_project_dir()
        if project_dir:
            dirs.add(project_dir)
        # Check mfix home
        dirs.add(get_mfix_home())

        # Now look for mfix/pymfix in these dirs
        for dir_ in dirs:
            for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                exe = os.path.abspath(os.path.join(dir_, name))
                if os.path.isfile(exe):
                    log = logging.getLogger(__name__)
                    log.debug("found %s executable in %s", name, dir_)
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
