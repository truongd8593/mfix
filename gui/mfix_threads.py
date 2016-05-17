"""Threads supporting MFIX-GUI"""

import os
import subprocess
import glob
import time

# import qt
from qtpy import QtCore, QtWidgets, QtGui, PYQT4, PYQT5
from qtpy.QtCore import QThread, pyqtSignal

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen


import logging
log = logging.getLogger(__name__)

from tools.general import get_mfix_home


# --- Threads ---
class MfixThread(QThread):

    line_printed = pyqtSignal(str, str, str)
    update_run_options = pyqtSignal()

    def __init__(self, parent):
        super(MfixThread, self).__init__(parent)
        self.cmd = None
        self.cwd = None
        self.mfixproc = None

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        mfixpid = self.mfixproc.pid
        try:
            self.mfixproc.terminate()
        except OSError as err:
            log.error("Error terminating process: %s", err)
        self.line_printed.emit("Terminating MFIX process (pid %s)" %
                                mfixpid, 'blue', '')

        # python >= 3.3 has subprocess.wait(timeout), which would be good to loop wait
        # os.waitpid has a nohang option, but it's not available on Windows
        if self.mfixproc:
            self.mfixproc.wait()
        self.mfixproc = None

        self.update_run_options.emit()

    def start_command(self, cmd, cwd, env):
        """Initialize local logging object, set local command and
        working directory to those provided.
        :param cmd: List to be passed as first parameter to subprocess.Popen()
        :param cwd: The working directory for the command (str)
        :param env: Variables to be set child process environment (dict)
        :return: None"""
        self.cmd = cmd
        self.cwd = cwd
        self.env = env
        self.mfixproc = subprocess.Popen(self.cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines = True,
                                        shell=False, cwd=self.cwd,
                                        env=self.env)
        self.start()

    def run(self):
        """Run a subprocess, with executable and arguments obtained
        from class-local self.cmd set in start_command()

        :return: None"""

        if self.cmd:
            mfixproc_pid = self.mfixproc.pid
            self.line_printed.emit(
                                "MFIX (pid %d) is running" % mfixproc_pid,
                                'blue', '')
            log.debug("Full MFIX startup parameters: %s" % ' '.join(self.cmd))
            log.debug("starting mfix output monitor threads")

            stdout_thread = MfixOutput(
                                name='stdout',
                                pipe=self.mfixproc.stdout,
                                signal=self.line_printed,
                                font='Courier')

            stderr_thread = MfixOutput(
                                name='stderr',
                                pipe=self.mfixproc.stderr,
                                signal=self.line_printed,
                                color='red',
                                font='Courier')

            stdout_thread.start()
            stderr_thread.start()

            self.update_run_options.emit()

            if self.mfixproc:
                self.mfixproc.wait()
            self.mfixproc = None

            self.line_printed.emit(
                                "MFIX (pid %s) has stopped" % mfixproc_pid,
                                'blue', '')
            self.update_run_options.emit()

            stderr_thread.quit()
            stdout_thread.quit()



class MfixOutput(QThread):
    """Generic class to handle streaming from a pipe and emiting read
        lines into a signal handler.

        :param name: Name of thread
        :type name: str
        :param pipe: Iterable object with readline method, assumed to be subprocess.PIPE
        :type pipe: subprocess.PIPE
        :param signal: Signal to be used to pass lines read from pipe
        :type signal: QtCore.pyqtSignal object
    #TODO: let's define message types which map to color/font settings
        :param color: Color to set when emitting signals
        :type color: str """

    def __init__(self, name, pipe, signal, color=None, font=None):
        super(MfixOutput, self).__init__()
        log.debug("Started thread %s" % name)
        self.name = name
        self.signal = signal
        self.pipe = pipe
        self.color = color
        self.font = font

    def __del__(self):
        #    # I suspect this is the source of QThread::wait errors
        # why needed?  - cgw
        self.wait()

    def run(self):
        lines_iterator = iter(self.pipe.readline, b"")
        for line in lines_iterator:
            lower = line.lower()
            if any(x in lower for x in ('error', 'warn', 'fail')):
                color='red'
            else:
                color = self.color
            self.signal.emit(str(line), color, self.font)



class MonitorThread(QThread):

    sig = pyqtSignal()

    def __init__(self, parent):
        QThread.__init__(self)
        self.parent = parent
        self.mfix_home = get_mfix_home()
        self.cache = {}
        self.executables = self.get_executables()
        self.outputs = self.get_outputs()

    def get_executables(self):
        """returns a dict mapping full [mfix|pymfix] paths
        to configuration options."""

        def mfix_print_flags(mfix_exe, cache=self.cache):
            """Determine mfix configuration by running mfix --print-flags.  Cache results"""
            try: # Possible race, file may have been deleted/renamed since isfile check!
                stat = os.stat(mfix_exe)
            except OSError:
                return ''

            cached_stat, cached_flags = cache.get(mfix_exe, (None, None))
            if cached_stat and cached_stat == stat:
                return cached_flags

            popen = subprocess.Popen(mfix_exe + " --print-flags",
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
        build_dir = os.path.join(self.mfix_home,'build')
        if os.path.exists(build_dir):
            for subdir in os.listdir(build_dir):
                dirs.add(os.path.join(build_dir, subdir))

        # Check run_dir
        project_dir = self.parent.get_project_dir()
        if project_dir:
            dirs.add(project_dir)
        # Check mfix home
        dirs.add(self.mfix_home)

        # Now look for mfix/pymfix in these dirs
        for dir_ in dirs:
            for name in 'mfix', 'mfix.exe', 'pymfix', 'pymfix.exe':
                exe = os.path.abspath(os.path.join(dir_, name))
                if os.path.isfile(exe):
                    log.debug("found %s executable in %s" % (name, dir_))
                    config_options[exe] = str(mfix_print_flags(exe))

        return config_options

    def get_res(self):
        if not self.parent.get_project_dir():
            return
        globb = os.path.join(self.parent.get_project_dir(),'*.RES')
        return glob.glob(globb)

    def get_outputs(self, patterns=[]):
        project_dir = self.parent.get_project_dir()
        if project_dir is None:
            return
        if len(patterns) == 0:
            patterns = [
                '*.LOG', '*.OUT', '*.RES', '*.SP?',
                '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        output_paths = [glob.glob(os.path.join(project_dir, n)) for n in patterns]
        outputs = []
        for path in output_paths:
            outputs += path
        return outputs

    def run(self):
        self.sig.emit()
        while True:  # FIXME - should be two different signals so we can
                     # determine what changed
            tmp = self.get_outputs()
            if tmp != self.outputs:
                self.outputs = tmp
                self.sig.emit()
            tmp = self.get_executables()
            if tmp != self.executables:
                self.executables = tmp
                self.sig.emit()
            time.sleep(1)


class UpdateResidualsThread(QThread):

    sig = pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            try:
                self.residuals = urlopen('http://localhost:5000/residuals').read()
            except Exception:
                log.debug("cannot retrieve residuals; pymfix process must have terminated.")
                self.job_done = True
                return
            time.sleep(1)
            self.sig.emit('update')
