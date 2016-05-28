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

# Message types, emitted with 'line_printed'

message_normal = 0
message_hi_vis = 1
message_error = 2
message_stdout = 3
message_stderr = 4

# --- Threads ---
class MfixThread(QThread):

    line_printed = pyqtSignal(str, int) # message, message type
    mfix_running = pyqtSignal(bool)  # running/not running

    def __init__(self, parent, name="MfixThread"):
        super(MfixThread, self).__init__(parent)
        self.stopped = False  # TODO use this to break out of loop (quit)
        self.cmd = self.cwd = None
        self.mfixproc = None
        self.name = name

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        if not self.mfixproc:
            log.debug("stop_mfix: mfixproc==None")
            return

        mfixpid = self.mfixproc.pid
        self.line_printed.emit("Terminating MFIX process (pid %s)" % mfixpid, message_hi_vis)

        try:
            self.mfixproc.terminate()
        except OSError as err:
            msg = "Error terminating process: %s" % err
            self.line_printed.emit(msg, message_error)
            return

        # python >= 3.3 has subprocess.wait(timeout), which would be good to loop wait
        # os.waitpid has a nohang option, but it's not available on Windows
        if self.mfixproc:
            self.mfixproc.wait() # FIXME timeout
        self.mfixproc = None
        self.mfix_running.emit(False)

    def start_command(self, cmd, cwd, env):
        """Start a command in subprocess, and start the thread to monitor it"""

        # TODO:  maybe instead of pipes, we should just write to files.  Then we could
        # leave mfix running after exiting the gui, without getting SIGPIPE

        self.cmd, self.cwd, self.env = cmd, cwd, env
        self.mfixproc = subprocess.Popen(self.cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines = True,
                                        shell=False, cwd=self.cwd,
                                        env=self.env)
        self.start() # calls 'run'


    def run(self):
        """Collect stdout/stderr from running command, using two helper threads.
        Waits for process to end"""

        mfixproc_pid = self.mfixproc.pid
        self.line_printed.emit("MFIX (pid %d) is running" % mfixproc_pid, message_hi_vis)
        self.mfix_running.emit(True)
        log.debug("MFIX command: %s" % ' '.join(self.cmd))
        log.debug("starting mfix output monitor threads")

        stdout_thread = OutputHelper(name='stdout',
                                     pipe=self.mfixproc.stdout,
                                     signal=self.line_printed,
                                     message_type=message_stdout)

        stderr_thread = OutputHelper(name='stderr',
                                     pipe=self.mfixproc.stderr,
                                     signal=self.line_printed,
                                     message_type=message_stderr)

        stdout_thread.start()
        stderr_thread.start()

        #self.update_run_options.emit()

        if self.mfixproc:
            self.mfixproc.wait() # TODO: detect/handle hung MFIX
        self.mfixproc = None

        self.line_printed.emit("MFIX (pid %s) has stopped" % mfixproc_pid, message_hi_vis)
        self.mfix_running.emit(False)

        # Allow remaining output to be collected - OutputHelper threads exit on reaching EOF
        stderr_thread.wait()
        stdout_thread.wait()
        #stderr_thread.terminate()
        #stdout_thread.terminate()



class OutputHelper(QThread):
    """Handle streaming from a pipe and emiting read lines into a signal handler.
    emit line_printed for each line, with a default message type
    params:
        name: Name of thread (str)
        pipe: pipe to child process
        message_type: message type to emit in line_printed (int) """

    def __init__(self, name, pipe, signal, message_type=message_normal):
        super(OutputHelper, self).__init__()
        self.message_type = message_type
        log.debug("Started thread %s" % name)
        self.stopped = False # See comment above
        self.name = name
        self.pipe = pipe
        self.name = name
        self.signal = signal # Share signal with parent class

    #def __del__(self):
    #    self.terminate()
    #    #    # I suspect this is the source of QThread::wait errors
    #    #self.wait()

    def run(self):
        while not self.stopped: # TODO: allow .quit() to set stopped=T
            line = str(self.pipe.readline()).lower()
            if not line:
                break # Process has exited - emit signal?
            self.signal.emit(line, self.message_type)


class MonitorThread(QThread):

    executables_changed = pyqtSignal()
    outputs_changed = pyqtSignal()

    def __init__(self, parent, name="MonitorThread"):
        QThread.__init__(self)
        self.parent = parent
        self.stopped = False # see comments above
        self.mfix_home = get_mfix_home()
        self.cache = {}
        self.executables = self.get_executables()
        self.outputs = self.get_outputs()
        self.name = name

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
                dirs.add(os.path.join(build_dir, subdir, 'build-aux'))

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
        pattern = os.path.join(self.parent.get_project_dir(),'*.RES')
        return glob.glob(pattern)

    def get_outputs(self, patterns=[]):
        project_dir = self.parent.get_project_dir()
        if project_dir is None:
            return
        if len(patterns) == 0:
            patterns = [
                '*.LOG', '*.OUT', '*.RES', '*.SP?',
                '*.pvd', '*.vtp', 'VTU_FRAME_INDEX.TXT']
        outputs = []
        for pat in patterns:
            outputs.extend(glob.glob(os.path.join(project_dir, pat)))
        return outputs

    def run(self):
        self.outputs_changed.emit()
        self.executables_changed.emit()
        while True:
            tmp = self.get_outputs()
            if tmp != self.outputs:
                self.outputs = tmp
                self.outputs_changed.emit()
            tmp = self.get_executables()
            if tmp != self.executables:
                self.executables = tmp
                self.executables_changed.emit()
            self.sleep(1) # don't use time.sleep in qthread


class UpdateResidualsThread(QThread):

    residuals_changed = pyqtSignal(object)

    def run(self):
        while True:
            self.job_done = False
            try:
                self.residuals = urlopen('http://localhost:5000/residuals').read()
            except Exception:
                log.debug("cannot retrieve residuals; is pymfix running")
                self.job_done = True
                return
            self.sleep(1)
            self.residuals_changed.emit('update')
