"""Threads supporting MFIX-GUI"""

import glob
import logging
import os
import subprocess
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

from tools.general import get_mfix_home

# Message types, emitted with 'line_printed'

message_normal = 0
message_hi_vis = 1
message_error = 2
message_stdout = 3
message_stderr = 4

class MfixJobManager(QThread):

    def __init__(self, parent):
        super(MfixJobManager, self).__init__(parent)
        self.parent = parent
        self.cmd = None
        self.cwd = None
        self.mfixproc = None
        self.name = name

    def is_running(self):
        """indicate whether an MFIX job is running"""
        return self.mfixproc is not None

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
            log = logging.getLogger(__name__)
            log.error("Error terminating process: %s", err)
            self.line_printed.emit("Terminating MFIX process (pid %s)" %
                                mfixpid, 'blue', '')

        # python >= 3.3 has subprocess.wait(timeout), which would be good to loop wait
        # os.waitpid has a nohang option, but it's not available on Windows
        if self.mfixproc:
            self.mfixproc.wait() # FIXME timeout
        self.mfixproc = None
        self.parent.update_run_options_signal.emit()

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
        """Run a subprocess, with executable and arguments obtained
        from class-local self.cmd set in start_command()

        :return: None"""

        if self.cmd:
            mfixproc_pid = self.mfixproc.pid
            self.line_printed.emit(
                                "MFIX (pid %d) is running" % mfixproc_pid,
                                'blue', '')
            log = logging.getLogger(__name__)
            log.debug("Full MFIX startup parameters: %s" % ' '.join(self.cmd))
            log.debug("starting mfix output monitor threads")

            stdout_thread = MfixOutput(
                                name='stdout',
                                pipe=self.mfixproc.stdout,
                                signal=self.line_printed,
                                font='Courier') # TODO find good cross-platform fixed width font

            stderr_thread = MfixOutput(
                                name='stderr',
                                pipe=self.mfixproc.stderr,
                                signal=self.line_printed,
                                color='red',
                                font='Courier')

            stdout_thread.start()
            stderr_thread.start()

            self.parent.update_run_options_signal.emit()

            if self.mfixproc:
                self.mfixproc.wait()
            self.mfixproc = None

            self.line_printed.emit(
                                "MFIX (pid %s) has stopped" % mfixproc_pid,
                                'blue', '')
            self.parent.update_run_options_signal.emit()

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
        log = logging.getLogger(__name__)
        log.debug("Started thread %s" % name)
        self.stopped = False # See comment above
        self.name = name
        self.pipe = pipe
        self.name = name
        self.signal = signal # Share signal with parent class

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
                    log = logging.getLogger(__name__)
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
                log = logging.getLogger(__name__)
                log.debug("cannot retrieve residuals; pymfix process must have terminated.")
                self.job_done = True
                return
            self.sleep(1)
            self.residuals_changed.emit('update')
