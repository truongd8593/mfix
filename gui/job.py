"""class to manage external MFIX process"""

import os
import time
import json
import logging

#
DEFAULT_TIMEOUT = 2 # seconds, for all socket ops
import socket
socket.setdefaulttimeout(DEFAULT_TIMEOUT)

log = logging.getLogger(__name__)

from qtpy.QtCore import QProcess, QTimer, QUrl
from qtpy.QtNetwork import QNetworkRequest, QNetworkAccessManager

class Job(object):
    """class for monitoring an MFIX job"""

    def __init__(self, parent):
        self.parent = parent
        self.cmd = None
        self.cwd = None
        self.env = None
        self.is_pymfix = False
        self.mfixproc = None
        self.mfix_pid = None # Keep track of pid separately
        self.status = {}
        self.timer = None
        self.pymfix_url = 'http://localhost:5000'
        self.control_manager = QNetworkAccessManager()
        self.status_manager = QNetworkAccessManager()

        def slot_update_status(reply):
            """update self.status when request finishes"""
            response_str = str(reply.readAll())
            try:
                self.status = json.loads(response_str)
                log.debug("status is %s", self.status)
                # FIXME: print JSON for now, plot it later
                import pprint
                status_str = pprint.PrettyPrinter(indent=4, width=50).pformat(self.status)
                self.parent.ui.residuals.setText(status_str)
            except ValueError:
                self.status.clear()
                if response_str:
                    log.error("could not decode JSON: %s", response_str)
            self.parent.update_run_options()

        def slot_control(reply):
            self.parent.update_run_options()
        self.control_manager.finished.connect(slot_control)
        self.status_manager.finished.connect(slot_update_status)

    def is_running(self):
        """indicate whether an MFIX/pymfix job is running
        Noted, pymfix will be "running" even if paused"""
        ret = (self.mfixproc is not None
               and self.mfixproc.state() == QProcess.Running)
        return ret

    def is_paused(self):
        """indicate whether pymfix is paused"""
        if not self.is_pymfix:
            return False
        return self.status.get('paused', False)

    def is_pausable(self):
        """indicate whether pymfix is pausable.
        Pymfix starts by reading the model, and
        will not be pausable at first"""
        if not self.is_pymfix:
            return False
        # Wait until least one update occurred
        return self.status.get('paused', None) != None

    def pause(self):
        "pause pymfix job"
        if not self.is_pymfix:
            log.error("pause() called on non-pymfix job")
            return
        qurl = QUrl('%s/pause' % self.pymfix_url)
        req = QNetworkRequest(qurl)
        self.control_manager.put(req, b"")

    def unpause(self):
        "unpause pymfix job"
        if not self.is_pymfix:
            log.error("unpause() called on non-pymfix job")
            return
        qurl = QUrl('%s/unpause' % self.pymfix_url)
        req = QNetworkRequest(qurl)
        self.control_manager.put(req, b"")

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        if self.is_pymfix:
            self.terminate_pymfix()
        else:
            mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
            try:
                open(mfix_stop_file, "ab").close()
            except OSError:
                pass

        while self.is_running():
            t0 = time.time()
            self.mfixproc.waitForFinished(2000)
            t1 = time.time()
            if self.is_running():
                log.warn("mfix still running after %.2f ms", (1000*(t1-t0)))
            else:
                break
            if self.parent.message(text="MFIX is not responding. Force kill?",
                                   buttons=['ok', 'cancel'],
                                   default='cancel')  != 'ok':
                log.warn("not killing mfix process %d at user request" % self.mfix_pid)
                return
            try:
                self.status.clear()
                log.warn("killing mfix process %d", self.mfix_pid)
                self.mfixproc.terminate()
            except OSError as err:
                log.error("Error terminating process: %s", err)

        self.mfixproc = None


    def start_command(self, is_pymfix, cmd, cwd, env):
        """Start MFIX in QProcess"""

        self.is_pymfix, self.cmd, self.cwd, self.env = is_pymfix, cmd, cwd, env
        if self.mfixproc:
            log.warn("Cannot start process, already have an mfix process")
            return

        mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        if os.path.exists(mfix_stop_file):
            try:
                os.remove(mfix_stop_file)
            except OSError:
                log.error("Cannot remove", mfix_stop_file)
                return

        self.cmdline = ' '.join(self.cmd) # cmd is a list
        self.mfixproc = QProcess()
        if not self.mfixproc:
            log.warn("QProcess creation failed")
            return
        self.mfixproc.setWorkingDirectory(self.cwd)

        def slot_start():
            self.mfix_pid = self.mfixproc.pid() # Keep a copy because it gets reset
            msg = "MFIX process %d is running" % self.mfix_pid
            self.parent.update_run_options_signal.emit(msg)
            log.debug("Full MFIX startup parameters: %s", self.cmdline)
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
            self.status.clear()
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
            self.status.clear()
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
            self.mfixproc = None
            self.parent.stderr_signal.emit(msg) # make the message print in red
            self.parent.update_run_options_signal.emit('')


        self.mfixproc.error.connect(slot_error)

        self.mfixproc.start(self.cmd[0], self.cmd[1:])

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        arg = 'enable' if state else 'disable'
        qurl = QUrl('%s/logging/%s' % (self.pymfix_url, arg))
        req = QNetworkRequest(qurl)
        self.control_manager.post(req, b"timeout=1")

    def terminate_pymfix(self):
        """update the status of  the pymfix monitor"""

        qurl = QUrl('%s/exit' % self.pymfix_url)
        req = QNetworkRequest(qurl)
        self.control_manager.post(req, b"timeout=1")

    def update_status(self):
        """update the status of  the pymfix monitor"""
        qurl = QUrl('%s/status' % self.pymfix_url)
        req = QNetworkRequest(qurl)
        self.status_manager.get(req)
