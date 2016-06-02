"""class to manage external MFIX process"""

import os
import time
import glob
import json
import logging

log = logging.getLogger(__name__)

try:
    # For Python 3.0 and later
    from urllib.request import urlopen, Request
    from urllib.parse import urlencode
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib import urlencode
    from urllib2 import urlopen, Request

from qtpy.QtCore import QProcess, QTimer

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

    def is_running(self):
        """indicate whether an MFIX job is running"""
        ret = (self.mfixproc is not None
               and self.mfixproc.state() == QProcess.Running)
        return ret

    def is_paused(self):
        """indicate whether pymfix is paused"""
        if not self.is_pymfix:
            return False
        return self.status.get('paused', False)

    def pause(self):
        "pause MFIX"
        if not self.is_pymfix:
            return
        # log = logging.getLogger(__name__)
        req = Request(url='%s/pause' % self.pymfix_url, data='')
        req.get_method = lambda: 'PUT'
        resp = urlopen(req)
        # log.info("response status: %s", resp.status)
        # log.info("response reason: %s", resp.reason)
        resp.close()
        self.parent.update_run_options()

    def unpause(self):
        "pause MFIX"
        if not self.is_pymfix:
            return
        # log = logging.getLogger(__name__)
        req = Request(url='%s/unpause' % self.pymfix_url, data='')
        req.get_method = lambda: 'PUT'
        resp = urlopen(req)
        # log.info("response status: %s", resp.status)
        # log.info("response reason: %s", resp.reason)
        resp.close()
        self.parent.update_run_options()

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

        def force_kill():
            """attempt to kill job; return True if further force_kill() calls may be necessary"""
            if not 'ok' in self.parent.message(text="MFIX is not responding. Force kill?",
                                               buttons=['ok', 'cancel'],
                                               default='cancel'):
                return

            if not self.is_running():
                return
            try:
                self.mfixproc.terminate()
            except OSError as err:
                log = logging.getLogger(__name__)
                log.error("Error terminating process: %s", err)

        while self.is_running():
            t0 = time.time()
            self.mfixproc.waitForFinished(1000)
            t1 = time.time()
            if self.is_running():
                log.warn("mfix still running after %.2f ms forcing kill" % (1000*(t1-t0)))
                force_kill()
        self.mfixproc = None


    def start_command(self, is_pymfix, cmd, cwd, env):
        """Start MFIX in QProcess"""

        self.is_pymfix, self.cmd, self.cwd, self.env = is_pymfix, cmd, cwd, env
        if self.mfixproc:
            log.warn("Cannot start process, already have an mfix process")
            return
        self.mfixproc = QProcess()
        if not self.mfixproc:
            log.warn("QProcess creation failed")
            return
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
        url = '%s/exit' % self.pymfix_url
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
        status_str = urlopen('%s/status' % self.pymfix_url).read()
        log = logging.getLogger(__name__)
        log.warning("status_str is %s", status_str)
        try:
            self.status = json.loads(status_str)
            log.debug("status is %s", self.status)
            # FIXME: print JSON for now, plot it later
            import pprint
            status_str = pprint.PrettyPrinter(indent=4, width=40).pformat(self.status)
            self.parent.ui.residuals.setText(status_str)
        except ValueError:
            log.error("could not decode JSON: %s", status_str)
        self.parent.update_run_options()
