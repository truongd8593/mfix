"""class to manage external MFIX process"""
import traceback

import json
import logging
from functools import wraps
import os
import tempfile
import time

DEFAULT_TIMEOUT = 2 # seconds, for all socket ops
import socket
socket.setdefaulttimeout(DEFAULT_TIMEOUT)

log = logging.getLogger(__name__)

from qtpy.QtCore import QProcess, QProcessEnvironment, QTimer, QUrl
from qtpy.QtNetwork import QNetworkRequest, QNetworkAccessManager


SUPPORTED_PYMFIXPID_FIELDS = ['url', 'pid', 'token']

class PymfixPID(object):
    """ Class to obtain Pymfix connection details"""

    def __init__(self, job_pid_file):

        job_pid_file = '%s.pid' % job_pid_file

        with open(job_pid_file) as pidfile:
            log.debug('opened pid file %s' % os.path.basename(job_pid_file))
            for line in pidfile.read().split('\n'):
                try:
                    (k, v) = line.split('=')
                    if k in SUPPORTED_PYMFIXPID_FIELDS:
                        log.debug('reading %s from pidfile: %s' % (k, v))
                        setattr(self, k, v)
                except ValueError:
                    continue
                except:
                    log.error('error processing %s' % job_name)
                    traceback.print_exc()


class PymfixAPI(QNetworkAccessManager):
    """ Class to extend QNetworkAccessManager with pymfix connection details
        set transparently.
    """

    def __init__(self, job_name, project_dir, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        self.connected = False
        self.pymfix = None
        self.pymfix_url = None
        self.job_name = job_name
        self.project_dir = project_dir
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}

    def connected(f):
        @wraps(f)
        def func(*args, **kwargs):
            self = args[0]
            self.connect()
            return f(*args, **kwargs)
        return func

    def connect(self):
        if bool(self.pymfix):
            return
        try:
            log.debug('connecting to API for job %s' % self.job_name)
            self.pymfix = PymfixPID(os.path.join(self.project_dir, self.job_name))
            # TODO: verify API is listening -- self.get('status')
            log.info('API connected: %s' % self.pymfix.url)
        except:
            log.warn("API connection to job %s failed" % self.job_name)
            traceback.print_exc()

    @connected
    def api_request(self, method, endpoint, data):
        log.debug("request sent: method=%s endpoint=%s" % (method, endpoint))
        method = str(method).lower()
        method = self.methods[method]
        req = QNetworkRequest(QUrl('%s/%s' % (self.pymfix.url, endpoint)))
        if self.pymfix.token:
            (name, value) = self.pymfix.token.split(':')
            req.setRawHeader(name, value)
        if data is not None:
            return method(req, data)
        else:
            return method(req)

    def get(self, endpoint, data=None):
        return self.api_request('get', endpoint, data)

    def put(self, endpoint, data=b""):
        return self.api_request('put', endpoint, data)

    def post(self, endpoint, data=b""):
        return self.api_request('post', endpoint, data)

class JobManager(object):

    """class for managing and monitoring an MFIX job"""

    def __init__(self, job_name, parent):
        self.parent = parent
        self.mfixproc = None
        self.mfix_pid = None # Keep track of pid separately
        self.status = {}
        self.timer = None
        self.job_name = job_name
        self.project_dir = self.parent.get_project_dir()
        self.api = PymfixAPI(self.job_name, self.project_dir, ignore_ssl_errors=True)
        self.pymfix_url = self.api.pymfix_url

        def slot_update_status(reply):
            """update self.status when request finishes"""
            print("slot_update_status")
            response_str = reply.readAll().data().decode('utf-8')
            try:
                self.status = json.loads(response_str)
                log.debug("status is %s", self.status)
                # FIXME: print JSON for now, plot it later
                import pprint
                status_str = pprint.PrettyPrinter(indent=4,
                                                width=50).pformat(self.status)
                self.parent.ui.residuals.setText(status_str)
            except ValueError:
                self.status.clear()
                log.error("could not decode JSON: %s", response_str)
            self.parent.update_run_options()

        def slot_control(reply):
            self.parent.update_run_options()

        def slot_ssl_error(reply):
            """ Handler for SSL connection errors. Check self.ignore_ssl_errors
                and continue as appropriate"""
            if self.api.ignore_ssl_errors:
                log.debug('call to %s:%s completed with ignored SSL errors' % \
                         (self.job_name, self.endpoint))
                reply.ignoreSsl()
            else:
                log.warn('call to %s:%s aborted due to SSL errors' % \
                         (self.job_name, self.endpoint))
                raise # find appropriate Qt exception or make this meaningful

        def slot_finished(reply):
            log.debug('call to %s:%s completed' % (self.job_name,
                                                   self.api.endpoint))
            self.api.finished.emit()

        self.api.sslErrors.connect(slot_ssl_error)

        self.api.finished.connect(slot_control)
        self.api.finished.connect(slot_update_status)

    def is_running(self):
        """indicate whether an MFIX/pymfix job is running
        Noted, pymfix will be "running" even if paused"""
        ret = self.pymfix_url or (self.mfixproc is not None and
                                  self.mfixproc.state() == QProcess.Running)
        return ret

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def pause(self):
        "pause pymfix job"
        if not self.pymfix_url:
            log.error("pause() called on non-pymfix job")
            return
        self.api.put('pause')

    def unpause(self):
        "unpause pymfix job"
        if not self.pymfix_url:
            log.error("unpause() called on non-pymfix job")
            return
        self.api.put('unpause')

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        if self.pymfix_url:
            self.terminate_pymfix()
            self.disconnect()
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
                log.warn("not killing mfix process %s at user request" % self.mfix_pid)
                return
            try:
                self.status.clear()
                log.warn("killing mfix process %d", self.mfix_pid)
                if self.mfixproc:
                    self.mfixproc.terminate()
            except OSError as err:
                log.error("Error terminating process: %s", err)

        self.mfixproc = None

    def connect(self, pymfix_url=None):
        """Connect to existing pymfix process"""

        # TODO: refactor instantiation to set whether this is pymfix or mfix
        if not bool(self.api.pymfix):
            self.api.connect()

        self.pymfix_url = self.api.pymfix_url

        """
        if not pymfix_url:
            log.warn("Cannot connect to process, invalid URL")
            return

        if self.pymfix_url:
            log.warn("Cannot connect to process, already connected to mfix process")
            return
        """

        self.timer = QTimer()
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.update_status)
        self.timer.start()
        self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)

    def disconnect(self):
        """stop pymfix update timer"""
        if self.timer:
            self.timer.stop()
            self.timer = None
        self.pymfix_url = None

    def transform_template(self, text, cmd):
        cores = self.parent.run_dialog.spinbox_cores_requested.value()
        mpirun = 'mpirun -np %d' % cores if self.parent.dmp_enabled() else ''
        return text.replace("${JOB_NAME}", self.parent.run_dialog.lineedit_job_name.text()) \
                   .replace("${CORES}", str(cores)) \
                   .replace("${QUEUE}", self.parent.run_dialog.combobox_queue_name.currentText()) \
                   .replace("${MODULES}", self.parent.run_dialog.lineedit_queue_modules.text()) \
                   .replace("${MPIRUN}", mpirun) \
                   .replace("${COMMAND}", ' '.join(cmd))

    def submit_command(self, cmd, port):
        with open(os.path.join(os.path.dirname(__file__), 'run_hpcee')) as qsub_template:
            template_text = qsub_template.read()

        qsub_script = tempfile.NamedTemporaryFile()
        qsub_script.write(self.transform_template(template_text, cmd))
        qsub_script.flush() # Keep tmpfile open

        cwd = os.getcwd()
        os.chdir(self.parent.get_project_dir())
        os.system('qsub %s' % qsub_script.name)
        os.chdir(cwd)

        if port: # port is not None iff using pymfix
            self.connect('http://%s:%s' % (socket.gethostname(), port))

        qsub_script.close() # deletes tmpfile

    def start_command(self, cmd, cwd, env):
        """Start MFIX in QProcess"""

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

        cmdline = ' '.join(cmd) # cmd is a list

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
            self.mfix_pid = self.mfixproc.pid() # Keep a copy because it gets reset
            msg = "MFIX process %d is running" % self.mfix_pid
            self.parent.update_run_options_signal.emit(msg)
            log.debug("Full MFIX startup parameters: %s", cmdline)
            QTimer.singleShot(1000, self.connect)

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
            self.mfix_pid = None
            #self.parent.stdout_signal.emit("MFIX (pid %s) has stopped" % self.mfixproc.pid())
            self.mfixproc = None
            self.disconnect()
            self.parent.update_run_options_signal.emit(msg)
        self.mfixproc.finished.connect(slot_finish)

        def slot_error(error):
            self.status.clear()
            if error == QProcess.FailedToStart:
                msg = "Process failed to start "+cmdline
            elif error == QProcess.Crashed:
                msg = "Process exit "+cmdline
            elif error == QProcess.Timedout:
                msg = "Process timeout "+cmdline
            elif error in (QProcess.WriteError, QProcess.ReadError):
                msg = "Process communication error "+cmdline
            else:
                msg = "Unknown error "+cmdline
            log.warn(msg)
            self.mfixproc = None
            self.parent.stderr_signal.emit(msg) # make the message print in red
            self.parent.update_run_options_signal.emit('')

        self.mfixproc.error.connect(slot_error)

        self.mfixproc.start(cmd[0], cmd[1:])

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.api.post(arg, data=b"")

    def terminate_pymfix(self):
        """ Request clean exit from MFIX """
        self.api.post('exit', b"timeout=1")

    def update_status(self):
        """update the status of  the pymfix monitor"""
        if not self.api.connected:
            return
        self.api.get('status')
