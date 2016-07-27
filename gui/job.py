"""class to manage external MFIX process"""

import json
import logging
import os
import pprint
import signal
import socket
import tempfile
import uuid

from subprocess import Popen

DEFAULT_TIMEOUT = 2 # seconds, for all socket ops
socket.setdefaulttimeout(DEFAULT_TIMEOUT)

log = logging.getLogger(__name__)

from qtpy.QtCore import QObject, QTimer, QUrl, Signal, pyqtSignal
from qtpy.QtNetwork import  QNetworkAccessManager, QNetworkReply, QNetworkRequest

from tools.general import make_callback
from tools.general import get_mfix_home

SUPPORTED_PYMFIXPID_FIELDS = ['url', 'pid', 'token', 'qjobid']

def get_dict_from_pidfile(pid_filename):

    try:
        pid_dict = {}
        pidfile = open(pid_filename)
        log.debug('opened pid file %s', os.path.basename(pid_filename))
        for line in pidfile.readlines():
            try:
                key, value = line.strip().split('=')
                if key in SUPPORTED_PYMFIXPID_FIELDS:
                    pid_dict[key] = value
            except ValueError:
                continue
        return pid_dict
    except: # TODO: make specific - file not found
        log.debug('PID file not found: %s', pid_filename)
        return {}


class PymfixAPI(QNetworkAccessManager):
    """ Class to extend QNetworkAccessManager with pymfix connection details
        set transparently.
    """


    def __init__(self, pidfile, response_handler, error_handler, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        self.pidfile = pidfile
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}
        self.pid_contents = get_dict_from_pidfile(self.pidfile)
        self.requests = set()
        self.default_response_handler = response_handler
        self.default_error_handler = error_handler
        self.api_available = False

        log.info('API connection for job %s', self.pidfile)
        log.debug(self)

    def test_connection(self, response, error):
        # use request/reply methods to call 'status', Job will provide the
        # signal handler
        log.debug('test_connection called %s', self)
        self.get('status', handlers={'response': response, 'error': error})

    def api_request(self, method, endpoint, data=None, handlers={}):
        req = QNetworkRequest(QUrl('%s/%s' % (self.pid_contents['url'], endpoint)))
        method = str(method).lower()
        method = self.methods[method]
        request_id = str(uuid.uuid4())
        self.requests.add(request_id)
        log.debug("API request: method=%s endpoint=%s data=%s id=%s; %s" % \
          (method, endpoint, data, request_id, req))
        # authentication header
        if self.pid_contents['token']:
            (name, value) = self.pid_contents['token'].split(':')
            req.setRawHeader(name, value)
        if data is not None:
            request_object = method(req, data)
        else:
            request_object = method(req)
        # connect response and error handlers
        response_handler = handlers.get('response', self.default_response_handler)
        error_handler = handlers.get('error', self.default_error_handler)
        request_object.finished.connect(
            make_callback(
              self.slot_api_reply, request_id, request_object, response_handler))
        request_object.error.connect(
            make_callback(
              self.slot_protocol_error, request_id, request_object, error_handler))
        return request_id

    def slot_api_reply(self, request_id, reply_object, signal):
        """Process QNetworkReply content then send a signal containing reply
        data.
        :arg reply: API call reply object
        :type reply: QtNetwork.QNetworkReply
        :arg signal: Signal object to be signalled
        :type signal: QtCore.pyqtSignal
        """
        # FIXME handle empty replies or Pymfix methods that don't return JSON
        try:
            reply = reply_object.readAll()
            response_json = reply.data()
            json.loads(response_json)
        except (TypeError, ValueError) as e:

            try:
                if reply_object.NetworkError == 1:
                    error_desc = 'Network connection failed'
                    error_code = 1
                else:
                    error_desc = 'API reponse could not be parsed as JSON'
                    error_code = None
            except:
                log.exception("reply_object.Network error test")
            #log.exception(error_desc)
            response_json = json.dumps({"pymfix_api_error": {
                                          "error_code": error_code,
                                          "error_desc": error_desc,
                                          "raw_api_message":  reply.data()}})
            log.debug(response_json)
            log.debug('raw headers: %s' % reply_object.rawHeaderList())
        finally:
            self.requests.discard(request_id)
            signal.emit(request_id, response_json)
            reply_object.close()

    def slot_protocol_error(self, request_id, reply_object, handler):
        log.debug('API protocol error')
        handler.emit(request_id, reply_object)

    def slot_ssl_error(self, reply):
        """ Handler for SSL connection errors. Check self.ignore_ssl_errors
            and continue as appropriate"""
        if self.api.ignore_ssl_errors:
            log.debug('call to %s:%s completed with ignored SSL errors' % \
                  (self.runname_pid, self.endpoint))
            reply.ignoreSsl()
        else:
            log.warn('call to %s:%s aborted due to SSL errors' % \
                  (self.runname_pid, self.endpoint))
            raise Exception # find appropriate Qt exception or make this meaningful

    def get(self, endpoint, handlers={}):
        request_id = self.api_request('get', endpoint, data=None, handlers=handlers)
        return request_id

    def put(self, endpoint, data=b'', handlers={}):
        request_id = self.api_request('put', endpoint, data=data, handlers=handlers)
        return request_id

    def post(self, endpoint, data=b'', handlers={}):
        request_id = self.api_request('post', endpoint, data=data, handlers=handlers)
        return request_id

    def get_job_status(self):
        # return (
        #       state=([ local | on_queue ], [ running | paused | error ]),
        #       queue=[(qid, qstate) | None ]
        #        )
        # use scheduler commands (qstat, etc)
        # TODO: abstract scheduler interface into QueueManager or somesuch
        #       to support grid environment plugins
        pass



class JobManager(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent):
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        # self.jobs = [] # in Justin's parametric branch?

    def try_to_connect(self, pidfile):
        if self.job:
            log.debug('try_to_connect called and we already have a Job object')
        elif os.path.isfile(pidfile):
            self.job = Job(pidfile)
            self.parent.signal_update_runbuttons.emit('')
            self.job.sig_job_exit.connect(self.slot_teardown_job)

            # TODO? handle with signal so this can be managed in gui
            self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)
            self.job.sig_update_job_status.connect(self.parent.slot_update_residuals)
            self.job.sig_update_job_status.connect(self.parent.slot_update_runbuttons)
            self.job.sig_change_run_state.connect(self.sig_change_run_state)
        return bool(self.job)

    def record(self,job_id):
        self.job_id = job_id
        self.parent.project.saveKeyword('job_id', job_id)
        self.parent.save_project()

    def slot_teardown_job(self):
        log.info('Job ended or connection to running job failed %s', self)
        log.debug('removing %s', self.job.pidfile)
        os.remove(self.job.api.pidfile)
        # update GUI
        self.job = None
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()

    def is_job_pending(self):
        return self.job and not self.job.api_available

    def submit_command(self, cmd, dmp_enabled, smp_enabled):

        with open(os.path.join(get_mfix_home(), 'gui', 'run_hpcee')) as qsub_template:
            template_text = qsub_template.read()

        cores = self.parent.run_dialog.spinbox_cores_requested.value()
        threads = self.parent.run_dialog.spinbox_threads.value()
        do_smp = 'env OMP_NUM_THREADS=%d' % threads if smp_enabled else ''
        mpirun = 'mpirun -np %d' % cores if dmp_enabled else ''
        qscript = template_text.replace("${JOB_NAME}", self.parent.run_dialog.lineedit_job_name.text()) \
                .replace("${CORES}", str(cores)) \
                .replace("${QUEUE}", self.parent.run_dialog.combobox_queue_name.currentText()) \
                .replace("${MODULES}", self.parent.run_dialog.lineedit_queue_modules.text()) \
                .replace("${MPIRUN}", mpirun) \
                .replace("${COMMAND}", ' '.join(cmd))

        qsub_script = tempfile.NamedTemporaryFile()
        print(qscript)
        qsub_script.write(qscript)
        qsub_script.flush() # Keep tmpfile open

        # FIXME: for testing without qsub installed, use:
        # Popen('qsub %s' % qsub_script.name, cwd=self.parent.get_project_dir())
        print('/bin/csh %s &' % qsub_script.name, self.parent.get_project_dir())
        proc = Popen('qsub %s' % qsub_script.name, shell=True, cwd=self.parent.get_project_dir())
        out, err = proc.communicate()
        job_id = out.split(' ')[2] # qsub stdout:  You job JOB_ID has been submitted

        qsub_script.close() # deletes tmpfile
        self.parent.job_manager.save_job_id(job_id)

    def stop_mfix(self):
        """Send stop request"""
        log.debug('stop_mfix')
        req_id = self.job.stop_mfix()
        open(os.path.join(self.parent.get_project_dir(), 'MFIX.STOP'), 'w').close()
        def force_stop():
            log.debug('force stop')
            pid = get_dict_from_pidfile(self.job.pidfile).get('pid', None)
            job_id = get_dict_from_pidfile(self.job.pidfile).get('qjobid', None)
            if pid and not job_id:
                os.kill(int(pid), signal.SIGKILL)
        QTimer.singleShot(1000, force_stop)
        #self.register_request(req_id, self.handle_stop)

class Job(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_job_exit = Signal()
    sig_read_status = pyqtSignal(str, str)
    sig_api_response = pyqtSignal(str, str)
    sig_api_error = pyqtSignal(str, QObject)

    sig_handle_api_test = pyqtSignal(str, str)
    sig_handle_api_test_error = pyqtSignal(str, QObject)

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()
    API_TIMEOUT = 3 # seconds, based on API timer interval

    def __init__(self, pidfile):

        super(Job, self).__init__()
        self.status = {}
        self.cached_status = ''
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.api_available = False
        self.test_api_timeout = 0
        self.mfix_pid = None
        self.sig_api_response.connect(self.slot_api_response)
        self.sig_api_error.connect(self.slot_api_error)

        self.sig_handle_api_test.connect(self.slot_handle_api_test)
        self.sig_handle_api_test_error.connect(self.slot_handle_api_test_error)

        self.api = PymfixAPI(self.pidfile,
                        self.sig_api_response,
                        self.sig_api_error,
                        ignore_ssl_errors=True)

        # TODO: more complete API testing at instantiation
        log.debug('testing API %s' % self.api)
        # test API before starting status loop
        self.test_api_timer = QTimer()
        self.test_api_timer.setInterval(1000)
        # sig_handle_api_test will fire slot_handle_api_test - API liveness check
        # Timer will be stopped by slot_handle_api_test
        self.test_api_connection()

    def test_api_connection(self):
        """Test API connection.
        Start API test timer if it is not already running. Check current test
        count and signal job exit if timeout has been exceeded."""
        log.debug('api_test_connection called: %s', self)
        if self.test_api_timeout > self.API_TIMEOUT:
            log.debug('API test timeout exceeded (test_api_connection)')
            self.api_test_timer.stop()
            self.sig_job_exit.emit()
            return
        else:
            # increment API test timeout
            self.test_api_timeout += 1
        self.api.test_connection(self.sig_handle_api_test,self.sig_handle_api_test_error)

    def slot_handle_api_test_error(self, request_id, response_object):
        # this can be hit multiple times: until self.api_test_timeout is reached
        # (incremented in self.test_api_connection)
        log.debug('API network error (slot_handle_api_test_error)')
        log.debug(response_object)
        #self.sig_job_exit.emit()

    def slot_handle_api_test(self, request_id, response_data):
        """Parse response data from API test call.
        If API works, stop API test timer, start API status timer.
        If API response includes a permenant failure or is not well formed,
        stop test timer and signal job exit."""

        # this can be hit multiple times: until self.api_test_timeout is reached
        # (incremented in self.test_api_connection)

        log.debug('RETURN FROM API (slot_handle_api_test)')
        try:
            response_json = json.loads(response_data)
            if response_json.get('pymfix_api_error'):
                log.debug('API response error (slot_handle_api_test):')
                log.debug("API error: %s" % json.dumps(response_json))
                #self.sig_job_exit.emit()
        except:
            # response was not parsable, API is not functional
            log.debug('API response format error (slot_handle_api_test)')
            #self.sig_job_exit.emit()
            return

        # API is available, stop test timer and start status timer
        log.debug("API test complete, stopping test timer (slot_handle_api_test)")
        #self.api_test_timer.stop()
        self.api_available = True
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()
        self.api_timer = QTimer()
        self.api_timer.setInterval(1000)
        self.api_timer.timeout.connect(self.update_job_status)
        self.api_timer.start()

    def register_request(self, request_id, handler):
        """bind handler to request_id"""
        registered_requests = self.requests.get(request_id, set([]))
        registered_requests.add(handler)
        self.requests[request_id] = registered_requests

    def slot_api_response(self, request_id, response_data=None):
        """Evaluate self.requests[request_id] and call registered handlers"""
        log.debug('processing API response for request %s' % request_id)
        handlers = self.requests.get(request_id, set())
        log.debug('Registered handlers: %s' % ', '.join([str(h) for h in handlers]))
        for slot in handlers:
            slot(response_data=response_data)
        try:
            self.requests.pop(request_id)
        except KeyError as e:
            # do we care?
            pass

    def slot_api_error(self, request_id, response_object):
        log.error("API error signal in Job, %s" % self)
        try:
            response_error_code = response_object.NetworkError
            if response_error_code == 1:
                log.debug('API network connection error (slot_api_error)')
            else:
                log.debug('Unknown network error, calling teardown (slot_api_error)')
        except:
            log.exception('unhandled API error (slot_api_error)')
        finally:
            # unconditionally exit; TODO: retry logic
            self.sig_job_exit.emit()

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def update_job_status(self):
        log.debug('Job:update_job_status()')
        req_id = self.api.get_job_status()
        self.register_request(req_id, self.sig_get_job_state)

    def pause(self):
        req_id = self.api.put('pause')
        self.register_request(req_id, self.handle_pause)

    def handle_pause(self, response_data=None):
        self.update_job_status()
        self.sig_change_run_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, response_data=None):
        self.update_job_status()
        self.sig_change_run_state.emit()

    def update_job_status(self):
        req_id = self.api.get('status')
        self.register_request(req_id, self.handle_status)

    def handle_status(self, response_data=None):
        response_json = json.loads(response_data)
        pretty_status = pprint.PrettyPrinter(indent=4,
            width=50).pformat(response_json)
        log.debug(pretty_status)
        if response_json.get('pymfix_api_error'):
            log.error("API error: %s" % json.dumps(response_json))
            self.status.clear()
        self.status = response_json
        self.cached_status = pretty_status
        self.sig_update_job_status.emit()

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.register_request(self.api.post(arg, data=b''), self.handle_output)

    def stop_mfix(self):
        """Send stop request"""
        self.api.post('exit')
        self.api_timer.stop()

    def slot_get_job_state(self, available):
        if available:
            self.mfix_pid = self.api.pid
            log.info('Job API available')
            self.sig_change_run_state.emit()
        else:
            log.error('Job API unavailable')
        # (this was in update_job_status. Refactor to here)
        # only relevant for local jobs
        #pid = int(self.api.pymfix['pid'])
        #try:
        #    os.kill(pid, 0)
        #except OSError:
        #    log.info("job with pid %d is no longer running", pid)
        #    self.disconnect()
        #    self.sig_change_run_state.emit()
        #    return
