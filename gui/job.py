"""class to manage external MFIX process"""

import json
import logging
import os
import pprint
import signal
import tempfile
import uuid

from subprocess import Popen

log = logging.getLogger(__name__)

from qtpy.QtCore import QObject, QTimer, QUrl, Signal, pyqtSignal
from qtpy.QtNetwork import  QNetworkAccessManager, QNetworkReply, QNetworkRequest

from tools.general import make_callback
from tools.general import get_mfix_home

SUPPORTED_PYMFIXPID_FIELDS = ['url', 'pid', 'token', 'qjobid']

def get_dict_from_pidfile(pid_filename):

    try:
        pid_dict = {}
        with open(pid_filename) as pidfile:
            log.debug('opened pid file %s', os.path.basename(pid_filename))
            for line in pidfile.readlines():
                try:
                    key, value = line.strip().split('=')
                    if key in SUPPORTED_PYMFIXPID_FIELDS:
                        pid_dict[key] = value
                        log.debug('PIDFILE %s = %s' % (key, value))
                except ValueError:
                    continue
            return pid_dict
    except (FileNotFoundError, IOError, OSError):
        log.exception('PID could not be opened: %s', pid_filename)
    return {}


class PymfixAPI(QNetworkAccessManager):
    """ Class to extend QNetworkAccessManager with pymfix connection details
        set transparently.
    """

    def __init__(self, pidfile, response_handler, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        log.debug('New PymfixAPI object %s', self)
        self.pidfile = pidfile
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}
        self.pid_contents = get_dict_from_pidfile(self.pidfile)
        for k, v in self.pid_contents.items():
            setattr(self, k, v)
        self.requests = set()
        self.default_response_handler = response_handler
        self.api_available = False
        log.debug('API object created for job %s', self.pidfile)

    def test_connection(self, response, error):
        # use request/reply methods to call 'status', Job will provide the
        # signal handler
        log.debug('test_connection called %s', self)
        if not self.pid_contents.get('url'):
            log.error('pidfile does not contain API url (test_connection)')
            raise Exception
        self.get('status', handlers={'response': response, 'error': error})

    def api_request(self, method, endpoint, data=None, handlers={}):
        url = self.pid_contents.get('url')
        token = self.pid_contents.get('token')
        req = QNetworkRequest(QUrl('%s/%s' % (url, endpoint)))
        method = str(method).lower()
        method = self.methods[method]
        request_id = str(uuid.uuid4())
        self.requests.add(request_id)
        for param in (('url', url),
                      ('endpoint', endpoint),
                      ('token', token),
                      ('data', data),
                      ('method', method)):
            log.debug("API request %s: %s=%s" % \
              (request_id, param[0], param[1]))
        # authentication header
        if token:
            (name, value) = token.split(':')
            req.setRawHeader(bytearray(name, 'utf-8'), bytearray(value, 'utf-8'))
        request_object = method(req) if data == None else method(req, data)
        # connect response and error handlers
        response_handler = handlers.get('response', self.default_response_handler)
        request_object.finished.connect(
            lambda slot=self.slot_api_response, ri=request_id, ro=request_object, rh=response_handler: slot(ri, ro, rh))
        request_object.error.connect(lambda network_error: log.error("network error: %s" % network_error))
        return request_id

    def slot_api_response(self, request_id, response_object, signal):
        """Process QNetworkReply content then send a signal containing response
        data.
        :arg response_object: API call response object
        :type response_object: QtNetwork.QNetworkReply
        :arg signal: Signal to emit
        :type signal: QtCore.pyqtSignal
        """
        if response_object.error != 0 and False:
            log.debug('API request error, no more processing of this request %s', request_id)
            response_object.close()
            self.requests.discard(request_id)
            return
        try:
            response = response_object.readAll()
            response_json = response.data().decode('utf-8')
            json.loads(response_json)
        except (TypeError, ValueError) as e:
            log.debug("API response parsing error")
            log.debug('response headers: %s' % response_object.rawHeaderList())
            error_code = response_object.error()
            error_desc = "API response could not be parsed"
            response_json = json.dumps(
                              {"mfix_status": {},
                               "command_output": "",
                               "internal_api_error": {
                                 "error_code": error_code,
                                 "error_desc": error_desc,
                                 "raw_api_message":  str(response.data())}})
            log.debug("processed response:\n%s", json.loads(response_json))
        finally:
            self.requests.discard(request_id)
            response_object.close()
            signal.emit(request_id, response_json)

    def slot_protocol_error(self, request_id, response_object):
        log.debug('API protocol error %s', request_id)
        try:
            response_error_code = response_object.error()
            if response_error_code == 1:
                log.error('API network connection error %s', self)
            else:
                log.debug('Unknown network error %s', self)
        except:
            log.exception('unhandled API error %s', self)
        finally:
            try:
                response_data = response_object.readAll()
            except:
                response_data = None
            error_desc = response_data.errorString() if response_data else ""
            raw_api_response = response_data.data() if response_data else ""
            error_code = response_object.error()
            response_string = {"mfix_status": {},
                               "command_output": "",
                               "internal_api_error": {
                                  "error_code": error_code,
                                  "error_desc": error_desc,
                                  "raw_api_response": raw_api_response}}
            try:
                response_object.writeData(json.dumps(response_string))
                #RuntimeError: no access to protected functions or signals for objects not created from Python
            except:
                pass

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


class JobManager(QObject):

    """class for managing and monitoring an MFIX job"""

    # TODO: state detection:

        #   no pid file -> state is new run |
        #   pid file exists -> (
        #       contains queue info -> (
        #           qstat state is running -> (
        #               API url in pid file -> (
        #                   API connection works -> state is running remote, may be paused |
        #                   API connection fails -> state is stale pid, enqueued remote and job error
        #               ) |
        #               API url is not in pid file -> state is enqueued, waiting for API
        #           ) |
        #           qstat state is error -> state is stale pid, queue error |
        #           qstat state is job finished -> state is stale pid, job error, mfix didn't write to pid
        #       ) |
        #       no queue info -> (
        #           contains API url -> (
        #               API connection works - state is running locally (may be paused) |
        #               API connection fails - state is stale pid |
        #           ) |
        #           no API url -> should not happen
        #       )
        #   )

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent):
        log.debug('New JobManager %s', self)
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        self.api_error_count = 0
        self.API_ERROR_SOFT_LIMIT = 5 # API request failures before pid check
        self.API_ERROR_HARD_LIMIT = 8 # trigger pid check

    def try_to_connect(self, pidfile):
        if self.job:
            log.debug('JobManager reusing %s', self.job)
        elif os.path.isfile(pidfile):
            self.reset_api_error_count()
            self.job = Job(pidfile)
            log.debug('JobManager created %s', self.job)
            self.job.sig_job_exit.connect(self.teardown_job)

            # TODO? handle with signal so this can be managed in gui
            self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)
            self.job.sig_update_job_status.connect(self.sig_update_job_status)
            self.job.sig_change_run_state.connect(self.sig_change_run_state)
            self.job.sig_api_error.connect(self.increment_api_error_count)
            self.job.sig_api_success.connect(self.reset_api_error_count)

    def mfix_proc_is_alive(self):
        """Handles process existence checks"""
        # TODO: handle queued job polling
        pid_contents = get_dict_from_pidfile(self.job.pidfile)
        # if local process
        pid = pid_contents.get('pid')
        if pid:
            try:
                return os.kill(int(pid), 0)
            except:
                pass
        return False
        # else check queue process (TODO)

    def reset_api_error_count(self):
        self.api_error_count = 0
        log.debug('Reset API error count %s', self)

    def increment_api_error_count(self, message=None):
        """Increment API error count and signal job exit if limit
        has been exceeded"""
        self.api_error_count += 1
        count = self.api_error_count
        log.debug('API error count incremented: %s', self.api_error_count)
        if count >= self.API_ERROR_SOFT_LIMIT:
            # check that the mfix process is still running
            if self.mfix_proc_is_alive() and count < self.API_ERROR_HARD_LIMIT:
                log.error('MFIX process is unresponsive, retry %s (of %s)' %\
                    ( self.api_error_count - self.API_ERROR_SOFT_LIMIT,
                      self.API_ERROR_HARD_LIMIT - self.API_ERROR_SOFT_LIMIT))
                return
            else:
                log.error('MFIX process has died or retry timeout reached')
                self.teardown_job()

    def record(self,job_id):
        self.job_id = job_id
        self.parent.project.saveKeyword('job_id', job_id)
        self.parent.save_project()

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
        self.job.stop_mfix()

    def teardown_job(self):
        """Job ended or exited. Destory Job object and remove pidfile"""
        log.debug('teardown_job')
        if not self.job:
            log.error('Job is not running')
            return
        pidfile = self.job.pidfile
        pid_contents = get_dict_from_pidfile(pidfile)
        pid = pid_contents.get('pid', None)
        job_id = pid_contents.get('qjobid', None)
        # needed?
        #open(os.path.join(self.parent.get_project_dir(), 'MFIX.STOP'), 'w').close()

        log.info('Job ended or connection to running job failed %s', self)
        # update GUI
        if self.job.api_test_timer and self.job.api_test_timer.isActive():
            self.job.api_test_timer.stop()
        if self.job.api_status_timer and self.job.api_status_timer.isActive():
            self.job.api_status_timer.stop()
        self.job = None
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()

        def force_stop():
            log.debug('force stop')
            # if mfix exited cleanly, there should be no pidfile
            if pid and not job_id:
                try:
                    os.kill(int(pid), 0)
                except:
                    log.debug("MFIX process %s does not exist", pid)
                    return
                try:
                    log.debug('MFIX process %s is running after exit request', pid)
                    os.kill(int(pid), signal.SIGKILL)
                except: pass
            if os.path.exists(pidfile):
                try:
                    os.remove(pidfile)
                    log.debug('removed %s', pidfile)
                except OSError:
                     log.debug("could not remove %s", pidfile)
        QTimer.singleShot(3000, force_stop)


class Job(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_job_exit = Signal()
    sig_api_response = pyqtSignal(str, str)
    sig_api_error = Signal()
    sig_api_success = Signal()

    sig_handle_api_test = pyqtSignal(str, str)
    sig_handle_api_test_error = pyqtSignal(str, QObject)

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, pidfile):

        log.debug('New JobManager.Job %s', self)
        super(Job, self).__init__()
        self.status = {}
        self.pretty_status = ""
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.api_available = False
        self.mfix_pid = None
        self.sig_api_response.connect(self.slot_api_response)

        self.sig_handle_api_test.connect(self.slot_handle_api_test)
        self.sig_handle_api_test_error.connect(self.slot_handle_api_test_error)

        self.api = PymfixAPI(self.pidfile,
                        self.sig_api_response,
                        ignore_ssl_errors=True)

        # test API before starting status loop
        # Timer will be stopped by slot_handle_api_test
        self.api_test_timer = QTimer()
        self.api_test_timer.setInterval(1000)
        self.api_test_timer.timeout.connect(self.test_api_connection)
        self.api_test_timer.start()
        # don't start status timer here, it is started once API is responsive
        self.api_status_timer = QTimer()
        self.api_status_timer.setInterval(1000)
        self.api_status_timer.timeout.connect(self.update_job_status)

        log.debug('testing API %s' % self.api)
        self.test_api_connection()

    def cleanup_and_exit(self):
        log.debug('cleanup_and_exit')
        if self.api_test_timer and self.api_test_timer.isActive():
            self.api_test_timer.stop()
        if self.api_status_timer and self.api_status_timer.isActive():
            self.api_status_timer.stop()
        self.api_available = False
        self.sig_job_exit.emit()

    def test_api_connection(self):
        """Test API connection.
        Start API test timer if it is not already running. Check current test
        count and signal job exit if timeout has been exceeded."""
        log.debug('api_test_connection called %s', self)
        try:
            self.api.test_connection(
              self.sig_handle_api_test,
              self.sig_handle_api_test_error)
        except Exception:
            log.exception('API connection test failed (test_api_connection)')
            self.sig_api_error.emit()

    def slot_handle_api_test_error(self, request_id, response_object):
        # this can be hit multiple times: until JobManager error limit is reached
        log.debug('API network error - req id %s' % request_id)
        self.sig_api_error.emit()

    def slot_handle_api_test(self, request_id, response_string):
        """Parse response data from API test call.
        If API works, stop API test timer, start API status timer.
        If API response includes a permenant failure or is not well formed,
        stop test timer and signal job exit."""
        # this can be hit multiple times: until JobManager error limit is reached
        log.debug('RETURN FROM API (slot_handle_api_test)')
        try:
            response_json = json.loads(response_string)
            if response_json.get('internal_api_error'):
                log.debug('API response error (slot_handle_api_test):')
                log.debug("API error: %s" % json.dumps(response_json))
                self.sig_api_error.emit()
                return
        except:
            # response was not parsable, API is not functional
            log.debug(response_string)
            log.exception('API response format error (slot_handle_api_test)')
            self.sig_api_error.emit()
            return

        # reset error count and mark API as available
        self.api_available = True
        # stop test timer and start status timer
        self.api_test_timer.stop()
        self.api_status_timer.start()
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()
        log.debug("API test complete, stopped test timer, started status timer (slot_handle_api_test)")

    def register_request(self, request_id, handler):
        """bind handler to request_id"""
        registered_requests = self.requests.get(request_id, set([]))
        registered_requests.add(handler)
        self.requests[request_id] = registered_requests

    def slot_api_response(self, request_id, response_string):
        """Evaluate self.requests[request_id] and call registered handlers"""
        log.debug('processing API response %s' % request_id)
        try:
            response_json = json.loads(response_string)
        except:
            log.error("API format error")
            self.sig_api_error.emit()
            return
        if response_json.get('internal_api_error'):
            log.error("API format error")
            log.debug("API response data: %s" % \
              json.dumps(response_json))
            self.sig_api_error.emit()
            #self.status.clear()
            return
        self.sig_api_success.emit()
        handlers = self.requests.get(request_id, set())
        log.debug('Registered handlers: %s' % ', '.join([str(h) for h in handlers]))
        for slot in handlers:
            slot(request_id=request_id, response_string=response_string)
        try:
            self.requests.pop(request_id)
        except KeyError as e:
            # do we care?
            pass

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def pause(self):
        req_id = self.api.put('pause')
        self.register_request(req_id, self.handle_pause)
        #req_id = self.api.get('foo')
        #self.register_request(req_id, self.handle_status)

    def handle_pause(self, request_id, response_string):
        self.handle_status(request_id, response_string)
        self.sig_change_run_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, request_id, response_string):
        self.handle_status(request_id, response_string)
        self.sig_change_run_state.emit()

    def update_job_status(self):
        req_id = self.api.get('status')
        self.register_request(req_id, self.handle_status)

    def handle_status(self, request_id, response_string):
        log.debug('Job:handle_status for %s', request_id)
        response_json = json.loads(response_string)
        self.status = json.loads(response_json.get('mfix_status'))
        if not self.status:
            log.debug("API status response was empty")
            self.status = {}
        self.pretty_status = pprint.PrettyPrinter(indent=4,
            width=50).pformat(self.status)
        log.debug('status:\n%s' % self.pretty_status)
        self.sig_update_job_status.emit()
        # remove once state handlers process API response directly:
        self.sig_change_run_state.emit()

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.register_request(self.api.post(arg, data=b''), self.handle_set_api_output)

    def handle_set_api_output(self, request_id, response_string):
        self.handle_status(request_id, response_string)
        self.sig_change_run_state.emit()

    def stop_mfix(self):
        """Send stop request"""
        req_id = self.api.post('exit')
        self.register_request(req_id, self.handle_stop_mfix)

    def handle_stop_mfix(self, request_id, response_string):
        log.debug('handle_stop_mfix')
        self.api_available = False
        self.api_status_timer.stop()
        self.handle_status(request_id, response_string)
        self.sig_change_run_state.emit()
        self.cleanup_and_exit()
