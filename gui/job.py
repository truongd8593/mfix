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
        log.debug('New PymfixAPI object %s', self)
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
            log.debug("API response parsing error")
            log.debug('raw headers: %s' % reply_object.rawHeaderList())
            log.debug(response_json)
            error_code = None
            error_desc = "API response could not be parsed"
            response_json = json.dumps(
                                {"pymfix_api_error": {
                                  "error_code": error_code,
                                  "error_desc": error_desc,
                                  "raw_api_message":  reply.data()}})
        finally:
            self.requests.discard(request_id)
            signal.emit(request_id, response_json)
            reply_object.close()

    def slot_protocol_error(self, request_id, reply_object, handler):
        log.debug('API protocol error %s', self)
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
    sig_teardown_job = Signal()

    def __init__(self, parent):
        log.debug('New JobManager %s', self)
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        # self.jobs = [] # in Justin's parametric branch?

    def try_to_connect(self, pidfile):
        if self.job:
            log.debug('try_to_connect reusing %s', self.job)
        elif os.path.isfile(pidfile):
            self.job = Job(pidfile)
            log.debug('try_to_connect created %s', self.job)
            self.parent.signal_update_runbuttons.emit('')
            self.job.sig_job_exit.connect(self.sig_teardown_job)

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
        """Job ended or exited. Destory Job object and remove pidfile"""
        if not self.job:
            return
        log.info('Job ended or connection to running job failed %s', self)
        try:
            log.debug('removing %s', self.job.pidfile)
            os.remove(self.job.pidfile)
        except OSError:
            log.debug("could not remove %s", self.job.pidfile)
        # update GUI
        if self.job.api_test_timer and self.job.api_test_timer.isActive():
            self.job.api_test_timer.stop()
        if self.job.api_status_timer and self.job.api_status_timer.isActive():
            self.job.api_status_timer.stop()
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
        if not self.job:
            log.debug('Job is not running')
        log.debug(self.job)
        open(os.path.join(self.parent.get_project_dir(), 'MFIX.STOP'), 'w').close()
        pid = get_dict_from_pidfile(self.job.pidfile).get('pid', None)
        job_id = get_dict_from_pidfile(self.job.pidfile).get('qjobid', None)
        # destroy Job object (stop timers, signal GUI to update)
        self.sig_teardown_job.emit()
        def force_stop():
            log.debug('force stop')
            if pid and not job_id:
                try:
                    os.kill(int(pid), signal.SIGKILL)
                except:
                    log.debug("MFIX process %s does not exist", pid)
        QTimer.singleShot(1000, force_stop)


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
    API_ERROR_LIMIT = 10

    def __init__(self, pidfile):

        log.debug('New JobManager.Job %s', self)
        super(Job, self).__init__()
        self.status = {}
        self.cached_status = {}
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.api_available = False
        self.api_error_count = 0 # consecutive errors allowed before job failure
        self.mfix_pid = None
        self.sig_api_response.connect(self.slot_api_response)
        self.sig_api_error.connect(self.slot_api_error)

        self.sig_handle_api_test.connect(self.slot_handle_api_test)
        self.sig_handle_api_test_error.connect(self.slot_handle_api_test_error)

        self.api = PymfixAPI(self.pidfile,
                        self.sig_api_response,
                        self.sig_api_error,
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
        if self.api_test_timer and self.api_test_timer.isActive():
            self.api_test_timer.stop()
        if self.api_status_timer and self.api_status_timer.isActive():
            self.api_status_timer.stop()
        self.sig_job_exit.emit()

    def increment_error_count(self, message=None):
        """Increment API error count and signal job exit if limit
        has been exceeded"""
        if self.api_error_count > self.API_ERROR_LIMIT:
            if message: log.debug(message)
            self.cleanup_and_exit()
            return
        log.debug('API error count incremented')
        self.api_error_count += 1

    def test_api_connection(self):
        """Test API connection.
        Start API test timer if it is not already running. Check current test
        count and signal job exit if timeout has been exceeded."""
        log.debug('api_test_connection called %s', self)
        try:
            self.api.test_connection(self.sig_handle_api_test,self.sig_handle_api_test_error)
        except Exception:
            log.debug('API connection test failed (test_api_connection)')
        finally:
            self.increment_error_count()

    def slot_handle_api_test_error(self, request_id, response_object):
        # this can be hit multiple times: until self.API_ERROR_LIMIT is reached
        log.debug('API network error (slot_handle_api_test_error)')
        self.increment_error_count()

    def slot_handle_api_test(self, request_id, response_data):
        """Parse response data from API test call.
        If API works, stop API test timer, start API status timer.
        If API response includes a permenant failure or is not well formed,
        stop test timer and signal job exit."""
        # this can be hit multiple times: until self.API_ERROR_LIMIT is reached
        log.debug('RETURN FROM API (slot_handle_api_test)')
        try:
            response_json = json.loads(response_data)
            response = json.loads(response_json.get('response', {}))
            if response.get('pymfix_api_error'):
                log.debug('API response error (slot_handle_api_test):')
                log.debug("API error: %s" % json.dumps(response))
                self.increment_error_count()
                return
        except:
            # response was not parsable, API is not functional
            log.debug(response_data)
            log.exception('API response format error (slot_handle_api_test)')
            self.increment_error_count()
            return

        # reset error count and mark API as available
        self.api_error_count = 0
        self.api_available = True
        # stop test timer and start status timer
        self.api_test_timer.stop()
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()
        self.api_status_timer.start()
        log.debug("API test complete, stopped test timer, started status timer (slot_handle_api_test)")

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
        log.error("API error signal %s" % self)
        self.api_error_count += 1
        try:
            response_error_code = response_object.NetworkError
            if response_error_code == 1:
                log.debug('API network connection error (slot_api_error)')
            else:
                log.debug('Unknown network error (slot_api_error)')
        except:
            log.exception('unhandled API error (slot_api_error)')
        finally:
            self.increment_error_count()

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def pause(self):
        req_id = self.api.put('pause')
        self.register_request(req_id, self.handle_pause)

    def handle_pause(self, response_data=None):
        self.handle_status(response_data)
        #self.update_job_status()
        self.sig_change_run_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, response_data=None):
        self.handle_status(response_data)
        #self.update_job_status()
        self.sig_change_run_state.emit()

    def update_job_status(self):
        req_id = self.api.get('status')
        self.register_request(req_id, self.handle_status)

    def handle_status(self, response_data=None):
        log.debug('response=%s, %s' % (type(response_data), response_data))
        response_json = json.loads(response_data)
        if response_json.get('pymfix_api_error'):
            log.error("API error: %s" % json.dumps(response_json))
            self.increment_error_count()
            self.status.clear()
            return
        status = json.loads(response_json.get('response', {}))
        pretty_status = pprint.PrettyPrinter(indent=4,
            width=50).pformat(status)
        log.debug(pretty_status)
        self.status = status
        self.cached_status = pretty_status
        self.sig_update_job_status.emit()
        # reset error count
        self.api_error_count = 0

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.register_request(self.api.post(arg, data=b''), self.handle_output)

    def stop_mfix(self):
        """Send stop request"""
        self.api.post('exit')
        # job cleanup deferred to JobManager

