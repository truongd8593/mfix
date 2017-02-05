"""
job.py
======

This module defines classes used to facilitate MFiX job management.

:platform: Unix, Windows
:license: Public Domain
"""

import json
import logging
import os
import pprint
import signal
import uuid
from functools import partial
import re

from subprocess import Popen, PIPE

log = logging.getLogger(__name__)

from qtpy.QtCore import QObject, QTimer, QUrl, Signal
from qtpy.QtNetwork import QNetworkAccessManager,  QNetworkRequest


from mfixgui.tools.general import replace_with_dict

#: List of valid keys to read from PID file
SUPPORTED_PYMFIXPID_FIELDS = ['url', 'pid', 'token', 'qjobid']

def get_dict_from_pidfile(pid_filename):
    """Read contents of provided MFiX job pid file and set supported
    key-value pairs as dictionary members

    :param pid_file: Name of file that contains MFiX process connection
                     information. This must be an absolute filename.
    """

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
    except Exception:
        log.debug('PID could not be opened: %s', pid_filename)
    return {}


class PymfixAPI(QNetworkAccessManager):
    """This class is used by :class:`job.Job` and it should not be necessary to
       use directly.

       The PymfixAPI class extends :mod:`QtNetwork.QNetworkAccessManager` and
       sets the MFiX API connection details transparently in API calls. This
       class is used by instances of :mod:`job.Job`
    """

    def __init__(self, pidfile, def_response_handler, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        log.debug('New PymfixAPI object %s', self)
        self.pidfile = pidfile
        self.pid_contents = get_dict_from_pidfile(self.pidfile)
        for k, v in self.pid_contents.items():
            setattr(self, k, v)
        self.requests = set()
        self.def_response_handler = def_response_handler
        self.api_available = False
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}
        log.debug('API object created for job %s', self.pidfile)

    def test_connection(self, response, error=None):
        # use request/reply methods to call 'status', Job will provide the
        # signal handler
        log.debug('test_connection called %s', self)
        if not self.pid_contents.get('url'):
            log.error('pidfile does not contain API url (test_connection)')
            raise Exception
        self.get('status', handlers={'response': response, 'error': error})

    def api_request(self, method, endpoint, data=None, headers=None, handlers=None):
        """This method abstracts HTTP requests. A unique request ID will be
        returned to the caller. Callbacks are set on the request to response
        and error handlers within this class.

        :param method: HTTP verb to be used in API request
        :type method: str
        :param endpoint: API endpoint (URL portion after the domain)
        :type endpoint: str
        :param data: Data to be included in request (only used in PUT and POST)
        :type data: str
        :param headers: HTTP headers to include in request
        :type headers: dict
        :param handlers: Response and error handling signals
        :type handlers: dict
        :return: Request ID
        :rtype: string

        The "handlers" dictionary must contain one or both of the keys 'error'
        and 'response'. These must be `QtCore.Signal` objects, and the internal
        handlers in this class will emit these signals after local response
        processing.
        """
        if not headers:
            headers = {}
        if not handlers:
            handlers = {}
        url = self.pid_contents.get('url')
        token = self.pid_contents.get('token')
        req = QNetworkRequest(QUrl('%s/%s' % (url, endpoint)))
        api_method = self.methods[str(method).lower()]
        req_id = str(uuid.uuid4())
        self.requests.add(req_id)
        # assemble headers
        if 'content-type' not in [x.lower for x in headers.keys() if type(x) == str]:
            headers['content-type'] = 'text/plain'
        headers['x-pymfix-requestid'] = req_id
        headers['content-length'] = str(len(data if data else ''))
        if token:
            (name, value) = token.split(':')
            headers[name] = value
        for k,v in headers.items():
            kb = k.encode('utf-8')
            vb = v.encode('utf-8')
            log.debug("setting header %s to %s" % (kb,vb))
            req.setRawHeader(kb,vb)
        for param in (('request id', req_id),
                      ('url', url),
                      ('endpoint', endpoint),
                      ('token', token),
                      ('data', data),
                      ('request handers', handlers)):
            log.debug("API %s request: %s=%s" % \
              (str(method).upper(), param[0], param[1]))
        if data is not None:
            log.debug('"data" arg in request %s is type %s' % (req_id, type(data)))
            log.debug("calling api_method %s: %s %s" % (api_method, req, data))
            request_object = api_method(req, data)
        else:
            log.debug("calling api_method %s", api_method)
            request_object = api_method(req)
        response_handler = handlers.get('response', self.def_response_handler)
        error_handler = handlers.get('error')
        request_object.finished.connect(
          partial(
            self.slot_api_response, req_id, request_object, response_handler))
        request_object.error.connect(
          partial(
            self.slot_protocol_error, req_id, request_object, error_handler))
        return req_id

    def slot_api_response(self, req_id, response_object, signal):
        """Process QNetworkReply content then emit a signal containing response
        content to the signal object provided. If reponse content is valid JSON
        the string content is emitted via the specified signal. If the response
        is not parsable as JSON, an error JSON string is assembled and emitted
        via the specified signal.

        :param response_object: API call response object
        :type response_object: QtNetwork.QNetworkReply
        :param signal: Signal to emit
        :type signal: QtCore.Signal
        """
        log.debug("processing %s for %s" % (req_id, signal))
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
            self.requests.discard(req_id)
            signal.emit(req_id, response_json)
            response_object.deleteLater()

    def slot_protocol_error(self, req_id, response_object, signal):
        """Process API response errors then emit a signal containing the
        request id and `QtNetwork.QNetworkReply` object to the signal object provided.

        :param req_id: API request ID
        :param response_object: API call response object
        :type response_object: QtNetwork.QNetworkReply
        :param signal: Signal to emit
        :type signal: QtCore.Signal
        """
        log.debug('API protocol error %s', req_id)
        try:
            response_error_code = response_object.error()
            if response_error_code == 1:
                log.warning('API network connection error %s', self)
            else:
                log.warning('Unknown network error %s', self)
        except:
            log.exception('unhandled API error %s', self)
        finally:
            response_object.deleteLater()

    def slot_ssl_error(self, reply):
        """Handler for SSL connection errors. Check self.ignore_ssl_errors
           and continue as appropriate"""
        if self.api.ignore_ssl_errors:
            log.debug('call to %s:%s completed with ignored SSL errors' % \
                  (self.runname_pid, self.endpoint))
            reply.ignoreSsl()
        else:
            log.warn('call to %s:%s aborted due to SSL errors' % \
                  (self.runname_pid, self.endpoint))
            # find appropriate Qt exception or make this meaningful
            raise Exception

    def get(self, endpoint, handlers=None):
        """API request via HTTP GET

        :param endpoint: API endpoint (URL suffix spit at and excluding domain)
        :param handlers: Dictionary containing ::QtCore.Signal:: objects to
                         connect to API response and error signals"""
        req_id = self.api_request(
                    'get', endpoint, handlers=handlers)
        return req_id

    def put(self, endpoint, data=b'', handlers=None):
        """API request via HTTP PUT

        :param endpoint: API endpoint (URL suffix spit at and excluding domain)
        :param data: Data to include in request body
        :param handlers: Dictionary containing ::QtCore.Signal:: objects to
                         connect to API response and error signals"""
        req_id = self.api_request(
                    'put', endpoint, data=data, handlers=handlers)
        return req_id

    def post(self, endpoint, data=b'', headers=None, handlers=None):
        """API request via HTTP POST

        :param endpoint: API endpoint (URL suffix spit at and excluding domain)
        :param data: Data to include in request body
        :param handlers: Dictionary containing ::QtCore.Signal:: objects to
                         connect to API response and error signals"""
        req_id = self.api_request(
                    'post', endpoint, data=data, headers=headers, handlers=handlers)
        return req_id


class JobManager(QObject):

    """class for managing and monitoring MFiX jobs"""

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

    #: `QtCore.Signal` to emit when job status changes
    sig_update_job_status = Signal()
    #: `QtCore.Signal` to emit when running job state changes
    sig_change_run_state = Signal()

    def __init__(self, parent):
        log.debug('New JobManager %s', self)
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        self.api_error_count = 0
        self.API_ERROR_SOFT_LIMIT = 3 # API request failures before pid check
        self.API_ERROR_HARD_LIMIT = 6 # trigger pid check

    def try_to_connect(self, pidfile):
        """Create Job object if one does not yet exist.

        :param pidfile: Name of file containing MFiX process and connection
                        details. The file name must be prefixed with an
                        absolute path.
        """

        if self.job:
            log.debug('JobManager reusing %s', self.job)
        elif os.path.isfile(pidfile):
            self.reset_api_error_count()
            self.job = Job(pidfile)
            self.job.connect()
            log.debug('JobManager created %s', self.job)

            # connect Job signals
            self.job.sig_job_exit.connect(self.teardown_job)
            self.job.sig_update_job_status.connect(self.sig_update_job_status)
            self.job.sig_change_run_state.connect(self.sig_change_run_state)
            self.job.sig_api_error.connect(self.increment_api_error_count)
            self.job.sig_api_success.connect(self.reset_api_error_count)
            # TODO? handle with signal so this can be managed in gui
        else:
            log.debug("pidfile is not available, can't start new job")

    def mfix_proc_is_alive(self):
        """Handles process existence checks"""
        # TODO: handle queued job polling
        if not self.job:
            return False
        pid_contents = get_dict_from_pidfile(self.job.pidfile)
        # if local process
        pid = pid_contents.get('pid')
        if pid:
            try:
                os.kill(int(pid), 0)
                return True
            except OSError:
                return False
        else:
            return False
        # else check queue process (TODO)

    def reset_api_error_count(self):
        """Set API connection error count to 0"""
        self.api_error_count = 0
        log.debug('Reset API error count %s', self)

    def increment_api_error_count(self):
        """Increment API error count. Signal job exit if limit
        :class:`JobManager.API_ERROR_HARD_LIMIT` has been exceeded.
        """
        self.api_error_count += 1
        count = self.api_error_count
        log.debug('API error count incremented: %s', self.api_error_count)
        if count >= self.API_ERROR_SOFT_LIMIT:
            if count == self.API_ERROR_SOFT_LIMIT:
                log.error("Soft error limit reached")
            if count < self.API_ERROR_HARD_LIMIT:
                log.error('MFiX process is unresponsive, retry %s (of %s)' %\
                    (self.api_error_count, self.API_ERROR_HARD_LIMIT))
                return
            if self.mfix_proc_is_alive():
               log.error("MFiX process is running and unresponsive. Giving up.")
            log.error('MFiX process has died or retry timeout reached')
            self.teardown_job()

    def record(self,job_id):
        self.job_id = job_id
        self.parent.project.saveKeyword('job_id', job_id)
        self.parent.save_project()

    def is_job_pending(self):
        return self.job and not self.job.api_available

    def is_job_ready(self):
        return self.job and self.job.api_available

    def submit_command(self, qscript, submit_cmd, delete_cmd, status_cmd, job_id_regex, replace_dict):

        project_dir = self.parent.get_project_dir()
        script_name = '.qsubmit_script'

        # write the script
        with open(os.path.join(project_dir, script_name), 'w') as f:
            f.write(qscript)

        replace_dict['SCRIPT'] = script_name

        submit_cmd = replace_with_dict(submit_cmd, replace_dict)

        # submit the job
        self.parent.print_internal("Job submit CMD: {}".format(submit_cmd),
                                   color='blue')

        proc = Popen(submit_cmd, shell=True, stdout=PIPE, stderr=PIPE,
                     cwd=self.parent.get_project_dir())
        out, err = proc.communicate()
        if job_id_regex is not None and out:
            job_id = re.findall(job_id_regex, out)
        else:
            job_id = []
        if job_id:
            job_id = job_id[0]
            self.parent.print_internal("Job successfully submitted with job id: {}".format(job_id),
                                       color='blue')
        else:
            self.parent.error('Could not determine job id')
            job_id = None
        if err:
            self.parent.error('Error with submission:\n{}'.format(err))

        #TODO:
        # use status_cmd to see if it queued, running, etc.

    def stop_mfix(self):
        if self.job is not None:
            self.job.stop_mfix()

    def teardown_job(self):
        """Job ended or exited. Destroy Job object and remove pidfile"""
        # TODO: this handles local only, not queued jobs
        log.debug('teardown_job')
        if not self.job:
            log.error('Job is not running')
            return
        log.info('Job ended or connection to running job failed %s', self)
        pidfile = self.job.pidfile
        try:
            pid_contents = get_dict_from_pidfile(pidfile)
            pid = pid_contents.get('pid', None)
            job_id = pid_contents.get('qjobid', None)
        except Exception:
            log.debug('Could not clean up job %s', pidfile)
            log.debug('pidfile was already removed, no pid available')
        finally:
            # update GUI
            if (self.job.api_test_timer
                and self.job.api_test_timer.isActive()):
                self.job.api_test_timer.stop()
            if (self.job.api_status_timer
                and self.job.api_status_timer.isActive()):
                self.job.api_status_timer.stop()
            self.job = None
            self.sig_update_job_status.emit()
            self.sig_change_run_state.emit()

        def force_stop():
            log.debug('force stop')
            # if mfix exited cleanly, there should be no pidfile
            if not os.path.isfile(pidfile):
                return
            if pid and not job_id:
                try:
                    os.kill(int(pid), 0)
                except:
                    log.debug("MFiX process %s does not exist", pid)
                    return
                try:
                    log.debug('MFiX process %s is still running', pid)
                    os.kill(int(pid), signal.SIGKILL)
                except: pass
            if os.path.exists(pidfile):
                try:
                    os.remove(pidfile)
                    log.debug('removed %s', pidfile)
                except OSError:
                     log.debug("could not remove %s", pidfile)
        # NOTE: this timeout needs to exceed the retry limit
        QTimer.singleShot(1000, force_stop)


class Job(QObject):

    """Class for managing and monitoring an MFiX job. This class contains
    methods for issuing requests to and handling responses from the MFiX API.

    :param pidfile: Name of file containing MFiX API connection information.
                    This string must contain the absolute path to the file.
    """

    #: Signal to emit at job exit. This is used for clean exit and fatal errors
    sig_job_exit = Signal()
    #: Signal to be bound to a local response handler
    sig_api_response = Signal(str, str)
    #: Signal to be emitted when an API error is encountered
    sig_api_error = Signal()
    #: Signal to be emitted after an API response is successfully parsed
    sig_api_success = Signal()

    #: Signal to be bound to API test response handler
    sig_handle_api_test = Signal(str, str)
    #: Signal to be bound to API test error handler
    sig_handle_api_test_error = Signal(str, QObject)

    #: Signal emitted when job status content has changed
    sig_update_job_status = Signal()
    #: Signal emitted when running job state has changed
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

        self.mfix_pid = self.api.pid

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

    def connect(self):
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
            self.api.get('status', handlers={'response':self.sig_handle_api_test})
        except Exception:
            log.exception('API connection test failed (test_api_connection)')
            self.sig_api_error.emit()

    def slot_handle_api_test_error(self, req_id, response_object):
        # this can be hit multiple times: until JobManager
        # error limit is reached
        log.debug('API network error - req id %s' % req_id)

    def slot_handle_api_test(self, req_id, response_string):
        """Parse response data from API test call. If response is in JSON
        format and does not contain an error message, then stop API test timer
        and start API status timer.

        If API response includes a permenant failure or is not well formed,
        emit :class:`Job.sig_api_error` signal.

        :param req_id: Request ID obtained from API call
        :param response_string: API response in JSON string format
        """
        # this can be hit multiple times: until we hit JobManager error limit
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
        log.debug("API test complete")
        self.api_test_timer.stop()
        log.debug("Stopped test timer")
        self.api_status_timer.start()
        log.debug("Started status timer")
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()

    def register_request(self, req_id, handler):
        """Bind handler to req_id.

        :param req_id: Request ID obtained from API call
        :param handler: :class:`Job` method to call when a response is recieved
        """
        registered_requests = self.requests.get(req_id, set([]))
        registered_requests.add(handler)
        self.requests[req_id] = registered_requests

    def slot_api_response(self, req_id, response_string):
        """Handler connected to :class:`Job.sig_api_response.

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        try:
            response_json = json.loads(response_string)
        except:
            log.warning("API response parsing error")
            self.sig_api_error.emit()
            return
        if response_json.get('internal_api_error'):
            log.warning("API format error")
            log.debug("API response data: %s" % \
              json.dumps(response_json))
            self.sig_api_error.emit()
            return
        self.sig_api_success.emit()
        handlers = self.requests.get(req_id, set())
        log.debug('Registered handlers: %s',
                     ', '.join([str(h) for h in handlers]))
        for slot in handlers:
            slot(req_id=req_id, response_string=response_string)
        try:
            self.requests.pop(req_id)
        except KeyError as e:
            log.debug("No known handler for request %s" % req_id)

    def reinit(self, project_file_contents):
        """reinitialize job. Sanity checks (paused, gui.unsaved_flag, etc)
        must be handled by caller"""
        req_id = self.api.post(
                    'reinitialize',
                    data=json.dumps({"project_file": project_file_contents}),
                    headers={'content-type': 'application/json'})
        self.register_request(req_id, self.handle_reinit)
        log.debug('reinitialize request made %s' % req_id)

    def handle_reinit(self, req_id, response_string):
        log.debug("reinitialize response recieved %s: %s" % \
                    (req_id, response_string))
        self.handle_status(req_id, response_string)

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def pause(self):
        req_id = self.api.put('pause')
        self.register_request(req_id, self.handle_pause)

    def handle_pause(self, req_id, response_string):
        """Handler for responses to `pause` API requests

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, req_id, response_string):
        """Handler for responses to `unpause` API requests

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def update_job_status(self):
        req_id = self.api.get('status')
        self.register_request(req_id, self.handle_status)

    def handle_status(self, req_id, response_string):
        """Handler for responses to `update_job_status`

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        log.debug('Job:handle_status for %s', req_id)
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
        """Toggle verbose Flask messages.

        :param state: State to set Flask's logging. This must be one of
                      'enable' or 'disable'."""

        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        req_id = self.api.post(arg, data=b'')
        self.register_request(req_id, self.handle_set_api_output)

    def handle_set_api_output(self, req_id, response_string):
        """Handler for responses to `set_pymfix_output` API requests

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def stop_mfix(self):
        """Send stop request"""
        req_id = self.api.post('exit')
        self.register_request(req_id, self.handle_stop_mfix)

    def handle_stop_mfix(self, req_id, response_string):
        """Handler for responses to `stop_mfix` API requests

        :param req_id: Request ID obtained from API call
        :param response_string: API call response in JSON string format
        """
        log.debug('processing API response %s' % req_id)
        log.debug('handle_stop_mfix')
        self.api_available = False
        self.api_status_timer.stop()
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()
        self.cleanup_and_exit()
