# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import errno
import json
import logging
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

def get_process_info(filename):
    """Read contents of provided MFiX job pid file and set supported
    key-value pairs as dictionary members
    """
    d = {}
    try:
        with open(filename) as f:
            for line in f:
                tok = line.strip().split('=', 1)
                if len(tok) == 2 and tok[0] in SUPPORTED_PYMFIXPID_FIELDS:
                    k, v = tok
                    d[k] = v
                else:
                    raise ValueError(line)

    except OSError as e:
        if e.errno != errno.ENOENT: # Return empty dict if no file
            raise
    return d


class PymfixAPI(QNetworkAccessManager):
    """This class is used by :class:`job.Job` and it should not be necessary to
       use directly.

       The PymfixAPI class extends :mod:`QtNetwork.QNetworkAccessManager` and
       sets the MFiX API connection details transparently in API calls. This
       class is used by instances of :mod:`job.Job`
    """

    def __init__(self, parent, def_response_handler, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        self.parent = parent
        self.warning = parent.warning
        self.error = parent.error
        self.pidfile = parent.pidfile
        self.process_info = get_process_info(self.pidfile)
        self.mfix_pid = self.process_info.get('pid')
        self.requests = set()
        self.def_response_handler = def_response_handler
        self.api_available = False
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}

    def test_connection(self, response, error=None):
        # use request/reply methods to call 'status', Job will provide the
        # signal handler
        if not self.process_info.get('url'):
            self.error('pidfile does not contain API url (test_connection)')
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
        url = self.process_info.get('url')
        token = self.process_info.get('token')
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
            req.setRawHeader(kb,vb)

        if data is not None:
            request_object = api_method(req, data)
        else:
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
        try:
            response = response_object.readAll()
            response_json = response.data().decode('utf-8')
            json.loads(response_json)
        except (TypeError, ValueError) as e:
            error_code = response_object.error()
            error_desc = "API response could not be parsed"
            response_json = json.dumps(
                              {"mfix_status": {},
                               "command_output": "",
                               "internal_api_error": {
                                 "error_code": error_code,
                                 "error_desc": error_desc,
                                 "raw_api_message":  str(response.data())}})
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
        try:
            response_error_code = response_object.error()
            if response_error_code == 1:
                self.warning("API network connection error %s" % self)
            else:
                self.warning("Unknown network error %s" % self)
        except Exception as e:
            self.error("Unhandled API error %s" % e)
        finally:
            response_object.deleteLater()

    def slot_ssl_error(self, reply):
        """Handler for SSL connection errors. Check self.ignore_ssl_errors
           and continue as appropriate"""
        if self.api.ignore_ssl_errors:
            reply.ignoreSsl()
        else:
            self.warning('call to %s:%s aborted due to SSL errors' %
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
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        self.warning = parent.warning
        self.error = parent.error
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
            pass
        elif os.path.isfile(pidfile):
            self.reset_api_error_count()
            self.job = Job(self, pidfile)
            self.job.connect()

            # connect Job signals
            self.job.sig_job_exit.connect(self.teardown_job)
            self.job.sig_update_job_status.connect(self.sig_update_job_status)
            self.job.sig_change_run_state.connect(self.sig_change_run_state)
            self.job.sig_api_error.connect(self.increment_api_error_count)
            self.job.sig_api_success.connect(self.reset_api_error_count)
            # TODO? handle with signal so this can be managed in gui
        else:
            self.warning("Can't start new job")

    def mfix_proc_is_alive(self):
        """Handles process existence checks"""
        # TODO: handle queued job polling
        if not self.job:
            return False
        process_info = get_process_info(self.job.pidfile)
        # if local process
        pid = process_info.get('pid')
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

    def increment_api_error_count(self):
        """Increment API error count. Signal job exit if limit
        :class:`JobManager.API_ERROR_HARD_LIMIT` has been exceeded.
        """
        self.api_error_count += 1
        count = self.api_error_count
        if count >= self.API_ERROR_SOFT_LIMIT:
            if count == self.API_ERROR_SOFT_LIMIT:
                self.error("Soft error limit reached")
            if count < self.API_ERROR_HARD_LIMIT:
                self.error('MFiX process is unresponsive, retry %s (of %s)' %
                    (self.api_error_count, self.API_ERROR_HARD_LIMIT))
                return
            if self.mfix_proc_is_alive():
               self.error("MFiX process is running and unresponsive. Giving up.")
            self.error('MFiX process has died or retry timeout reached')
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
                     cwd=self.parent.get_project_dir(), env=dict(os.environ, LD_PRELOAD=""))
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
        if not self.job:
            self.error('Job is not running')
            return
        self.warning('Job ended or connection to running job failed')
        pidfile = self.job.pidfile
        try:
            process_info = get_process_info(pidfile)
            pid = process_info.get('pid', None)
            job_id = process_info.get('qjobid', None)

        except Exception as e:
            self.error("Cannot get PID: %s" % e)
            pid = job_id = None

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
            # if mfix exited cleanly, there should be no pidfile
            if not os.path.isfile(pidfile):
                return
            if pid and not job_id:
                try:
                    os.kill(int(pid), 0)
                except:
                    self.warning('MFiX process %s does not exist' % pid)
                    return
                try:
                    self.warning('MFiX process %s is still running' % pid)
                    os.kill(int(pid), signal.SIGKILL)
                except Exception as e:
                    self.warning("MFIX process %s: %s" % (pid, e))

            if os.path.exists(pidfile):
                try:
                    os.remove(pidfile)
                except OSError as e:
                    self.warning("could not remove %s: %s" % (pidfile, e))
        # NOTE: this timeout needs to exceed the retry limit
        QTimer.singleShot(1000, force_stop) # why?


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

    def __init__(self, parent, pidfile):
        super(Job, self).__init__()
        self.warning = parent.warning
        self.error = parent.error
        self.mfix_pid = None
        self.status = {}
        self.pretty_status = ""
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.api_available = False
        self.sig_api_response.connect(self.slot_api_response)

        self.sig_handle_api_test.connect(self.slot_handle_api_test)
        self.sig_handle_api_test_error.connect(self.slot_handle_api_test_error)

        self.api = PymfixAPI(self,
                        self.sig_api_response,
                        ignore_ssl_errors=True)

        self.mfix_pid = self.api.mfix_pid # ?

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
        self.test_api_connection()

    def cleanup_and_exit(self):
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
        try:
            self.api.get('status', handlers={'response':self.sig_handle_api_test})
        except Exception as e:
            self.error("test_api_connection: %s" % e)
            self.sig_api_error.emit()

    def slot_handle_api_test_error(self, req_id, response_object):
        # this can be hit multiple times: until JobManager
        # error limit is reached
        self.warning('API network error - req id %s' % req_id)

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
        try:
            response_json = json.loads(response_string)
            if response_json.get('internal_api_error'):
                self.warning('Internal error %s' % response_string)
                self.sig_api_error.emit()
                return
        except Exception as e:
            # response was not parsable, API is not functional
            self.warning('API response format error: %s, response=%s' %
                         (e, response_string))
            self.sig_api_error.emit()
            return

        # reset error count and mark API as available
        self.api_available = True
        # stop test timer and start status timer
        self.api_test_timer.stop()
        self.api_status_timer.start()
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
        try:
            response_json = json.loads(response_string)
        except Exception as e:
            self.warning("API response parsing error: %s" % e)
            self.sig_api_error.emit()
            return
        if response_json.get('internal_api_error'):
            self.warning("Internal error: response=%s" % response_string)
            self.sig_api_error.emit()
            return
        self.sig_api_success.emit()
        handlers = self.requests.get(req_id, set())
        for slot in handlers:
            slot(req_id=req_id, response_string=response_string)
        try:
            self.requests.pop(req_id)
        except KeyError as e:
            self.warning("No known handler for request %s" % req_id)

    def reinit(self, project_file_contents):
        """reinitialize job. Sanity checks (paused, gui.unsaved_flag, etc)
        must be handled by caller"""
        req_id = self.api.post(
                    'reinitialize',
                    data=json.dumps({"project_file": project_file_contents}),
                    headers={'content-type': 'application/json'})
        self.register_request(req_id, self.handle_reinit)

    def handle_reinit(self, req_id, response_string):
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
        response_json = json.loads(response_string)
        self.status = json.loads(response_json.get('mfix_status'))
        if not self.status:
            self.status = {}
        self.pretty_status = pprint.PrettyPrinter(indent=4,
            width=50).pformat(self.status)
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
        self.api_available = False
        self.api_status_timer.stop()
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()
        self.cleanup_and_exit()
