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
from qtpy.QtNetwork import (QNetworkAccessManager,
                            QNetworkRequest,
                            QNetworkReply)


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
       sets the MFiX connection details transparently in API calls. This
       class is used by instances of :mod:`job.Job`
    """

    def __init__(self, parent, default_response_handler, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        self.parent = parent
        self.warning = parent.warning
        self.error = parent.error
        self.pidfile = parent.pidfile
        self.process_info = get_process_info(self.pidfile)
        self.mfix_pid = self.process_info.get('pid')
        self.requests = set()
        self.default_response_handler = default_response_handler
        self.api_available = False
        self.ignore_ssl_errors = ignore_ssl_errors
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}


    def test_connection(self, response, error=None):
        # use request/reply methods to call 'status', Job will provide the
        # signal handler
        if not self.process_info.get('url'):
            self.error('pidfile does not contain URL')
            raise Exception
        self.get('status', handlers={'response': response, 'error': error})


    def request(self, method, endpoint, data=None, headers=None, handlers=None):
        """This method abstracts HTTP requests. A unique request ID will be
        returned to the caller. Callbacks are set on the request to response
        and error handlers within this class.

        The "handlers" dictionary may contain one or both of the keys 'error'
        and 'response'. These must be `QtCore.Signal` objects, and the internal
        handlers in this class will emit these signals after local response
        processing.
        """
        #print("REQUEST", method)
        if not headers:
            headers = {}
        if not handlers:
            handlers = {}
        url = self.process_info.get('url')
        token = self.process_info.get('token')
        req = QNetworkRequest(QUrl('%s/%s' % (url, endpoint)))
        method = self.methods[str(method).lower()]
        req_id = str(uuid.uuid4())
        self.requests.add(req_id)
        # assemble headers
        if 'content-type' not in [x.lower() for x in headers.keys() if type(x) == str]:
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
            request_object = method(req, data)
        else:
            request_object = method(req)
        response_handler = handlers.get('response', self.default_response_handler)
        error_handler = handlers.get('error')
        request_object.finished.connect(
          partial(
            self.slot_response, req_id, request_object, response_handler))
        request_object.error.connect(
          partial(
            self.slot_protocol_error, req_id, request_object, error_handler))
        return req_id


    def slot_response(self, req_id, response_object, signal):
        """Process QNetworkReply content then emit a signal containing response
        content to the signal object provided. If reponse content is valid JSON
        the string content is emitted via the specified signal. If the response
        is not parsable as JSON, an error JSON string is assembled and emitted
        via the specified signal.
        """
        #print("SLOT_RESPONSE", req_id)
        try:
            response = response_object.readAll()
            if not response:  # Nothing read. We should not have gotten here
                response_json = None # for "finally" clause
                return
            response_json = response.data().decode('utf-8')
            json.loads(response_json)
        except (TypeError, ValueError, json.decoder.JSONDecodeError) as e:
            error_code = response_object.error()
            error_desc = str(e)
            response_json = json.dumps(
                              {"mfix_status": {},
                               "command_output": "",
                               "internal_error": {
                                 "error_code": error_code,
                                 "error_desc": error_desc,
                                 "raw_message":  str(response.data())}})
        finally:
            self.requests.discard(req_id)
            if signal and response_json:
                signal.emit(req_id, response_json)
            response_object.deleteLater()

    def slot_protocol_error(self, req_id, response_object, signal):
        """Process response errors then emit a signal containing the
        request id and `QtNetwork.QNetworkReply` object to the signal object provided.
        """
        #print("SLOT_PROTOCOL_ERROR", req_id)
        try:
            error_code = response_object.error()
            error_string = response_object.errorString()
            self.error(error_string)
            if error_code == QNetworkReply.ConnectionRefusedError:
                self.parent.parent.increment_error_count() # XXX
        finally:
            if signal:
                signal.emit(req_id, response_object)
            response_object.deleteLater()


    def slot_ssl_error(self, reply):
        #print("SLOT_SSL_ERROR", reply)
        if self.ignore_ssl_errors:
            reply.ignoreSsl()
        else:
            raise RuntimeError("SSL connection error")


    def get(self, endpoint, handlers=None):
        """HTTP GET"""
        #print("GET", endpoint)
        req_id = self.request(
                    'get', endpoint, handlers=handlers)
        return req_id

    def put(self, endpoint, data=b'', handlers=None):
        """HTTP PUT"""
        #print("PUT", endpoint, data)
        req_id = self.request(
                    'put', endpoint, data=data, handlers=handlers)
        return req_id

    def post(self, endpoint, data=b'', headers=None, handlers=None):
        """HTTP POST"""
        #print("POST", endpoint, data)
        req_id = self.request(
                    'post', endpoint, data=data, headers=headers, handlers=handlers)
        return req_id


class JobManager(QObject):

    """class for managing and monitoring MFiX jobs"""

    # TODO: state detection:

        #   no pid file -> state is new run |
        #   pid file exists -> (
        #       contains queue info -> (
        #           qstat state is running -> (
        #               url in pid file -> (
        #                   http connection works -> state is running remote, may be paused |
        #                   http connection fails -> state is stale pid, enqueued remote and job error
        #               ) |
        #               url is not in pid file -> state is enqueued, waiting for API
        #           ) |
        #           qstat state is error -> state is stale pid, queue error |
        #           qstat state is job finished -> state is stale pid, job error, mfix didn't write to pid
        #       ) |
        #       no queue info -> (
        #           contains url -> (
        #               http connection works - state is running locally (may be paused) |
        #               http connection fails - state is stale pid |
        #           ) |
        #           no url -> should not happen
        #       )
        #   )

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent):
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        self.pidfile = None
        self.warning = parent.warning
        self.error = parent.error
        self.error_count = 0
        self.ERROR_SOFT_LIMIT = 3 # API request failures before pid check
        self.ERROR_HARD_LIMIT = 6 # trigger pid check


    def cleanup_stale_pidfile(self, pidfile):
        """Delete pidfile if it refers to a job that is no longer running.
        Return True if file was stale/deleted"""
        #  TODO extend for queued jobs
        try:
            process_info = get_process_info(pidfile)
        except Exception as e:
            self.error("Removing file %s due to exception %s" % (pidfile, e))
            return True
        pid = process_info.get('pid')
        job_id = process_info.get('qjobid')
        if job_id is not None: # TODO check for remote job
            return False
        if pid is not None:
            try:
                os.kill(int(pid), 0) # Succeeds if process is running
                return False
            except Exception as e:
                self.error('MFiX process %s: %s' % (pid, e))
        try:
            self.warning("Removing stale file %s" % pidfile)
            os.unlink(pidfile)
        except Exception as e:
            self.error('%s: %s' % (pidfile, e))
        return True


    def try_to_connect(self, pidfile):
        """Create Job object if one does not yet exist."""
        if self.job: # Already connected
            pass
        elif os.path.isfile(pidfile):
            self.reset_error_count()
            if self.cleanup_stale_pidfile(pidfile):
                return
            self.job = Job(self, pidfile)
            self.pidfile = pidfile
            self.job.connect()
            # connect Job signals
            self.job.sig_job_exit.connect(self.teardown_job)
            self.job.sig_update_job_status.connect(self.sig_update_job_status)
            self.job.sig_change_run_state.connect(self.sig_change_run_state)
            self.job.sig_error.connect(self.increment_error_count)
            self.job.sig_success.connect(self.reset_error_count)
            # TODO? handle with signal so this can be managed in gui


    def mfix_proc_is_alive(self):
        """Handles process existence checks"""
        # TODO: handle queued job polling
        if not self.job:
            return False
        try:
            process_info = get_process_info(self.pidfile)
        except Exception as e:
            self.error("%s: %s" % (self.pidfile, e))
            return False
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

    def reset_error_count(self):
        self.error_count = 0

    def increment_error_count(self):
        """Increment API error count. Signal job exit if limit has been exceeded.        """
        self.error_count += 1
        count = self.error_count
        if count >= self.ERROR_SOFT_LIMIT:
            if count == self.ERROR_SOFT_LIMIT:
                self.error("Soft error limit reached")
            if count <= self.ERROR_HARD_LIMIT:
                self.error('MFiX process is unresponsive, retry %s (of %s)' %
                    (self.error_count, self.ERROR_HARD_LIMIT))
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

    def submit_command(self,
                       qscript,
                       submit_cmd,
                       delete_cmd, # XXX UNUSED
                       status_cmd,
                       job_id_regex,
                       replace_dict):

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

    def force_stop_mfix(self):
        if self.job is not None:
            self.job.force_stop_mfix()


    def teardown_job(self):
        """Job ended or exited. Destroy Job object and remove pidfile"""
        # TODO: this handles local only, not queued jobs
        if not self.job:
            self.warning('Job is not running')
        else:
            self.warning('Job ended or connection to running job failed')
        pidfile = self.pidfile
        if pidfile:
            try:
                process_info = get_process_info(pidfile)
                pid = process_info.get('pid')
                job_id = process_info.get('qjobid')

            except Exception as e:
                self.error("Cannot get PID: %s" % e)
                pid = job_id = None

            finally:
                # update GUI
                if self.job:
                    self.job.stop_timers() # cleanup_and_exit recurses
                self.job = None
                self.sig_update_job_status.emit()
                self.sig_change_run_state.emit()


        def force_stop(): # XXX UNUSED
            # if mfix exited cleanly, there should be no pidfile
            if not os.path.isfile(pidfile):
                return
            if pid and not job_id:
                try:
                    os.kill(int(pid), 0)
                except:
                    self.error('MFiX process %s does not exist' % pid)
                try:
                    self.error('MFiX process %s is still running' % pid)
                    os.kill(int(pid), signal.SIGKILL)
                except Exception as e:
                    self.error("MFIX process %s: %s" % (pid, e))
                    return

                if os.path.exists(pidfile):
                    try:
                        os.remove(pidfile)
                    except OSError as e:
                        self.error("could not remove %s: %s" % (pidfile, e))


class Job(QObject):
    """Class for managing and monitoring an MFiX job. This class contains
    methods for issuing requests to and handling responses from the MFiX API.

    pidfile: Name of file containing MFiX API connection information.
    """

    # Signals
    sig_job_exit = Signal()
    sig_response = Signal(str, str)
    sig_error = Signal()
    sig_success = Signal()
    sig_handle_test = Signal(str, str)
    sig_handle_test_error = Signal(str, QObject)
    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent, pidfile):
        super(Job, self).__init__()
        self.parent = parent
        self.warning = parent.warning
        self.error = parent.error
        self.mfix_pid = None
        self.status = {}
        self.pretty_status = ""
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.api_available = False
        self.sig_response.connect(self.slot_response)

        self.sig_handle_test.connect(self.slot_handle_test)
        self.sig_handle_test_error.connect(self.slot_handle_test_error)

        self.api = PymfixAPI(self,
                        self.sig_response,
                        ignore_ssl_errors=True)

        self.mfix_pid = self.api.mfix_pid # ?

        # test API before starting status loop
        # Timer will be stopped by slot_handle_test
        self.api_test_timer = QTimer()
        self.api_test_timer.setInterval(1000)
        self.api_test_timer.timeout.connect(self.test_connection)
        self.api_test_timer.start()
        # don't start status timer here, it is started once API is responsive
        self.api_status_timer = QTimer()
        self.api_status_timer.setInterval(1000)
        self.api_status_timer.timeout.connect(self.update_job_status)

    def connect(self):
        self.test_connection()

    def stop_timers(self):
        if self.api_test_timer and self.api_test_timer.isActive():
            self.api_test_timer.stop()
        if self.api_status_timer and self.api_status_timer.isActive():
            self.api_status_timer.stop()

    def cleanup_and_exit(self):
        self.stop_timers()
        self.api_available = False
        self.sig_job_exit.emit()

    def test_connection(self):
        """Test API connection.
        Start API test timer if it is not already running. Check current test
        count and signal job exit if timeout has been exceeded."""
        try:
            self.api.get('status', handlers={'response':self.sig_handle_test})
        except Exception as e: ## FIXME, too generic
            #print("SIG_ERROR", e) #increments error count
            self.sig_error.emit()

    def slot_handle_test_error(self, req_id, response_object):
        # this can be hit multiple times: until JobManager
        # error limit is reached
        self.error('Network error: request id %s' % req_id)


    def slot_handle_test(self, req_id, response_string):
        """Parse response data from API test call. If response is in JSON
        format and does not contain an error message, then stop API test timer
        and start API status timer.

        If API response includes a permenant failure or is not well formed,
        emit :class:`Job.sig_error` signal.
        """
        # this can be hit multiple times: until we hit JobManager error limit
        try:
            response_json = json.loads(response_string)
            if response_json.get('internal_error'):
                self.error('Internal error %s' % response_string)
                self.sig_error.emit()
                return
        except Exception as e:
            # response was not parsable, API is not functional
            self.error('Response format error: %s, response=%s' %
                         (e, response_string))
            self.sig_error.emit()
            return

        # reset error count and mark API as available
        self.api_available = True
        # stop test timer and start status timer
        self.api_test_timer.stop()
        self.api_status_timer.start()
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()

    def register_request(self, req_id, handler):
        """Bind handler to req_id."""
        #print("REGISTER", req_id, handler)
        registered_requests = self.requests.get(req_id, set([]))
        registered_requests.add(handler)
        self.requests[req_id] = registered_requests

    def slot_response(self, req_id, response_string):
        #print("RESPONSE", req_id, response_string)
        try:
            response_json = json.loads(response_string)
        except Exception as e:
            self.error("Response parsing error: %s" % e)
            self.sig_error.emit()
            return
        if response_json.get('internal_error'):
            self.error("Internal error: response=%s" % response_string)
            self.sig_error.emit()
            return
        self.sig_success.emit()
        handlers = self.requests.get(req_id, set())
        for slot in handlers:
            slot(req_id=req_id, response_string=response_string)
        try:
            self.requests.pop(req_id)
        except KeyError as e:
            self.error("No handler for request %s" % req_id)

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
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, req_id, response_string):
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def update_job_status(self):
        req_id = self.api.get('status')
        self.register_request(req_id, self.handle_status)

    def handle_status(self, req_id, response_string):
        response_json = json.loads(response_string)
        self.status = json.loads(response_json.get('mfix_status'))
        if not self.status:
            self.status = {}
        self.pretty_status = pprint.PrettyPrinter(indent=4,
            width=50).pformat(self.status)
        self.sig_update_job_status.emit()
        # remove once state handlers process response directly:
        self.sig_change_run_state.emit()

    def set_pymfix_output(self, state):
        """Toggle verbose Flask messages."""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        req_id = self.api.post(arg, data=b'')
        self.register_request(req_id, self.handle_set_output)

    def handle_set_output(self, req_id, response_string):
        """Handler for responses to `set_pymfix_output` requests"""
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()

    def stop_mfix(self):
        """Send stop request"""
        req_id = self.api.post('exit')
        self.register_request(req_id, self.handle_stop_mfix)

    def force_stop_mfix(self):
        """Send stop request"""
        req_id = self.api.post('die')
        self.register_request(req_id, self.handle_stop_mfix)


    def handle_stop_mfix(self, req_id, response_string):
        """Handler for responses to stop requests"""
        self.api_available = False
        self.api_status_timer.stop()
        self.handle_status(req_id, response_string)
        self.sig_change_run_state.emit()
        self.cleanup_and_exit()
