"""class to manage external MFIX process"""

from functools import wraps
import json
import logging
import os
import pprint
import socket
import sys
import time
import traceback
import uuid

DEFAULT_TIMEOUT = 2 # seconds, for all socket ops
socket.setdefaulttimeout(DEFAULT_TIMEOUT)

log = logging.getLogger(__name__)

from qtpy.QtCore import QObject, QTimer, QUrl, Signal, pyqtSignal
from qtpy.QtNetwork import  QNetworkAccessManager, QNetworkReply, QNetworkRequest

from tools.general import make_callback

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
        log.debug('PID file not found: %s' % pid_filename)
        return {}


class PymfixAPI(QNetworkAccessManager):
    """ Class to extend QNetworkAccessManager with pymfix connection details
        set transparently.
    """

    sig_api_available = pyqtSignal(str)
    # sig_api_available_public = pyqtSignal(bool)
    # sig_api_call_finished = pyqtSignal(Signal)

    def __init__(self, pidfile, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        # self.project_dir = project_dir
        # self.project_name = project_name
        # self.pidfile = os.path.join(project_dir, project_name + '.pid')
        self.pidfile = pidfile
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}
        self.pid_contents = get_dict_from_pidfile(self.pidfile)
        self.requests = set()

        # test API calls
        log.info('connecting to API for job %s' % self.pidfile)
        self.sig_api_available.connect(self.slot_api_available)
        # if self.pid_contents.get('url'):
        #     self.sig_api_available_public.emit(False)

    def slot_api_available(self, response):
        log.info('slot_api_available')
        try:
            response = json.loads(response)
            log.info('API available: %s' % self.pid_contents['url'])
            self.sig_api_available_public.emit(True)
            self.api_available = True
        except Exception:
            log.exception('API connection failed: %s' % self.pid_contents['url'])
            self.api_available = False
            self.sig_api_available_public.emit(False)

    def api_request(self, method, endpoint, data):
        log.warning("API request: method=%s endpoint=%s data=%s" % (method, endpoint, data))
        request_id = str(uuid.uuid4())
        method = str(method).lower()
        method = self.methods[method]

        req = QNetworkRequest(QUrl('http://%s/%s' % (self.pid_contents['url'], endpoint)))
        if self.pid_contents['token']:
            (name, value) = self.pid_contents['token'].split(':')
            req.setRawHeader(name, value)
        if data is not None:
            request_object = method(req, data)
        else:
            request_object = method(req)
        log.debug("request sent, request_id %s: method=%s endpoint=%s" % (request_id, method, endpoint))
        self.requests.add(request_id)
        return (request_id, request_object)

    def slot_read_api_reply(self, request_id, reply_object, signal):
        """Process QNetworkReply content then send a signal containing reply
        data.
        :arg reply: API call reply object
        :type reply: QtNetwork.QNetworkReply
        :arg signal: Signal object to be signalled
        :type signal: QtCore.pyqtSignal
        """
        def emit_error(message, obj=None):
            log.error(message)
            json_msg = {'error':message}
            if obj:
                obj = str(obj)
                json_msg['raw_object'] = obj
            json_str = json.dumps(json_msg)
            signal.emit(json_str)

        # FIXME handle empty replies or Pymfix methods that don't return JSON
        # try:
        if True:
            reply = reply_object.readAll()
            response_json = reply.data() if bool(reply.data()) else json.dumps({})
            reply_object.close()
            json.loads(response_json)
        # JSON parsing errors:
        # except (TypeError, ValueError) as e:
        #     emit_error('API reponse could not be parsed as JSON')
        #     self.parent_api_signal.emit(False)
        #     log.exception('API reponse could not be parsed as JSON')
        #     response_json = json.dumps('{"raw_api_message": %s"}' % reply.data())
        # QNetwork* errors: (TODO)
        #except (QNetworkReply.NetworkError) as e:
        #    emit_error('API reponse could not be parsed as JSON')
        #    self.parent_api_signal.emit(False)

        # All other errors:
        # except Exception as e:
        #     emit_error('API response could not be read')
        #     self.parent_api_signal.emit(False)
        #     traceback.print_exception(*sys.exc_info())
        #     return
        self.requests.discard(request_id)
        signal.emit(request_id, response_json)

    def get_job_status(self):
        # return (
        #       state=([ local | on_queue ], [ running | paused | error ]),
        #       queue=[(qid, qstate) | None ]
        #        )
        # use scheduler commands (qstat, etc)
        # TODO: abstract scheduler interface into QueueManager or somesuch
        #       to support grid environment plugins
        pass

    def get(self, endpoint, handler):
        (request_id, request) = self.api_request('get', endpoint, data=None)
        if handler:
            request.finished.connect(
              make_callback(self.slot_read_api_reply,request_id,request,handler))
        return request_id

    def put(self, endpoint, handler, data=b""):
        (request_id, request) = self.api_request('put', endpoint, data)
        if handler:
            request.finished.connect(
              make_callback(self.slot_read_api_reply,request_id,request,handler))
        return request_id

    def post(self, endpoint, handler, data=b""):
        (request_id, request) = self.api_request('post', endpoint, data)
        if handler:
            request.finished.connect(
              make_callback(self.slot_read_api_reply,request_id,request,handler))
        return request_id



class JobManager(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_change_job_state = Signal()

    def __init__(self, parent):
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent
        # self.jobs = [] # in Justin's parametric branch?

    def try_to_connect(self, pidfile):
        if os.path.isfile(pidfile):
            self.job = Job(pidfile)
            self.parent.signal_update_runbuttons.emit('')
            self.job.sig_job_exit.connect(self.slot_teardown_job)

            # TODO? handle with signal so this can be managed in gui
            self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)
            self.job.sig_update_status.connect(self.parent.slot_update_residuals)

    def record(self,job_id):
        self.job_id = job_id
        self.parent.project.saveKeyword('job_id', job_id)
        self.parent.save_project()

    def slot_teardown_job(self):
        log.info('Job ended or connection to running job failed')
        self.job = None

    def submit_command(self, cmd):

        with open(os.path.join(get_mfix_home(), 'gui', 'run_hpcee')) as qsub_template:
            template_text = qsub_template.read()

        qsub_script = tempfile.NamedTemporaryFile()
        qsub_script.write(self.transform_template(template_text, cmd))
        qsub_script.flush() # Keep tmpfile open

        # FIXME: for testing without qsub installed, use:
        # Popen('qsub %s' % qsub_script.name, cwd=self.parent.get_project_dir())
        print('/bin/csh %s &' % qsub_script.name, self.parent.get_project_dir())
        proc = Popen('qsub %s' % qsub_script.name, shell=True, cwd=self.parent.get_project_dir())
        proc.wait()
        job_id = proc.output().split(' ')[1]

        qsub_script.close() # deletes tmpfile
        self.parent.job_manager.save_job_id(job_id)

class Job(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_job_exit = Signal()
    sig_read_status = pyqtSignal(str, str)
    sig_api_available = pyqtSignal(bool)
    sig_api_response = pyqtSignal(str, str)

    sig_trytohandlestatusfornow = pyqtSignal(str, str)
    sig_update_status = Signal()

    def __init__(self, pidfile):

        super(Job, self).__init__()
        self.api_timer = QTimer()
        self.api_timer.setInterval(1000)
        self.api_timer.timeout.connect(self.update_status)
        self.api_timer.start()
        self.status = {'paused': False}
        self.cached_status = ''
        self.pidfile = pidfile
        self.requests = {}
        self.api = None
        self.mfix_pid = None
        self.sig_trytohandlestatusfornow.connect(self.handle_status)

        #self.sig_read_status.connect(self.slot_read_status)
        self._connect()

    def register_request(self, request_id, handler):
        """Store bind handler to request_id"""
        registered_requests = self.requests.get(request_id, set([]))
        registered_requests.add(handler)
        self.requests[request_id] = registered_requests

    def slot_api_response(self, request_id, response_data=None):
        """Evaluate self.requests[request_id] and call registered handlers"""
        log.debug('processing API response for request %s' % request_id)
        handlers = self.requests.get(request_id, set())
        log.debug('Registered handlers: %s' % ', '.join([h for h in handlers]))
        for slot in handlers:
            slot(response_data=response_data)
        try:
            self.requests.pop(request_id)
        except KeyError as e:
            pass

    def is_paused(self):
        """indicate whether pymfix is paused"""
        log.info('is_paused: %s' % str(self.status.get('paused')))
        return self.status.get('paused')

    def update_job_state(self):
        req_id = self.api.get_job_status()
        self.register_request(req_id, self.sig_get_job_state)

    def handle_update_job_state(self, **kwargs):
        pass

    def pause(self):
        req_id = self.api.put('pause')
        self.register_request(req_id, self.handle_pause)

    def handle_pause(self, response_data=None):
        self.update_status()
        self.sig_change_job_state.emit()

    def unpause(self):
        req_id = self.api.put('unpause')
        self.register_request(req_id, self.handle_unpause)

    def handle_unpause(self, response_data=None):
        self.update_status()
        self.sig_change_job_state.emit()

    def update_status(self):
        req_id = self.api.get('status', self.sig_trytohandlestatusfornow)
        self.register_request(req_id, self.handle_status)

    def handle_status(self, response):
        try:
            self.status = json.loads(response)
            self.cached_status = pprint.PrettyPrinter(indent=4,
                                    width=50).pformat(self.status)
        except ValueError:
            self.status.clear()
            log.error("could not decode JSON: %s", response)
        except TypeError as e:
            self.status.clear()
            log.error("could not decode JSON")
        finally:
            self.sig_update_status.emit()

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.register_request(self.api.post(arg, data=b''), self.handle_output)

    def stop_mfix(self):
        """Send stop request"""
        self.register_request(self.api.put('stop'), self.handle_stop)

    def slot_get_job_state(self, available):
        # TODO? move this into PymfixAPI?
        #
        # TODO: support these states:
        #
        # state check -> (
        #
        #   no pid file -> state is new run |
        #   pid file without valid URL -> pid file is corrupt
        #
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
        #           ) |
        #       ) |
        #       no queue info -> (
        #           contains API url -> (
        #               API connection works - state is running locally (may be paused) |
        #               API connection fails - state is stale pid |
        #           no API url -> should not happen
        #       )
        #   )


        if available:
            self.mfix_pid = self.api.pid
            log.info('Job API available')
            self.sig_change_job_state.emit()
        else:
            log.error('Job API unavailable')
            #self.disconnect()
        # (this was in update_status. Refactor to here)
        # only relevant for local jobs
        #pid = int(self.api.pymfix['pid'])
        #try:
        #    os.kill(pid, 0)
        #except OSError:
        #    log.info("job with pid %d is no longer running", pid)
        #    self.disconnect()
        #    self.sig_change_job_state.emit()
        #    return

    def _connect(self):
        """Connect to existing pymfix process"""

        # refactor to be more specific
        def slot_control(reply):
            log.info('Job._connect.slot_control')
            pass
            #self.parent.signal_update_runbuttons.emit('')

        # move to API class
        def slot_ssl_error(reply):
            """ Handler for SSL connection errors. Check self.ignore_ssl_errors
                and continue as appropriate"""
            if self.api.ignore_ssl_errors:
                log.debug('call to %s:%s completed with ignored SSL errors' % \
                         (self.runname_pid, self.endpoint))
                reply.ignoreSsl()
            else:
                log.warn('call to %s:%s aborted due to SSL errors' % \
                         (self.runname_pid, self.endpoint))
                raise # find appropriate Qt exception or make this meaningful

        self.api = PymfixAPI(self.pidfile, ignore_ssl_errors=True)
        self.api.sslErrors.connect(slot_ssl_error)
        self.api.finished.connect(slot_control)
