"""class to manage external MFIX process"""

from functools import wraps
import json
import logging
import os
import pprint
import socket

DEFAULT_TIMEOUT = 2 # seconds, for all socket ops
socket.setdefaulttimeout(DEFAULT_TIMEOUT)

log = logging.getLogger(__name__)

from qtpy.QtCore import (QObject, QTimer, QUrl, Signal)
from qtpy.QtNetwork import QNetworkRequest, QNetworkAccessManager


SUPPORTED_PYMFIXPID_FIELDS = ['url', 'pid', 'token']

def get_dict_from_pidfile(pid_filename):

    pid_dict = {}

    with open(pid_filename) as pidfile:
        log.debug('opened pid file %s', os.path.basename(pid_filename))
        for line in pidfile.readlines():
            try:
                key, value = line.strip().split('=')
                if key in SUPPORTED_PYMFIXPID_FIELDS:
                    pid_dict[key] = value
            except ValueError:
                continue
    return pid_dict


class PymfixAPI(QNetworkAccessManager):
    """ Class to extend QNetworkAccessManager with pymfix connection details
        set transparently.
    """

    def __init__(self, pidfile, ignore_ssl_errors=False):

        super(PymfixAPI, self).__init__()
        self.connected = False
        self.pymfix = None
        self.pidfile = pidfile
        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}
        log.debug('connecting to API for job %s' % self.pidfile)
        self.pymfix = get_dict_from_pidfile(self.pidfile)
        # TODO: verify API is listening -- self.get('status')
        log.info('API connected: %s' % self.pymfix['url'])

    def connected(f):
        @wraps(f)
        def func(*args, **kwargs):
            return f(*args, **kwargs)
        return func

    @connected
    def api_request(self, method, endpoint, data):
        log.debug("request sent: method=%s endpoint=%s" % (method, endpoint))
        method = str(method).lower()
        method = self.methods[method]

        import requests
        response = requests.get('%s/status' % self.api.pymfix['url']).text

        req = QNetworkRequest(QUrl('http://%s/%s' % (self.pymfix['url'], endpoint)))
        if self.pymfix['token']:
            (name, value) = self.pymfix['token'].split(':')
            req.setRawHeader(name, value)
        if data is not None:
            return method(req, data)
        else:
            return method(req)

    def get(self, endpoint):
        import requests
        response = requests.get('%s/%s' % (self.pymfix['url'], endpoint))
        return response.text
        # return self.api_request('get', endpoint, data)

    def put(self, endpoint, data=b""):
        import requests
        response = requests.put('%s/%s' % (self.pymfix['url'], endpoint), data=data)
        return response.text
        # return self.api_request('put', endpoint, data)

    def post(self, endpoint, data=b""):
        import requests
        response = requests.post('%s/%s' % (self.pymfix['url'], endpoint), data=data)
        return response.text
        # return self.api_request('post', endpoint, data)

class JobManager(QObject):

    """class for managing and monitoring an MFIX job"""

    sig_update_parent = Signal()
    sig_job_exit = Signal()

    def __init__(self, runname_pid, parent):

        super(JobManager, self).__init__()

        self.parent = parent
        self.mfix_pid = None # Keep track of pid separately
        self.timer = None
        self.status = {}
        self.cached_status = ''
        self.runname_pid = runname_pid
        self.api = None

        self.sig_update_parent.connect(self.parent.slot_update_runbuttons)
        self.sig_update_parent.connect(self.parent.update_residuals)
        self._connect()

    def is_paused(self):
        """indicate whether pymfix is paused"""
        return self.status.get('paused')

    def pause(self):
        "pause pymfix job"
        self.api.put('pause')

    def unpause(self):
        "unpause pymfix job"
        self.api.put('unpause')

    def stop_mfix(self):
        """Terminate a locally running instance of mfix"""
        self.terminate_pymfix()
        self.disconnect()

        mfix_stop_file = os.path.join(self.parent.get_project_dir(), 'MFIX.STOP')
        try:
            open(mfix_stop_file, "ab").close()
        except OSError:
            pass

    def _connect(self):
        """Connect to existing pymfix process"""

        def slot_control(reply):
            self.parent.signal_update_runbuttons.emit('')

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

        self.api = PymfixAPI(pidfile=self.runname_pid, ignore_ssl_errors=True)
        self.api.sslErrors.connect(slot_ssl_error)
        self.api.finished.connect(slot_control)

        self.timer = QTimer()
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self.update_status)
        self.timer.start()
        # TODO? handle with signal so this can be managed in gui
        self.parent.ui.tabWidgetGraphics.setCurrentWidget(self.parent.ui.plot)

    def disconnect(self):
        """stop pymfix update timer"""
        if self.timer:
            self.timer.stop()
            self.timer = None
        self.api = None
        self.parent.job_manager = None
        self.parent.signal_update_runbuttons.emit('')

    def set_pymfix_output(self, state):
        """toggle Flask messages"""
        state = 'enable' if state else 'disable'
        arg = 'logging/%s' % state
        self.api.post(arg, data=b"")

    def terminate_pymfix(self):
        """ Request clean exit from MFIX """
        self.api.post('exit', {"timeout":"1"})

    def update_status(self):
        pid = int(self.api.pymfix['pid'])
        try:
            os.kill(pid, 0)
        except OSError:
            log.info("job with pid %d is no longer running", pid)
            self.disconnect()
            self.sig_update_parent.emit()
            return
        try:
            response = self.api.get('status')
            self.status = json.loads(response)
            self.cached_status = pprint.PrettyPrinter(indent=4,
                                    width=50).pformat(self.status)
        except ValueError:
                self.status.clear()
                log.error("could not decode JSON: %s", response)
        except TypeError:
                self.status.clear()
                log.error("could not decode JSON: %s", response)
        finally:
            self.sig_update_parent.emit()
