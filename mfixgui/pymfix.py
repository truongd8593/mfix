#!/usr/bin/env python

"""The pymfix script starts mfix from Python, with a web server running for
interactive control of the run."""

__version_str__ = "2017.1" #TODO link to mfix version

import argparse
import copy
import json
import logging
import numpy.core # used?
import os
try:
    import packaging
    import packaging.requirements
    import packaging.specifiers
    import packaging.version
except ImportError:
    print("warning: can't import the module packaging")
import random
import socket
import string
import sys
import tempfile
import threading
import time
import traceback

from timeit import default_timer as timer
from functools import wraps
from flask import Flask, jsonify, make_response, render_template, request, redirect, url_for

sys.path.append(os.getcwd())
pidfilename = None

# append the path of the symlink __file__ (not its realpath)
sys.path.append(os.path.dirname(__file__))

# Fortran modules are in uppercase since Fortran uses uppercase (even though it's
# conventional to only use uppercase for constants)
from mfixsolver import compar as COMPAR
from mfixsolver import des_time_march as DES_TIME_MARCH
from mfixsolver import discretelement as DEM
from mfixsolver import iterate as ITERATE
from mfixsolver import main as MAIN
from mfixsolver import debug as DEBUG
from mfixsolver import parallel_mpi as PARALLEL_MPI
from mfixsolver import residual as RESIDUAL
from mfixsolver import run as RUN
from mfixsolver import step as STEP

PYMFIX_DIR = os.path.dirname(os.path.realpath(__file__))

FLASK_APP = Flask(__name__)
FLASK_APP.config['SECRET_KEY'] = \
            ''.join([random.choice(string.digits + string.ascii_letters)
            for x in range(0,64)])
FLASK_APP.config['TOKEN_NAME'] = 'x-pymfix-auth'
DEBUG_FLAG = False

log = logging.getLogger('werkzeug')
logformat = logging.Formatter('pymfix: %(levelname)s %(message)s')
loghandler = logging.StreamHandler()
loghandler.setFormatter(logformat)
log.setLevel(logging.INFO)
log.addHandler(loghandler)
log.disabled = True

mfix_thread = None

class WSGICopyBody(object):
    """Copy wsgi request body into environment variable for
    logging within Flask"""
    def __init__(self, application):
        self.application = application

    def __call__(self, environ, start_response):
        from cStringIO import StringIO
        length = environ.get('CONTENT_LENGTH', '0')
        length = 0 if length == '' else int(length)
        body = environ['wsgi.input'].read(length)
        meta = "WSGI INPUT\nLENGTH\t%s\nBODY\n%s" % (length, body)
        environ['wsgi_req_copy'] = meta
        environ['wsgi.input'] = StringIO(body)
        # Call the wrapped application
        app_iter = self.application(environ,
                                    self._sr_callback(start_response))
        # Return modified response
        return app_iter

    def _sr_callback(self, start_response):
        def callback(status, headers, exc_info=None):
            # Call upstream start_response
            start_response(status, headers, exc_info)
        return callback

FLASK_APP.wsgi_app = WSGICopyBody(FLASK_APP.wsgi_app)

def find_free_port():
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("",0))
    sock.listen(1)
    port = sock.getsockname()[1]
    sock.close()
    return port

def get_run_name():
    return RUN.run_name.tobytes().decode('utf-8').replace('\x00', '').strip()

def main():
    """The main function starts MFIX on a separate thread, then
       start the Flask server. """

    mfix_dat, paused, port, keyword_args = parse_command_line_arguments()

    global mfix_thread
    mfix_thread = Mfix(mfix_dat, paused, keyword_args, port=port)
    mfix_thread.start()

    def setup_ssl():
        # not too critical until remote host connections are supported
        # Use default cert/key, read commandline options, etc
        # disabled for now

        # try:
        #     import ssl
        #     sslcert = sslkey = None
        #     sslcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
        #     sslcontext.load_cert_chain(sslcert, sslkey)
        #     return True
        # except Exception as e:
        #     sslcontext = None
        #     return False
        return False

    protocol = 'https' if setup_ssl() else 'http'

    def _start_flask(host, port, debug, use_reloader):
        try:
            global pidfilename
            while not get_run_name():
                # wait for mfix thread to initialize RUN_NAME
                pass
            pidfilename = '%s.pid' % get_run_name()
            with open(pidfilename, 'w') as pid:
                pid.write('pid=%s\n' % (os.getpid(),))
                pid.write('url=%s://%s:%s\n' % (
                                    protocol, socket.gethostname(), port))
                pid.write('token=%s:%s\n' % (
                                    FLASK_APP.config['TOKEN_NAME'],
                                    FLASK_APP.config['SECRET_KEY']))
            log.debug('flask starting on port %d' % port)
            FLASK_APP.run(host=host,
                          port=port, debug=debug,
                          use_reloader=use_reloader)

        except OSError:
            os.remove(pidfilename)
            log.exception('cannot bind to port %d' % port)
            _start_flask(host=host,
                         port=find_free_port(), debug=DEBUG_FLAG,
                         use_reloader=use_reloader)

    # start the Flask server on rank 0
    if COMPAR.mype == 0:
        try:
            _start_flask(host='0.0.0.0',
                         port=port, debug=DEBUG_FLAG,
                         use_reloader=False)
        except KeyboardInterrupt:
            traceback.print_exc()
            # If we get here, the user hit Ctrl-C to shutdown the server,
            # so we call _exit() to kill the run_mfix thread.
            os.remove(pidfilename)
            os._exit(0)

    else:
        # nothing else for rank>0 to do
        mfix_thread.thread.join()


# FIXME: it would be make sense to subclass Thread
class Mfix(object):
    requests = {}
    responses = {}

    def __init__(self, mfix_dat, paused, keyword_args, port):
        self.keyword_args = keyword_args
        self.thread = None
        self.status = json.dumps({})
        self.stopped = False
        self.paused = paused
        self.mfix_dat = mfix_dat
        self.port = port
        self.t0 = 0.
        self.walltime = 0.
        self.walltime_remaining = 'unknown'
        self.time_step_init_walltime = 'unknown'
        self.time_step_end_walltime = 'unknown'
        self.do_iteration_walltime = 'unknown'
        self.des_time_init_walltime = 'unknown'
        self.des_time_steps_walltime = 'unknown'
        self.des_time_end_walltime = 'unknown'

    def start(self):
        " start the MFIX thread"
        self.thread = threading.Thread(target=self.run_mfix, kwargs={"keyword_args":self.keyword_args})
        self.thread.start()

    def run_mfix(self, keyword_args=None):
        "Main thread for running MFIX itself"

        RUN.interactive = True

        for arg in keyword_args:
            MAIN.add_command_line_keyword(arg)

        PARALLEL_MPI.parallel_init()
        DEBUG.is_pymfix = True

        # Read input data, check data, do computations for IC and BC locations
        # and flows, and set geometry parameters such as X, X_E, DToDX, etc.
        MAIN.get_data(self.mfix_dat)
        # DEBUG.good_config is not instantaneously set correctly
        time.sleep(0.25)

        while not DEBUG.good_config:
            # loaded project file can't be run, enter loop to wait for good
            # config to be uploaded
            self.paused = True
            self.check_requests()
            print('reading data from %s' % self.mfix_dat)
            MAIN.get_data(self.mfix_dat)
            time.sleep(0.1)

        MAIN.initialize(self.mfix_dat)

        self.t0 = float(RUN.time)

        self.update_status()

        if DEM.discrete_element and not DEM.des_continuum_coupled:
            DES_TIME_MARCH.des_time_init()
            if RUN.dem_solids:
                self.dem_time_march()
            if RUN.pic_solids:
                DEM.pic_time_march()
        else:

            while not self.stopped:
                self.check_requests()
                iteration_start = timer()
                self.do_step()
                self.walltime += float(timer() - iteration_start)
                if RUN.tstop <= RUN.time + 0.1*RUN.dt:
                    self.walltime_remaining = str(0.)
                    self.paused = True
                elif RUN.time > self.t0:
                    self.walltime_remaining = str(self.walltime * (RUN.tstop/(RUN.time-self.t0) - 1))

                if RUN.steady_state:
                    break

        MAIN.finalize()

    def do_step(self):
        """Run MFIX for a single timestep"""
        start = timer()
        STEP.time_step_init(self.mfix_dat)
        self.time_step_init_walltime = float(timer() - start)

        step_incomplete = True
        while step_incomplete:
            ITERATE.iterate_init()
            while ITERATE.nit < ITERATE.max_nit and not (ITERATE.converged or ITERATE.diverged):
                ITERATE.nit = ITERATE.nit + 1
                start = timer()
                ITERATE.do_iteration(self.mfix_dat)
                self.do_iteration_walltime = float(timer() - start)
                self.check_requests()

            ITERATE.post_iterate()

            step_incomplete = ITERATE.adjustdt(self.mfix_dat) and not RUN.steady_state

        STEP.check_low_dt()
        STEP.chem_mass()

        if RUN.dem_solids:
            start = timer()
            DES_TIME_MARCH.des_time_init()
            self.des_time_init_walltime = float(timer() - start)

            self.dem_time_march()

            start = timer()
            DES_TIME_MARCH.des_time_end()
            self.des_time_end_walltime = float(timer() - start)

        start = timer()
        STEP.time_step_end()
        self.time_step_end_walltime = float(timer() - start)

    def dem_time_march(self):
        """Run DEM timesteps"""
        start = timer()
        for ii in range(DES_TIME_MARCH.factor):
            print("DEM timestep %d / %d" % (ii, DES_TIME_MARCH.factor))
            DES_TIME_MARCH.des_time_step(ii)
            self.check_requests()
        self.des_time_steps_walltime = float(timer() - start)

    def update_status(self):
        """save status as JSON """

        output = {}
        output['paused'] = self.paused
        output['time'] = float(RUN.time)
        output['tstop'] = float(RUN.tstop)
        output['dt'] = float(RUN.dt)
        output['walltime_elapsed'] = self.walltime
        output['walltime_remaining'] = self.walltime_remaining
        output['profiling'] = (('time_step_init', self.time_step_init_walltime),
                               ('do_iteration', int(ITERATE.nit), self.do_iteration_walltime),
                               ('time_step_end', self.time_step_end_walltime),
                               ('des_time_init', self.des_time_init_walltime),
                               ('des_time_steps', int(DES_TIME_MARCH.factor), self.des_time_steps_walltime),
                               ('des_time_end', self.des_time_end_walltime))
        output['nit'] = int(ITERATE.nit)
        output['residuals'] = []
        if RESIDUAL.group_resid:
            for res_id in range(len(RESIDUAL.resid_grp_string)):
                output['residuals'].append((str(RESIDUAL.get_resid_grp_string(res_id)),
                                            str(RESIDUAL.get_resid_grp(res_id))))
        else:
            for res_id in range(len(RESIDUAL.resid_string)):
                output['residuals'].append((str(RESIDUAL.get_resid_string(res_id)),
                                            str(RESIDUAL.get_resid(res_id))))

        try:
            self.status = json.dumps(output)
        except UnicodeDecodeError:
            log.exception("exception when decoding:  %s" % output)

    def check_pidfile(self):
        if 0 != COMPAR.mype:
            return True

        global pidfilename
        try:
            pidfile = open(pidfilename)
            pid = pidfile.readline()[4:]
            url = pidfile.readline()[4:]
            token = pidfile.readline()[6:]
            pidfile.close()
        except IOError:
            log.exception("could not find PID file ", pidfilename)
            return False

        if int(pid) != os.getpid():
            log.error('did not find pidfile containing PID %d', os.getpid())
            return False
        return True

    def check_requests(self):
        "check for requests sent by the Flask thread"

        # exit if pidfile is missing
        self.check_pidfile()

        while True:
            if self.requests:
                # requests would only arrive at rank 0
                req_id, cmd_args = self.requests.popitem()
                if not self.check_pidfile():
                    cmd_args = ('EXIT', None)
            else:
                # command is empty for rank>0, or when rank 0 hasn't received anything
                req_id = None
                cmd_args = (None, None)

            json_cmd_args = json.dumps(cmd_args)
            # broadcast command from rank 0 to all ranks
            json_cmd_args = MAIN.do_mpi_bcast(json_cmd_args)
            json_cmd_args = json_cmd_args.tostring().rstrip()
            json_cmd_args = json_cmd_args.decode('utf-8')
            command, args = json.loads(json_cmd_args)

            if command:
                cmd = command.split(' ')[0].lower().strip()

                if hasattr(self, cmd):
                    self.responses[req_id] = getattr(self, cmd)(args)
                else:
                    self.responses[req_id] = 500, 'UNRECOGNIZED COMMAND\n'

            self.update_status()
            if not self.paused:
                return
            if not DEBUG.good_config:
                return
            # loop until unpaused
            time.sleep(0.1)

    def unpause(self, _):
        " unpause "
        if DEBUG.good_config:
            self.paused = False
            return 200, "UNPAUSING MFIX"
        else:
            return 200, "UNABLE TO UNPAUSE MFIX"

    def pause(self, _):
        " paused "
        self.paused = True
        return 200, "PAUSING MFIX"

    def write_dbg_vt(self, _):
        " call write_dbg_vtu_and_vtp_files "
        MAIN.do_write_dbg_vtu_and_vtp_files()
        return 200, 'Calling WRITE_DBG_VTU_AND_VTP_FILES\n'

    def backupres(self, _):
        " backup resource files"
        MAIN.do_backupres()
        return 200, 'BACKING UP RESOURCE FILES\n'

    def reinit(self, _):
        " reinitialize "
        self.mfix_dat = _.get('mfix_dat')
        MAIN.do_reinit(self.mfix_dat)
        return 200, 'REINITIALIZING MFIX\n'

    def exit(self, _):
        " run_mfix thread should exit cleanly "
        self.stopped = True
        self.paused = False
        return 200, 'EXITING MFIX\n'

    def step(self, args):
        " take one or more timesteps "
        stepcount = int(args.get('stepcount', None)[0])
        if RUN.tstop <= RUN.time:
            RUN.tstop = RUN.tstop + stepcount*RUN.dt
        for _ in range(stepcount):
            self.do_step()
        return 200, 'DOING %s TIMESTEP(S)\n' % stepcount

    def do_command(self, cmd, args=None):
        "Puts a command that was received over the web interface on the queue"
        req_id = threading.current_thread().ident
        self.requests[req_id] = (cmd, args)
        while req_id not in self.responses:
            time.sleep(0.1)
        resp = self.responses[req_id]
        del self.responses[req_id]

        return resp


# Flask request routing

def token_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        client_token = request.headers.get(FLASK_APP.config['TOKEN_NAME'])
        server_token = FLASK_APP.config['SECRET_KEY']
        client_token = server_token
        if client_token == server_token:
            return f(*args, **kwargs)
        else:
            return "Authentication required", 401, \
                  {'Content-Type': 'text/plain; charset=utf-8'}
    return decorated_function

def api_response(status_code, command_output):
    return make_response(
      jsonify(mfix_status=mfix_thread.status,
              command_output=command_output),
      status_code,
      {'Content-Type': 'application/json; charset=utf-8'})

@FLASK_APP.before_request
def debug_before():
    if DEBUG_FLAG:
        print('REQUEST HEADERS')
        print(request.headers)
        print('REQUEST BODY')
        print(request.get_data())
        print('WSGI REQUEST')
        print(request.environ['wsgi_req_copy'])

@FLASK_APP.route('/reinitialize', methods=['POST'])
@token_required
def reinitialize():
    """Save submitted project file and pass it to mfix_thread.do_reinit"""

    # request content type should be "application/json",
    # request body must be a JSON formatted string with one
    # key, "project_file", containing the project file text

    project_file = request.get_json(force=True)
    project_str = project_file.get('project_file')
    try:
        prefix = '%s.' % get_run_name()
        with tempfile.NamedTemporaryFile(prefix=prefix, delete=DEBUG_FLAG, dir=os.getcwd()) as tmp:
            tmp.write(project_str)
            tmp.flush()
            # MFIX truncates path to 80 characters, so try to keep it short.
            # split tmp name (defaults to absolute) on cwd, remove leading slash
            relative_name = tmp.name.split(os.getcwdu())[-1].lstrip('/')
            status_code, command_output = \
              mfix_thread.do_command("REINIT", args={'mfix_dat': relative_name})
    except Exception as e:
        status_code = 500
        command_output = "Error saving submitted project file"
    return api_response(status_code, command_output)

@FLASK_APP.route('/set/<modname>/<varname>', methods=['POST'])
@token_required
def set_variable(modname, varname):
    "sets a variable"
    args = dict(request.form)
    args['modname'] = modname
    args['varname'] = varname
    status_code, command_output = mfix_thread.do_command("SET", args=args)
    return api_response(status_code, command_output)

@FLASK_APP.route('/set/<modname>/<varname>/<elem>', methods=['POST'])
@token_required
def set_variable_array(modname, varname, elem):
    "sets a variable"
    args = dict(request.form)
    args['modname'] = modname
    args['varname'] = varname
    args['elem'] = elem
    status_code, command_output = mfix_thread.do_command("SET", args=args)
    return api_response(status_code, command_output)

@FLASK_APP.route('/get/<modname>/<varname>', methods=['GET'])
@token_required
def get_variable(modname, varname):
    "retrieves a variable"
    args = dict(request.args)
    args['modname'] = modname
    args['varname'] = varname
    status_code, command_output = mfix_thread.do_command("GET", args=args)
    return api_response(status_code, command_output)

@FLASK_APP.route('/get/<modname>/<varname>/<elem>', methods=['GET'])
@token_required
def get_variable_array(modname, varname, elem):
    "retrieves a variable"
    args = dict(request.args)
    args['modname'] = modname
    args['varname'] = varname
    args['elem'] = elem
    status_code, command_output = mfix_thread.do_command("GET", args=args)
    return api_response(status_code, command_output)


@FLASK_APP.route('/write_dbg_vt', methods=['POST'])
@token_required
def write_dbg_vt():
    "calls WRITE_DBG_VTU_AND_VTP_FILES"
    status_code, command_output = mfix_thread.do_command("WRITE_DBG_VT")
    return api_response(status_code, command_output)


@FLASK_APP.route('/backupres', methods=['POST'])
@token_required
def backupres():
    status_code, command_output = mfix_thread.do_command("BACKUPRES")
    return api_response(status_code, command_output)


@FLASK_APP.route('/exit', methods=['POST'])
@token_required
def exit_mfix():
    status_code, command_output = mfix_thread.do_command("EXIT")
    mfix_thread.thread.join()
    global pidfilename
    os.remove(pidfilename)
    shutdown = request.environ.get('werkzeug.server.shutdown')
    if shutdown is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    shutdown() # will finish current request, then shutdown
    return api_response(status_code, command_output)


@FLASK_APP.route('/step', methods=['POST'])
@token_required
def step():
    """runs mfix for one timestep, regardless of TIME and TSTOP"""
    args = dict(request.form)
    status_code, command_output = mfix_thread.do_command("STEP", args=args)
    return api_response(status_code, command_output)


@FLASK_APP.route('/pause', methods=['PUT'])
@token_required
def pause():
    "pauses MFIX if unpaused"
    status_code, command_output = mfix_thread.do_command("PAUSE")
    return api_response(status_code, command_output)

@FLASK_APP.route('/unpause', methods=['PUT'])
@token_required
def unpause():
    "unpause MFIX if paused"
    status_code, command_output = mfix_thread.do_command("UNPAUSE")
    return api_response(status_code, command_output)

@FLASK_APP.route('/status', methods=['GET'])
@token_required
def get_status():
    "returns current status: paused state, current time, and residuals"
    # status is included in api_response
    return api_response(200, '')

@FLASK_APP.route('/logging/<state>', methods=['POST'])
@FLASK_APP.route('/logging/<state>/<level>', methods=['POST'])
@token_required
def set_logging(state=None, level=None):
    """toggle Flask logging"""
    if state == 'enable':
        log.disabled = False
    if state == 'disable':
        log.disabled = True
    # maybe handle loglevel someday
    return api_response(200, 'OK')


def conv_bool(s):
    s = s.lower().replace('.','')
    if s == 'true' or s == 't':
        return '.T.'
    elif s == 'false' or s == 'f':
        return '.F.'
    else:
        raise ValueError

# from http://stackoverflow.com/questions/18800328/python-read-in-multiple-key-value-dict-from-command-line-into-a-variable
class DictAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        #Implementation is from argparse._AppendAction
        items = copy.copy(argparse._ensure_value(namespace, self.dest, {}))  # Default mutables, use copy!
        for value in values:
            try:
                k, v = value.split("=", 1)

                # all values are strings, try infering type form the string so
                # that mfix handles it
                converted_v = None
                for conv in [int, float, conv_bool]:
                    try:
                        converted_v = str(conv(v))
                        break
                    except ValueError:
                        continue

                if converted_v is None:
                    converted_v = "'"+v+"'"

                if not v:
                    raise argparse.ArgumentError(self, "Value not provided for: %s" % k)
            except ValueError:
                raise argparse.ArgumentError(self, "Format must be key=value")
            items[k] = converted_v
        setattr(namespace, self.dest, items)

class PrintFlagsAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            MAIN.print_flags()
            sys.exit(1)

def check_port(value):
    ivalue = int(value)
    if ivalue > 1025 or ivalue < 65536:
         raise argparse.ArgumentTypeError("%s is an invalid port number, must be 1025 < port < 65536 " % value)
    return ivalue

def parse_command_line_arguments():
    "handle command line arguments"
    parser = argparse.ArgumentParser(description='Welcome to PYMFIX')
    parser.add_argument('MFIX_KEY=VALUE', action=DictAction, nargs='*',
                        help='Series of MFIX_KEY=VALUE to over-ride values in the project file. Does not support indices.')
    parser.add_argument('-f', '--file',  metavar='FILE', action='store',
                        help='specify an input file (*.mfx or *.dat)',
                        required=True)
    parser.add_argument('-p', '--print-flags', action=PrintFlagsAction, nargs=0,
                        help='return the compile flags and exit')
    parser.add_argument('-P', '--port', metavar='PORT', action='store', default=random.randint(1025, 65536),
                        type=check_port,
                        help='specify a port number to use')
    parser.add_argument('-s', '--start', action='store_false',
                        help='do not wait for api connection to run')
    parser.add_argument('-v', '--version', action='version', version=__version_str__)

    args = parser.parse_args()

    passed_kwargs = ['='.join([k,v]) for k,v in vars(args)['MFIX_KEY=VALUE'].items()]
    return args.file.ljust(80), args.start, args.port, passed_kwargs

if __name__ == '__main__':
    main()
