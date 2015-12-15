import os, sys
sys.path.append(os.getcwd())

import thread
import time

from twisted.internet import reactor, protocol

import pymfix

requests = {}
responses = {}

mfix_stopped = False

def run_mfix():
    global mfix_stopped

    pymfix.main.setup()
    pymfix.main.start()
    while pymfix.run.time < pymfix.run.tstop:
        if requests:
            # requests would only arrive at rank 0
            req_id,command = requests.popitem()
        else:
            # command is empty for rank>0
            req_id,command = '',''

        # broadcast command from rank 0 to all ranks
        command = pymfix.main.do_mpi_bcast(command)

        command = ''.join(command).strip()

        if command:
            handle(req_id,command)

        if mfix_stopped:
            time.sleep(1)
        else:
            pymfix.main.step()
    pymfix.main.end()

def handle(req_id,command):
    global mfix_stopped

    cmd = command.split(' ')[0].lower().strip()
    print "THE COMMAND IS",cmd
    if cmd=='help':
        responses[req_id] = '''Usage:
        stop          - stop mfix
        go            - (re)start mfix
        hello         - say hello from each rank
        step          - execute one timestep
        set VAR=VALUE - set variable to value
        \n'''

    elif cmd=='stop':
        responses[req_id] = 'STOPPING MFIX\n'
        mfix_stopped = True

    elif cmd=='go':
        responses[req_id] = 'STARTING MFIX\n'
        mfix_stopped = False

    elif cmd=='hello':
        print("HELLO FROM MPI RANK %d" % (pymfix.compar.mype))
        responses[req_id] = 'PRINTING MPI RANK INFO TO STDOUT\n'

    elif cmd=='step':
        pymfix.main.step()
        responses[req_id] = 'DOING ONE TIMESTEP\n'

    else:
        responses[req_id] = 'UNRECOGNIZED COMMAND\n'

class Echo(protocol.Protocol):
    """This is just about the simplest possible protocol"""

    def dataReceived(self, data):
        "Called with data sent from client."
        # Should make sure this routine is synchronized

        req_id = thread.get_ident()

        requests[req_id] = data

        while req_id not in responses:
            pass

        self.transport.write(responses[req_id])
        del responses[req_id]

        # change any variable...
        # pymfix.run.time  = 1000
        # pymfix.run.tstop = pymfix.run.tstop/2
        # pymfix.run.dt    = 10

        # response = 'the pressure is %s \n tstop is %s' % (str(pymfix.fldvar.p_g[:100]),pymfix.run.tstop)

def main():
    """This runs the protocol on port 8000"""
    thread.start_new_thread(run_mfix, ())
    time.sleep(1)

    if 0==pymfix.compar.mype:
        factory = protocol.ServerFactory()
        factory.protocol = Echo
        reactor.listenTCP(8000,factory)
        reactor.run()
    else:
        while(True): time.sleep(1)

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
