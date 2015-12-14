import os, sys
sys.path.append(os.getcwd())

import thread
import time

from twisted.internet import reactor, protocol

import pymfix

def run_mfix():
    pymfix.main.setup()
    pymfix.main.start()
    while pymfix.run.time < pymfix.run.tstop:
        if 1 == pymfix.main.request_pending:
            pymfix.main.mfix_stopped = 1
        while 1 == pymfix.main.mfix_stopped:
            time.sleep(1)
        pymfix.main.step()
    pymfix.main.end()


def handle():
    req = queue.get()
    q.task_done()

class Echo(protocol.Protocol):
    """This is just about the simplest possible protocol"""

    def dataReceived(self, data):
        "Called with data sent from client."
        # Should make sure this routine is synchronized
        pymfix.main.request_pending = 1

        # wait for mfix to stop
        while 1 != pymfix.main.mfix_stopped:
            pass

        # change any variable...
        # pymfix.run.time  = 1000
        # pymfix.run.tstop = pymfix.run.tstop/2
        # pymfix.run.dt    = 10

        cmd = data.split(' ')[0].lower().strip()
        print "THE DATA IS",data
        print "THE COMMAND IS",cmd
        if cmd=='stop':
            response = 'STOPPING MFIX\n'
            pymfix.main.mfix_stopped = 1
        elif cmd=='go':
            response = 'STARTING MFIX\n'
            pymfix.main.mfix_stopped = 0
        elif cmd=='hello':
            print("HELLO FROM MPI RANK %d" % (pymfix.compar.mype))
            response = 'PRINTING MPI RANK INFO TO STDOUT\n'
        elif cmd=='step':
            pymfix.main.step()
            response = 'DOING ONE TIMESTEP\n'
        else:
            response = 'UNRECOGNIZED COMMAND\n'

        # response = 'the pressure is %s \n tstop is %s' % (str(pymfix.fldvar.p_g[:100]),pymfix.run.tstop)
        self.transport.write(response)

        # resume mfix
        pymfix.main.request_pending = 0

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
