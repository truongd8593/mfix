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
        go                                  - (re)start mfix
        list                                - list variables that can be set
        set VAR=VALUE [ i1 j1 k1 i2 j2 k2 ] - set variable to value [ for a certain index range ]
        step                                - execute one timestep
        stop                                - stop mfix
        \n'''

    elif cmd=='list':
        resp = ''
        for module in ('compar','discretelement','fldvar','main','param','run','run_dp'):
            attr = getattr(pymfix,module)
            resp = resp + 'in module %s you can set attributes %s \n\n' % (module,dir(attr))
        resp = resp + " Except for those of the above that are really subroutines/functions instead of variables. Also, the above list is incomplete, and is missing most of the interesting variables. \n\n  Don't forget to prefix with pymfix; e.g. set pymfix.run.dt 0.001\n\n"
        responses[req_id] = resp
        mfix_stopped = True

    elif cmd=='stop':
        responses[req_id] = 'STOPPING MFIX\n'
        mfix_stopped = True

    elif cmd=='go':
        responses[req_id] = 'STARTING MFIX\n'
        mfix_stopped = False

    elif cmd=='step':
        pymfix.main.step()
        responses[req_id] = 'DOING ONE TIMESTEP\n'

    elif cmd=='set':
        pymfix.main.step()
        var = command.split(' ')[1].lower().strip()
        val = command.split(' ')[2].lower().strip()

        # slow things down for development
        time.sleep(pymfix.compar.mype)

        if len(command.split(' ')) > 3:
            ii,jj,kk,i2,j2,k2 = command.split(' ')[3:9]
            for i in range(int(ii),int(i2)+1):
                for j in range(int(jj),int(j2)+1):
                    for k in range(int(kk),int(k2)+1):
                        if is_on_mype_owns(i,j,k):
                            exec_string = ('%s[%d] = %s' % (var,funijk(i,j,k),val))
                            print "GOING TO EXECUTE: ",exec_string
                            exec(exec_string)
                            print "rank ",pymfix.compar.mype," set value for",i,j,k,funijk(i,j,k)
                        else:
                            print "rank ",pymfix.compar.mype,"does not own",i,j,k,funijk(i,j,k)
        else:
            exec('%s = %s' % (var,val))
        responses[req_id] = 'ok, I set %s to %s\n' % (var,val)

    else:
        responses[req_id] = 'UNRECOGNIZED COMMAND\n'

def funijk(i,j,k):
    return pymfix.compar.ijk_array_of[i,j,k]

def is_on_mype_owns(li, lj, lk):
      return li >= pymfix.compar.istart and li <= pymfix.compar.iend and lj >= pymfix.compar.jstart and lj <= pymfix.compar.jend and lk >= pymfix.compar.kstart and lk <= pymfix.compar.kend

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
