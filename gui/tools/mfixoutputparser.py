# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:35:42 2014

@author: jweber
"""
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import re
import numpy as np

from mfixgui.tools.general import numToTime

class OutputParser(object):
    def __init__(self):
        self.clear()

        self.timeDtRegex = re.compile('time[\ ]+=[\ ]+ ([\d.e\+-]+)[\ ]+dt[\ ]+=[\ ]+([\d.e-]+)', re.IGNORECASE|re.MULTILINE)
        self.smRegex = re.compile('t=[\ ]+([\d.e-]+).*sm=[\ ]*([\d.e-]+)', re.IGNORECASE|re.MULTILINE)
        self.wtrRegex = re.compile('wall time remaining[\ ]+=[\ ]+ ([\d.e-]+)[\ ]+([a-z]+)', re.IGNORECASE|re.MULTILINE)
        self.wtrRegex2015 = re.compile('Est. Remaining:[\ ]+([\d.e-]+)([a-z]+)', re.IGNORECASE|re.MULTILINE)
        self.twtRegex = re.compile('total wall time used[\ ]+=[\ ]+([\d.e-]+)[\ ]+([a-z]+)', re.IGNORECASE|re.MULTILINE)
        self.cpuRegex = re.compile('elapsed cpu time[\ ]+=[\ ]+([\d\.e\-\+]+)[\ ]+sec', re.IGNORECASE|re.MULTILINE)
        self.waltimeelapsed2015 = re.compile('Wall Time - Elapsed:[\ ]+([\d.e-]+)([a-z]+)', re.IGNORECASE|re.MULTILINE)
        self.errorRegex = re.compile(r'error[ ]+?([\d]+)\:([a-z _\(\)\d:\n\,\.\*]+?)\n \*+?\n', re.IGNORECASE|re.MULTILINE)
        self.termRegex = re.compile(r'program[ ]+terminated', re.IGNORECASE|re.MULTILINE)
        self.mfixStopRegex = re.compile(r'[= ]+(MFIX STOP SIGNAL DETECTED)[= ]+', re.IGNORECASE|re.MULTILINE)

        self.setResidualList(['nit', 'p0', 'p1', 'u0', 'v0', 'u1', 'v1','maxres'])

    def setResidualList(self, resList):
        self.residualList = resList
        self.residualRegex = re.compile(('([\d]+)'+
                                         '[\ ]+([\d\.e\-\+\*]+)'*(len(self.residualList)-2)+
                                         '[\ ]+([\dpvu]+)'
                                         ), re.IGNORECASE|re.MULTILINE)

    def parse(self, text):

        self.tail = text
        text = ''.join(text)

        timeDt = self.timeDtRegex.findall(text)
        if timeDt:
            self.time = float(timeDt[-1][0])

            self.dtlist+=[float(dt[1]) for dt in timeDt]

        wtr = self.wtrRegex.findall(text)
        if wtr:
            self.wallTimeRemaining = [float(wtr[-1][0]), wtr[-1][1]]
        # 2015-1 release
        wtr = self.wtrRegex2015.findall(text)
        if wtr:
            self.wallTimeRemaining = [float(wtr[-1][0]), wtr[-1][1]]

        twt = self.twtRegex.findall(text)
        if twt:
            self.totalWallTimeUsed = [float(twt[-1][0]), twt[-1][1]]

        sm = self.smRegex.findall(text)
        if sm:
            self.smtimelist+=[float(m[0]) for m in sm]
            self.smlist+=[float(m[1]) for m in sm]

        cpu = self.cpuRegex.findall(text)
        if cpu:
            self.cputime=float(cpu[-1])

        # 2015-1 release
        cpu = self.waltimeelapsed2015.findall(text)
        if cpu:
            self.cputime=(float(cpu[-1][0]), cpu[-1][1])

        residuals = self.residualRegex.findall(text)
        if residuals:
            if len(residuals[0]) == 8:
                for res in residuals:

                    for var, num in zip(self.residualList, res):
                        if var in ['nit', 'maxres']:
                            self.residualDict[var].append(num)
                        else:
                            try:
                                f = float(num)
                            except ValueError:
                                f = float('nan')

                            self.residualDict[var].append(f)

        errors = self.errorRegex.findall(text)
        if errors:
            self.errorList+=errors

        if self.termRegex.findall(text):
            self.terminated = True

        if self.mfixStopRegex.findall(text):
            self.mfixStop = True

        # Calculate Stats
        if self.dtlist:
            self.dtmean = np.mean(self.dtlist)
            self.dtstd = np.std(self.dtlist)

        try:
            if isinstance(self.cputime, tuple):
                self.sperday = self.time/numToTime(self.cputime[0], self.cputime[1], outunit='s')*60*60*24 #s/day
            else:
                self.sperday = self.time/self.cputime*60*60*24 #s/day
        except:
            self.sperday = 0

        # time in dd:hh:mm:ss
        #self.wallTimeRemainingTime = reduce(lambda ll,b : divmod(ll[0],b) + ll[1:], [(self.wallTimeRemaining,),60,60,24])

    def clear(self):

        # Variables
        self.terminated = False
        self.mfixStop = False
        self.errorList = []
        self.tail = []
        self.time = 0
        self.dtlist = []
        self.dtmean = 0
        self.dtstd = 0
        self.smlist = []
        self.smtimelist = []
        self.cputime = 0
        self.sperday = 0
        self.wallTimeRemaining = [float('inf'), 's']
        self.totalWallTimeUsed = None
        self.residualDict = {'nit':[],
                             'p0':[],
                             'p1':[],
                             'u0':[],
                             'v0':[],
                             'u1':[],
                             'v1':[],
                             'maxres':[],
                             }

def follow(thefile):
    thefile.seek(0, 2) # Go to the end of the file
    while True:
        line = thefile.readline()
        if not line:
            time.sleep(0.1) # Sleep briefly
            continue
        yield line

if __name__ == "__main__":
    import time

    runlogpath = 'C:/Users/weberjm/Desktop/test2/run.log'

    outputParser = OutputParser()
    updatetime = 1
    nsamples=1

    prevtime = time.time()
    sample=0
    with open(runlogpath, 'r') as logfile:
        while True:
            sample+=1
            lines = logfile.readlines()
            if lines:
                print('parsing')
                outputParser.parse(''.join(lines))
                print('done')

            print(outputParser.time)

            if sample>=nsamples:
                break

            time.sleep(updatetime)

    print(outputParser.time)
    print(outputParser.wallTimeRemaining)
    print(outputParser.totalWallTimeUsed)
    print(len(outputParser.dtlist))
    print(len(outputParser.smlist))
    print(outputParser.terminated)
    print(outputParser.errorList)