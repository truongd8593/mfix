# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 16:11:22 2015

@author: Weberjm
"""

import re
import os

class FortranWalker(object):
    def __init__(self):

        self.subroutine_re = re.compile('[ ]*SUBROUTINE[ ]+(.*)', re.IGNORECASE)

        self.use_re = re.compile('[ ]*use[ ]+(.+)', re.IGNORECASE)
        self.useonly_re = re.compile(' only[ :]*(.*)', re.IGNORECASE)

        self.integer_re = re.compile('[ ]*integer[^:]*::\s+(.+)', re.IGNORECASE)
        self.real_re = re.compile('[ ]*real[^:]*::\s+(.*)', re.IGNORECASE)
        self.double_re = re.compile('[ ]*DOUBLE PRECISION[^:]*::\s+(.*)', re.IGNORECASE)
        self.character_re = re.compile('[ ]*character[^:]*::\s+(.*)', re.IGNORECASE)
        self.logical_re = re.compile('[ ]*logical[^:]*::\s+(.*)', re.IGNORECASE)

    def walk(self, fname):

        print(self.parse(fname))

    def parse(self, fname):

        base = os.path.basename(fname)
        dataDict = {'modules':{},
                    'ints':[],
                    'reals':[],
                    'doubles':[],
                    'bools':[],
                    'chars':[],
                    }

        with open(fname) as f:
            for line in f:

                # chop comments
                line = line.split('!')[0]

                # subroutine
                subroutine = self.subroutine_re.findall(line)
                if subroutine:
                    print(subroutine)

                # Modules (use)
                use = self.use_re.findall(line)
                if use:
                    only = self.useonly_re.findall(use[0])
                    if only:
                        dataDict['modules'][use[0].split(',')[0]]=only[0].split(',')
                    else:
                        dataDict['modules'][use[0]]='*'

                # look for integers
                ints = self.integer_re.findall(line)
                if ints:
                    dataDict['ints']+=ints[0].split(',')

                # real
                reals = self.real_re.findall(line)
                if reals:
                    print(reals)

                # double
                doubles = self.double_re.findall(line)
                if doubles:
                    dataDict['doubles']+=doubles[0].split(',')

                # boolean
                bools = self.logical_re.findall(line)
                if bools:
                    dataDict['bools']+=[b.split('=')[0].strip() for b in bools[0].split(',')]

                # character
                chars = self.character_re.findall(line)
                if chars:
                    dataDict['chars']+=chars[0].split(',')

        return base, dataDict

if __name__ == '__main__':
    walker = FortranWalker()
    walker.walk('N:/MFIX/mfix/model/mfix.f')
