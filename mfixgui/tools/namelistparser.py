# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
MFIX [Multiphase Flow with Interphase eXchanges] is a general-purpose
computer code developed at the National Energy Technology Laboratory
[NETL] for describing the hydrodynamics, heat transfer and chemical
reactions in fluid-solid systems.

Please visit: https://mfix.netl.doe.gov/

This python file contains a function to parse the MFiX init_namelist
files to extract documentation of the keywords.

@author: Justin Weber
"""
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

import re, os, glob
import json
import codecs
from collections import OrderedDict

def buildKeywordDoc(mfixSourcePath):

    searchPath = os.path.join(mfixSourcePath, 'model')

    if not os.path.exists(searchPath):
        searchPath = mfixSourcePath

    mfixKeywordDict = {'keywordlist':[], 'categories':[]}

    searchPathList = []
    for root, directory, filename in os.walk(searchPath):
        searchPathList += findNamePath(os.path.join(searchPath,root))

    sortedPathListTemp = searchPathList

    for i in range(2):
        sortedPathList = []
        for fname in sortedPathListTemp:

            base = os.path.basename(fname)
            if base.lower() == 'init_namelist.f':
                sortedPathList.insert(0, fname)
            elif base.lower() == 'cartesian_grid_init_namelist.f':
                sortedPathList.insert(1, fname)
            elif base.lower() == 'des_init_namelist.f':
                sortedPathList.insert(2, fname)
            else:
                sortedPathList.append(fname)
        sortedPathListTemp = sortedPathList

    for fname in sortedPathList:
        parsedNameList = parse(fname = fname)

        keywordlist = []
        if 'keywordlist' in mfixKeywordDict:
            keywordlist = mfixKeywordDict['keywordlist']

        catlist = []
        if 'categories' in mfixKeywordDict:
            catlist = mfixKeywordDict['categories']

        mfixKeywordDict.update(parsedNameList)
        mfixKeywordDict['keywordlist'] = keywordlist+parsedNameList['keywordlist']
        mfixKeywordDict['categories'] = catlist+parsedNameList['categories']
    cgs_re = re.compile(r' *\[.*CGS.*\] *', flags=re.IGNORECASE)
    def redact(v):
        if not isinstance(v, dict): # 'categories' and 'keywordlist'
            return v
        d = v.get('description')
        if d:
            d = cgs_re.sub('', d)
            d = d.replace("Youngs", "Young's")
            d = d.replace("Poissons", "Poisson's")
            d = d.replace("Poisson ratio", "Poisson's ratio")
            d = d.replace(" (in units of g/cm^2.s)", "")
            d = d.replace(" By default, the gravity force acts in the negative y-direction.", "")
            v['description'] = d
        return v

    return dict((k.lower(),redact(v)) for k, v in mfixKeywordDict.items())

def findNamePath(path):
    pattern = os.path.join(path, insensitiveGlobPattern('*namelist*.f'))
    pattern = os.path.realpath(pattern)
    return glob.glob(pattern)

def insensitiveGlobPattern(pattern):
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return ''.join(map(either,pattern))

def cleanString(string):
    string = string.strip()
    if string:
        strings = string.split('!')
        strings = [s.strip() for s in strings]
        return ' '.join(strings)

    return string


def parse(fname = None, string = None):
    '''
    Read mfix namelists to generate documentation.

    returns dictionary

    !<keyword category="category name" required="true/false"
    !    tfm="true/false" dem="true/false" pic="true/false"
    !                                    legacy="true/false">
    !  <description></description>
    !  <arg index="" id="" max="" min=""/>
    !  <dependent keyword="" value="DEFINED"/>
    !  <conflict keyword="" value="DEFINED"/>
    !  <valid value="" note="" alias=""/>
    !  <range min="" max="" />
      MFIX_KEYWORD = INIT_VALUE
    !</keyword>


    UNDEFINED <- double
    UNDEFINED_I <-integer
    UNDEFINED_C <- character
    ZERO <- double
    .true./.false. <- logical
    \d* <- integer
    \d*\.\d*[Dd]?\d? <- double

    '''

    keyvalueMatch = re.compile('''([^=^ ]+)=["'](.+?)["']''', re.DOTALL|re.M)          # key = "value" matching
    keywordAttMatch = re.compile('<keyword(.*?)>', re.DOTALL|re.M)                     # <keyword ... > matching
    keywordMatch = re.compile('(! *?<keyword.*?)</keyword>', re.DOTALL|re.M)           # <keyword ... </keyword> matching
    descriptionMatch = re.compile('<description>(.*?)</description>', re.DOTALL|re.M)  # <description> ... </description> matching
    argMatch = re.compile('<arg(.*?)/>', re.DOTALL|re.M)                               # <arg ... /> matching
    dependentMatch = re.compile('<dependent(.*?)/>',  re.DOTALL|re.M)                  # <dependent ... /> matching
    conflictMatch = re.compile('<conflict(.*?)/>',  re.DOTALL|re.M)                    # <conflict ... /> matching
    validMatch = re.compile('<valid(.*?)/>',  re.DOTALL|re.M)                          # <valid ... /> matching
    rangeMatch = re.compile('<range(.*?)/>',  re.DOTALL|re.M)                          # <range ... /> matching

    parMatch = re.compile('(\(.*\))')

    doubleMatch = re.compile('^[+-]?\d*\.?\d*[Dd]?[+-]?\d?$')
    intMatch = re.compile('^[+-]?[0-9]+$')
    strMatch = re.compile('''^['"]?([a-zA-Z]*)['"]?''')
    dtypeDict = {'undefined': 'DP',
                 'undefined_i': 'I',
                 'undefined_c': 'C',
                 'undefined_l': 'L',
                 'zero': 'DP',
                 'one': 'DP',
                 '-one': 'DP',
                 '.true.': 'L',
                 '.false.': 'L',
                 '.t.': 'L',
                 '.f.': 'L',
                 'LARGE_NUMBER'.lower(): 'DP',
                 '-LARGE_NUMBER'.lower(): 'DP',
                 }

    keywordDocDict = {}
    keywordDocList = []
    keywordDocCat = set()

    if fname:
        with codecs.open(fname, "r", "utf-8") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if 'subroutine' in line.lower() and line[0] != '!':
                    break

        keywordBlocks = keywordMatch.findall(''.join(lines[i:]))

    elif string:
        keywordBlocks = keywordMatch.findall(string)

    else:
        return {}

    for keywordBlockCount, keywordBlock in enumerate(keywordBlocks):

        keywordList = []
        for line in keywordBlock.split('\n'):

            splitLine = line.split('!')

            if len(splitLine[0]) > 2:
                keyword = splitLine[0].split('=')
                if len(keyword[0]) > 0:
                    keywordList.append(keyword)

        if len(keywordList) > 1:
#            mylogger.warning('Multiple keywords found in keyword comment block number: {}'.format(keywordBlockCount))
            pass
        elif len(keywordList) <= 0:
#            mylogger.warning('No keyword found in keyword comment block number: {}, line: {}'.format(keywordBlockCount, line))
            continue

        init_value = keywordList[0][-1].strip()
        initPythonValue = None
        keyword = keywordList[0][0].strip()

        #  keyword
        for parin in parMatch.findall(keyword):
            keyword = keyword.replace(parin, '')

        keyword = keyword.strip()

        # get dtype from init_value
        dtype = None
        init_value = init_value.lower().replace(' ', '')
        if init_value in dtypeDict.keys():
            dtype = dtypeDict[init_value]
        elif intMatch.findall(init_value):
            dtype = 'I'
        elif doubleMatch.findall(init_value):
            dtype = 'DP'
        elif strMatch.findall(init_value):
            dtype = 'C'
        else:
#            mylogger.warning('Cold not determine init value: %s' % init_value)
            pass

        # Parse initial value, convert to python type
        if 'undefined' not in init_value:
            if dtype == 'I':
                try:
                    initPythonValue = int(init_value)
                except:
                    pass
            elif dtype == 'DP':
                initPythonValue = init_value.replace('d', 'e').replace('one', '1').replace('zero', '0')

                try:
                    initPythonValue = float(initPythonValue)
                except ValueError:
                    pass

            elif dtype == 'L':
                if init_value in ['.f.', '.false.']:
                    initPythonValue = False
                elif init_value in ['.t.', '.true.']:
                    initPythonValue = True

        # find keyword attritbutes: category="" required="" legacy=""
        keywordAtts = None
        for match in keywordAttMatch.findall(keywordBlock):
            keywordAtts = keyvalueMatch.findall(match)

        # find description
        description = descriptionMatch.findall(keywordBlock)
        if not description:
#            mylogger.warning('No description defined for: %s' % keyword)
            description = ['None']

        # find arguments and argument attritbutes: index="" id="" max="" min=""
        args = []
        for match in argMatch.findall(keywordBlock):
            args.append(keyvalueMatch.findall(match))

        # find dependents and dependent attritbutes: keyword="" value=""
        dependents = []
        for match in dependentMatch.findall(keywordBlock):
            dependents.append(keyvalueMatch.findall(match))

        # find conflicts and conflict attritbutes: keyword="" value=""
        conflicts = []
        for match in conflictMatch.findall(keywordBlock):
            conflicts.append(keyvalueMatch.findall(match))

        # find valids and valid attritbutes: value="" note="" alias=""
        valids = []
        for match in validMatch.findall(keywordBlock):
            valids.append(keyvalueMatch.findall(match))

        # find ranges and range attritbutes: min="" max=""
        def try_float(kv):
            k,v = kv
            try:
                return k, float(v)
            except ValueError:
                return k, v

        validrange = None
        for match in rangeMatch.findall(keywordBlock):
            validrange = map(try_float, keyvalueMatch.findall(match))
            break

        # keyword list (to preserve order)
        keywordDocList.append(keyword.lower())

        # Build dictionary
        keywordDocDict[keyword] = {'init': init_value,
                                   'initpython': initPythonValue,
                                   'dtype': dtype,
                                   'description': cleanString(description[0]),
                                   'category': 'other',
                                   'required': False,
                                   'tfm': False,
                                   'dem': False,
                                   'pic': False,
                                   'legacy': False,
                                   'args': OrderedDict(),
                                   'dependents': OrderedDict(),
                                   'conflicts': OrderedDict(),
                                   'valids': OrderedDict(),
                                   'validrange': OrderedDict(),
                                   }

        def tryBool(s):
            return True if s=='true' else False if s=='false' else s

        for keywordAtt in keywordAtts:
            key = keywordAtt[0].strip().lower()
            val = keywordAtt[-1].lower()
            keywordDocDict[keyword][key] = tryBool(cleanString(val))

        keywordDocCat.add(keywordDocDict[keyword]['category'])

        for arg in args:
            argDict = dict(arg)
            for key in ['id', 'min', 'max', 'index']:
                if key not in argDict.keys():
                    argDict[key] = None
            if argDict['index']:
                keywordDocDict[keyword]['args'][argDict['index']]= {'id': argDict['id'],
                                                                    'min': argDict['min'],
                                                                    'max': argDict['max'],
                                                                    }
        for dependent in dependents:
            dependentDict = dict(dependent)
            for key in ['keyword', 'value']:
                if key not in dependentDict.keys():
                    dependentDict[key] = None
                else:
                    dependentDict[key] = cleanString(dependentDict[key])
            if dependentDict['keyword']:
                keywordDocDict[keyword]['dependents'][dependentDict['keyword']]=dependentDict

        for conflict in conflicts:
            conflictDict = dict(conflict)
            for key in ['keyword', 'value']:
                if key not in conflictDict.keys():
                    conflictDict[key] = None
                else:
                    conflictDict[key] = cleanString(conflictDict[key])
            if conflictDict['keyword']:
                keywordDocDict[keyword]['conflicts'][conflictDict['keyword']]=conflictDict

        for valid in valids:
            validDict = dict(valid)
            for key in ['value', 'note', 'alias']:
                if key not in validDict.keys():
                    validDict[key] = None
                else:
                    validDict[key] = cleanString(validDict[key])
            if validDict['value']:
                keywordDocDict[keyword]['valids'][validDict['value']]=validDict

        if validrange:
            keywordDocDict[keyword]['validrange'] = dict(validrange)




    keywordDocDict['keywordlist'] = keywordDocList
    keywordDocDict['categories'] = list(keywordDocCat)

    return keywordDocDict

def writeFiles(d):
    json.dump(d, open('./keywordDoc.json','w'))
    with open('./keywordList.txt','w') as fb:
        fb.write('\n'.join(d['keywordlist']))

def test_1():
    testString = '''!<keyword category="category name" required="true/false"
                    !                                     legacy="true/false">
                     !    <description>Test description
                    !               over multiple lines
                      !   </description>
                    !  <arg index="1" id="argID"
                     !        max="+Inf" min="-inf"/>
                     !  <arg index="" id="" max="" min=""/>
                    !   <dependent keyword="dependentKeyword"
                    !                        value="DEFINED"/>
                     !  <dependent keyword="" value="DEFINED"/>
                    !  <conflict keyword="conflictKeyword" value="DEFINED"/>
                    !   comment inside block
                    !  <conflict keyword="" value="DEFINED"/>
                    !  <valid value="VALIDvALUE" note="this is a
                    !      valid value" alias=""/>
                    !  <valid value="" note="" alias=""/>
                       !  <range min="" max="" />
                    !  <range min="0" max="2D0" />
                        MFIX_KEYWORD = INIT_VALUE !comment
                      !   </keyword>
                '''
    return parse(string = testString)

def test_2():
    testString = '''
                    !<keyword>
                    !    <description>character</description>
                        CHAR1 = UNDEFINED_C
                    !</keyword>
                    !<keyword>
                    !    <description>character</description>
                        CHAR2 = 'LEBOWITZ'
                    !</keyword>
                    !<keyword>
                    !    <description>character</description>
                        CHAR3 = 'SYAM_OBRIEN'
                    !</keyword>
                    !<keyword>
                    !    <description>integer</description>
                        INT1 = UNDEFINED_I
                    !</keyword>
                    !<keyword>
                    !    <description>integer</description>
                        INT2 = 10
                    !</keyword>
                    !<keyword>
                    !    <description>bool</description>
                        BOOL1 = .true.
                    !</keyword>
                    !<keyword>
                    !    <description>bool</description>
                        BOOL2 = .t.
                    !</keyword>
                    !<keyword>
                    !    <description>bool</description>
                        BOOL3 = .false.
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE1 = 1.
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE2 = 1.0D0
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE3 = 1.0
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE4 = UNDEFINED
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE5 = 0d0
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE6 = LARGE_NUMBER
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE7 = -LARGE_NUMBER
                    !</keyword>
                    !<keyword>
                    !    <description>double</description>
                        DOUBLE8 = - LARGE_NUMBER
                    !</keyword>
                    !<keyword>
                    !    <description>zero</description>
                        DOUBLE9 = ZERO
                    !</keyword>
                    !<keyword>
                    !    <description>zero</description>
                        DOUBLE10 = ONE
                    !</keyword>
                    !<keyword>
                    !    <description>zero</description>
                        DOUBLE11 = - one
                    !</keyword>
                 '''
    docDict = parse(string = testString)
    for key in sorted(docDict.keys()):
        if 'dtype' in docDict[key]:
            print(key, '\t',  docDict[key]['dtype'])
    return {}

def test_3():
    return parse(fname = './init_namelist_test.f')

def test_4():
    return buildKeywordDoc('./')

def test_5():
    testString = '''!<keyword category="category name" required="true/false"
                    !                                     legacy="true/false">
                     !    <description>Test description
                    !               over multiple lines
                      !   </description>
                    !  <arg index="1" id="argID"
                     !        max="+Inf" min="-inf"/>
                     !  <arg index="" id="" max="" min=""/>
                    !   <dependent keyword="dependentKeyword"
                    !                        value="DEFINED"/>
                     !  <dependent keyword="" value="DEFINED"/>
                    !  <conflict keyword="conflictKeyword" value="DEFINED"/>
                    !   comment inside block
                    !  <conflict keyword="" value="DEFINED"/>
                    !  <valid value="VALIDvALUE" note="this is a
                    !      valid value" alias="dude"/>
                    !  <valid value="" note="" alias="cool"/>
                       !  <range min="" max="" />
                    !  <range min="0" max="2D0" />
                        MFIX_KEYWORD = INIT_VALUE !comment
                      !   </keyword>
                '''
    return parse(string = testString)

if __name__ == '__main__':

    docDict = test_5()
    print(docDict)

    #writeFiles(docDict)

    #from pprint import pprint
    #pprint(docDict)
