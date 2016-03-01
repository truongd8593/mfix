# -*- coding: utf-8 -*-
"""
MFIX [Multiphase Flow with Interphase eXchanges] is a general-purpose
computer code developed at the National Energy Technology Laboratory
[NETL] for describing the hydrodynamics, heat transfer and chemical
reactions in fluid-solid systems.

Please visit: https://mfix.netl.doe.gov/

This python file contains code for analyzing mfix.dat input files

Last update: 07/09/2013

@author: Justin Weber
"""

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

import re, shlex, logging
mylogger = logging.getLogger(__name__)
from mfixgui.guilib import config
CONFIG = config.guiConfig()
mylogger.setLevel(CONFIG.get('debugging', 'level').upper())

# Custom modules

def any(name, alternates):
    "Return a named group pattern matching list of alternates."
    return "(?P<%s>" % name + "|".join(alternates) + ")"

class MfixCodeAnalyzer(object):
    def __init__(self):
        # Checking Regular Expressions
        self.re_allvariables = re.compile("(\w+(?:\(\d+\))?)\s*")
        self.regX_keyValue = re.compile(r'(\w+(?:\(([ \d,]+)\))?)\s*=\s*(.*?)(?=(!|$|\w+(\(\d+\))?\s*=))')
        self.regX_float_exp = re.compile("^[+-]?[0-9]+\.([0-9]+)?(?:[eEdD][+-]?[0-9]+)?$")

        sqstring = r"(\b[rRuU])?'[^'\\\n]*(\\.[^'\\\n]*)*'?"
        dqstring = r'(\b[rRuU])?"[^"\\\n]*(\\.[^"\\\n]*)*"?'
        boolstr = [r'\.true\.', r'\.false\.', r'\.t\.', r'\.f\.']
        string = any("string", [sqstring, dqstring]+boolstr)
        #number = any("number",
        #             [r"\b[+-]?[0-9]+[lL]?\b",
        #              r"\b[+-]?0[xX][0-9A-Fa-f]+[lL]?\b",
        #              r"\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eEdD][+-]?[0-9]+)?\b",
        #              r"(@\([0-9+-/\*pi]+\))",
        #              ])
        self.re_remove = re.compile('|'.join([string]), re.S|re.I)

        self.mfixKeyWordDoc={}

    def setKeywordDoc(self, keywordDoc):
        self.mfixKeyWordDoc = keywordDoc

    def check_mfix_source(self, source_code):
        """Finds mfix.dat file errors"""
        results = []
        chemSection = [False, None, None]
        thermoSection = False
        keyval = {}

        mylogger.debug('Checking mfix source')

        for line, text in enumerate(source_code.splitlines()):
            # check for line > 80 characters
            if len(text)> 80:
                results.append(('Syntax Error:\nLine exceeds 80 characters.', line+1))

            # ignore comments
            if '#' in text or '!' in text:
                code = text.split('#')[0]
                if len(code)>0:
                    if '!' in code:
                        code = code.split('!')[0]
                    if len(code)<=0:
                        continue
                else:
                    continue
            else:
                code = text

            # look for Chemistry
            if '@(RXNS)'.lower() in text.lower():
                chemSection[0] = True
                chemSection[1] = line
            elif '@(END)'.lower() in text.lower() and chemSection[0]:
                chemSection[0] = False
                chemSection[2] = line

            # Check thermodata
            if 'THERMO DATA'.lower() in text.lower():
                thermoSection = True

            # Check code
            # Look for valid keywords
            if not chemSection[0] and not thermoSection and line>chemSection[2]:

                # Look for key, value pair
                for match in self.regX_keyValue.findall(code):
                    try:
                        vals = shlex.split( match[2].strip())
                    except ValueError:
                        vals=[]
                    cleanVals = []
                    for val in vals:
                        if '*' in val:
                            cleanVals = cleanVals + self.expandShortHand(val)
                        else:
                            val = self.cleanString(val)
                            if val is not None:
                                cleanVals.append(val)

                    if re.sub(r'\([^)]*\)', '', match[0].lower()) in self.mfixKeyWordDoc.keys():
                        if match[0].lower() not in keyval:
                            keyval[match[0].lower()] = {'line':line+1,
                                                        'vals':cleanVals,
                                                        'args':[i for i in match[1].split(',') if len(i)>0]
                                                        }
                        else:
                            results.append(("Repeated Keyword, '"+str(match[0].lower()) +
                                            "'\nFrom line: "+str(keyval[match[0].lower()]['line']), line+1))
                    else:
                        results.append(("Unknown keyword, '"+str(re.sub(r'\([^)]*\)', '', match[0].lower()))+"'", line+1))


                    code = code.replace(match[2],'')
                    code = code.replace(match[0].replace(match[2],''),'')
                code =code.replace('=','')

                if len(code.strip())>0:
                    results.append(("Syntax Error, '"+code.strip()+"'", line+1))

        # Check individual keyword values, etc.
        for keywordArgs in keyval.keys():
            keyword = re.sub(r'\([^)]*\)', '', keywordArgs).lower()

            # Check Values
            if keyval[keywordArgs]['vals']:
                for val in keyval[keywordArgs]['vals']:
                    # Check data type
                    if self.getType(val)!=self.mfixKeyWordDoc[keyword]['dtype']:
                        if not self.getType(val)=='I' and not self.mfixKeyWordDoc[keyword]['dtype']=='DP':
                            results.append(("Type Error: '"+str(val)+"' \nNeeds to be a "+self.mfixKeyWordDoc[keyword]['dtype'], keyval[keywordArgs]['line']))

                    # Check Boolean
                    elif self.mfixKeyWordDoc[keyword]['dtype']=='L':
                        if type(val)!=type(True):
                            results.append(("Value Error: '"+str(val)+"' \nIs not a valid value.", keyval[keywordArgs]['line']))

                    # Checl valid values
                    elif self.mfixKeyWordDoc[keyword]['valids']:
                        validList = [valid.lower() for valid in self.mfixKeyWordDoc[keyword]['valids'].keys()]

                        for valid in self.mfixKeyWordDoc[keyword]['valids'].values():
                            if valid['alias']:
                                validList.append(valid['alias'].lower())

                        if val.lower() not in validList:
                            results.append(("Value Error: '"+str(val)+"' \nIs not a valid value.", keyval[keywordArgs]['line']))

                    # Check valid range
                    elif self.mfixKeyWordDoc[keyword]['validrange']:
                        if float(val)< float(self.mfixKeyWordDoc[keyword]['validrange']['min']) or float(val)>float(self.mfixKeyWordDoc[keyword]['validrange']['max']):
                            results.append(("Value Error: '"+str(val)+"' \nIs outside the valid range ("+
                            self.mfixKeyWordDoc[keyword]['validrange']['min']+", "+self.mfixKeyWordDoc[keyword]['validrange']['max']+").", keyval[keywordArgs]['line']))
            else:
                results.append(("Value Error: '"+keyword+"' \nNeeds value(s).", keyval[keywordArgs]['line']))

            # Check arguments
            if keyval[keywordArgs]['args']:
                if len(keyval[keywordArgs]['args'])!=len(self.mfixKeyWordDoc[keyword]['args'].keys()):
                    results.append(("Argument Error: '"+keywordArgs+"' \nNeeds to have "+str(len(self.mfixKeyWordDoc[keyword]['args']))+' argument(s).', keyval[keywordArgs]['line']))

            # Check for legacy keywords
            if self.mfixKeyWordDoc[keyword]['legacy']=='true':
                results.append(("Legacy Keyword: '"+keyword+"' \nDo not use.", keyval[keywordArgs]['line']))

        mylogger.debug('Checking mfix source ... Done')
        return results

    def getType(self,string):
        if type(string)==type(True):
            return 'L'

        try:
            int(string)
            return 'I'
        except ValueError:
            pass

        try:
            float(string.replace('d','e'))
            return 'DP'
        except ValueError:
            return 'C'

    def cleanString(self, string):
        cleanVal = None

        string =  string.replace("'",'').replace('"','').replace(',','')
        if len(string)>0:
            if '.t.' == string.lower() or '.true.' == string.lower():
                cleanVal = True
            elif '.f.' == string.lower() or '.false.' == string.lower():
                cleanVal = False
            else:
                cleanVal = string.lower()

        return cleanVal

    def expandShortHand(self, string):
        splitString = string.split('*')
        expandList = [string]
        if len(splitString)==2:
            if splitString[0].isdigit():
                expandList = int(splitString[0])*[self.cleanString(splitString[1])]
        return expandList

    def get_completion_list(self, text, source_code, offset, filename):
        proposals = self.mfixKeyWordDoc.keys()
        return proposals

    def get_calltip_and_docs(self, text, source_code, offset, filename):
        objectName = text.upper()
        if objectName.lower() in self.mfixKeyWordDoc.keys():
            doc_text = self.mfixKeyWordDoc[objectName.lower()]['description']
            rtn = re.split('([.!?] *)', doc_text)
            doc_text =''.join([each.capitalize() for each in rtn])
            argspec = 'args:'
            note = 'notes:'
            objectName = objectName #+' [dtype: '+mfixkwdata[objectName.lower()]['dtype'].upper()+']'
        else:
            doc_text = "Unknown keyword: '"+ objectName +"'"
            argspec = 'args:'
            note = 'notes:'
        return [objectName, doc_text, argspec, note]

    def get_definition_location(self, source_code, offset, filename):
        filename = None
        lineno = None
        return filename, lineno


mfixCodeAnalyzerObject = None
def get_mfixDoc():
    """Create a single"""
    global mfixCodeAnalyzerObject
    if mfixCodeAnalyzerObject is None:
        mylogger.debug('Creating Single mfixCodeAnalyzerObject')
        mfixCodeAnalyzerObject = MfixCodeAnalyzer()
    return mfixCodeAnalyzerObject

def checkMfixDat(source_code):
    return get_mfixDoc().check_mfix_source(source_code)

if __name__ == "__main__":
    mfixDoc = get_mfixDoc()
    mfixDoc.setKeywordDoc({'description':{'legacy':'true','dtype':'C','valids':{},'validrange':{}}})
    source = """dude \ndescription = 'hello world'\n description=false"""
    print(mfixDoc.check_mfix_source(source))
    print(mfixDoc.get_completion_list('test', source, 10, 'none'))
    print(mfixDoc.get_calltip_and_docs('test', source, 10, 'none'))
