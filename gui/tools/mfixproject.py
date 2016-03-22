# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
MFIX [Multiphase Flow with Interphase eXchanges] is a general-purpose
computer code developed at the National Energy Technology Laboratory
[NETL] for describing the hydrodynamics, heat transfer and chemical
reactions in fluid-solid systems.

Please visit: https://mfix.netl.doe.gov/

This python file contains code to manage a mfix project.

Last update: 07/03/2014

@author: Justin Weber
"""
# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

# Python core imports
import shlex
import re
import math
import copy
try:
    # Python 2.X
    from StringIO import StringIO
except ImportError:
    # Python 3.X
    from io import StringIO

if __name__ == '__main__' and __package__ is None:
    import os, sys
    SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir, os.pardir))
    sys.path.append(SCRIPT_DIRECTORY)

# local imports
from .simpleeval import simple_eval
# from mfixgui.guilib.spydereditor.sourcecode.encoding import to_unicode_from_fs
# from mfixgui.guilib.spydereditor.sourcecode.p3compat import is_string, is_unicode, to_text_string

class FloatExp(float):
    fmt = '4'
    def __repr__(self):
      return '{:.{}e}'.format(self, self.fmt)

class Equation(object):
    def __init__(self, eq):

        # Check to make sure eq is not already an equation object
        if isinstance(eq, Equation):
            self.eq = eq.eq
        else:
            self.eq = str(eq)

    def _eval(self):
        if len(self.eq)==0:
            return 0
        elif self.eq is None or self.eq == 'None':
            return float('nan')
        else:
            try:
                return simple_eval(self.eq.lower(), names={"pi": math.pi, "e": math.e, "p": 1})
            except SyntaxError:
                return 0

    def __nonzero__(self):
        # Python 2
        if float(self)==float('nan'):
            return False
        else:
            return True

    def __bool__(self):
        # Python 3
        if float(self)==float('nan'):
            return False
        else:
            return True

    def __float__(self):
        return float(self._eval())

    def __int__(self):
        return int(self._eval())

    def __repr__(self):
        return ''.join(['@(', str(self.eq), ')'])

class KeyWord(object):
    def __init__(self, key, val, comment=None, row=None, col=None, dtype=None):

        self.key = key
        self.value = val
        self.dtype = dtype
        self.min = None
        self.max = None
        self.valids = None
        self.fmt = None
        self.col = col
        self.row = row
        self.comment = comment

        self.regX_expression = re.compile('([eEpiPI\+\-/*\^\(\)]+)')

        if dtype is None:
            self._checkdtype()

    def __deepcopy__(self, memo):
        return KeyWord(copy.copy(self.key),
                       copy.copy(self.value),
                       comment=copy.copy(self.comment),
                       row=copy.copy(self.row),
                       col=copy.copy(self.col),
                       dtype=copy.copy(self.dtype),
                      )

    def __float__(self):
        try:
            return float(self.value)
        except ValueError:
            return float('nan')
        except TypeError:
            return float('nan')

    def __int__(self):
        return int(self.value)

    def __cmp__(self, other):
        if self.value == other:
            return 0
        elif self.value > other:
            return 1
        elif self.value < other:
            return -1

    def __nonzero__(self):
        # python 2
        if self.dtype == bool:
            return bool(self.value)
        elif self.value is None:
            return False
        elif self.value == '':
            return False
        elif float(self) == float('nan'):
            return False
        else:
            return True

    def __bool__(self):
        # python 3
        if self.dtype == bool:
            return bool(self.value)
        elif self.value is None:
            return False
        elif self.value == '':
            return False
        elif float(self) == float('nan'):
            return False
        else:
            return True

    def __str__(self):
        if self.dtype == FloatExp:
            return str(self.value)
        elif self.dtype == Equation:
            return str(self.value)
        elif self.dtype == float:
            return str(self.value)
        elif self.dtype == int:
            return str(self.value)
        elif self.dtype == str and self.value is not None:
            if isinstance(self.value, str):
                self.value = self.value.replace('"','').replace("'",'')
            else:
                self.value = str(self.value)
            return ''.join(["'", str(self.value), "'"])
        elif self.dtype == bool:
            return ''.join([".", str(self.value), "."])
        else:
            return str(self.value)

    def _checkdtype(self):
        # Check data type
        for dtype in [float, int, bool, FloatExp, Equation]:
            if isinstance(self.value, dtype):
                self.dtype = dtype

        # If still None, assume string
        if self.dtype is None:
            self.dtype = str

    def updateValue(self, value):

        if value is None:
            self.value = None
        elif self.dtype == Equation and isinstance(value, str):
            self.value.eq = value
        elif self.regX_expression.findall(str(value)) and self.dtype == float:
            self.value = Equation(value)
        elif self.dtype == FloatExp and isinstance(value, float):
            self.value = FloatExp(value)
        else:
            self.value = value

        self._checkdtype()

    def lower(self):
        if isinstance(self.value, str):
            return self.value.lower()
        else:
            raise TypeError('Keyword is not a str')

class Base(object):
    def __init__(self, ind):
        self.ind = ind
        self._keywordDict = {}
        self.delete = False

    def __getattr__(self, name):
        if name in self._keywordDict:
            return self._keywordDict[name].value
        else:
            # Default behaviour
            raise AttributeError(name)

    def __getitem__(self, name):
        return self._keywordDict[name]

    def __setitem__(self, name, value):
        self._keywordDict[name] = value

    def __contains__(self, item):
        return item in self._keywordDict

    def addKeyword(self, key, value):
        self._keywordDict[key] = KeyWord(key, value)

    def deleteDict(self):
        self.delete = True
        return dict([(key, None) for key, val in self._keywordDict.items()])


class CondBase(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

        self.gasSpecies = SpeciesCollection()
        self.solids = SolidsCollection()

class BC(CondBase):
    def __init__(self, ind):
        CondBase.__init__(self, ind)
        self.name = 'BC'

    def __str__(self):
        try:
            bctype = self.bc_type
        except AttributeError:
            bctype = 'UNDEFINED_C'

        if self.gasSpecies:
            gasSpec = ['  Gas Species:'] + [''.join(['    ', val]) for val in
                                            self.gasSpecies._prettyPrintList()]
        else:
            gasSpec = []

        if self.solids:
            solids = ['  Solids:'] + self.solids._prettyPrintList()
        else:
            solids = []

        return '\n'.join(["Boundary Condition {}: {}".format(self.ind, bctype)]
                         + [''.join(['  ', str(key), ': ', str(value)]) for key, value in
                            self._keywordDict.items()] +
                            gasSpec + solids
                         )


class IC(CondBase):
    def __init__(self, ind):
        CondBase.__init__(self, ind)
        self.name = 'IC'

    def __str__(self):
        try:
            ictype = self.ic_type
        except AttributeError:
            ictype = 'UNDEFINED_C'
        return "Initial Condition {}: {}".format(self.ind, ictype)


class PS(CondBase):
    def __init__(self, ind):
        CondBase.__init__(self, ind)
        self.name = 'PS'

    def __str__(self):
        return "Point Source {}:".format(self.ind)


class IS(CondBase):
    def __init__(self, ind):
        CondBase.__init__(self, ind)
        self.name = 'IS'


class ConditionCollection(list):
    def __init__(self, condtype=None):
        list.__init__(self)
        self.condtype = condtype

    def __contains__(self, item):
        if item in [itm.ind for itm in self]:
            return True
        else:
            return False

    def __getitem__(self, item):
        for itm in self:
            if itm.ind == item:
                return itm

        if self.condtype == 'ic':
            self.append(IC(item))
        elif self.condtype == 'bc':
            self.append(BC(item))
        elif self.condtype == 'ps':
            self.append(PS(item))
        elif self.condtype == 'is':
            self.append(IS(item))
        else:
            raise Exception("Collection Type Not Set.")

        return self[item]

    def delete(self, itm):
        self.pop(self.index(itm))


class Species(Base):
    """
    Class to hold properties and functions for gas and solid species.

    Usage: Species(index, phase)
        index = int
        phase = 'g' for gas, or 's' for solid, or l for liquid.
    """
    def __init__(self, ind, phase='g'):
        Base.__init__(self, ind)
        self.phase = phase

    def __str__(self):
        if self.phase == 'g':
            p = 'Gas'
        else:
            p = "Solid"
        return '\n'.join(["{} Species: {}".format(p, self.ind), ] +
                         [''.join(['  ', str(key), ': ', str(value)])
                         for key, value in self._keywordDict.items()])


class Solid(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

        self.species = SpeciesCollection()
        self.name = 'Solid {}'.format(self.ind)

    def addSpecies(self, ind):
        return self.species.new(ind, phase='s')

    def _prettyPrintList(self):
        solidAttr = [''.join(['  ', str(key), ': ', str(value)])
                     for key, value in self._keywordDict.items()]

        lst = [self.name] + solidAttr +\
              ['  Species:'] + [''.join(['    ', val]) for val in
                                self.species._prettyPrintList()]
        return lst

    def __str__(self):

        return '\n'.join(self._prettyPrintList())


class Collection(list):
    def __init__(self):
        list.__init__(self)
        self.indStart = 1

    def __getattr__(self, name):
        if name in self:
            for item in [itm.name for itm in self]:
                return item
        else:
            # Default behaviour
            raise AttributeError(name)

    def __contains__(self, item):
        if item in [itm.ind for itm in self]:
            return True
        else:
            return False

    def __getitem__(self, item):
        for itm in self:
            if itm.ind == item:
                return itm

    def _checkind(self, ind):
        currentSet = [itm.ind for itm in self]
        if ind in currentSet:
            raise Exception("An index of {} already exsits".format(ind))

        if ind < self.indStart:
            raise Exception("An index of {} not allowed. \
                             Range starts at {}".format(ind, self.indStart))

        if ind is None:
            if len(currentSet) < 1:
                ind = 1
            else:
                full_set = set(xrange(1, max(currentSet) + 1))
                ind = sorted(full_set - currentSet)[0]

        return ind

    def clean(self):
        for itm in self:
            if hasattr(itm, '_keywordDict') and len(itm._keywordDict.keys())<1:
                self.pop(self.index(itm))


    def delete(self, itm):
        self.pop(self.index(itm))


class SpeciesCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None, phase='g'):
        ind = self._checkind(ind)
        self.append(Species(ind, phase))
        return self[ind]

    def _prettyPrintList(self):
        printDict = {}
        for spec in self:
            if hasattr(spec, 'species_alias_s'):
                printDict[spec.ind] = spec.species_alias_s
            elif hasattr(spec, 'species_alias_g'):
                printDict[spec.ind] = spec.species_alias_g
            else:
                printDict[spec.ind] = None

        return [' '.join([str(key), str(val)])
                for key, val in printDict.items()]

    def __str__(self):
        return '\n'.join(self._prettyPrintList())

class SpeciesEq(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

class SpeciesEqCollection(Collection):
    def __init__(self):
        Collection.__init__(self)
        self.indStart = 0

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(SpeciesEq(ind))
        return self[ind]

class SolidsCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(Solid(ind))
        return self[ind]

    def _prettyPrintList(self):
        printList = []
        for solid in self:
            printList += solid._prettyPrintList()

        return printList

    def __str__(self):
        return '\n'.join(self._prettyPrintList())

class LinearEquationSolver(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

class LinearEquationSolverCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(LinearEquationSolver(ind))
        return self[ind]

class SPX(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

class SPXCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(SPX(ind))
        return self[ind]

class VtkVar(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

class VtkVarCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(VtkVar(ind))
        return self[ind]

class VariableGridVar(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

class  VariableGridCollection(Collection):
    def __init__(self):
        Collection.__init__(self)
        self.indStart = 0

    def new(self, ind=None):
        ind = self._checkind(ind)
        self.append(VariableGridVar(ind))
        return self[ind]

class Project(object):
    '''
    Class for managing a mfix project
    '''
    def __init__(self, mfixDatFile=None):

        self.mfixDatFile = mfixDatFile
        self._keywordDict = {}
        self.lines = []
        self.line_numbers = {}

        # Regular Expressions
        self.regX_keyValue = re.compile(r'(\w+)(?:\(([\d, ]+)\))?\s*=\s*(.*?)(?=(!|$|\w+(\([\d, ]+\))?\s*=))')
        self.regX_float_exp = re.compile("^([+-]?[0-9]*\.?[0-9]+?(?:[eEdD][+-]?[0-9]+))$")
        self.regX_float = re.compile("^[+-]?[0-9]+\.([0-9]+)?$")
        self.regX_expression = re.compile('@\(([0-9.eEpiPI+-/*\(\))]+)\)')
        #self.regX_stringShortHand = re.compile("""(\d+\*["'].+?["'])""")
        self.regX_stringShortHand = re.compile("""([\d\.]+\*((["']+?.+?["']+?)|([\d\.]+)))""")

        self.__initDataStructure__()

        if self.mfixDatFile is not None:
            self.parsemfixdat()

    def __initDataStructure__(self):
        # Conditions
        self.ics = ConditionCollection('ic')
        self.bcs = ConditionCollection('bc')
        self.pss = ConditionCollection('ps')
        self.iss = ConditionCollection('is')

        # Species/solids
        self.gasSpecies = SpeciesCollection()
        self.solids = SolidsCollection()
        self.speciesEq = SpeciesEqCollection()

        # LEQ
        self.linearEq = LinearEquationSolverCollection()

        # SPX
        self.spx = SPXCollection()

        # VTK_VAR
        self.vtkvar = VtkVarCollection()

        # variablegrid
        self.variablegrid = VariableGridCollection()

    def __getattr__(self, name):
        return self._keywordDict[name]

    def __getitem__(self, key):
        return self._keywordDict.get(key, None)

    def __contains__(self, item):
        return item in self._keywordDict

    def save(self, fname):
        for key, value in self._keywordDict.iteritems():
            if type(value) == type(''):
                self.lines[self.line_numbers[key]] = '{} = {}\n'.format(key.upper(), value)
            elif type(value) == type(0.0):
                self.lines[self.line_numbers[key]] = '{} = {:f}\n'.format(key.upper(), value)
            elif type(value) == type(False):
                self.lines[self.line_numbers[key]] = '{} = {}\n'.format(key.upper(), '.T.' if value else '.F.')

        save_file = open(fname, mode='w')
        for line in self.lines:
            save_file.write(line)
        save_file.close()

    def parsemfixdat(self, fname=None):
        """
        Parse the mfix.dat file with regular expressions.

        fname can be a StringIO instance, path, or a plain string.

        Saves the results in self.mfixDatKeyDict dictionary.
        """
        if fname:
            self.mfixDatFile = fname

        assert self.mfixDatFile is not None

        # check to see if the file is a StringIO object
        if isinstance(self.mfixDatFile, StringIO):
            self._parsemfixdat(self.mfixDatFile)
            return

        # Try to open the file
        try:
            with open(self.mfixDatFile) as mfixDatFile:
                self._parsemfixdat(mfixDatFile)
            return
        except (IOError, OSError):
            pass
        # maybe it is just a string?
        if is_string(self.mfixDatFile) or is_unicode(self.mfixDatFile):
            self._parsemfixdat(StringIO(self.mfixDatFile))
            return

    def _parsemfixdat(self, fobject):
        """
        This does the actual parsing.
        """
        reactionSection = False
        self.comments = {}
        self._keywordDict = {}
        self.__initDataStructure__()

        for i, line in enumerate(fobject):
            self.lines.append(line)
            line = line.strip()
#            line = to_text_string(line, 'utf8').strip()

            if '@(RXNS)' in line:
                reactionSection = True
            elif '@(END)' in line and reactionSection:
                reactionSection = False
            elif not reactionSection:
                # remove comments
                if line.startswith('#') or line.startswith('!'):
                    self.comments[i] = line
                    line = ''
                elif '#' in line or '!' in line:
                    line, keywordComment = re.split('[#!]', line)[:2]

                for match in self.regX_keyValue.findall(line):
                    # match chould be: [keyword, args, value,
                    #                   nextKeywordInLine, something]

                    # keyword
                    key = match[0].lower().strip()
                    self.line_numbers[key] = i

                    # values
                    valString = match[2].strip()

                    # look for shorthand [count]*[value] and expand.
                    # Replace short hand string
                    shortHandList = self.regX_stringShortHand.findall(valString)
                    if shortHandList:
                        for shortHand in shortHandList:
                            # check for expression
                            if not re.findall('@\('+shortHand[0].replace('*','\*')+'\)',valString):
                                valString = valString.replace(shortHand[0],
                                                              self.expandshorthand(
                                                              shortHand[0])
                                                              )

                    # split values using shlex, it will keep quoted strings
                    # together.
                    try:
                        vals = shlex.split(valString.strip())
                    except ValueError:
                        if 'description' in key:
                            vals = [shlex.split(line.strip())[-1]]
                        else:
                            vals = []

                    # A bug with shlex and unicode adds \x00?
                    for i, val in enumerate(vals):
                        if '\x00' in val:
                            vals[i] = val.replace('\x00','')


                    # clean the values converting to python types
                    cleanVals = []
                    for val in vals:
                        val = self.cleanstring(val)
                        if val is not None:
                            cleanVals.append(val)

                    # I don't know if this is needed?
                    if not cleanVals:
                        cleanVals = [None]

                    # add args to legacy keys if no keys
                    if not match[1] and len(cleanVals) == 1 and \
                            key in ['d_p0', 'ro_s0']:
                        match = list(match)
                        match[1] = '1'

                    # If args or multiple values
                    if match[1] or len(cleanVals) > 1:

                        keyWordArgList = []
                        if match[1]:
                            args = [int(arg) for arg in match[1].split(',')]
                        else:
                            args = []

                        numVals = len(cleanVals)

                        if numVals > 1:
                            if args:
                                for i in range(0, numVals):
                                    keyWordArgList.append([i+args[0]]+args[1:])
                            else:
                                for i in range(1, numVals+1):
                                    keyWordArgList.append([i]+args[1:])
                        else:
                            keyWordArgList.append(args)

                        # Find condition keywords and separate
                        if key.startswith('ic_'):
                            cond = self.ics
                        elif key.startswith('bc_') and not key.endswith('_q'):
                            cond = self.bcs
                        elif key.startswith('ps_'):
                            cond = self.pss
                        elif key.startswith('is_'):
                            cond = self.iss
                        else:
                            cond = None

                        for keyWordArgs, cleanVal in zip(keyWordArgList, cleanVals):

                            # Save conditions
                            if cond is not None:
                                if len(keyWordArgs) == 1:
                                    cond[keyWordArgs[0]][key] = \
                                        KeyWord(key, cleanVal)

                                # Gas Keys
                                elif len(keyWordArgs) == 2 and\
                                        key.endswith('_g'):

                                    condItm = cond[keyWordArgs[0]]

                                    if keyWordArgs[1] not in \
                                            condItm.gasSpecies:
                                        spec = condItm.gasSpecies.new(
                                            keyWordArgs[1])
                                    else:
                                        spec = condItm.gasSpecies[
                                            keyWordArgs[1]]

                                    spec[key] = KeyWord(key, cleanVal)

                                # Solid Keys
                                elif key.endswith('_s') or key in ['ic_theta_m']:

                                    condItm = cond[keyWordArgs[0]]

                                    if keyWordArgs[1] not in condItm.solids:
                                        solid = condItm.solids.new(
                                            keyWordArgs[1])
                                    else:
                                        solid = condItm.solids[keyWordArgs[1]]


                                    if len(keyWordArgs) == 2:
                                        solid[key] = KeyWord(key, cleanVal)
                                    elif len(keyWordArgs) == 3:
                                        if keyWordArgs[2] not in solid.species:
                                            spec = solid.addSpecies(
                                                keyWordArgs[2])
                                        else:
                                            spec = solid.species[
                                                keyWordArgs[2]]

                                        spec[key] = KeyWord(key, cleanVal)
#                                    else:
#                                        print(key, keyWordArgs, cleanVals)
#
#                                else:
#                                    print(key, keyWordArgs, cleanVals)

                            # Solid Species
                            elif key in ['species_s', 'species_alias_s',
                                         'mw_s', 'd_p0', 'ro_s', 'nmax_s',
                                         'c_ps0', 'k_s0', 'x_s0', 'ro_xs0',
                                         'solids_model', 'close_packed',]:

                                if keyWordArgs[0] not in self.solids:
                                    solid = self.solids.new(keyWordArgs[0])
                                else:
                                    solid = self.solids[keyWordArgs[0]]

                                if len(keyWordArgs) == 1:
                                    solid[key] = KeyWord(key, cleanVal)
                                else:
                                    if keyWordArgs[1] not in solid.species:
                                        spec = solid.addSpecies(keyWordArgs[1])
                                    else:
                                        spec = solid.species[keyWordArgs[1]]
                                    spec[key] = KeyWord(key, cleanVal)

                            # Gas Species
                            elif key in ['species_g', 'species_alias_g',
                                         'mw_g']:
                                if keyWordArgs[0] not in self.gasSpecies:
                                    spec = self.gasSpecies.new(keyWordArgs[0])
                                else:
                                    spec = self.gasSpecies[keyWordArgs[0]]

                                spec[key] = KeyWord(key, cleanVal)

                            # Species_eq
                            elif key in ['species_eq']:
                                if keyWordArgs[0] not in self.speciesEq:
                                    leq = self.speciesEq.new(keyWordArgs[0])
                                else:
                                    leq = self.speciesEq[keyWordArgs[0]]

                                leq[key] = KeyWord(key, cleanVal)

                            # LEQ
                            elif key in ['leq_method', 'leq_tol', 'leq_it',
                                         'leq_sweep', 'leq_pc', 'ur_fac',
                                         ]:
                                if keyWordArgs[0] not in self.linearEq:
                                    leq = self.linearEq.new(keyWordArgs[0])
                                else:
                                    leq = self.linearEq[keyWordArgs[0]]

                                leq[key] = KeyWord(key, cleanVal)

                            # SPX
                            elif key in ['spx_dt']:
                                if keyWordArgs[0] not in self.spx:
                                    spx = self.spx.new(keyWordArgs[0])
                                else:
                                    spx = self.spx[keyWordArgs[0]]

                                spx[key] = KeyWord(key, cleanVal)

                            # VTK
                            elif key in ['vtk_var']:
                                if keyWordArgs[0] not in self.vtkvar:
                                    vtkvar = self.vtkvar.new(keyWordArgs[0])
                                else:
                                    vtkvar = self.vtkvar[keyWordArgs[0]]

                                vtkvar[key] = KeyWord(key, cleanVal)

                            # variable grid
                            elif key in ['cpx', 'ncx', 'erx', 'first_dx',
                                         'last_dx', 'cpy', 'ncy', 'ery',
                                         'first_dy', 'last_dy', 'cpz', 'ncz',
                                         'erz', 'first_dz', 'last_dz']:

                                if keyWordArgs[0] not in self.variablegrid:
                                    variablegrid = self.variablegrid.new(keyWordArgs[0])
                                else:
                                    variablegrid = self.variablegrid[keyWordArgs[0]]

                                variablegrid[key] = KeyWord(key, cleanVal)

                            # Save other keywords with one arg to dict
                            # {arg:value}
                            elif len(keyWordArgs) == 1:

                                if key not in self._keywordDict:
                                    self._keywordDict[key] = {}

                                if keyWordArgs[0] not in self._keywordDict[key]:
                                    self._keywordDict[key][keyWordArgs[0]] = \
                                        KeyWord(key, cleanVal)

                                elif isinstance(self._keywordDict[key][keyWordArgs[0]], KeyWord):
                                    self._keywordDict[key][keyWordArgs[0]].updateValue(cleanVal)
                                else:
                                    self._keywordDict[key][keyWordArgs[0]] = \
                                        KeyWord(key, cleanVal)
                            else:
                                print(key, keyWordArgs, cleanVals)
                    else:
                        self._keywordDict[key] = KeyWord(key, cleanVals[0])

    def addKeyword(self, key, value, dtype = None):
        self._keywordDict[key] = KeyWord(key, value, dtype=dtype)

    def cleanstring(self, string):
        """
        Attempt to clean strings of '," and convert .t., .f., .true., .false.
        to booleans, and catch math expressions i.e. @(3*3) and remove the @()
        """
        cleanVal = None
        string = string.replace("'", '').replace('"', '')
        if len(string) > 0:
            if '.t.' == string.lower() or '.true.' == string.lower():
                cleanVal = True
            elif '.f.' == string.lower() or '.false.' == string.lower():
                cleanVal = False
            elif self.regX_expression.findall(string):
                cleanVal = Equation(self.regX_expression.findall(string)[0])
            elif self.regX_float_exp.findall(string):
                cleanVal = FloatExp(string.lower().replace('d', 'e'))
            elif any([val.isdigit() for val in string]):
                try:
                    cleanVal = int(string)
                except ValueError:
                    try:
                        cleanVal = float(string)
                    except ValueError:
                        cleanVal = string
            else:
                cleanVal = string
        return cleanVal

    def expandshorthand(self, string):
        """
        Expand mfix.dat short hands:

        expandShortHand("key = 6*'ISIS'") = ['ISIS','ISIS','ISIS','ISIS',
                                             'ISIS','ISIS']
        """
        splitString = string.split('*')
        expandList = [string]
        if len(splitString) == 2:
            if splitString[0].isdigit():
                expandList = int(splitString[0])*[str(self.cleanstring(
                    splitString[1]))]
        return ' '.join(expandList)

    def cleanDeletedItems(self):
        '''
        Purge objects marked with self.delete==True.
        '''

        for condType in [self.ics, self.bcs, self.pss, self.iss]:
            for cond in condType:
                for gas in cond.gasSpecies:
                    if gas.delete:
                        cond.gasSpecies.delete(gas)
                for solid in cond.solids:
                    for species in solid.species:
                        if species.delete:
                            solid.species.delete(species)
                    if solid.delete:
                        cond.solids.delete(solid)

                if cond.delete:
                    condType.delete(cond)

        # Species/solids
        for species in self.gasSpecies:
            if species.delete:
                self.gasSpecies.delete(species)

        for solid in self.solids:
            for species in solid.species:
                if species.delete:
                    solid.species.delete(species)
            if solid.delete:
                self.solids.delete(solid)

        for specieseq in self.speciesEq:
            if specieseq.delete:
                self.speciesEq.delete(specieseq)

        # LEQ
        for lineq in self.linearEq:
            if lineq.delete:
                self.linearEq.delete(lineq)

        # SPX
        for spx in self.spx:
            if spx.delete:
                self.spx.delete(spx)

        # VTK_VAR
        for vtk in self.vtkvar:
            if vtk.delete:
                self.vtkvar.delete(vtk)

        # variablegrid
        for vargrid in self.variablegrid:
            if vargrid.delete:
                self.variablegrid.delete(vargrid)


def readMfixDat():
    mfixproj = Project('./mfixJordan.dat')

    print(mfixproj._keywordDict.keys())
    print(mfixproj.run_name)
    print(mfixproj['run_name'])

    print('coordinates' in mfixproj)

    for bc in mfixproj.bcs:
        print(bc)

    for ic in mfixproj.ics:
        print(ic)

    for ps in mfixproj.pss:
        print(ps)

    for solid in mfixproj.solids:
        print(solid)

    for species in mfixproj.gasSpecies:
        print(species)


def readStringAsFile():

    strFile = StringIO("""
                       key = 6*'value'
                       run_name = "test"
                       coordinates = 'cartesian'
                       """
                       )

    mfixproj = Project(strFile)
    print(mfixproj._keywordDict.keys())
    print(mfixproj.run_name)
    print(mfixproj['run_name'])
    print(mfixproj['coordinates'] == 'cylindrical')
    mfixproj.run_name = 10
    print(mfixproj.run_name)

def testEmpty():
    mfixproj = Project()

    'test' in mfixproj

def test_keyword():
    keyword = KeyWord('cool', True)

    if keyword:
        print(True)


if __name__ == '__main__':
    #readMfixDat()
    #readStringAsFile()
    #testEmpty()
    #test_keyword()

    test = """
           ic_ep_g = .4      1.0
           ic_x_w		=  2*1.0  @(10*2)
           IC_ROP_s(1,1)	=  @(0.6*2484.0)	0.0
#           ic_x_s(1,1,1) = 0.3
#           ic_t_s(2,1) = 300
#           ic_x_g(2,1) = 0.7
#           ic_x_g(2,2) = 0.3
#           species_eq(0) = .True.
#           species_eq(1) = .False.
#           spx_dt = 10*.01
           """

    proj = Project(test)

    keyB = KeyWord('test', True)
    keyS = KeyWord('test', 'test')
    keyF = KeyWord('test', 0.0)

    keyCopy = copy.deepcopy(keyF)

    if keyF:
        print('keyF has a value')


    keyEq = KeyWord('test', None, dtype = Equation)

    if keyEq:
        print('keyEq has a value')

    for ic in proj.ics:
#        if 'ic_ep_g' in ic:
#            print(ic.ind, ic['ic_ep_g'])
        if 'ic_x_w' in ic:
            print(ic.ind, ic['ic_x_w'])

    print(Equation('5*2'))
    print(float(Equation('2*5')))
