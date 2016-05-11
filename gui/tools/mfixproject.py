# -*- coding: utf-8 -*-
"""
This file is part of the pymfix library

Licence
-------
As a work of the United States Government, this project is in the public domain
within the United States. As such, this code is licensed under
CC0 1.0 Universal public domain.

Please see the LICENSE.md for more information.
"""

# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

# Python core imports
import shlex
import re
import math
import copy
import warnings
try:
    # Python 2.X
    from StringIO import StringIO
except ImportError:
    # Python 3.X
    from io import StringIO

import logging
log = logging.getLogger()

# local imports
from tools.simpleeval import simple_eval
from tools.general import (recurse_dict, recurse_dict_empty, get_from_dict,
                           to_unicode_from_fs, is_string, is_unicode)


NaN = float('NaN')

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
        if len(self.eq) == 0:
            return 0
        elif self.eq is None or self.eq == 'None':
            return NaN
        else:
            try:
                return simple_eval(self.eq.lower(),
                                   names={"pi": math.pi, "e": math.e, "p": 1})
            except SyntaxError:
                return 0

    def __nonzero__(self):  # Python 2
        return not math.isnan(self._eval())

    __bool__ = __nonzero__  # Python 3

    def __float__(self):
        return float(self._eval())

    def __int__(self):
        # Will raise ValueError if equation evaluates to NaN
        return int(self._eval())

    def __repr__(self):
        return ''.join(['@(', str(self.eq), ')'])

    def __add__(self, value):
        return float(self._eval()) + float(value)

    def __sub__(self, value):
        return float(self._eval()) - float(value)

    def __mul__(self, value):
        return float(self._eval()) * float(value)

    def __pow__(self, value):
        return float(self._eval()) ** float(value)


def format_key_with_args(key, args=None):
    if args:
        return "%s(%s)" % (key, ','.join(str(a) for a in args))
    else:
        return str(key)

class Keyword(object):
    def __init__(self, key, val, comment='', dtype=None, args=None):

        self.key = key
        self.value = val
        self.dtype = dtype
        self.min = None
        self.max = None
        self.valids = None
        self.fmt = None
        self.comment = comment
        if args is None:
            args = []
        self.args = args

        # FIXME:  this is not just matching 'pi' but any jumble of p's and i's
        # see also def in Project.__init__
        self.re_expression = re.compile('([eEpiPI\+\-/*\^\(\)]+)')

        if dtype is None:
            self._checkdtype()

    def __deepcopy__(self, memo):
        return Keyword(copy.copy(self.key),
                       copy.copy(self.value),
                       comment=copy.copy(self.comment),
                       dtype=copy.copy(self.dtype),
                       args=copy.copy(self.args),
                       )

    def __float__(self):
        try:
            return float(self.value)
        except (ValueError, TypeError, ZeroDivisionError):
            return NaN

    def __int__(self):
        return int(self.value)

    def __cmp__(self, other):
        if self.value == other:
            return 0
        elif self.value > other:
            return 1
        elif self.value < other:
            return -1

    def __nonzero__(self): # python 2
        if self.dtype == bool:
            return bool(self.value)
        elif self.value is None:
            return False
        elif self.value == '':
            return False
        elif math.isnan(float(self)):
            return False
        else:
            return True

    __bool__ = __nonzero__ # python 3

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
                self.value = self.value.replace('"', '').replace("'", '')
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

    def line(self):
        if len(self.args) == 0:
            line = '  {} = {}'.format(self.key, str(self))
        else:
            line = '  {}({}) = {}'.format(self.key,
                                          ','.join([str(x) for x in self.args]),
                                          str(self))

        if len(self.comment) > 0:
            line = '    !'.join([line, self.comment])

        return line

    def updateValue(self, value):
        if value is None:
            self.value = None
        elif self.dtype == Equation and isinstance(value, str):
            self.value.eq = value
        elif self.re_expression.findall(str(value)) and self.dtype == float:
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
        self.ind = ind # index
        self._keyword_dict = {}
        self.delete = False

    def __getattr__(self, name):
        if name in self._keyword_dict:
            return self._keyword_dict[name].value
        else:
            # Default behaviour
            raise AttributeError(name)

    def __getitem__(self, name):
        return self._keyword_dict[name]

    def __setitem__(self, name, value):
        if name in self._keyword_dict and isinstance(self._keyword_dict[name],
                                                    Keyword):
            if isinstance(value, Keyword):
                self._keyword_dict[name].updateValue(value.value)
            else:
                self._keyword_dict[name].updateValue(value)
        else:
            self._keyword_dict[name] = value

    def __contains__(self, item):
        return item in self._keyword_dict

    def __len__(self):
        return len(self._keyword_dict)

#    def updateKeyword(self, key, value, args=[]):
#        self._keyword_dict[key] = Keyword(key, value, args=args)

#    def deleteDict(self):
#        self.delete = True
#        return dict([(key, None) for key, val in self._keyword_dict.items()])


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
                            self._keyword_dict.items()] +
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
                         for key, value in self._keyword_dict.items()])


class Solid(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)

        self.species = SpeciesCollection()
        self.name = 'Solid {}'.format(self.ind)

    def addSpecies(self, ind):
        return self.species.new(ind, phase='s')

    def _prettyPrintList(self):
        solidAttr = [''.join(['  ', str(key), ': ', str(value)])
                     for key, value in self._keyword_dict.items()]

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
        for itm in self: # XX FIXME - O(n^2)
            if itm.ind == item:
                return itm

    def _checkind(self, ind):
        currentSet = [itm.ind for itm in self]
        if ind in currentSet:
            raise Exception("An index of {} already exists".format(ind))

        if ind < self.indStart:
            raise Exception("An index of {} not allowed. "
                            "Range starts at {}".format(ind, self.indStart))

        if ind is None:
            if len(currentSet) < 1:
                ind = 1
            else:
                full_set = set(range(1, max(currentSet) + 1))
                ind = sorted(full_set - currentSet)[0]

        return ind

    def clean(self):
        for itm in self:
            if hasattr(itm, '_keyword_dict') and \
                    len(itm._keyword_dict.keys()) < 1:
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
    def __init__(self, dat_file=None, keyword_doc=None):

        self.dat_file = dat_file
        self.keyword_doc = keyword_doc
        self._keyword_dict = {}
        self.dat_file_list = []
        self.thermoindex = None

        # Regular Expressions
        self.re_keyValue = re.compile(r"""
            (\w+)                           # Alphanumeric, key name
            (?:\(([\d, ]+)\))?              # Indices, :(NUM,)
            \s*                             # Possible whitespace
            =                               # Equal sign
            \s*                             # Possible whitespace
            (.*?)                           # Value
            (?=(!|$|\w+(\([\d, ]+\))?\s*=)) # comment ?  further keywords ?
        """, re.VERBOSE)

        self.re_float_exp = re.compile(r"""
            ^([+-]?[0-9]*\.?[0-9]+?(?:[eEdD][+-]?[0-9]+))$
        """, re.VERBOSE)

        self.re_float = re.compile(r"""
            ^[+-]?[0-9]+\.([0-9]+)?$
        """, re.VERBOSE)

        self.re_expression = re.compile(r"""
            @\(([ 0-9.eEpiPI+-/*\(\))]+)\)
        """, re.VERBOSE)

        self.re_stringShortHand = re.compile(r"""
            [\d]+                                      # unsigned integer
            \*                                         # literal *
            (?:                                        # non-capturing group
              '[\w]+'                                  # single-quoted word
            | "[\w]+"                                  # double-quoted word
            | \.[\w]+\.                                # .TOKEN.
            | [-+]?[\d]+\.?[\d]*(?:[dDeE][-+]?[\d]+)?  # number
            )                                          # end group
        """, re.VERBOSE)

        self.__initDataStructure__()

        if self.dat_file is not None:
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
        return self._keyword_dict[name]

    def __getitem__(self, key):
        if not isinstance(key, list) and not isinstance(key, tuple):
            key = [key]
        return get_from_dict(self._keyword_dict, key)

    def __contains__(self, item):
        try:
            if not isinstance(item, list) and not isinstance(item, tuple):
                item = [item]
            get_from_dict(self._keyword_dict, item)
            return True
        except KeyError:
            return False

    def __deepcopy__(self, memo):
        # TODO: this is not efficient
        return Project(''.join(self.convertToString()))

    def parsemfixdat(self, fname=None):
        """Parse the mfix.dat file with regular expressions.
        fname can be a StringIO instance, path, or a plain string.
        Saves the results in self.mfixDatKeyDict dictionary.
        """
        if fname:
            self.dat_file = fname

        assert self.dat_file is not None

        # check to see if the file is a StringIO object
        if isinstance(self.dat_file, StringIO):
            self._parsemfixdat(self.dat_file)
            return

        # Try to open the file
        try:
            with open(self.dat_file) as dat_file:
                self._parsemfixdat(dat_file)
            return
        except (IOError, OSError):
            pass
        # maybe it is just a string?
        if is_string(self.dat_file) or is_unicode(self.dat_file):
            self._parsemfixdat(StringIO(self.dat_file))
            return

    def parseKeywordLine(self, line):
        matches = self.re_keyValue.findall(line)
        single_key = False
        if matches:
            for match in matches:
                # match could be: [keyword, args, value,
                #                   nextKeywordInLine, something]

                # convert to list
                match = list(match)

                # keyword
                key = match[0].lower().strip()

                # values
                val_string = match[2].strip()

                # remove spaces from equations: @( 2*pi)
                exps = self.re_expression.findall(val_string)
                for exp in exps:
                    val_string = val_string.replace(exp, exp.replace(' ',''))

                # look for shorthand [count]*[value] and expand.
                # Don't expand shorthand inside $( or 'description' field
                if not (val_string.startswith('@(') or 'description' in key):
                    for shorthand in self.re_stringShortHand.findall(val_string):
                        val_string = val_string.replace(shorthand,
                                                        self.expandshorthand(shorthand))
                # split values using shlex, it will keep quoted strings
                # together.
                if 'description' in key:
                    vals = [shlex.split(line.strip())[-1]]
                    single_key = True
                else:
                    try:
                        vals = shlex.split(val_string.strip())
                    except ValueError:
                        vals = []

                # clean the values converting to python types
                cleanVals = []
                for val in vals:
                    val = self.cleanstring(val)
                    if val is not None:
                        cleanVals.append(val)

                # check keyword_doc to see if the keyword is an array
                if (self.keyword_doc is not None and key in self.keyword_doc
                        and not match[1] and self.keyword_doc[key]['args']):
                    match[1] = ','.join(str(arg['min']) for arg in
                                        self.keyword_doc[key]['args'].values())

                # clean up arguments
                if match[1]:
                    args = [int(arg) for arg in match[1].split(',')]
                else:
                    args = []

                # If multiple values, split apart into multiple key, args, values
                if len(cleanVals) > 1:
                    keywordArgList = []

                    numVals = len(cleanVals)

                    if numVals > 1:
                        if args:
                            for val in range(0, numVals):
                                keywordArgList.append([val+args[0]]+args[1:])
                        else:
                            # hack for species eq (why?)
                            start = 0 if key == 'species_eq' else 1
                            for val in range(start, numVals+1):
                                keywordArgList.append([val]+args[1:])
                    else:
                        keywordArgList.append(args)

                    for keywordArgs, cleanVal in zip(keywordArgList, cleanVals):
                        yield (key, keywordArgs, cleanVal)

                else:
                    yield (key, args, cleanVals[0])

                if single_key:
                    break # no more keywords on this line
        else:
            yield (None, None, None)

    def _parsemfixdat(self, fobject):
        """
        This does the actual parsing.
        """
        reactionSection = False
        self.comments = {}
        self._keyword_dict = {}
        self.__initDataStructure__()
        self.dat_file_list = []

        for i, line in enumerate(fobject):
            line = to_unicode_from_fs(line).strip('\n')
            if '@(RXNS)' in line:
                reactionSection = True
            elif '@(END)' in line and reactionSection:
                reactionSection = False
            elif 'thermo data' in line.lower():
                self.thermoindex = i
                self.dat_file_list.append(line)
            elif not reactionSection:
                # remove comments
                commentedline = ''
                if line.startswith('#') or line.startswith('!'):
                    commentedline = line
                    line = ''
                elif '#' in line or '!' in line:
                    line, keywordComment = re.split('[#!]', line, maxsplit=1)
                else:
                    keywordComment = ''

                # loop through all keywords in the line
                try:
                    for (key, args, value) in self.parseKeywordLine(line):
                        if key is None:
                            self.dat_file_list.append(line+commentedline)
                        else:
                            try:
                                self.updateKeyword(key, value, args, keywordComment)
                            except ValueError:
                                # error at line i
                                warnings.warn("Cannot set %s=%s" % (format_key_with_args(key, args), value))
                except Exception as e:
                    warnings.warn("Parse error: %s: line %d, %s" % (e, i, line))

    def updateKeyword(self, key, value, args=None,  keywordComment=''):
        '''
        Update or add a keyword to the project.  Raises ValueError if there is a
        problem with the key or value.

        Parameters
        ----------
        key (str):
            the keyword to be added
        value (int, float, bool, str):
            the value of the keyword
        args (list):
            list of arguments for the keyword, or None
        keywordComment (str):
            a comment to be included with the keyword (default: '')
        '''
        # TODO:  refactor
        if args is None:
            args = []
        # check to see if the keyword already exists
        if [key]+args in self:
            self[[key]+args].updateValue(value)
            if keywordComment:
                self[[key]+args].comment = keywordComment
            return self[[key]+args]

        keywordobject = None

        if args:
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
            # Save conditions
            if cond is not None:
                if len(args) == 1: # 1-dimensional
                    keywordobject = Keyword(key, value, args=args,
                                            comment=keywordComment)
                    cond[args[0]][key] = keywordobject
                # Gas Keys
                elif len(args) == 2 and key.endswith('_g'):
                    condItm = cond[args[0]]
                    if args[1] not in condItm.gasSpecies:
                        spec = condItm.gasSpecies.new(args[1])
                    else:
                        spec = condItm.gasSpecies[args[1]]
                    keywordobject = Keyword(key, value, args=args,
                                            comment=keywordComment)
                    spec[key] = keywordobject
                # Solid Keys
                elif key.endswith('_s'):
                    condItm = cond[args[0]]
                    if args[1] not in condItm.solids:
                        solid = condItm.solids.new(
                            args[1])
                    else:
                        solid = condItm.solids[args[1]]
                    if len(args) == 2:
                        keywordobject = Keyword(key, value, args=args,
                                                comment=keywordComment)
                        solid[key] = keywordobject
                    elif len(args) == 3:
                        if args[2] not in solid.species:
                            spec = solid.addSpecies(args[2])
                        else:
                            spec = solid.species[args[2]]
                        keywordobject = Keyword(key, value, args=args,
                                                comment=keywordComment)
                        spec[key] = keywordobject
                else:
                    # Create keyword objects for other keys not handled above,
                    # eg IC_THETA_M(IC, Phase),
                    # (see 'save everything else' comment below)
                    # TODO - is there more to do here?
                    keywordobject = Keyword(key, value, args=args,
                                            comment=keywordComment)
                    log.debug("Created keyword object for",
                              format_key_with_args(key, args))
            # Solid Species
            elif key in ['species_s', 'species_alias_s', 'mw_s', 'd_p0',
                         'ro_s', 'nmax_s', 'c_ps0', 'k_s0', 'x_s0', 'ro_xs0',
                         'solids_model', 'close_packed', ]:
                if args[0] not in self.solids:
                    solid = self.solids.new(args[0])
                else:
                    solid = self.solids[args[0]]
                if len(args) == 1:
                    keywordobject = Keyword(key, value, args=args,
                                            comment=keywordComment)
                    solid[key] = keywordobject
                else:
                    if args[1] not in solid.species:
                        spec = solid.addSpecies(args[1])
                    else:
                        spec = solid.species[args[1]]
                    keywordobject = Keyword(key, value, args=args,
                                            comment=keywordComment)
                    spec[key] = keywordobject
            # Gas Species
            elif key in ['species_g', 'species_alias_g', 'mw_g']:
                if args[0] not in self.gasSpecies:
                    spec = self.gasSpecies.new(args[0])
                else:
                    spec = self.gasSpecies[args[0]]
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                spec[key] = keywordobject
            # Species_eq
            elif key in ['species_eq']:
                if args[0] not in self.speciesEq:
                    leq = self.speciesEq.new(args[0])
                else:
                    leq = self.speciesEq[args[0]]
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                leq[key] = keywordobject
            # LEQ
            elif key in ['leq_method', 'leq_tol', 'leq_it', 'leq_sweep',
                         'leq_pc', 'ur_fac', ]:
                if args[0] not in self.linearEq:
                    leq = self.linearEq.new(args[0])
                else:
                    leq = self.linearEq[args[0]]
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                leq[key] = keywordobject
            # SPX
            elif key in ['spx_dt']:
                if args[0] not in self.spx:
                    spx = self.spx.new(args[0])
                else:
                    spx = self.spx[args[0]]
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                spx[key] = keywordobject
            # VTK
            elif key in ['vtk_var']:
                if args[0] not in self.vtkvar:
                    vtkvar = self.vtkvar.new(args[0])
                else:
                    vtkvar = self.vtkvar[args[0]]

                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                vtkvar[key] = keywordobject
            # variable grid
            elif key in ['cpx', 'ncx', 'erx', 'first_dx', 'last_dx', 'cpy',
                         'ncy', 'ery', 'first_dy', 'last_dy', 'cpz', 'ncz',
                         'erz', 'first_dz', 'last_dz']:
                if args[0] not in self.variablegrid:
                    variablegrid = self.variablegrid.new(args[0])
                else:
                    variablegrid = self.variablegrid[args[0]]
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
                variablegrid[key] = keywordobject
            # Save everything else
            else:
                keywordobject = Keyword(key, value, args=args,
                                        comment=keywordComment)
        else: # no args
            keywordobject = Keyword(key, value, args=None,
                                    comment=keywordComment)

        # add keyword to other data structures
        if keywordobject is not None:
            if self.thermoindex is not None:
                self.dat_file_list.insert(self.thermoindex, keywordobject)
                self.thermoindex += 1
            else:
                self.dat_file_list.append(keywordobject)
            self._recursiveAddKeyToKeywordDict(keywordobject,
                                               [key]+args)
        else:
            raise ValueError(format_key_with_args(key, args))

        return keywordobject

    def _recursiveAddKeyToKeywordDict(self, keyword, keys=()):
        keywordDict = self._keyword_dict
        for key in keys[:-1]:
            keywordDict = keywordDict.setdefault(key, {})
        keywordDict[keys[-1]] = keyword

    def _recursiveRemoveKeyToKeywordDict(self, keys=(), warn=True):
        keywordDict = self._keyword_dict
        for key in keys[:-1]: # ? Why ?
            keywordDict = keywordDict[key]
        keywordDict.pop(keys[-1])

    def _purgeemptydicts(self):
        for keys, value in list(recurse_dict_empty(self._keyword_dict)):
            if isinstance(value, dict) and not value:
                keywordDict = self._keyword_dict
                for key in keys[:-1]:
                    keywordDict = keywordDict[key]
                keywordDict.pop(keys[-1])

    def keywordItems(self):
        for keys, value in recurse_dict_empty(self._keyword_dict):
            yield value

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
            elif self.re_expression.findall(string):
                cleanVal = Equation(self.re_expression.findall(string)[0])
            elif self.re_float_exp.findall(string):
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

    def expandshorthand(self, shorthand):
        """Expand mfix.dat shorthand:
        expandShortHand("5*'IKE'") = "'IKE' 'IKE' 'IKE' 'IKE' 'IKE'"
        """
        count, word = shorthand.split('*', 1)
        count = int(count)
        ret =' '.join(count * [word])
        return ret

    def removeKeyword(self, key, args=None, warn=True):
        '''Remove a keyword from the project.  Returns True if item deleted.'''

        if args is None:
            args = []
        keyword = self.keywordLookup(key, args, warn=False)
        if not keyword:
            if warn:
                raise KeyError(key)
            else:
                return False

        # pop from dat_file_list
        self.dat_file_list.remove(keyword)
        if self.thermoindex is not None:
            self.thermoindex -= 1

        # remove from dict
        self._recursiveRemoveKeyToKeywordDict([key]+args, warn)
        for i in range(len(args)):
            self._purgeemptydicts()

        # purge
        keyword.delete = True
        self._cleanDeletedItems()
        return True

    def _cleanDeletedItems(self):
        '''Purge objects marked with self.delete==True.'''

        for condType in [self.ics, self.bcs, self.pss, self.iss]:
            for cond in condType:
                for gas in cond.gasSpecies:
                    if gas.delete:
                        cond.gasSpecies.delete(gas)
                        self.dat_file_list.pop(gas)
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

    def keywordLookup(self, keyword, args=None, warn=True):
        '''
        Search the project for a keyword and return the Keyword object, else
        raise an exception.

        Parameters
        ----------
        keyword (str):
            a keyword to search for
        args (list):
            a list of args to search for

        Returns
        -------
        keyword (pymfix.Keyword):
            the Keyword object
        '''

        if args is None:
            args = []

        keyword = keyword.lower()
        keywordobject = None

        for key, value in recurse_dict(self._keyword_dict):
            if key[0] == keyword:
                if args and value.args == args:
                    keywordobject = value
                    break
                elif not args:
                    keywordobject = value
                    break
                else:
                    continue

        if keywordobject is None:
            if warn:
                raise KeyError(keyword)

        return keywordobject

    def changekeywordvalue(self, key, value, args=None):
        '''
        Change a value of a keyword.

        Parameters
        ----------
        key (str):
            keyword to change
        value (int, float, bool, str):
            the value to set the keyword to
        args (list):
            a list of arguments for that keyword
        '''

        if args is None:
            args = []

        keywordobject = self.keywordLookup(key, args)

        keywordobject.updateValue(value)

        return keywordobject

    def convertToString(self):
        for line in self.dat_file_list:
            if hasattr(line, 'line'):
                yield line.line()+'\n'
            else:
                yield str(line)+'\n'

    def writeDatFile(self, fname):
        '''
        Write the project to a text file.

        Parameters
        ----------
        fname (str):
            the file name to write the project to
        '''

        with open(fname, 'w') as dat_file:
            for line in self.convertToString():
                dat_file.write(line)


if __name__ == '__main__':
    proj = Project()
    print(list(proj.parseKeywordLine('key = @( 2* 10)')))
