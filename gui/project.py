"""

This module defines Project and Keyword classes.

A Project holds Keyword settings representing an mfix project.
It can read/write files in 'mfix.dat' format

Keyword values are coerced to lower-case.

Keywords can take arguments, which are typically indices representing
phase number, eg   momentum_x_eq(0)  for phase 0.   Some keywords take
multiple arguments.

Keyword objects have a 'key' (the string) and a 'value'.

There are three paths to the Keyword objects in the Project:
 dict, list, and Collections

 list for preserving the order for writing
 dict for programmatically changing keywords
 collections for user interaction with the keywords


History
-------

This file was derived from the pymfix library
License
-------
As a work of the United States Government, this project is in the public domain
within the United States. As such, this code is licensed under
CC0 1.0 Universal public domain.

Please see the LICENSE.md for more information.
"""


# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

# Python core imports
import sys
import math
import warnings
import traceback
from collections import OrderedDict
try:
    # Python 2.X
    from StringIO import StringIO
except ImportError:
    # Python 3.X
    from io import StringIO

import logging
log = logging.getLogger(__name__)

# local imports
from tools.simpleeval import simple_eval
from tools.comparable import Comparable
from tools.general import (recurse_dict, recurse_dict_empty, get_from_dict,
                           to_unicode_from_fs, to_fs_from_unicode,
                           is_text_string, to_text_string,
                           safe_shlex_split, format_key_with_args)

from regexes import *
from constants import *

NaN = float('NaN')

class FloatExp(float):
    fmt = '4'
    def __repr__(self):
        return '{:.{}e}'.format(self, self.fmt)

    __str__ = __repr__

def make_FloatExp(val): # Handle Fortran formats
    # (we can't do this in FloatExp.__init__)
    try:
        return FloatExp(val)
    except ValueError:
        val = val.lower().replace('d', 'e')
        if val.endswith('e'):
            val += "+0"
        return FloatExp(val)


def clean_string(string):
    """Attempt to clean strings of '," and convert .t., .f., .true., .false.
    to booleans, and catch math expressions i.e. @(3*3) and remove the @()"""

    if not string:
        return ''

    # remove all quotes
    string = string.replace("'", '').replace('"', '')
    # remove any leading or trailing space, after removing quotes
    string = string.strip()
    # lower-case version of string
    s_low = string.lower()

    if s_low in ('.t.',  '.true.'):
        return True
    elif s_low in ('.f.', '.false.'):
        return False

    # Look for @() expression
    if string.startswith('@(') and string.endswith(')'):
        return Equation(string)

    # Look for exponential-notation
    match = re_float_exp.match(string)
    if match:
        return make_FloatExp(string)

    # Maybe it's a number (would a regex be better?)
    if any(char.isdigit() for char in string):
        try:
            return int(string)
        except ValueError:
            try:
                return float(string)
            except ValueError:
                pass

    # default - return string unchanged
    return string

def expand_shorthand(string):
    """Expand mfix.dat shorthand:
    expand_shorthand("bill 5*'IKE' fred") = "bill 'IKE' 'IKE' 'IKE' 'IKE' 'IKE' fred"
    But don't do it inside parenthesized expressions.
    """
    # TODO:  we could do better about handling quoted strings
    # and escapes
    def _expand(s): # Inner function to do expansion, once we've handled quoting and parens
        for shorthand in re_shorthand.findall(s):
            count, word = shorthand.split('*', 1)
            count = int(count)
            expansion =' '.join(count * [word])
            s = s.replace(shorthand, expansion)
        return s
    paren_level = 0
    ret = []
    part = ''
    for c in string:
        if c == '(':
            paren_level += 1
            if paren_level == 1: # Just entered paren. context
                ret.append(_expand(part))
                part = ''
        if c == ')':
            paren_level -= 1
            if paren_level == 0:
                ret.append(part+c)
                part = ''
        else:
            part += c
    if part: # leftover at end
        ret.append(_expand(part))
    ret = ''.join(part for part in ret if part)
    return ret


def remove_spaces_from_equations(string):
    # Can't do this with regex, here's a simple lexical scanner.
    quote = None
    escape = False
    paren_level = 0
    in_eq = False
    ret = ''
    for c in string:
        if c == '\\':
            escape = not escape
        else:
            if quote and c==quote and not escape:
                quote = False
            elif c == '"' or c == "'": # Start of quoted text
                                # We're not doing triple-strings!
                quote = c
            elif not (quote or escape):
                if c == '@':
                    in_eq = True
                if c == '(':
                    paren_level += 1
                if c == ')':
                    paren_level -= 1
                    if paren_level == 0:
                        in_eq = False

        if not(in_eq and c.isspace()):
            ret += c

    return ret


class Equation(object):
    """represents a simple arithmetic expression, which can be evaluated
    by simple_eval, in a namespace which includes the math constants 'e'
    and 'pi'.  Calling "float" causes evaluation, which will result in
    ValueError on any failures.  Evaluating an empty string returns 0.0,
    while evaluating None returns NaN"""

    def __init__(self, eq, dtype=float):
        eq = str(eq)
        while eq.startswith("@(") and eq.endswith(")"):
            eq = eq[2:-1]
        self.eq = eq
        self.dtype = dtype

    def get_used_parameters(self):
        av_params = PARAMETER_DICT.keys()
        eq = re.split('[\*\/\-\+ \(\)]', self.eq)
        return [p for p in av_params if p in eq]

    def _eval(self):
        if len(self.eq) == 0:
            return 0
        elif self.eq is None or self.eq == 'None':
            return NaN
        else:
            try:
                return float(simple_eval(self.eq.lower(),
                                         names=PARAMETER_DICT
                                         ))
            except:
                raise ValueError(self.eq)

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

    def dumps(self):
        return '%s #!MFIX-GUI eq{%s}' % (self.dtype(self._eval()), self.eq)

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


class Keyword(Comparable):
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

        if dtype is None:
            self._update_dtype()

    def __float__(self):
        try:
            return float(self.value)
        except (ValueError, TypeError, ZeroDivisionError):
            return NaN

    def __int__(self):
        return int(self.value)

    def _cmpkey(self):
        return self.value

    def __nonzero__(self): # python 2
        if self.dtype == bool:
            return bool(self.value)
        elif self.value is None:
            return False
        elif self.value == '':
            return False
        elif math.isnan(float(self)): # Failed expression evaluation.
            return False
        else:
            return True

    __bool__ = __nonzero__ # python 3

    def __str__(self):
        sval = to_text_string(self.value)
        if self.dtype == FloatExp:
            return sval
        elif self.dtype == Equation:
            return sval
        elif self.dtype == float:
            return sval
        elif self.dtype == int:
            return sval
        elif self.dtype == bool:
            return '.%s.' % sval
        elif self.dtype == str and self.value is not None:
            # this is weird, we should not be modifying self.value in a __str__ method
            #if is_text_string(self.value):
            #    self.value = self.value.replace('"', '').replace("'", '')
            #else:
            #    self.value = to_text_string(self.value)
            sval = sval.replace('"','').replace("'",'')
            return "'%s'" % sval
        else:
            return sval

    def _update_dtype(self):
        # Check data type & modify if needed
        for dtype in (bool, float, int, FloatExp, Equation):
            if isinstance(self.value, dtype):
                self.dtype = dtype
                break

        # If still None, assume string
        if self.dtype is None:
            self.dtype = str

    def line(self):
        if self.dtype == Equation:
            val = self.value.dumps()
        else:
            val = to_text_string(self)

        if len(self.args) == 0:
            line = '  {} = {}'.format(self.key, val)
        else:
            line = '  {}({}) = {}'.format(self.key,
                                          ','.join([to_text_string(x) for x in self.args]),
                                          val)

        if len(self.comment) > 0:
            line = '    !'.join([line, self.comment])

        return line

    def updateValue(self, value):
        sval = to_text_string(value)
        if value is None or sval=='': # is this the right place to check this?
            self.value = None
        elif self.dtype == Equation and is_text_string(value):
            self.value.eq = value
        elif (self.dtype == float and not re_float.match(sval) and not re_int.match(sval)):
            if re_float_exp.match(sval):
                self.value = make_FloatExp(value)
            else:
                eq = Equation(value)
                try:
                    f = float(eq)
                    self.value = eq
                except ValueError:
                    pass # Don't update invalid equation

        elif self.dtype == FloatExp and (isinstance(value, (int, float)
                                                    and not isinstance(value, bool))):
            self.value = make_FloatExp(value)
        else:
            self.value = value
        # TODO> make sure we don't change to invalid type
        self._update_dtype()

    def lower(self):
        return self.value.lower()


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

    def get(self, key, default=None):
        # Note, this only works with dynamic attributes, not static ones defined
        # in subclasses of Base (eg Solid.name)
        d =  self._keyword_dict.get(key)
        return default if d is None else d.value

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

        return '\n'.join(
            ["Boundary Condition %s: %s"%(self.ind, bctype)]
            + [''.join(['  ', to_text_string(key), ': ', to_text_string(value)])
               for key, value in self._keyword_dict.items()]
            + gasSpec
            + solids)



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

    #def get(self, key, default=None):
    #    return self._keyword_dict.get(key, default)

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

    #def get(self, key, default=None):
    #    return self._keyword_dict.get(key, default)

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

    def _check_ind(self, ind):
        # Check index for valid range
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
        ind = self._check_ind(ind)
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
        ind = self._check_ind(ind)
        self.append(SpeciesEq(ind))
        return self[ind]


class SolidsCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._check_ind(ind)
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
        ind = self._check_ind(ind)
        self.append(LinearEquationSolver(ind))
        return self[ind]


class SPX(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)


class SPXCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._check_ind(ind)
        self.append(SPX(ind))
        return self[ind]


class VtkVar(Base):
    def __init__(self, ind):
        Base.__init__(self, ind)


class VtkVarCollection(Collection):
    def __init__(self):
        Collection.__init__(self)

    def new(self, ind=None):
        ind = self._check_ind(ind)
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
        ind = self._check_ind(ind)
        self.append(VariableGridVar(ind))
        return self[ind]


class Project(object):
    """holds keywords and thermodynamic data for an MFIX project,
    reads and writes project file"""

    def __init__(self, dat_file=None, keyword_doc=None):

        self.dat_file = dat_file
        self.keyword_doc = keyword_doc
        self.comments = {}
        self._keyword_dict = {}
        self.dat_file_list = [] # contains the project file, lines are replaced with
                            # keywords as parsed
        self.thermo_data =  []
        self.mfix_gui_comments = OrderedDict() # lines starting with #!MFIX-GUI
        self.parameter_key_map = {}  #data structure to hold parameter->keyword mapping
        # See also 'reset'

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

    def get_value(self, key, default=None, args=None):
        # move to gui.py since it's just a convenience func?
        if not isinstance(key, list) and not isinstance(key, tuple):
            key = [key]
        if args:
            if isinstance(args, int):
                args = [args]
            key += list(args)
        try:
            r =  get_from_dict(self._keyword_dict, key)
            return r.value
        except KeyError:
            return default

    def __deepcopy__(self, memo):
        # TODO: this is not efficient
        return Project(''.join(self.convertToString()))

    def parsemfixdat(self, fname=None):
        """Parse the mfix.dat file with regular expressions.
        fname can be a StringIO instance, path, or a plain string.
        Saves the results in self.mfixDatKeyDict dictionary.
        """
        # TODO. move this to another module
        #  maybe write a real tokenizer and grammar.
        if fname:
            self.dat_file = fname
        assert self.dat_file is not None

        # check to see if the file is a StringIO object
        if isinstance(self.dat_file, StringIO):
            self._parsemfixdat(self.dat_file)
            return

        # Try to open the file
        try:
            with open(self.dat_file, 'rb') as dat_file:
                self._parsemfixdat(dat_file)
            return
        except (IOError, OSError):
            pass
        # maybe it is just a string?
        if is_text_string(self.dat_file):
            self._parsemfixdat(StringIO(self.dat_file))
            return

    def parseKeywordLine(self, line, equation_str=None):
        if not line.strip():
            yield(None, None, None)
            return
        matches = re_keyValue.findall(line)
        single_key = False
        if matches:
            for match in matches:
                # match could be: [keyword, args, value,
                #                   nextKeywordInLine, something]

                # convert to list so we can reassign
                match = list(match)

                # keyword
                key = match[0].lower().strip()

                # values
                val_string = match[2].strip()

                # remove spaces from equations: @( 2*pi)
                # (doing this with regex won't work because it can't detect
                # balanced parens)
                val_string = remove_spaces_from_equations(val_string)

                # look for shorthand [count]*[value] and expand.
                val_string = expand_shorthand(val_string)

                # split values using shlex, it will keep quoted strings together.
                if 'description' in key: # Don't do anything with description line, it may contain *, etc.
                    vals = [safe_shlex_split(line.strip())[-1]]
                    single_key = True
                else:
                    try:
                        vals = safe_shlex_split(val_string.strip())
                    except ValueError:
                        vals = []

                # if equation in comment, replace last value
                if equation_str is not None and vals:
                    vals[-1] = '@(' + equation_str +')'

                # clean the values converting to python types
                cleanVals = []
                for val in vals:
                    val = clean_string(val)
                    if val is not None:
                        cleanVals.append(val)

                # check keyword_doc to see if the keyword is an array
                if (self.keyword_doc is not None and key in self.keyword_doc
                        and not match[1] and self.keyword_doc[key]['args']):
                    match[1] = ','.join(str(arg['min']) for arg in
                                        self.keyword_doc[key]['args'].values())

                # clean up arguments
                args = []
                colon_arg = None # if one of the args is of the form lo:hi
                colon_lo = None
                colon_hi = None
                if match[1]:
                    for (i,arg) in enumerate(match[1].split(',')):
                        if ':' in arg:
                            if colon_arg is not None: # Only one index can have :
                                raise ValueError(match[1])
                            colon_arg = i
                            colon_lo, colon_hi = map(int,
                                                     arg.split(':'))
                            args.append(colon_lo)
                        else:
                            args.append(int(arg))

                numVals = len(cleanVals)

                if colon_arg is not None:
                    expect = 1 + colon_hi - colon_lo
                    if numVals != expect:
                        raise ValueError('Expected %s values, got %s' % (expect, numVals))

                # If multiple values, split apart into multiple key, args, values

                if numVals > 1:
                    keywordArgList = []

                    if numVals > 1:
                        if colon_arg is not None:
                            for n in range(colon_lo, colon_hi+1):
                                keywordArgList.append(
                                    args[:colon_arg] + [n] + args[colon_arg+1:])

                        elif args: # This distributes over the first index - is that
                                    # correct?
                                    # a(3,4) = 11*5 sets a(3,4) through a(13,4) = 5
                                    # instead of a(3,4) through a(3,14)

                            for n in range(args[0], args[0]+numVals):
                                keywordArgList.append([n] + args[1:])
                        else:
                            # hack for species eq
                            # FIXME - do this for more keywords which start
                            # at 0, like momentum_eq
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

    def parse_mfix_gui_comments(self, fobject):
        """read through the file looking for #!MFIX-GUI"""
        self.mfix_gui_comments.clear()
        for i, line in enumerate(fobject):
            line = to_unicode_from_fs(line).strip()
            if line.startswith("#!MFIX-GUI"):
                if "=" in line:
                    key, val = line[10:].split('=', 1)
                    key, val = key.strip(), val.strip()
                    self.mfix_gui_comments[key] = val
                    # For now, just ignore any other lines in this block  - it's just
                    # comments, and an experimental feature.  Don't want to create
                    # tight dependencies on GUI version - treat them like HTML tags,
                    # ignore the ones you can't handle

    def _parsemfixdat(self, fobject):
        """This does the actual parsing."""
        self.comments.clear()
        self._keyword_dict.clear()
        self.__initDataStructure__()
        self.dat_file_list = []
        self.thermo_data = []

        reactionSection = False
        thermoSection = False
        for i, line in enumerate(fobject):
            line = to_unicode_from_fs(line).strip()
            if line.startswith("#!MFIX-GUI"):
                # these should be already parsed
                continue

            if '@(RXNS)' in line:
                reactionSection = True
            elif '@(END)' in line and reactionSection:
                reactionSection = False
            elif 'thermo data' in line.lower():
                thermoSection = True
                # Don't save 'THERMO SECTION' line - we'll regenerate it.
            elif thermoSection:
                self.thermo_data.append(line)
            elif not reactionSection and not thermoSection:
                equation_str = None
                # remove comments
                commentedline = ''
                if line.startswith('#') or line.startswith('!'):
                    commentedline = line
                    line = ''
                elif '#' in line or '!' in line:
                    line, keywordComment = re.split('[#!]', line, maxsplit=1)

                    # look for equation
                    if 'MFIX-GUI' in keywordComment and 'eq{' in keywordComment:
                        start = keywordComment.find('eq{')
                        equation_str = keywordComment[start+3:keywordComment.find('}', start)]
                        keywordComment = '' # clears the comment so that it is not saved again
                else:
                    keywordComment = ''

                # loop through all keywords in the line
                try:
                    for (key, args, value) in self.parseKeywordLine(line, equation_str):
                        if key is None:
                            self.dat_file_list.append(line+commentedline)
                        else:
                            try:
                                self.updateKeyword(key, value, args, keywordComment)
                            except ValueError:
                                # error at line i
                                warnings.warn("Cannot set %s=%s" % (format_key_with_args(key, args), value))
                except Exception as e:
                    traceback.print_exception(*sys.exc_info())
                    warnings.warn("Parse error: %s: line %d, %s" % (e, i, line))

    def updateKeyword(self, key, value, args=None,  keywordComment=''):
        """Update or add a keyword to the project.  Raises ValueError if there is a
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
        """
        # TODO:  refactor
        if args is None:
            args = []

        # check to see if the keyword already exists
        if [key]+args in self:
            # if previous value is equation, update the parameter dict
            if isinstance(self[[key]+args].value, Equation) or isinstance(value, Equation):
                self.update_parameter_map(value, key, args)
            self[[key]+args].updateValue(value)
            if keywordComment:
                self[[key]+args].comment = keywordComment
            return self[[key]+args]

        # if equation, update the parameter dict
        if isinstance(value, Equation):
            self.update_parameter_map(value, key, args)

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

    def _recursiveRemoveKeyFromKeywordDict(self, keys=(), warn=True):
        keywordDict = self._keyword_dict
        for key in keys[:-1]:
            keywordDict = keywordDict[key]
        keywordDict.pop(keys[-1])

    def _purgeemptydicts(self):
        # TODO: write a test for this!
        while True:
            changed = False
            for (path,val) in list(recurse_dict_empty(self._keyword_dict)):
                if val == {}:
                    keywordDict = self._keyword_dict
                    for key in path[:-1]:
                        keywordDict = keywordDict[key]
                    keywordDict.pop(path[-1])
                    changed = True
            if not changed:
                break

    def keywordItems(self):
        for keys, value in recurse_dict_empty(self._keyword_dict):
            yield value

    def removeKeyword(self, key, args=None, warn=True):
        """Remove a keyword from the project.
        return True if item deleted.
        if warn=True, raise KeyError if item not present, else return False."""
        if args is None:
            args = []
        keyword = self.keywordLookup(key, args, warn=False)
        if keyword is None:
            if warn:
                raise KeyError(key)
            else:
                return False
        # remove from dat_file_list
        # Note. list.remove does not work because of the way keyword comparison is implemented
        #  (remove will remove the first key with matching value)
        #  Do we really need comparison of keyword object by value to work?
        #try:
        #    self.dat_file_list.remove(keyword) #
        #except ValueError:
        #    pass
        self.dat_file_list = [k for k in self.dat_file_list if k is not keyword]

        # remove from dict
        self._recursiveRemoveKeyFromKeywordDict([key]+args, warn)
        for i in range(len(args)):
            self._purgeemptydicts()

        # purge
        keyword.delete = True
        self._cleanDeletedItems()
        return True

    def _cleanDeletedItems(self):
        """Purge objects marked with self.delete==True"""

        for condType in [self.ics, self.bcs, self.pss, self.iss]:
            for cond in condType:
                for gas in cond.gasSpecies:
                    if gas.delete:
                        cond.gasSpecies.delete(gas)
                        #self.dat_file_list.remove(gas) FIXME this won't work
                for solid in cond.solids:
                    for species in solid.species:
                        if species.delete:
                            solid.species.delete(species)
                            # dat file list?
                    if solid.delete:
                        cond.solids.delete(solid)

                if cond.delete:
                    condType.delete(cond)

        for collection in (self.gasSpecies, self.solids, self.speciesEq,
                     self.linearEq, self.speciesEq, self.spx,
                     self.vtkvar, self.variablegrid):
            if hasattr(collection, 'species'):
                for item in collection.species:
                    if item.delete:
                        collection.species.delete(item)
            for item in collection:
                collection.delete(item)


    def keywordLookup(self, key, args=None, warn=True):
        """Search the project for a key and return the Keyword object
        If not found, raise KeyError if warn, else return None"""

        if args is None:
            args = []

        key = key.lower()
        keywordobject = None

        # FIXME
        # This is unfortunate - defeats the efficiency of dictionaries
        # because all lookups iterate through the whole dictionary
        for (k, v) in recurse_dict(self._keyword_dict):
            if k[0] == key:
                if args and args == v.args:
                    keywordobject = v
                    break
                elif not args: #? sloppy matching - why needed?
                    keywordobject = v
                    break
                else:
                    continue

        if warn and keywordobject is None:
            raise KeyError(key)
        return keywordobject


    def changekeywordvalue(self, key, value, args=None):
        """Change a value of a keyword.

        Parameters
        ----------
        key (str):
            keyword to change
        value (int, float, bool, str):
            the value to set the keyword to
        args (list or None):
            a list of arguments for that keyword"""

        if args is None:
            args = []
        keywordobject = self.keywordLookup(key, args)
        keywordobject.updateValue(value)
        return keywordobject

    def convertToString(self):
        for line in self.dat_file_list:
            if hasattr(line, 'line'):
                yield to_text_string(line.line() + '\n')
            else:
                yield to_text_string(line + '\n')

        if self.mfix_gui_comments: # Special comment block to hold non-keyword gui params
            yield '\n'
            yield '#!MFIX-GUI SECTION\n'
            for (key, val) in self.mfix_gui_comments.items():
                yield '#!MFIX-GUI %s = %s\n' % (key, val)

        if self.thermo_data:
            yield '\n'
            yield 'THERMO DATA\n'
        for line in self.thermo_data:
            yield line+'\n'

    def writeDatFile(self, fname):
        """ Write the project to specified text file"""
        ### TODO:  format species sections
        # delimit new additions from initial file contents (comment line)

        last_line = None
        with open(fname, 'wb') as dat_file:
            for line in self.convertToString():
                if line == last_line == '\n': # Avoid multiple blank lines
                    continue
                last_line = line
                dat_file.write(to_fs_from_unicode(line))

    def reset(self):
        self.dat_file = None
        self.dat_file_list = []
        self._keyword_dict.clear()
        self.comments.clear()
        self.thermo_data = []
        self.mfix_gui_comments.clear()
        self.parameter_key_map = {}

        for name in dir(self):
            attr = getattr(self, name)
            if isinstance(attr, Collection):
                Collection.__init__(attr)

    def update_parameter_map(self, new_value, key, args):
        """update the mapping of parameters and keywords"""
        key_args = format_key_with_args(key, args)

        # new params
        if isinstance(new_value, Equation):
            new_params = new_value.get_used_parameters()
        else:
            new_params = []

        # old params
        if [key]+args in self and isinstance(self[[key]+args].value, Equation):
            old_params = self[[key]+args].value.get_used_parameters()
        else:
            old_params = []

        add = set(new_params)-set(old_params)
        for param in add:
            if param not in self.parameter_key_map:
                self.parameter_key_map[param] = set()
            self.parameter_key_map[param].add(key_args)

        remove = set(old_params)-set(new_params)
        for param in remove:
            self.parameter_key_map[param].remove(key_args)
            if len(self.parameter_key_map[param]) == 0:
                self.parameter_key_map.pop(param)


if __name__ == '__main__':
    proj = Project()
    print(list(proj.parseKeywordLine('key = @( 2* 10)')))
