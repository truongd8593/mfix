# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

try:
    # Python 3
    from functools import reduce
except:
    pass

import re
import os
import sys
import locale
import logging
import shlex
import copy

log = logging.getLogger(__name__)


# TODO: factor out util funcs which don't require Qt
# import qt
from qtpy import QtGui, QtWidgets, QtCore

SCRIPT_DIRECTORY = './'
PY2 = sys.version_info.major == 2
PY3 = sys.version_info.major == 3


def set_script_directory(script):
    global SCRIPT_DIRECTORY
    SCRIPT_DIRECTORY = script


# Helper functions
def get_mfix_home():
    "return the top level directory"
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


def format_key_with_args(key, args=None):
    if args:
        return "%s(%s)" % (key, ','.join(str(a) for a in args))
    else:
        return str(key)


def unformat_key_with_args(string):
    """companion function to "undo" format_key_with_args"""
    if string.endswith(')'):
        key, args = string[:-1].split('(')
        args = [int(arg) for arg in args.split(',')]
    else:
        key = string
        args = []
    return key, args


def plural(n, word):
    fmt = "%d %s" if n == 1 else "%d %ss"
    return fmt % (n, word)


def set_item_noedit(item):
    item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)

def set_item_enabled(item, enabled):
    """Enable/disable items which do not have a setEnabled method, like menu items"""
    flags = item.flags()
    if enabled:
        flags |= QtCore.Qt.ItemIsEnabled
    else:
        flags &= ~QtCore.Qt.ItemIsEnabled
    item.setFlags(flags)

def get_combobox_item(combobox, n):
    """Return the n'th menu item from a combobox"""
    model = combobox.model()
    return model.item(n, 0)

def get_selected_row(table):
    """get index of selected row from a QTableWidget"""
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedIndexes())
    return None if not rows else rows.pop()


def num_to_time(time, unit='s', outunit='time'):
    """Convert time with a unit to another unit."""
    unit = unit.lower()
    time = float(time)

    # convert time to seconds
    if unit in ['d', 'days']:
        time *= 60 * 60 * 24
    elif unit in ['h', 'hr']:
        time *= 60 * 60
    elif unit in ['m', 'min']:
        time *= 60

    if outunit == 'time':
        time = reduce(
            lambda ll, b: divmod(ll[0], b) + ll[1:], [(time, ), 60, 60, 24])

        tlist = []
        for i, num in enumerate(time):
            if i == 3:
                tlist.append('{:.2f}'.format(num))
            elif num > 0 or len(tlist) > 0:
                tlist.append(int(num))

        return ':'.join([str(t) for t in tlist])

    elif outunit in ['d', 'days']:
        return time/(60.0*60.0*24.0)

    elif outunit in ['h', 'hr', 'hrs']:
        return time/(60.0*60.0)

    elif outunit in ['m', 'min', 'mins']:
        return time/(60.0)

    else:
        return time


def get_image_path(name):
    """"get path to images"""
    path = os.path.join(SCRIPT_DIRECTORY, 'icons', name)

    if os.name == 'nt':
        path = path.replace('\\', '//')

    return path


def make_callback(func, *args, **kwargs):
    """Helper function to make sure lambda functions are cached and not lost."""
    return lambda: func(*args, **kwargs)


icon_cache = {}
def get_icon(name, default=None, resample=False):
    """Return image inside a QIcon object
    default: default image name or icon
    resample: if True, manually resample icon pixmaps for usual sizes
    (16, 24, 32, 48, 96, 128, 256). This is recommended for QMainWindow icons
    created from SVG images on non-Windows platforms due to a Qt bug (see
    http://code.google.com/p/spyderlib/issues/detail?id=1314)."""

    icon = icon_cache.get((name,default,resample))
    if icon:
        return icon

    if default is None:
        icon = QtGui.QIcon(get_image_path(name))
    elif isinstance(default, QtGui.QIcon):
        icon_path = get_image_path(name)
        icon = default if icon_path is None else QtGui.QIcon(icon_path)
    else:
        icon = QtGui.QIcon(get_image_path(name, default))
    if resample:
        icon0 = QtGui.QIcon()
        for size in (16, 24, 32, 48, 96, 128, 256, 512):
            icon0.addPixmap(icon.pixmap(size, size))
        icon_cache[(name, default, resample)] = icon0
        ret = icon0
    else:
        ret= icon
    icon_cache[(name, default, resample)] = ret
    return ret


def get_unique_string(base, listofstrings):
    "uniquify a string"
    if base in listofstrings:
        # look for number at end
        nums = re.findall('[\d]+$', base)
        if nums:
            number = int(nums[-1]) + 1
            base = base.replace(nums[-1], '')
        else:
            number = 1
        base = get_unique_string(base + str(number), listofstrings)

    return base

def widget_iter(widget):
    for child in widget.children():
        if child.children():
            for child2 in widget_iter(child):
                yield child2
        yield child


def get_from_dict(data_dict, map_list):
    return reduce(lambda d, k: d[k], map_list, data_dict)


def set_in_dict(data_dict, map_list, value):
    get_from_dict(data_dict, map_list[:-1])[map_list[-1]] = value


def recurse_dict(d, path=()):
    """Depth-first iterator though nested dictionaries,
    Yields (path, value), eg if d['a']['b']['c]=3,
    one of the yielded values will be (('a','b','c'), 3).
    See test_recurse_dict for an example"""

    # Not quite depth-first b/c order of dictionary key is arbitrary
    for (k,v) in d.items():
        subpath = path + (k,)
        if isinstance(v, dict):
            for r in recurse_dict(v, subpath):
                yield r
        else:
            yield (subpath, v)

def test_recurse_dict():
    d = {1: {2:3,
             3: {},
             4: {5:6},
             5: {6:{}, 7:8}},
         2: {}}

    l = list(recurse_dict(d))
    l.sort()
    assert l == [((1,2), 3), ((1,4,5), 6), ((1,5,7), 8)]


def recurse_dict_empty(d, path=()):
    """Depth-first iterator though nested dictionaries
    Yields (keytuple, value), eg if d['a']['b']['c]=3,
    one of the yielded values will be (('a','b','c'), 3)
    Differs from recurse_dict in that an empty dictionary
    encountered as a value will be yielded, eg
    if d[1][2]={}, then ((1,2), {}) will be one of the
    yielded values.
    See test_recurse_dict_empty for an example"""

    # Not quite depth-first b/c order of dictionary key is arbitrary
    for (k,v) in d.items():
        subpath = path + (k,)
        if isinstance(v, dict):
            if v == {}:
                yield(subpath, v)
            else:
                for r in recurse_dict_empty(v, subpath):
                    yield r
        else:
            yield (subpath, v)

def test_recurse_dict_empty():
    d = {1: {2:3,
             3: {},
             4: {5:6},
             5: {6:{}, 7:8}},
         2: {}}

    l = list(recurse_dict_empty(d))
    l.sort()
    assert l == [((1,2), 3), ((1,3), {}), ((1,4,5), 6), ((1,5,6), {}), ((1,5,7), 8), ((2,), {})]

#http://stackoverflow.com/questions/14218992/shlex-split-still-not-supporting-unicode
#see also notes at https://pypi.python.org/pypi/ushlex/
def safe_shlex_split(string):
    """shlex.split is not unicode-safe"""
    if PY2:
        return [s.decode('utf-8') for s in shlex.split(string.encode('utf-8'))]
    else:
        return shlex.split(string)

# Debugging hooks
def debug_trace():
    """Set a tracepoint in the Python debugger that works with Qt"""
    from qtpy.QtCore import pyqtRemoveInputHook
    from pdb import set_trace
    pyqtRemoveInputHook()
    set_trace()


# These functions were extracted from spyder's p3compat.py code.
def is_text_string(obj):
    """Return True if `obj` is a text string, False if it is anything else,
    like binary data (Python 3) or QString (Python 2, PyQt API #1)"""
    if PY2:
        # Python 2
        return isinstance(obj, basestring)
    else:
        # Python 3
        return isinstance(obj, str)


def is_binary_string(obj):
    """Return True if `obj` is a binary string, False if it is anything else"""
    if PY2:
        # Python 2
        return isinstance(obj, str)
    else:
        # Python 3
        return isinstance(obj, bytes)


def is_string(obj):
    """Return True if `obj` is a text or binary Python string object,
    False if it is anything else, like a QString (Python 2, PyQt API #1)"""
    return is_text_string(obj) or is_binary_string(obj)


def is_unicode(obj):
    """Return True if `obj` is unicode"""
    if PY2:
        # Python 2
        return isinstance(obj, unicode)
    else:
        # Python 3
        return isinstance(obj, str)


def to_text_string(obj, encoding=None, errors=None):
    """Convert `obj` to (unicode) text string"""
    if PY2:
        # Python 2
        if encoding is None and errors is None:
            return unicode(obj)
        elif encoding is None and errors:
            return unicode(obj, errors=errors)
        elif encoding and errors is None:
            return unicode(obj, encoding)
        else:
            return unicode(obj, encoding, errors=errors)
    else:
        # Python 3
        if isinstance(obj, str):
            # In case this function is not used properly, this could happen
            return obj
        elif encoding is None and errors is None:
            return str(obj)
        elif encoding is None and errors:
            return str(obj, errors=errors)
        elif encoding and errors is None:
            return str(obj, encoding)
        else:
            return str(obj, encoding, errors=errors)


# These functions were extracted from spyder's econding.py code.
PREFERRED_ENCODING = locale.getpreferredencoding()


def getfilesystemencoding():
    """
    Query the filesystem for the encoding used to encode filenames
    and environment variables.
    """
    encoding = sys.getfilesystemencoding()
    if encoding is None:
        # Must be Linux or Unix and nl_langinfo(CODESET) failed.
        encoding = PREFERRED_ENCODING
    return encoding

FS_ENCODING = getfilesystemencoding()


def to_unicode_from_fs(string):
    """Return a unicode version of string decoded using the file system encoding."""
    if not is_string(string):  # string is a QString
        string = to_text_string(string.toUtf8(), 'utf-8')
    else:
        if is_binary_string(string):
            try:
                return string.decode(encoding=FS_ENCODING, errors='replace')
            except (UnicodeDecodeError, TypeError) as e:
                log.warn("%s: %s" % (e, string))

    return string

def to_fs_from_unicode(string):
    return string.encode(encoding=FS_ENCODING, errors='replace')

class CellColor(object):
    """
    A class to store color information and return '' if str or print is called
    on it. This is used to store colors in cells of a table.
    """
    def __init__(self, color=[1, 0, 0], text=''):

        self.color = color
        self.text = text

    @property
    def color_int(self):
        return [255*c for c in self.color]

    @property
    def color_float(self):
        return self.color

    @property
    def qcolor(self):
        return QtGui.QColor(*self.color_int)

    def __repr__(self):
        return self.text


def insert_append_action(menu, action, insert=None):
    if insert:
        menu.insertAction(insert, action)
    else:
        menu.addAction(action)


def insert_append_separator(menu, insert=None):
    if insert:
        menu.insertSeparator(insert)
    else:
        menu.addSeparator()

def topological_sort(dependency_dict):
    '''
    Sort the dependency tree.
    Inspired by: http://code.activestate.com/recipes/578272-topological-sort/
    '''
    
    data = copy.deepcopy(dependency_dict)

    # Ignore self dependencies.
    for k, v in data.items():
        v.discard(k)
    # Find all items that don't depend on anything.
    extra_items_in_deps = reduce(set.union, data.itervalues()) - set(data.iterkeys())
    # Add empty dependences where needed
    data.update({item: set() for item in extra_items_in_deps})
    while True:
        ordered = set(item for item, dep in data.iteritems() if not dep)
        if not ordered:
            break
        yield ordered
        data = {item: (dep - ordered)
                for item, dep in data.iteritems()
                if item not in ordered}
    assert not data, "Cyclic dependencies exist among these items:\n%s" % '\n'.join(repr(x) for x in data.iteritems())

if __name__ == '__main__':
    test_recurse_dict()
    test_recurse_dict_empty()
