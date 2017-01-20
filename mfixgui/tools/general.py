# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals, division

try:
    # Python 3
    from functools import reduce
except:
    pass

import copy
import locale
import logging
import operator
import os
import random
import re
import shlex
import site
import subprocess
import sys

from collections import OrderedDict

log = logging.getLogger(__name__)

# TODO: factor out util funcs which don't require Qt
# import qt
from qtpy import QtGui, QtWidgets, QtCore

PY2 = sys.version_info.major == 2
PY3 = sys.version_info.major == 3

SCRIPT_DIRECTORY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))

# Helper functions
def get_mfix_home():
    """return the top level MFIX directory"""
    # TODO:  This is confusing, add comment why both get_mfix_home and
    # SCRIPT_DIRECTORY are needed, and how they differ
    top_level_pkg_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    # if configure_mfix is present, we weren't installed via setup.py
    if os.path.exists(os.path.join(top_level_pkg_dir,'configure_mfix')):
        return top_level_pkg_dir
    elif os.path.exists(os.path.join(site.USER_BASE, 'share', 'mfix')):
        return os.path.join(site.USER_BASE, 'share', 'mfix')
    elif os.path.exists(os.path.join(sys.prefix, 'share', 'mfix')):
        return os.path.join(sys.prefix, 'share', 'mfix')
    elif os.path.exists(os.path.join('/usr', 'local', 'share', 'mfix')):
        return os.path.join('/usr', 'local', 'share', 'mfix')
    elif os.path.exists(os.path.join('/usr', 'share', 'mfix')):
        return os.path.join('/usr', 'share', 'mfix')
    else:
        for p in sys.path:
            mfix_path = os.path.join(p, 'share', 'mfix')
            if os.path.exists(mfix_path):
                return mfix_path
            mfix_path = os.path.join(p, os.pardir, os.pardir, os.pardir, 'share', 'mfix')
            if os.path.exists(mfix_path):
                return mfix_path
        raise Exception("Unable to find MFIX_HOME")


def format_key_with_args(key, args=None):
    if args is not None and args != []:
        if isinstance(args, int):
            args = [args]
        return "%s(%s)" % (key, ','.join(str(a) for a in args))
    else:
        return str(key)


def parse_key_with_args(string):
    # inverse of format_key_with_args
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


def item_enabled(item):
    flags = item.flags()
    return bool(flags & QtCore.Qt.ItemIsEnabled)


def get_combobox_item(combobox, n):
    """Return the n'th menu item from a combobox"""
    model = combobox.model()
    return model.item(n, 0)


def get_selected_row(table):
    """get index of selected row from a QTableWidget"""
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedIndexes())
    return None if not rows else rows.pop()


def get_selected_rows(table):
    """get index of selected row from a QTableWidget"""
    # note, currentRow can return  >0 even when there is no selection
    rows = set(i.row() for i in table.selectedIndexes())
    return sorted(list(rows))


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


pixmap_cache = {}
def get_pixmap(name, width, height):
    pixmap = pixmap_cache.get((name, width, height))
    if pixmap is None:
        pixmap = QtGui.QPixmap(get_image_path(name)).scaled(
            width, height, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        pixmap_cache[(name, width, height)] = pixmap
    return pixmap


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

def random_pastel_color():
    return [(random.randint(0, 128) + 100)/255.0 for i in range(3)]

class CellColor(object):
    """
    A class to store color information and return '' if str or print is called
    on it. This is used to store colors in cells of a table.
    """
    def __init__(self, color=[1, 0, 0], text=''):

        if isinstance(color, (list, tuple)):
            self.color = QtGui.QColor(*color)
        else:
            self.color = QtGui.QColor(color)
        self.text = text

    @property
    def color_int(self):
        return [255*c for c in self.color_float]

    @property
    def color_float(self):
        return self.color.getRgbF()[:3]

    @property
    def color_hex(self):
        return self.color.name()

    @property
    def qcolor(self):
        return self.color

    def __repr__(self):
        return self.text

    def __deepcopy__(self, memo):
        return CellColor(copy.deepcopy(self.color_hex), copy.deepcopy(self.text))

    def rand(self):
        self.color.setRgbF(*random_pastel_color())


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
    extra_items_in_deps = reduce(set.union, data.values()) - set(data.keys())
    # Add empty dependences where needed
    data.update({item: set() for item in extra_items_in_deps})
    while True:
        ordered = set(item for item, dep in data.items() if not dep)
        if not ordered:
            break
        yield ordered
        data = {item: (dep - ordered)
                for item, dep in data.items()
                if item not in ordered}
    assert not data, "Cyclic dependencies exist among these items:\n%s" % '\n'.join(repr(x) for x in data.items())

def drop_row_column_triangular(a, n, r):
    # Inputs:
    #  a :  upper-triangular n by n matrix represented as a list (n)(n+1)/2
    #  n :  size of input
    #  r :  index of row/column to drop, 1-based
    # Return:
    #  list representing matrix "a" with row/column "r" removed
    ret = []
    i = j = 1
    for x in a:
        if i != r and j != r:
            ret.append(x)
        j += 1
        if j > n:
            i += 1
            j = i
    return ret

def append_row_column_triangular(a, n, fill_value = None):
    # Append a row and column to an upper-triangular rank-n matrix
    ret = []
    i = j = 1
    for x in a:
        ret.append(x)
        if j == n:
            ret.append(fill_value)
        j += 1
        if j > n:
            i += 1
            j = i
    ret.append(fill_value)
    return ret

def sort_dict(dict_, key, start=0):
    """given an dict of dicts and a key, sort the outside dict based on the
    value of one of the the internal dict's keys and return the sorted
    OrderedDict"""
    return OrderedDict(
        [(k, dict_[old_k])
         for k, (old_k, v) in enumerate(sorted([(k, v[key])
         for k, v in dict_.items()], key=operator.itemgetter(1)), start)])

def clear_layout(layout):
    """given a layout, clear all widgets"""
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget:
            widget.deleteLater()

def extract_config(path):
    '''Given a path to a file, extract the section that starts with ## CONFIG
    and ends with ## END CONFIG'''
    config = []
    script = []
    record = False
    with open(path) as f:
        for line in f:
            l = line.rstrip()
            if '## CONFIG' in l:
                record = True
            elif '## END CONFIG' in l:
                record = False
            elif record:
                config.append(l)
            else:
                script.append(l)
    return '\n'.join(config), '\n'.join(script+[''])

def replace_with_dict(string, dict_):
    # Deprecate in favor of built-in Python string operators
    '''given a string and a dict, replace all dict.key found in the string
    with the corresponding dict.value'''

    for key, value in dict_.items():
        string = string.replace('${'+key+'}', str(value))
    return string

def is_vnc():
    """determine if the gui is running in vnc"""
    if os.name == 'nt':
        return False
    xdpyinfo = subprocess.Popen('xdpyinfo', stdout=subprocess.PIPE).communicate()[0]
    return 'vnc' in str(xdpyinfo)

def get_separator(vertical=True):
    """create a QFrame that looks like a separator"""
    f = QtWidgets.QFrame
    line = f()
    if vertical:
        line.setFrameShape(f.VLine)
    else:
        line.setFrameShape(f.HLine)
    line.setFrameShadow(f.Sunken)
    return line

def get_username():
    """attempt to retunr the currnet user name"""
    name = 'unknown'
    for e in ['USER', 'user', 'USERNAME', 'username']:
        name = os.environ.get(e)
        if name:
            break
    if not name:
        name = 'unknown'
    return name

def convert_string_to_python(string):
    """Attempt to convert strings to python types"""

    if not string:
        return ''

    # remove all quotes
    string = string.replace("'", '').replace('"', '')
    # remove any leading or trailing space, after removing quotes
    string = string.strip()
    # Remove comma separators if present
    if string.endswith(','):
        string = string[:-1]

    # lower-case version of string
    s_low = string.lower()

    if s_low in ('.t.',  '.true.'):
        return True
    elif s_low in ('.f.', '.false.'):
        return False

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


def safe_float(value, default=0.0):
    """try to convert the value to a float, if ValueError, return default"""
    try:
        return float(value)
    except ValueError:
        return default


def safe_int(value, default=0):
    """try to convert the value to a int, if ValueError, return default"""
    try:
        return int(value)
    except ValueError:
        return default

def deepcopy_dict(dirty_dict, qobjects=False):
    '''deep copy a dictionary that has Qt objects in it
    Note: python 3.6+ can't copy Qt objects like QPixmap
    setting qobjects=True will copy the qt objects'''

    clean_dict = OrderedDict() if isinstance(dirty_dict, OrderedDict) else {}

    for key, value in dirty_dict.items():
        if isinstance(value, (dict, OrderedDict)):
            clean_dict[key] = deepcopy_dict(value, qobjects)
        elif isinstance(value, QtGui.QPixmap):
            if qobjects:
                clean_dict[key] = QtGui.QPixmap(value)
        elif isinstance(value, QtGui.QColor):
            if qobjects:
                clean_dict[key] = QtGui.QColor(value)
        else:
            clean_dict[key] = copy.deepcopy(value)
    return clean_dict

if __name__ == '__main__':
    def test_recurse_dict():
        d = {1: {2:3,
                 3: {},
                 4: {5:6},
                 5: {6:{}, 7:8}},
             2: {}}

        l = list(recurse_dict(d))
        l.sort()
        assert l == [((1,2), 3), ((1,4,5), 6), ((1,5,7), 8)]

    test_recurse_dict()

    def test_recurse_dict_empty():
        d = {1: {2:3,
                 3: {},
                 4: {5:6},
                 5: {6:{}, 7:8}},
             2: {}}

        l = list(recurse_dict_empty(d))
        l.sort()
        assert l == [((1,2), 3), ((1,3), {}), ((1,4,5), 6), ((1,5,6), {}), ((1,5,7), 8), ((2,), {})]


    test_recurse_dict_empty()

    def test_drop_add_triangular():
        a = [ 1, 2, 3, 4,
                 5, 6, 7,
                    8, 9,
                       10]
        b = drop_row_column_triangular(a, 4, 2)
        assert b == [1, 3, 4,
                        8, 9,
                           10]

        c = drop_row_column_triangular(b, 3, 3)
        assert c == [1, 3,
                        8]

        d = append_row_column_triangular(c, 2, 999)
        assert d == [1, 3, 999,
                        8, 999,
                           999]

    test_drop_add_triangular()
