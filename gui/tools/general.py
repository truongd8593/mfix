# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

try:
    # Python 3
    from functools import reduce
except:
    pass

import re
import os
import sys
import locale

# import qt
from qtpy import QtGui

SCRIPT_DIRECTORY = './'
PY2 = sys.version[0] == '2'
PY3 = sys.version[0] == '3'


def set_script_directory(script):
    global SCRIPT_DIRECTORY
    SCRIPT_DIRECTORY = script


def num_to_time(time, unit='s', outunit='time'):
    '''
    Function to convert time with a unit to another unit.
    '''
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
    " get path to images "
    path = os.path.join(SCRIPT_DIRECTORY, 'icons', name)

    if os.name == 'nt':
        path = path.replace('\\', '//')

    return path


def make_callback(func, *args, **kwargs):
    '''
    Helper function to make sure lambda functions are cached and not lost.
    '''
    return lambda: func(*args, **kwargs)


def get_icon(name, default=None, resample=False):
    """Return image inside a QIcon object
    default: default image name or icon
    resample: if True, manually resample icon pixmaps for usual sizes
    (16, 24, 32, 48, 96, 128, 256). This is recommended for QMainWindow icons
    created from SVG images on non-Windows platforms due to a Qt bug (see
    http://code.google.com/p/spyderlib/issues/detail?id=1314)."""
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
        return icon0
    else:
        return icon


def get_unique_string(base, listofstrings):
    " uniquify a string "
    if base in listofstrings:
        # look for number at end
        nums = re.findall('[\d]+', base)
        if nums:
            number = int(nums[-1]) + 1
            base = base.replace(nums[-1], '')
        else:
            number = 1

        base = get_unique_string(''.join([base, str(number)]), listofstrings)

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


def recurse_dict(d, keys=()):
    '''
    Recursively loop through a dictionary of dictionaries.

    Parameters
    ----------
    d (dict):
        a dictionary to loop over

    Returns
    -------
    (keys (tuple), value):
        a tuple of keys (tuple) and the value.
    '''
    if type(d) == dict:
        for k in d:
            for rv in recurse_dict(d[k], keys + (k, )):
                yield rv
    else:
        yield (keys, d)


def recurse_dict_empty(d, keys=()):
    '''
    Recursively loop through a dictionary of dictionaries.

    Parameters
    ----------
    d (dict):
        a dictionary to loop over

    Returns
    -------
    (keys (tuple), value):
        a tuple of keys (tuple) and the value.
    '''
    if type(d) == dict and d:
        for k in d:
            for rv in recurse_dict_empty(d[k], keys + (k, )):
                yield rv
    else:
        yield (keys, d)


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
    """
    Return a unicode version of string decoded using the file system encoding.
    """
    if not is_string(string):  # string is a QString
        string = to_text_string(string.toUtf8(), 'utf-8')
    else:
        if is_binary_string(string):
            try:
                unic = string.decode(FS_ENCODING)
            except (UnicodeError, TypeError):
                pass
            else:
                return unic
    return string
