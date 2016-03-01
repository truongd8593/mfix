# Import from the future for Python 2 and 3 compatability!
from __future__ import print_function, absolute_import, unicode_literals

try:
    # Python 3
    from functools import reduce
except:
    pass

import re
import os

# import qt
from qtpy import QtGui

SCRIPT_DIRECTORY = './'

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


def make_callback(func, param):
    '''
    Helper function to make sure lambda functions are cached and not lost.
    '''
    return lambda: func(param)


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
