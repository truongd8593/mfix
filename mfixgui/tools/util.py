""" module for functions that don't depend on Qt """

from __future__ import print_function, absolute_import, unicode_literals, division


import os
import sys

def find_mfixgui_module_directory():
    """ return the directory corresponding to the module mfixgui/__init__.py """

    general_mod = os.path.realpath(__file__)
    tools_mod = os.path.dirname(general_mod)
    mfixgui_mod = os.path.dirname(tools_mod)
    return mfixgui_mod

SCRIPT_DIRECTORY = find_mfixgui_module_directory()

# Helper functions
def get_mfix_home():
    """return the top level MFiX directory containing directories defaults,model,tutorials"""
    mfix_src_root = os.path.dirname(find_mfixgui_module_directory())
    if os.path.isfile(os.path.join(mfix_src_root, 'configure_mfix')):
        # if configure_mfix is present, we are in a source directory
        return mfix_src_root
    else:
        # we are in an installed package, search PYPATH for directory named "mfix"
        for pypath in sys.path:
            mfix_path = os.path.join(pypath, 'mfix')
            if os.path.isfile(os.path.join(mfix_path, 'configure_mfix')):
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
