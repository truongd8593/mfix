""" import stuff into top level package module to simplify imports elsewhere in the code """

from mfixgui.tools.util import get_mfix_home, SCRIPT_DIRECTORY, format_key_with_args
from mfixgui.tools.general import (get_icon,
                                   get_pixmap,
                                   get_separator,
                                   safe_shlex_split,
                                   to_unicode_from_fs,
                                   to_fs_from_unicode)
# from mfixgui.tools.general import (get_icon, get_pixmap, get_separator)
