"""This module has the __version__ string imported by other code both at
buildtime and runtime. version numbering from PEP 440"""

import os
import subprocess
import pkg_resources


GIT_DESCRIBE_TAG = os.environ.get('GIT_DESCRIBE_TAG', None)
GIT_DESCRIBE_NUMBER = os.environ.get('GIT_DESCRIBE_NUMBER', None)

if GIT_DESCRIBE_TAG:
    if GIT_DESCRIBE_NUMBER:
        __version__ = '%s-%s' % (GIT_DESCRIBE_TAG, GIT_DESCRIBE_NUMBER)
    else:
        __version__ = '%s' % GIT_DESCRIBE_TAG
else:
    try:
        # if mfix is installed with pip, return installed version
        __version__ = pkg_resources.get_distribution("mfix").version
    except pkg_resources.DistributionNotFound:
        # running from development version, perhaps with 'python -m mfixgui'
        __version__ = get_git_revision_short_hash()


def get_git_revision_short_hash():
    """Try to get the current git hash"""
    default_version = u"unknown"
    try:
        git_hash = subprocess.check_output(['git', 'describe', '--always']).strip().decode('utf-8')
    except subprocess.CalledProcessError:
        git_hash = default_version
    except OSError:
        git_hash = default_version

    return git_hash
