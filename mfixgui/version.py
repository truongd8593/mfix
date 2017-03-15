"""This module has the __version__ string imported by other code both at
buildtime and runtime. version numbering from PEP 440"""

import os
import subprocess
import pkg_resources


GIT_DESCRIBE_TAG = os.environ.get('GIT_DESCRIBE_TAG', None)
GIT_DESCRIBE_NUMBER = os.environ.get('GIT_DESCRIBE_NUMBER', None)

if GIT_DESCRIBE_TAG and GIT_DESCRIBE_NUMBER:
    __version__ = '%s-%s' % (GIT_DESCRIBE_TAG, GIT_DESCRIBE_NUMBER)
else:
    try:
        # return the pip installed version if running from installed package, otherwise default
        __version__ = pkg_resources.get_distribution("mfix").version
    except pkg_resources.DistributionNotFound:
        __version__ = u"17.1dev0"


def get_git_revision_short_hash():
    """Try to get the current git hash"""
    try:
        git_hash = subprocess.check_output(['git', 'describe', '--always']).strip().decode('utf-8')
    except subprocess.CalledProcessError:
        git_hash = __version__
    except OSError:
        git_hash = __version__

    return git_hash
