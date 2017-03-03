""" version module stores three ways of getting the version: pip (setuptools), git describe, and fallback __version__ """


import subprocess
import pkg_resources


# version numbering from PEP 440
__version__ = u"17.1dev"


def get_version():
    """ return the pip installed version if running from installed package, otherwise default """
    try:
        return get_pkg_version()
    except pkg_resources.DistributionNotFound:
        return __version__


def get_pkg_version():
    """ return currently installed mfix version found by pip/setuptools """
    return pkg_resources.get_distribution("mfix").version


def get_git_revision_short_hash():
    """Try to get the current git hash"""
    try:
        git_hash = subprocess.check_output(['git', 'describe', '--always']).strip().decode('utf-8')
    except subprocess.CalledProcessError:
        git_hash = __version__

    return git_hash
