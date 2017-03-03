""" This module builds the mfixsolver as a Python extension library. It is used
by setup.py when building the mfix distribution package, and is run as a
standalone command to build the custom mfixsolvers. """

import atexit
import os
import shutil
import subprocess
import sys
import tempfile

# must import setuptools and cygwinccompiler before numpy.distutils
import setuptools
import distutils.cygwinccompiler

from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.core import Extension, setup

from mfixgui.tools.general import get_mfix_home

from mfixgui.version import get_version

HERE = os.path.abspath(os.path.dirname(__file__))
NAME = 'mfix'

CONFIGURE_ARGS = 'CC=gcc FC=gfortran FCFLAGS=-fPIC FFLAGS=-fPIC '
MODEL_DIR = os.path.join(HERE, 'model')


DEFAULT_CC = 'gcc'
DEFAULT_FC = 'gfortran'
DEFAULT_CFLAGS = '-O2 -fPIC'
DEFAULT_FCFLAGS = '-O2 -fPIC'


def get_configure_args(cc=DEFAULT_CC,
                       fc=DEFAULT_FC,
                       cflags=DEFAULT_CFLAGS,
                       fcflags=DEFAULT_FCFLAGS):
    return 'CC="%s" FC="%s" CFLAGS="%s" FCFLAGS="%s" FFLAGS="%s"' % (cc, fc, cflags, fcflags, fcflags)


def make_mfixsolver():
    configure_args_dirname = get_configure_args()
    configure_args_dirname = configure_args_dirname.replace(' ', '_')
    configure_args_dirname = configure_args_dirname.replace('=', '_')
    configure_args_dirname = configure_args_dirname.replace('"', '_')
    configure_args_dirname = configure_args_dirname.replace('__', '_').strip('_')
    build_dir = os.path.join('build', configure_args_dirname)

    return Extension(name='mfixsolver',
                     sources=get_pymfix_src(),
                     extra_f90_compile_args=['-cpp'],
                     module_dirs=[os.path.join(build_dir, 'model')],
                     extra_objects=[
                         '.build/read_database.o',
                         '.build/read_namelist.o',
                         os.path.join(build_dir, 'build-aux/libmfix.a'),
                     ])


def get_pymfix_src():
    """ copies those Fortran sources to be built with Python to .f90 extension """
    pymfix_src = [
        'param_mod.f',
        'param1_mod.f',
        'dmp_modules/compar_mod.f',
        'dmp_modules/debug_mod.f',
        'dmp_modules/parallel_mpi_mod.f',
        'fldvar_mod.f',
        'des/discretelement_mod.f',
        'des/des_time_march.f',
        'iterate.f',
        'residual_mod.f',
        'time_step.f',
        'main.f',
        'run_mod.f',
    ]

    f90_tmp = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, f90_tmp)

    os.makedirs(os.path.join(f90_tmp, 'dmp_modules'))
    os.makedirs(os.path.join(f90_tmp, 'des'))
    for src in pymfix_src:
        shutil.copyfile(os.path.join(get_mfix_home(), 'model', src),
                        os.path.join(f90_tmp, src)+'90')

    return [os.path.join(f90_tmp, s)+'90' for s in pymfix_src]


class BuildMfixCommand(setuptools.Command):
    """ builds libmfix (which does not depend on the Python version) """
    description = "build mfix (Fortran code)"
    user_options = [
        # The Format is (long option, short option, description)
        ('cc=', None, 'C compiler'),
        ('fc=', None, 'Fortran 90 compiler'),
        ('fcflags=', None, 'flags for Fortran 90 compiler'),
        ('cflags=', None, 'flags for C compiler'),
    ]

    def initialize_options(self):
        self.cc = DEFAULT_CC
        self.fc = DEFAULT_FC
        self.fcflags = DEFAULT_FCFLAGS
        self.cflags = DEFAULT_CFLAGS

    def finalize_options(self):
        if '-fPIC' not in self.fcflags:
            self.fcflags += ' -fPIC'
        if '-fPIC' not in self.cflags:
            self.cflags += ' -fPIC'

    def run(self):
        """ configure and build mfix; requires bash to be in PATH """

        configure_mfix = os.path.join(get_mfix_home(), 'configure_mfix')

        args = get_configure_args(self.cc, self.fc, self.cflags, self.fcflags)
        cmd = 'bash %s --buildrundir %s' % (configure_mfix, args)

        subprocess.check_call(cmd, shell=True)
        subprocess.check_call("make", shell=True)


def mfix_prereq(command_subclass):
    """A decorator for classes subclassing one of the setuptools commands.

    It modifies the run() method to run build_mfix before build_ext
    """
    orig_run = command_subclass.run

    def modified_run(self):
        self.run_command('build_mfix')
        orig_run(self)

    command_subclass.run = modified_run
    return command_subclass


@mfix_prereq
class BuildExtCommand(build_ext):
    pass


def main():
    """ Custom installation: set --prefix to RUNDIR
     https://docs.python.org/3/install/#custom-installation """

    rundir = os.getcwd()
    builddir = os.path.join(rundir, '.build')

    pypath = os.path.join(rundir,
                          'lib',
                          'python%d.%d' % sys.version_info[:2],
                          'site-packages')
    os.environ['PYTHONPATH'] = pypath
    for dirname in (builddir, pypath):
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

    sys.argv[1:] = ['install_lib',
                    '--install-dir', pypath,
                   ]

    os.chdir(builddir)

    setup(
        name='custom_mfixsolver',

        cmdclass={
            'build_ext': BuildExtCommand,
            'build_mfix': BuildMfixCommand,
        },

        version=get_version(),

        description='MFiX computational fluid dynamics solver',

        install_requires=[
            'mfix',
        ],

        ext_modules=[make_mfixsolver(),],
        entry_points={},
        scripts={},
    )

    mfixsolver = os.path.join(rundir, 'mfixsolver')
    with open(mfixsolver, 'w') as wrapper:
        wrapper.write('#!/bin/sh\n')
        wrapper.write('\n')
        wrapper.write('env PYTHONPATH=%s %s -m mfixgui.pymfix "$@"\n' % (pypath, sys.executable))

    # make executable
    permissions = os.stat(mfixsolver)
    os.chmod(mfixsolver, permissions.st_mode | 0o111)

if __name__ == '__main__':
    main()
