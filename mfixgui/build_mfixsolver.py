""" This module builds the mfixsolver as a Python extension library. It is used
by setup.py when building the mfix distribution package, and is run as a
standalone command to build the custom mfixsolvers. """

import atexit
import glob
import os
import platform
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

MODEL_DIR = os.path.join(HERE, 'model')
CONFIGURE_ARGS = []
MAKE_ARGS = []


def make_mfixsolver():
    configure_args_dirname = ' '.join(CONFIGURE_ARGS)
    configure_args_dirname = configure_args_dirname.replace('--dmp', '--enable-dmp')
    configure_args_dirname = configure_args_dirname.replace('--smp', '--enable-smp')
    configure_args_dirname = configure_args_dirname.replace('--python', '--enable-python')
    configure_args_dirname = configure_args_dirname.replace('--buildrundir', '')
    configure_args_dirname = configure_args_dirname.replace(' ', '_')
    configure_args_dirname = configure_args_dirname.replace('=', '_')
    configure_args_dirname = configure_args_dirname.replace('"', '_')
    configure_args_dirname = configure_args_dirname.replace('__', '_').strip('_')
    build_dir = os.path.join('build', configure_args_dirname)

    udfs = []
    if not os.path.isfile('configure_mfix'):
        for root, dirnames, filenames in os.walk('.'):
            for filename in glob.glob(os.path.join(root, '*.f')):
                udfs.append(filename)

        udfs = [ os.path.splitext(f)[0]+'.o' for f in udfs if not f.startswith('build') ]

    extra_objects = udfs + [
        '.build/read_database.o',
        '.build/read_namelist.o',
        os.path.join(build_dir, 'build-aux/libmfix.a') ,
    ]

    return Extension(name='mfixsolver',
                     sources=get_pymfix_src(),
                     extra_f90_compile_args=['-cpp'],
                     module_dirs=[os.path.join(build_dir, 'model')],
                     extra_objects=extra_objects)


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
        pass

    def finalize_options(self):
        pass

    def run(self):
        """ configure and build mfix; requires bash to be in PATH """

        configure_mfix = os.path.join(get_mfix_home(), 'configure_mfix')

        global CONFIGURE_ARGS
        global MAKE_ARGS

        cmd = 'bash %s --buildrundir %s' % (configure_mfix, ' '.join(CONFIGURE_ARGS))
        subprocess.check_call(cmd, shell=True)

        cmd = 'make libmfix %s' % ' '.join(MAKE_ARGS)
        subprocess.check_call(cmd, shell=True)


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
    """ Custom installation: run install_lib command with install-dir==RUNDIR/lib/python/site-packages
     https://docs.python.org/3/install/#custom-installation """

    rundir = os.getcwd()

    pypath = os.path.join(rundir,
                          'lib',
                          'python%d.%d' % sys.version_info[:2],
                          'site-packages')
    os.environ['PYTHONPATH'] = pypath
    if not os.path.isdir(pypath):
        os.makedirs(pypath)


    # copy sys.argv to BUILD_ARGS, route arguments to either configure_mfix or make
    BUILD_ARGS = list(sys.argv[1:])
    valid_make_args = ['-k', '-j']

    global CONFIGURE_ARGS
    CONFIGURE_ARGS = [arg for arg in BUILD_ARGS if arg not in valid_make_args]

    global MAKE_ARGS
    MAKE_ARGS = [arg for arg in BUILD_ARGS if arg in valid_make_args]

    sys.argv[1:] = [
        'install_lib',
        '--install-dir', pypath,
    ]

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

    if platform.system() == 'Windows':
        mfixsolver = os.path.join(rundir, 'mfixsolver.bat')
        template = BAT_TEMPLATE
    else:
        mfixsolver = os.path.join(rundir, 'mfixsolver')
        template = SH_TEMPLATE

    wrapper_contents = template.replace('{PYPATH}', pypath).replace('{PYEXE}', sys.executable)

    with open(mfixsolver, 'wt') as wrapper:
        wrapper.write(wrapper_contents)

    # set executable permission for wrapper
    permissions = os.stat(mfixsolver)
    os.chmod(mfixsolver, permissions.st_mode | 0o111)

SH_TEMPLATE = \
r"""#!/bin/sh

env PYTHONPATH={PYPATH}:${PYTHONPATH} {PYEXE} -m mfixgui.pymfix "$@"
"""

BAT_TEMPLATE = \
r"""@echo on

set PYTHONPATH={PYPATH};%PYTHONPATH%
call "{PYEXE}" -m mfixgui.pymfix %*
"""

if __name__ == '__main__':
    main()
