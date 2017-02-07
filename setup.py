"""A setuptools based setup module for the MFiX GUI.

See:
https://packaging.python.org/en/latest/distributing.html
http://mfix.netl.doe.gov/
"""

# Always prefer setuptools over distutils
from codecs import open
from distutils.command.build_ext import build_ext
from glob import glob
from os import makedirs, path, walk
from numpy.f2py import f2py2e

# must import setuptools before numpy.distutils
import setuptools
from numpy.distutils.core import Extension, setup

from zipfile import ZipFile
import numpy as np
import platform
import tempfile
from shutil import copyfile

import subprocess
import sys

from mfixgui.tools.namelistparser import buildKeywordDoc, writeFiles

exec(open('mfixgui/version.py').read())

HERE = path.abspath(path.dirname(__file__))
NAME = 'mfix'

# Get the long description from the README file
with open(path.join(HERE, 'README'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

MODEL_DIR = path.join(HERE, 'model')
writeFiles(buildKeywordDoc(MODEL_DIR))

data_files = []

for subdir in ['defaults', 'model', 'tutorials', 'benchmarks', 'tests']:
    for root,dirs,files in walk(subdir):
        dir_files = []
        for f in files:
            dir_files.append(path.join(root,f))
        data_files.append((path.join('share', NAME, root), dir_files))

if platform.system() == 'Windows':
    zf = ZipFile(os.path.join('build-aux','Win64','FORTRAN_DLLS.zip'))
    data_files += zf.namelist()
    zf.extractall()

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

f90 = tempfile.mkdtemp()
makedirs(path.join(f90, 'dmp_modules'))
makedirs(path.join(f90, 'des'))
for s in pymfix_src:
    copyfile(path.join('model', s), path.join(f90, s)+'90')

pymfix_src = [ path.join(f90, s)+'90' for s in pymfix_src ]

mfixsolver = Extension(name = 'mfixsolver',
                       sources = pymfix_src,
                       extra_f90_compile_args = ['-cpp'],
                       module_dirs = ['build/_/model'],
                       extra_objects = [
                           '.build/read_database.o',
                           '.build/read_namelist.o',
                           'build/_/build-aux/libmfix.a',
                       ]
)

setup(
    name='mfixgui',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='A GUI for the MFiX computational fluid dynamics solver',
    long_description=LONG_DESCRIPTION,

    # The project's main homepage.
    url='http://mfix.netl.doe.gov/',

    # Author details
    author='Multiflow Science Group at NETL',
    author_email='mfix-gui@mfix.netl.doe.gov',
    platforms=["any"],

    # Choose your license
    license='public domain',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Computational Fluid Dynamics :: GUI',

        # Pick your license as you wish (should match "license" above)
        'License :: public domain',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=[
        'mfixgui',
        'mfixgui.colormaps',
        'mfixgui.doc',
        'mfixgui.icons',
        'mfixgui.tests',
        'mfixgui.tools',
        'mfixgui.uifiles',
        'mfixgui.widgets',
    ],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        'flask',
        'numpy==1.11.3',
        'packaging',
        'qtpy',
    ],

    ext_modules=[mfixsolver,],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'mfixgui.widgets': ['burcat.pickle'],
        'mfixgui.tools': ['keyword_args.txt', 'keywordDoc.json', 'keywordList.txt'],
        'mfixgui.icons': ['*.png'],
        'mfixgui.uifiles': ['*'],
        'mfixgui.colormaps': ['*'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('lib/python2.7/site-packages', ['mfix.so']),
    #             ('tutorials', ['tutorials/fluidBed.pdf']),
    # ],
    data_files=data_files,

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'mfixgui=mfixgui.gui:main',
            'pymfix=mfixgui.pymfix:main',
        ],
    },
)

# clean tempdir
try:
    shutil.rmtree(f90)
except OSError as e:
    if e.errno != errno.ENOENT:
        raise
