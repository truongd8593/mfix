"""A setuptools based setup module for the MFiX GUI.

See:
https://packaging.python.org/en/latest/distributing.html
http://mfix.netl.doe.gov/
"""

# Always prefer setuptools over distutils
from codecs import open
from distutils.command.build_ext import build_ext
from glob import glob
from os import makedirs, path
from shutil import copyfile

from setuptools import setup, Extension

HERE = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(HERE, 'README'), encoding='utf-8') as f:
    long_description = f.read()

# setup_args = { ... }
# if platform.system() == 'Windows':

class MfixBuildExt(build_ext):
    """Override build_extension to copy the shared library file"""

    def build_extension(self, ext):
        ''' Copies the already-compiled pyd
        '''
        # subprocess.call(["./configure", "--python"])
        # subprocess.call(["make"])
        extpath = path.dirname(self.get_ext_fullpath(ext.name))
        if not path.exists(extpath):
            makedirs(extpath)

        mfixsolver_sharedlib = [p for p in glob('mfixsolver*') if path.isfile(p)]
        if not mfixsolver_sharedlib:
            raise EnvironmentError("setup requires mfixsolver shared library; run './configure --python; make' before 'python setup.py'")
        mfixsolver_sharedlib = mfixsolver_sharedlib[0]
        copyfile(path.join(HERE, mfixsolver_sharedlib), self.get_ext_fullpath(ext.name))

setup(
    name='mfixgui',
    cmdclass={'build_ext': MfixBuildExt},

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.0',

    description='A GUI for the MFiX computational fluid dynamics solver',
    long_description=long_description,

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
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=[
        'mfixgui',
        'mfixgui.doc',
        'mfixgui.icons',
        'mfixgui.tests',
        'mfixgui.tools',
        'mfixgui.widgets',
        'mfixgui.uifiles',
    ],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['qtpy'],

    ext_modules=[Extension('mfixsolver', sources=[])],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'mfixgui.widgets': ['burcat.pickle'],
        'mfixgui.tools': ['keyword_args.txt'],
        'mfixgui.icons': ['*.png'],
        'mfixgui.uifiles': ['*'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('lib/python2.7/site-packages', ['mfix.so']),
    #             ('tutorials', ['tutorials/fluidBed.pdf']),
    # ],
    data_files=[
        ('model', glob('model/*.f')),
        ('model/cartesian_grid', glob('model/cartesian_grid/*.f')),
        ('model/des', glob('model/des/*.f')),
    ],

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
