"""A setuptools based setup module for the MFiX GUI.

See:
https://packaging.python.org/en/latest/distributing.html
http://mfix.netl.doe.gov/
"""

import platform
import shutil
import io
import subprocess
import sys
import zipfile

from glob import glob
from os import makedirs, path, walk, environ

# must import setuptools and cygwinccompiler before numpy.distutils
import setuptools
import distutils.cygwinccompiler

from numpy.distutils.core import setup
from mfixgui.tools.namelistparser import buildKeywordDoc, writeFiles

exec(io.open('mfixgui/version.py').read())

from mfixgui.build_mfixsolver import BuildExtCommand, BuildMfixCommand, make_mfixsolver

HERE = path.abspath(path.dirname(__file__))
NAME = 'mfix'

# Get the long description from the README file
with io.open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

MODEL_DIR = path.join(HERE, 'model')
writeFiles(buildKeywordDoc(MODEL_DIR))

def get_data_files():
    """ walks subdirectories to generate a list of all files that get packaged as data_files """
    data_files = []

    # to run autoreconf and generate build-aux autotools files
    cmd = 'bash configure_mfix'
    subprocess.check_call(cmd, shell=True)

    data_files.append((NAME, ['configure_mfix']))

    subdirs = [
        'build-aux',
        'defaults',
        'model',
        'tutorials',
        'benchmarks',
        'tests',
        'queue_templates',
    ]

    for subdir in subdirs:
        for root, dirs, files in walk(subdir):
            dir_files = []
            for f in files:
                dir_files.append(path.join(root, f))
            data_files.append((path.join(NAME, root), dir_files))

    if platform.system() == 'Windows':
        fortran_dlls = zipfile.ZipFile(path.join('build-aux', 'Win64', 'FORTRAN_DLLS.zip'))
        data_files += fortran_dlls.namelist()
        fortran_dlls.extractall()

    return data_files


def build_doc():
    pandoc_args = ['-s', '--toc', '-N', '-m']
    docs = [('', 'README'),
            ('doc', 'SETUP_GUIDE'),
            ('doc', 'USER_GUIDE'),
            ('doc', 'TUTORIALS')]
    rendered_docs = []
    for src in docs:

        doc_src = path.join(*src) + '.md'

        # doc_dest is useful for previewing docs using build_doc directly
        doc_dest = path.join(*src) + '.html'

        # doc_pkg is the location for the docs distributed in the package
        doc_pkg = path.join('mfixgui', path.join(*src) + '.html')

        subprocess.call(['pandoc', doc_src, '-o', doc_dest] + pandoc_args)

        with io.open(doc_dest, encoding="utf8") as doc:
            data = doc.read()

        # fix links in README to the GUIDES
        data = data.replace('SETUP_GUIDE.md', 'SETUP_GUIDE.html')
        data = data.replace('USER_GUIDE.md', 'USER_GUIDE.html')
        data = data.replace('TUTORIALS.md', 'TUTORIALS.html')

        with io.open(doc_dest, 'w', encoding='utf8') as doc:
            doc.write(data)

        # fix links to images in packaged docs
        data = data.replace('../mfixgui/icons', '../icons')

        with io.open(doc_pkg, 'w', encoding='utf8') as doc:
            doc.write(data)

        rendered_docs.append(path.basename(doc_pkg))

    return rendered_docs

def build_doc_media():
    dest_dir = path.join('mfixgui', 'doc', 'media')
    if not path.isdir(dest_dir):
        makedirs(dest_dir)

    media = []
    for image in glob(path.join('doc', 'media', '*')):
        filename = path.basename(image)
        shutil.copy(image, dest_dir)
        media.append(filename)

    return media

class BuildDocCommand(setuptools.Command):
    """ renders setup guide and user guide from Markdown to HTML """

    description = "build mfix documentation (setup guide and user guide)"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        build_doc()

class TestLoadAllCommand(setuptools.Command):
    """ search for all mfix.dat files and open each in GUI with -t option """

    description = "load every mfix.dat file for testing"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cases = []
        for root, dirs, files in walk('.'):
            if 'mfix.dat' in files:
                cases.append(path.join(root, 'mfix.dat'))

        for case in cases:
            cmd = '%s -m mfixgui.gui -d -linfo -t %s' % (sys.executable, case)
            environ["MFIX_NO_VTK"] = "1"
            subprocess.check_call(cmd, shell=True)


setup(
    name=NAME,

    cmdclass={
        'build_doc': BuildDocCommand,
        'build_ext': BuildExtCommand,
        'build_mfix': BuildMfixCommand,
        'test_load_all': TestLoadAllCommand,
    },

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
        'mfixgui.doc.media',
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
        'qtpy',
    ],

    ext_modules=[make_mfixsolver(),],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'mfixgui.colormaps': ['*'],
        'mfixgui.icons': ['*.png'],
        'mfixgui.tools': ['keyword_args.txt', 'keywordDoc.json', 'keywordList.txt', 'template_data.json'],
        'mfixgui.uifiles': ['*'],
        'mfixgui.widgets': ['burcat.pickle'],
        'mfixgui.doc': build_doc(),
        'mfixgui.doc.media': build_doc_media(),
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('lib/python2.7/site-packages', ['mfix.so']),
    #             ('tutorials', ['tutorials/fluidBed.pdf']),
    # ],
    data_files=get_data_files(),

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'mfix=mfixgui.gui:main',
            'mfixsolver=mfixgui.pymfix:main',
            'build_mfixsolver=mfixgui.build_mfixsolver:main',
        ],
    },
)
