"""A setuptools based setup module for the MFiX GUI.

See:
https://packaging.python.org/en/latest/distributing.html
http://mfix.netl.doe.gov/
"""

import codecs
import errno
import platform
import shutil
import subprocess
import tempfile
import zipfile
from glob import glob
from os import makedirs, path, walk

# must import setuptools and cygwinccompiler before numpy.distutils
import setuptools
import distutils.cygwinccompiler

from numpy.distutils.command.build_ext import build_ext
from numpy.distutils.core import Extension, setup

from mfixgui.tools.namelistparser import buildKeywordDoc, writeFiles

exec(codecs.open('mfixgui/version.py').read())

HERE = path.abspath(path.dirname(__file__))
NAME = 'mfix'

# Get the long description from the README file
with codecs.open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

F90_TMP = tempfile.mkdtemp()
MODEL_DIR = path.join(HERE, 'model')
writeFiles(buildKeywordDoc(MODEL_DIR))

def get_data_files():
    """ walks subdirectories to generate a list of all files that get packaged as data_files """
    data_files = []

    for subdir in ['defaults', 'model', 'tutorials', 'benchmarks', 'tests', 'queue_templates']:
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

    makedirs(path.join(F90_TMP, 'dmp_modules'))
    makedirs(path.join(F90_TMP, 'des'))
    for src in pymfix_src:
        shutil.copyfile(path.join('model', src), path.join(F90_TMP, src)+'90')

    return [path.join(F90_TMP, s)+'90' for s in pymfix_src]

def cleanup_tmp():
    # clean tempdir
    try:
        shutil.rmtree(F90_TMP)
    except OSError as ex:
        if ex.errno != errno.ENOENT:
            raise


DEFAULT_CC = 'gcc'
DEFAULT_FC = 'gfortran'
DEFAULT_CFLAGS = '-O2 -fPIC'
DEFAULT_FCFLAGS = '-O2 -fPIC'


def make_mfixsolver():
    configure_args = get_configure_args(DEFAULT_CC, DEFAULT_FC, DEFAULT_CFLAGS, DEFAULT_FCFLAGS)
    build_dir = path.join('build', configure_args.replace(' ', '_').replace('=', '_')).replace('"', '_').replace('__', '_')

    return Extension(name='mfixsolver',
                     sources=get_pymfix_src(),
                     extra_f90_compile_args=['-cpp'],
                     module_dirs=[path.join(build_dir, 'model')],
                     extra_objects=[
                         '.build/read_database.o',
                         '.build/read_namelist.o',
                         path.join(build_dir, 'build-aux/libmfix.a'),
                     ])


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

        with codecs.open(doc_dest, encoding="utf8") as doc:
            data = doc.read()

        # fix links in README to the GUIDES
        data = data.replace('SETUP_GUIDE.md', 'SETUP_GUIDE.html')
        data = data.replace('USER_GUIDE.md', 'USER_GUIDE.html')
        data = data.replace('TUTORIALS.md', 'TUTORIALS.html')

        with codecs.open(doc_dest, 'w', encoding='utf8') as doc:
            doc.write(data)

        # fix links to images in packaged docs
        data = data.replace('../mfixgui/icons', '../icons')

        with codecs.open(doc_pkg, 'w', encoding='utf8') as doc:
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


def get_configure_args(cc, fc, cflags, fcflags):
    return 'CC="%s" FC="%s" CFLAGS="%s" FCFLAGS="%s"' % (cc, fc, cflags, fcflags)


class BuildMfixCommand(setuptools.Command):
    """ builds libmfix (Python version agnostic) """
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

        # should work Linux/Mac/Windows as long as bash is in PATH
        cmd = 'bash ./configure_mfix %s' % get_configure_args(self.cc, self.fc, self.cflags, self.fcflags)
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            raise EnvironmentError("Failed to configure_mfix correctly")

        returncode = subprocess.call("make", shell=True)
        if returncode != 0:
            raise EnvironmentError("Failed to build mfix correctly")



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

setup(
    name=NAME,

    cmdclass={
        'build_doc': BuildDocCommand,
        'build_ext': BuildExtCommand,
        'build_mfix': BuildMfixCommand,
    },

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=get_git_revision_short_hash(),

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
            'pymfix=mfixgui.pymfix:main',
        ],
    },
)

cleanup_tmp()
