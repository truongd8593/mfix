# MFIX Setup Guide

This document explains how to install or build binary packages for MFIX 17.1.

## About MFIX

MFIX is an open-source multiphase flow solver and is free to download and use. A
one-time no-cost registration is required prior to downloading the source code.
To register, go to https://mfix.netl.doe.gov/ and click on the "Register" button
in the upper right corner. Once you have read the notice, you can submit your
application by clicking on "REGISTER." After your application has been reviewed
and accepted, you will receive an email notification and instructions on how to
download the code. Please allow for 2-3 business days for your registration to
be processed.

Potential users may find reviewing the Frequently Asked Questions section of the
MFIX website useful before downloading the code.

For information on running MFIX, please see the user guide: [USER_GUIDE.html](USER_GUIDE.html)


To run MFIX simulations, you need to have the MFIX gui (referred to as mfixgui) and the MFIX flow solver (referred to as mfixsolver) installed. The MFIX flow solver can be the default solver (no UDF), or the customed solver (with UDF, for example for chemically reacting flows that required coded reaction rates).It is recommended to start installing MFIX from the provided binary packages on a given platform (Linux, MAC, or Windows OS). The binary packages include mfixgui and the default mfixsolver. A list of known operating systems where MFIX has been sucessfully installed is given below.

## Known operating systems where MFIX can be installed/run:

- Linux: Ubuntu (insert versions), OpenSuse (insert versions)
- macOS: insert OS versions
- Windows: Windows 7, Windows 10 (need to verify that)

Once users have gained experience with generic simulations (no UDFs), advanced simulations (with UDFs) can be run after building a custom solver, that take into account the UDFs. Building a custom solver require additional packages and is more involved than installing a binary package. Knowledge of entering commnands at the prompt is required for all platforms (Linux, MAC and Windows OS).

It is also possible to build and install MFIX from source, but this will typically be suitable for advanced users and developpers.

Finally, running the mfix solver without the GUI is still possible (Linux environment only). The solver must be built from source, and the simulation is set by editing a text file. This correspond to the way MFIX was run in previous versions (versions 2016-1 and earlier).

To recap:
- If you want to install MFIX from binary packages, see [Installing MFIX](#installing-mfix).
  This is the recommended starting point to learn and use MFIX for most users.
  The default solver will be installed which will allow basic simulations (no UDFs).
- If you already have MFIX installed and want to build a custom solver (with UDFs), see [build a custom interactive mfixsolver](#building-custom-mfixsolver).
- If you want to build the MFIX binary from source , see [Building MFIX](#building-mfix).
- If you want to run the command line version of MFIX without the GUI (Linux only, similar to previous versions), see



# Prerequisites

Using the GUI requires a Python environment, whether installing the package from binaries, or building from source. The recommended Python environment is the Anaconda Python distribution.

Download [Miniconda](https://conda.io/miniconda.html) for your platform
(Linux, macOS, or Windows). (Miniconda is a minimal distribution of Anaconda.)
The 64-bit version is required; Python 3 version is recommended. Install
Anaconda, and at the end of the install make sure to add the install location to
your PATH.

The `conda` command is the package manager for Anaconda, and will be used to
install MFIX dependencies.

At the time of this writing, the default anaconda package for vtk needs to be
installed separately. Also, MFIX only works with numpy 1.11.3 so a specific
version needs to be installed.


## Linux


Installation instructions for dependencies are distribution-specific. On Ubuntu
Linux, they can be installed with:

```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt qtpy pandoc
```

To install MFIX from binaries, proceed to [Installing MFIX](#installing-mfix)

To build MFIX from source or build the custom solver, install the following dependencies with:

```shell
> sudo apt install gcc gfortran autoconf automake make
```

and proceed to [Building MFIX](#building-mfix) or [build a custom interactive mfixsolver](#building-custom-mfixsolver).


## macOS


Homebrew is the easiest way to install MFIX build dependencies. Go
to [the Homebrew website](http://brew.sh) and follow the installation
instructions.

Once homebrew is installed, install MFIX  dependencies with the commands:

```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt qtpy pandoc
> pip install pypandoc
```

To install MFIX from binaries, proceed to [Installing MFIX](#installing-mfix)


To build MFIX from source or build the custom solver, install the following dependencies with:

```shell
> brew install gcc autoconf automake make gnu-sed
```

and proceed to [Building MFIX](#building-mfix) or [build a custom interactive mfixsolver](#building-custom-mfixsolver).


To install MFIX from binaries, proceed to [Installing MFIX](#installing-mfix)
To build MFIX from source, proceed to [Building MFIX](#building-mfix).


## Windows

Open  the Anaconda terminal, and enter the following commands:


```shell
C:\> conda install -c menpo vtk=7.0.0
C:\> conda install numpy==1.11.3 pyqt qtpy pandoc
```


To install MFIX from binaries, proceed to [Installing MFIX](#installing-mfix)


To build MFIX from source or build the custom solver:

### Install MSYS2 environment (a Unix-like environment for Windows based on Cygwin):


```shell
C:\> conda install m2-base m2-autoconf m2-automake-wrapper m2-make m2-tar m2w64-gcc m2w64-gcc-fortran
```


### Patch distutils for Python 3

If you are using Python 3, you need to make edits to the following file.

<!-- TODO: Explain what ANACONDA_HOME is (maybe give an example) -->

In `ANACONDA\_HOME/envs/py3/Lib/site-packages/numpy/distutils/misc_util.py:402` add the line

<!-- This patch is in numpy master branch on Github. We can probably remove this patch once numpy 1.13 is released. -->

```python
'1900': 'msvcr140',  # MVCS 14 (aka 2014)
```

to:

```python
        lib = {'1300': 'msvcr70',    # MSVC 7.0
               '1310': 'msvcr71',    # MSVC 7.1
               '1400': 'msvcr80',    # MSVC 8
               '1500': 'msvcr90',    # MSVC 9 (VS 2008)
               '1600': 'msvcr100',   # MSVC 10 (aka 2010)
               '1900': 'msvcr140',   # MSVC 14 (aka 2014)
              }.get(msc_ver, None)
```

### Set Windows PATH

- Open Windows Control Panel
- Open Environment Variables dialog
- Add the following paths to beginning of your PATH environment variable:
  - `ANACONDA_HOME\Library\mingw-w64\bin`
  - `ANACONDA_HOME\Library\usr\bin`

### Run bash

From the command prompt (CMD.exe), run `bash`. The rest of this document assumes
all commands are run from a bash prompt.

```shell
C:\bash
user@hostname /c
$
```

and proceed to [Building MFIX](#building-mfix)or [build a custom interactive mfixsolver](#building-custom-mfixsolver).


# Installing MFIX

If you built MFIX as described in [Building MFIX](#building-mfix), use those binaries.

Otherwise, download the latest binaries for your platform:

- [latest Linux binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build_linux_py3)
- [latest Mac binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build_windows_py3)
- [latest Windows binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build_windows_py3)

<!-- TODO: replace with links to website for actual release -->
<!-- TODO: Navigate to location where the downloaded binary is located -->

```shell
> pip install mfix-17.1-<PYTHON_VERSION>-<PLATFORM>.whl
```

You are now ready to set up and run MFIX simulation!


# Running MFIX simulations

## From the GUI with default interactive solver (preferred for beginers)

This will use the mfixsolver Python library installed with the package. Only the
location of pymfix needs to be defined in the GUI. You can pause, unpause, stop,
or get info from the solver.

The GUI is launched from the prompt with:
```shell
> mfix
```

## From the GUI, with user interactive solver

Here you use a source distribution of MFIX
to [build a custom interactive mfixsolver](#building-custom-mfixsolver) located
in the project directory. When running `mfix`, in
the [Run Dialog](USER_GUIDE.html#run-dialog) select the mfixsolver.so file you
have just built. You can pause, unpause, stop, or get info from the solver.

The GUI is launched from the prompt with:
```shell
> mfix
```

## From the GUI, with user non-interactive mfixsolver executable (should we keep this?)

 The solver is built from [Building MFIX Solver](#building-mfixsolver).
 You cannot pause, unpause, stop, or get info from the solver.

The GUI is launched from the prompt with:
```shell
> mfix
```

## From the command line with a custom mfixsolver executable (without the GUI)

This is the same as in previous MFIX releases. The solver is built
from [Building MFIX Solver](#building-mfixsolver). This option does not require
Anaconda Python.





## Uninstalling MFIX

If you ever want to uninstall MFIX:
```shell
> pip uninstall mfix
```









# Building custom mfixsolver

For some cases, you may want to use a custom mfixsolver. For instance, when
running cases with User Defined Files (UDFs), it is necessary to build a
separate mfixsolver extension module for that case.

This requires downloading a source distribution of MFIX ([Building MFIX](#building-mfix)).

Assume the source distribution tarball is extracted to `MFIX_HOME`.

On Linux and Mac, there will be a default system python command
`/usr/bin/python`. Make sure you are NOT using this Python by explicitly specifying
which python distribution you are installing it with using PYTHON_BIN.

```shell
> cd my_example_case
> ls *.f
usr0.f   write_usr0.f
> $MFIX_HOME/configure_mfix --python PYTHON_BIN=$HOME/miniconda3/python3.5 CC=gcc
> make mfixsolver.so
> ls *.so
mfixsolver.so
```

<!-- FIXME: mention mfixsolver.so vs mfixsolver.pyd -->

When running `mfixgui`, in the [Run Dialog](USER_GUIDE.html#run-dialog) select the mfixsolver.so file you have just built.


# Building mfixsolver

This option is currently limited to Linux environment.

If you are using the MFIX GUI, this section is not needed. However, you may want
to use the command-line version of MFIX (same as in previous versions) if you do
not have a Python distribution installed, or if you are using any of the
following features:

-	SMP/OpenMP
-	DMP/MPI
-	Compilers other than GCC

## Prerequisites

MFIX is built using GNU Autoconf, which is a general tool for producing
configure scripts for building and installing software on different computer
systems. First, run the shell script configure_mfix to create a Makefile, then
run GNU Make to build the MFIX executable. A step-by-step tutorial is presented
at the end of this section.

To build MFIX, the following must be installed on your system. Contact your system administrator for assistance if necessary.

-	Fortran compiler. Commonly available compilers include:
-	GCC (gfortran) version 4.3 and above
-	Intel (ifort) version 11.1 and above

GNU Autoconf version 2.69 or greater is required when building from source code
obtained from the MFIX git repository. This does not apply to the source code
tarball on the MFIX website.

## Extracting the MFIX directory

MFIX is distributed as a compressed source tar ball named mfix-17.1.tar.gz. To
decompress and extract the tar file:

```shell
> tar xzf mfix-17.1.tar.gz
```

Here it is assumed that you are in the directory containing the tar ball.
Extracting the tar ball creates a directory named mfix-17.1 containing the
MFIX and POSTMFIX source codes, tests and tutorials, as well as some additional
documentation and utilities.

## Configuring with configure_mfix

### Passing arguments to the build script

Arguments may be passed to the build script to specify various options:
compiler, optimization flags, SMP/DMP support, and other options. All
configuration options can be displayed with:

```shell
> ~/mfix-17.1/configure_mfix --help
```

The most common arguments are given in the following table.

| Argument          | Description                                    |
| --------          | -----------                                    |
| FC=NAME           | Specify the Fortran compiler command           |
| F77=NAME          | Specify the Fortran 77 compiler command        |
| `FCFLAGS='FLAGS'` | Specify compiler flags                         |
| --dmp             | Enable distributed memory support (MPI)        |
| --smp             | Enable shared memory parallel support (OpenMP) |

##	Building mfix with GNU make

If the configure script successfully created a Makefile (see above), then the
next step is to build MFIX by running GNU make command.

```shell
> make
```

The -j option may be used build in parallel which may decrease compile time.
```shell
> make -j
```

Note that on some systems parallel builds may fail due to file dependencies
(e.g., file1 depends of file2 which has not been compiled yet). Typically,
running the make command again will overcome these errors. If compiling and
linking are successful, an executable named mfixsolver along with a few intermediate
build files will be in the current directory.

## Building custom solver for DMP (MPI)

DMP support is only tested on Linux.

DMP support `--dmp` can be combined with interactive support `--python`, but has
only been tested for GCC.

### Prerequisites

```shell
> apt-get install libopenmpi-dev openmpi-bin
```

### Building solver with MPI support


```shell
> configure_mfix --dmp
> make
```

## Building custom solver with non-GCC compilers

MFIX can be built with other compilers, but interactive features `--python` are not supported at this time.

```shell
> configure_mfix --dmp
> make
```

Specifying the Fortran compiler (optional but recommended)

The Fortran compiler is specified by passing the FC argument to
`configure_mfix`. If the FC argument is not specified, the script will search
the environment for a Fortran compiler (usually gfortran). The `configure_mfix`
script will test the Fortran compiler with any specified flags (covered next).
If the compiler is not found or if the compiler gives an error, the Makefile
will not be generated and detailed error messages will be in config.log. When
seeking help on the mailing list for build issues, include the config.log file
for assistance.

The Fortran 77 compiler is specified by the F77 argument to the configure_mfix
script. By default, F77 is set to the same compiler specified by the FC
argument. However, on some systems F77 may need to be defined separately.

Common compilers are given in the following table.

| Compiler | Description                                        |
| -------- | -----------                                        |
| gfortran | GNU Fortran compiler (serial/smp)                  |
| ifort    | Intel Fortran compiler (serial/smp)                |
| mpif90   | Generic MPI Fortran wrapper (dmp/combined smp-dmp) |
| mpifort  | Generic MPI Fortran wrapper (dmp/combined smp-dmp) |
| mpiifort | Intel MPI Fortran wrapper (dmp/combined smp-dmp)   |

Older versions commonly use the mpif90 command.


## Specifying custom compiler flags

Compiler flags are specified by passing the FCFLAGS argument to
`configure_mfix`. If the FCFLAGS argument is not specified, the compiler
defaults are used in building the executable. The `configure_mfix` script will
test the compiler with the specified flags. If the compiler is not found or if
the compiler gives an error with the given flags, the Makefile will not be
generated and detailed error messages will be in config.log. When seeking help
on the mailing list for build issues, include the config.log file for
assistance.

Common compiler flags for GNU Fortran are given in the following table.

 GNU  Fortran

| Option                      | Description                                                                                                                                  |
| --------                    | -----------                                                                                                                                  |
| -g                          | Produce debugging information                                                                                                                |
| -O0, -O1, -O2, -O3          | Optimization level (refer to the compiler manual for specific information about each optimization level). Typically, the following are used: |
| -O0	debugging             |                                                                                                                                              |
| -O2	production executables|                                                                                                                                              |
| -fcheck=all                 | Generates runtime checks useful for debugging.                                                                                               |
| -fbacktrace                 | Proved a full backtrace when a runtime error is encountered for debugging.                                                                   |

Common compiler flags for Intel Fortran are given in the following table.

 Intel  Fortran

| Option             | Description                                                                                                                                  |
| --------           | -----------                                                                                                                                  |
| -g                 | Produce debugging information                                                                                                                |
| -O0, -O1, -O2, -O3 | Optimization level (refer to the compiler manual for specific information about each optimization level). Typically, the following are used: |
| -O0                | debugging                                                                                                                                    |
| -O2                | production executables                                                                                                                       |
| -check all         | Generates runtime checks useful for debugging.                                                                                               |
| -debug             | Generates complete debugging information                                                                                                     |


## Building MFIX: A step-by-step tutorial

The following example shows how to build and run the fluidbed1 tutorial. This
example assumes you have the GNU Fortran compiler (gfortran) installed.

To begin, go to the fluidbed1 tutorial directory and list its contents:

```shell
> cd ~/mfix-17.1/tutorials/fluidbed1
> ls
mfix.dat
```

As shown this folder contains the input file, mfix.dat, which specifies the
simulation setup. This file is only needed to compile MFIX when user-defined
reaction files are present (e.g., usr_rates.f). However, the mfix.dat file must
be present when running MFIX, the last step of this tutorial.

The directory where MFIX will be run is called the run directory (RUNDIR). Here
the run directory is fluidbed1, and it currently contains only the mfix.dat
file.

Running the `configure_mfix` script and specifying the compiler and optimization flags:

```shell
> ../../configure_mfix FC=gfortran FCFLAGS='-g -O2'
```

The above commands run `configure_mfix` by specifying the relative path (two
levels up). Alternatively, we could specify the absolute path:

```shell
> ~/mfix-17.1/configure_mfix FC=gfortran FCFLAGS='-g -O2'
```

or use the alias created in section 4.3.1:

```shell
> configure_mfix FC=gfortran FCFLAGS='-g -O2'
```

If the configure script successfully created a Makefile, build MFIX by running
make.

```shell
> make
```

The build may take several minutes to complete. If the build is successful, the
executable will be in the run directory.

```shell
> ls
mfix.dat    mfix
Finally, the simulation is started by entering:

```shell
> ./mfixsolver
```

# Appendix - Building MFIX

This appendix describes how to build MFIX packages for distribution. Most users will not need to do this.

## Building MFIX

This section describes how to build the mfix Python package. At the end of the build process, you will have a distribution wheel to the ones available for download.


Building from source requires the MFIX source tarball:

- [mfix-17.1.tar.gz](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build%3Asrc)
<!-- TODO: replace with link to website for actual release -->

Build MFIX with:

```shell
> tar xzf mfix-17.1.tar.gz
> cd mfix-17.1
> python setup.py bdist_wheel
```

You can now install with the `dist/mfix*.whl` file as described in [Installing MFIX](#installing-mfix).
