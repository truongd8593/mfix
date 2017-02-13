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

Potential users may find reviewing the Frequently Asked Questions section of the MFIX website useful before downloading the code.

For information on running MFIX, please see the user guide: [USER_GUIDE.md](USER_GUIDE.md)

- If you want to install MFIX with binary packages, see [Installing MFIX](#installing-mfix)
- If you want to build and install MFIX from source , see [Building MFIX](#building-mfix)
- If you want to run the command line version of MFIX from previous versions, see [Building MFIX Solver](#building-mfix-solver)


# Installing MFIX

Download [Anaconda](https://www.continuum.io/downloads)
or [Miniconda](https://conda.io/miniconda.html) for your platform (Windows,
macOS, or Linux). Miniconda is a smaller download because it is a minimal
distribution of Anaconda; either can be used to run MFIX. The 64-bit version is
required; Python 3 version is recommended. Install Anaconda, and at the end of
the install make sure to add the install location to your PATH.

The `conda` command is the package manager for Anaconda, and will be used to
install MFIX dependencies.

At the time of this writing, the default anaconda package for vtk needs to be
installed separately. Also, MFIX only works with numpy 1.11.3 so a specific
version needs to be installed. The pyqt library must be installed manually as well.

## Linux

### Prerequisites
```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt
```

### MFIX install
Download the [latest Linux binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_linux_py3) for MFIX.

<!-- TODO: replace with link to website for actual release -->
```shell
> pip install mfix-17.1-cp35-cp35m-linux_x86_64.whl
```

If you ever want to uninstall MFIX:
```shell
> pip uninstall mfixgui
```

## macOS
### Prerequisites
```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt
```

### MFIX install
Download the [latest Mac binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_windows_py2) for MFIX.
<!-- TODO: replace with link to website for actual release -->

```shell
> pip install mfix-17.1-cp35-cp35m-macosx_10_6_x86_64.whl
```

If you ever want to uninstall MFIX:
```shell
> pip uninstall mfixgui
```

## Windows
### Prerequisites

```shell
C:\> conda install -c menpo vtk=7.0.0
C:\> conda install numpy==1.11.3 pyqt
```

### MFIX install

Download the [latest Windows binaries](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_windows_py2) for MFIX.

<!-- TODO: replace with link to website for actual release -->


```shell
C:\> pip install mfix-17.1-cp27-cp27m-win_amd64.whl
```

If you ever want to uninstall MFIX:
```shell
C:\> pip uninstall mfixgui
```

# Building MFIX

This section describes how to build the mfix Python package.
For running command-line MFIX without the GUI, see the [next section](#building-and-running-mfix-command-line-solver).

You can use `python setup.py bdist_wheel` to build binaries packages to be
installed as described in [Installing MFIX](#installing-mfix), or you can use
`python setup.py install` to install the package in a single command. In either
case, you can uninstall with `pip uninstall mfix`.

## Prerequisites

Building also requires gcc, autoconf, automake, and GNU Make.
Building with other compilers is not yet supported.

Building the package requires numpy 1.11.3 is another build dependency. (The
other dependencies in [Installing MFIX](#installing-mfix) are run-time, but not
build-time depedencies.)

Building from source requires the MFIX source tarball [mfix-17.1.tar.gz](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build%3Asrc)
<!-- TODO: replace with link to website for actual release -->


## Linux

Download and install Anaconda (as described in [Installing MFIX](#installing-mfix).)

Installation instructions for depedencies are distribution-specific. On Ubuntu Linux, they can be installed with:

```shell
> sudo apt install gcc autoconf automake make
```

Install numpy with:
```shell
> conda install numpy==1.11.3
```

Now that the prerequisites are installed, build MFIX with:

```shell
> tar xzf mfix-17.1.tar.gz
> cd mfix-17.1
> python setup.py bdist_wheel
```


## macOS

Download and install Anaconda (as described in [Installing MFIX](#installing-mfix).)

Homebrew is the easiest way to install MFIX build dependencies.
 - Go to http://brew.sh and follow the installation instructions.
 - Once homebrew is installed, install MFIX build dependencies with the command:

```shell
> brew install gcc autoconf automake make gnu-sed
```

Install numpy with:
```shell
> conda install numpy==1.11.3
```

Now that the prerequisites are installed, build MFIX with the following command (make sure the python command used is from Anaconda, not /usr/bin/python).
```shell
> tar xzf mfix-17.1.tar.gz
> cd mfix-17.1
> python setup.py bdist_wheel
```


## Windows

Download and install Anaconda (as described in [Installing MFIX](#installing-mfix).)

MSYS2 (a Unix-like environment for Windows based on Cygwin) is the easiest way to install MFIX build dependencies. The MSYS2 environment can be installed from anaconda
```shell
> conda install numpy==1.11.3
> conda install m2-base m2-autoconf m2-automake-wrapper m2-make m2-tar m2w64-gcc
```

Now that the prerequisites are installed, build MFIX with the following command (make sure the python command used is from Anaconda).
```shell
C:\> tar xzf mfix-17.1.tar.gz
C:\> cd mfix-17.1
C:\> python setup.py bdist_wheel
```


# Building for UDFs

When running cases with User Defined Files (UDFs), it is necessary to build a separate mfixsolver extension module for that case.

This requires building from a source distribution of MFIX ([Building MFIX](#building-mfix)), not installed from a binary distribution.

Assume the source distribution is installed to `MFIX_HOME`.

```shell
> cd my_example_case
> ls *.f
usr0.f   write_usr0.f
> $MFIX_HOME/configure_mfix --python
> make mfixsolver.so
> ls *.so
mfixsolver.so
```

When running `pymfix` from `mfixgui`, the case-specific MFIX solver will override the default mfixsolver implementation.


#	Building and Running MFIX command-line solver

If you are using the MFIX GUI, this section is not needed. However, you may want
to use the command-line version of MFIX (same as in previous versions) if you do
not have a Python distribution installed, or if you are using any of the
following features:

-	SMP/OpenMP
-	DMP/MPI
-	Compilers other than GCC

##     Prerequisites
To build MFIX, the following must be installed on your system. Contact your system administrator for assistance if necessary.

-	Fortran compiler. Commonly available compilers include:
-	GCC (gfortran) version 4.3 and above
-	Intel (ifort) version 11.1 and above

GNU Autoconf version 2.69 or greater is required when building from source code
obtained from the MFIX git repository. This does not apply to the source code
tarball on the MFIX website.


##     Extracting the MFIX directory

MFIX is distributed as a compressed source tar ball named mfix-2017.1.tar.gz. To
decompress and extract the tar file:

```shell
> tar xzf mfixgui-2017.1.tar.gz
```

Here it is assumed that you are in the directory containing the tar ball.
Extracting the tar ball creates a directory named mfix-2017.1 containing the
MFIX and POSTMFIX source codes, tests and tutorials, as well as some additional
documentation and utilities.

MFIX is built using GNU Autoconf, which is a general tool for producing configure scripts for building and installing software on different computer systems.  First, run the shell script configure_mfix to create a Makefile, then run GNU Make to build the MFIX executable. A step-by-step tutorial is presented at the end of this section.

##	Configuring with configure_mfix

This section focuses on the configure_mfix script and the availability of several flags. configure_mfix is a wrapper for the usual GNU Autoconf configure script.
Alias creation (optional)
For convenience, an alias can be created to avoid specifying the path to the MFIX configure script. Assuming the MFIX source was extracted in the home directory, and you are using the C shell, an alias can be created by executing:

```shell
> echo "alias configure_mfix ~/mfix-2017.1/configure_mfix" >> ~/.cshrc
```
This appends the quoted text to the .cshrc file located in the home directory. The alias will take effect the on the next login or after sourcing the .cshrc file. To source the .cshrc file in the current terminal, enter at the prompt:
```shell
> source ~/.cshrc
```
Afterwards, MFIX can be configured from any directory by running the alias, configure_mfix. Users familiar with the creation of aliases may choose other various ways to define an alias (e.g., directly editing the .cshrc or .cshrc_aliases files). Users that use a different shell should create an alias with the respective shell resource file (e.g. bash shell users ~/.bashrc).
Passing arguments to the build script
Arguments may be passed to the build script to specify various options: compiler, optimization flags, SMP/DMP support, and other options. All configuration options can be displayed with:
```shell
> ~/mfix-2017.1/configure_mfix --help
```

The most common arguments are given in the following table.

| Argument          | Description                                    |
| --------          | -----------                                    |
| FC=NAME           | Specify the Fortran compiler command           |
| F77=NAME          | Specify the Fortran 77 compiler command        |
| `FCFLAGS='FLAGS'` | Specify compiler flags                         |
| --dmp             | Enable distributed memory support (MPI)        |
| --smp             | Enable shared memory parallel support (OpenMP) |

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

### Specifying compiler options (optional)

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


##	Building mfix with GNU make

If the configure script successfully created a Makefile (see above), then the next step is to build MFIX by running GNU make command.

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
linking are successful, an executable named mfix along with a few intermediate
build files will be in the current directory.


##	Building MFIX: A step-by-step tutorial

The following example shows how to build and run the fluidbed1 tutorial. This
example assumes you have the GNU Fortran compiler (gfortran) installed.

To begin, go to the fluidbed1 tutorial directory and list its contents:

```shell
> cd ~/mfix-2017.1/tutorials/fluidbed1
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
> ~/mfix-2017.1/configure_mfix FC=gfortran FCFLAGS='-g -O2'
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
