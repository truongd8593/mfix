# Introduction

This document explains how to install or build the binary packages for MFIX 17.1


# Installing MFIX

## Prerequisites

Download Anaconda (https://www.continuum.io/downloads) or Miniconda
(https://conda.io/miniconda.html) for your platform (Windows, macOS, or Linux).
Miniconda is a smaller download because it is a minimal distribution of
Anaconda; either can be used to run MFIX. The 64-bit version is required; Python
3 version is recommended. Install Anaconda, and at the end of the install make
sure to add the install location to your PATH.

The `conda` command is the package manager for Anaconda, and will be used to
install MFIX dependencies.

At the time of this writing, the default anaconda package for vtk needs to be
installed separately. Also, MFIX only works with numpy 1.11.3 so a specific
version needs to be installed. The pyqt library must be installed manually as well.

### Linux
```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt
```

### macOS
```shell
> conda install -c menpo vtk=7.0.0
> conda install numpy==1.11.3 pyqt
```

### Windows

```shell
C:\> conda install -c menpo vtk=7.0.0
C:\> conda install numpy==1.11.3 pyqt
```


## Install MFIX

### Linux

Download the latest version of the package for your platform:

TODO: replace with link to website for actual release
https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_linux_py3

```shell
> pip install mfix-17.1-cp35-cp35m-linux_x86_64.whl
```

### macOS
Download the latest version of the package for your platform:

TODO: replace with link to website for actual release
https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_mac_py3

```shell
> pip install mfix-17.1-cp35-cp35m-macosx_10_6_x86_64.whl
```

### Windows
Download the latest version of the package for your platform:

TODO: replace with link to website for actual release
https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build_windows_py2


```shell
C:\> pip install mfix-17.1-cp27-cp27m-win_amd64.whl
```

## Uninstalling MFIX

MFIX can be uninstalled with:

(Mac or Linux)
```shell
> pip uninstall mfixgui
```

or

(Windows)
```shell
C:\> pip uninstall mfixgui
```


# Building MFIX

## Prerequisites

Building also requires gcc, autoconf, automake, and GNU Make.
Building with other compilers is not yet supported. For running command-line MFIX without the GUI, see Section 6.

Building MFIX from source requires the source tarball mfix-17.1.tar.gz.
TODO: replace with link to website for actual release
https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/gui/download?job=build%3Asrc


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


# Running MFIX


## Linux
Start MFIX with:
```shell
> mfixgui
```


## macOS
Start MFIX with:

```shell
> mfixgui
```


## Windows
Start MFIX with:

```shell
C:\> mfixgui
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
