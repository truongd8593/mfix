# MFIX Setup Guide

This document explains how to install or build binary packages for MFIX 17.1.

## About MFIX

MFIX is an open-source multiphase flow solver and is free to download and use. A
one-time no-cost registration is required prior to downloading the source code.
To register, go to [https://mfix.netl.doe.gov/](https://mfix.netl.doe.gov/) and click on the "Register" button
in the upper right corner. Once you have read the notice, you can submit your
application by clicking on "REGISTER." After your application has been reviewed
and accepted, you will receive an email notification and instructions on how to
download the code. Please allow for 2-3 business days for your registration to
be processed.

Potential users may find reviewing the Frequently Asked Questions section of the
MFIX website useful before downloading the code.

For information on running MFIX, please see the user guide: [USER_GUIDE.html](USER_GUIDE.html)

To run MFIX simulations, you need to have the MFIX gui (referred to as mfix)
and the MFIX flow solver (referred to as mfixsolver) installed. The MFIX flow
solver can be the default solver (no UDF), or the custom solver (with UDF, for
example for chemically reacting flows that required coded reaction rates).It is
recommended to start installing MFIX from the provided binary packages on a
given platform (Linux, MacOS, or Microsoft Windows). The binary packages include mfix
and the default mfixsolver. A list of known operating systems where MFIX has
been successfully installed is given below.

## Known operating systems where MFIX can be installed/run:

- Linux: tested on Ubuntu 16.04, OpenSUSE 10.4
- macOS: tested on 10.12 Sierra
- Windows: tested on Windows 7 and Windows 10

Once users have gained experience with generic simulations (no UDFs), advanced
simulations (with UDFs) can be run after building a custom solver, that take
into account the UDFs. Building a custom solver require additional packages and
is more involved than installing a binary package. Knowledge of entering
commands at the prompt is required for all platforms (Linux, MAC and Windows
OS).

It is also possible to build and install MFIX from source, but this will
typically be suitable for advanced users and developers.

Finally, running mfixsolver without the GUI is still possible (Linux environment
recommended). The solver must be built from source, and the simulation is set by
editing a text file. This corresponds to the way MFIX has been run in previous
versions (versions 2016-1 and earlier).

To recap:

- If you want to install MFIX from binary packages, see [Installing MFIX](#installing-mfix) and skip [installing build dependencies](#install-build-dependencies).
  This is the recommended starting point to learn and use MFIX for most users.
  The default solver will be installed which will allow basic simulations (no UDFs).
- If you have MFIX installed and want to build a custom solver (with UDFs), you need to [install build dependencies](#install-build-dependencies) and then [build a custom interactive mfixsolver](#building-custom-mfixsolver).
- If you are a developer and want to build and install MFIX from source, see [Building MFIX](#building-mfix-for-developers).
- If you want to run the command line version of MFIX without Python and without the GUI (Linux recommended, similar to previous versions), see [building a non-interactive solver](#running-non-interactive-fortran-only-solver).

# Installing MFIX

Using the GUI requires a Python environment, whether installing the package from
binaries, or building from source. The recommended Python environment is the
Anaconda Python distribution.

## Install Miniconda or Anaconda

Download [Miniconda](https://conda.io/miniconda.html) for your platform
(Linux, macOS, or Windows). (Miniconda is a minimal distribution of Anaconda.)
The 64-bit version is required; Python 3 version is recommended.

Install Anaconda or Miniconda, and at the end of the install make sure to add
the install location to your PATH.

The `conda` command is the package manager for Anaconda, and will be used to
install MFIX dependencies. On Windows, use the "Anaconda Prompt" installed
under the Start Menu to make sure the `conda` command is in your PATH.

Windows
: Under the Start Menu, type "Anaconda Prompt", select the Anaconda Prompt.

Linux or macOS
: Exit your shell and start a new shell to ensure `conda` is in your `PATH`.

##  Install build dependencies

***(Optional but Recommended)***

To build MFIX from source or build the custom solver, you will also need to
install the build dependencies, primarily GNU Fortran. (Other Fortran compilers
may be supported in the future.)

### Linux

Installation instructions for build dependencies are distribution-specific. On
Ubuntu Linux, they can be installed with:

```shell
> sudo apt install gcc gfortran autoconf automake make
```

### macOS

Homebrew is the easiest way to install MFIX build dependencies. Go
to [the Homebrew website](http://brew.sh) and follow the installation
instructions.

Once homebrew is installed, install MFIX build dependencies with the command:

```shell
> brew install gcc autoconf automake make gnu-sed
```

### Windows

On Windows, `conda` can be used to install both MSYS2 (a Unix-like environment for
Windows based on Cygwin) and MinGW (a minimalist GCC compiler targeting
Windows).


```shell
### Anaconda Prompt ###
C:\> conda install m2-base m2-autoconf m2-automake-wrapper m2-make m2w64-toolchain
```

## Install MFiX

Install MFIX using the following command: **TODO replace with URL for actual release**

```shell
> conda install -c http://condatest:condatest12345@aeolustec.com/dist mfix
```



You are now ready to set up and run MFIX simulation!

# Uninstalling MFIX

If you ever want to uninstall MFIX:
```shell
> conda uninstall mfix
```

# Running MFIX

On Windows, click on the "MFIX" logo on the Desktop or Start Menu to launch the GUI.

On Linux or Mac, The GUI is launched from the prompt with:
```shell
> mfix
```

For information on running MFIX, please see the user guide: [USER_GUIDE.html](USER_GUIDE.html)

# Running MFIX simulations

## From the GUI with default interactive solver (recommended for beginners)

This will use the mfixsolver Python library installed with the package. Only the
location of mfixsolver needs to be defined in the GUI. You can pause, unpause, stop,
or get info from the solver.

The GUI is launched from the prompt with:
```shell
> mfix
```

## From the GUI, with custom interactive solver

Here you [build a custom interactive mfixsolver](#building-custom-mfixsolver) in
the project directory. When running `mfix`, in
the [Run Dialog](USER_GUIDE.html#run-dialog) select the mfixsolver executable
you have just built.

You can pause, unpause, stop, or get info from the solver.

The GUI is launched from the prompt with:
```shell
> mfix
```

## From the command line (without the GUI)

Here you [build a custom interactive mfixsolver](#building-custom-mfixsolver) in
the project directory.

The default solver is launched in the project directory with:
```shell
> mfixsolver --help
> mfixsolver -f <RUN_NAME>.mfx
```

Note that the default `mfixsolver` should be in your PATH.

The custom solver is launched in the project directory with:
```shell
> build_mfixsolver
> ./mfixsolver -f <RUN_NAME>.mfx
```

Note that the custom `mfixsolver` is in your project directory.

The solver is built from [Building MFIX Solver](#building-custom-mfixsolver). This option does not require
Python.

## From the command line, with non-interactive mfixsolver

The solver is built as in previous version of MFiX.
See [the Appendix](#running-non-interactive-fortran-only-solver). This does not
require Python, just a Fortran compiler.

You cannot pause, unpause, stop, or get info from the solver.

```shell
> ./mfixsolver -f <RUN_NAME>.mfx
```


# Building custom mfixsolver

For some cases, you may want to use a custom mfixsolver. For instance, when
running cases with User Defined Functions (UDFs), it is necessary to build a
separate mfixsolver extension module for that case.

For example, suppose you have a project `MY_EXAMPLE.mfx` with UDFs `usr0.f` and
`write_usr0.f`:

```shell
> cd my_example
> ls *
MY_EXAMPLE.mfx
usr0.f
write_usr0.f
> build_mfixsolver
...build output...
> ls  *
MY_EXAMPLE.mfx
usr0.f
write_usr0.f
mfixsolver
build/
lib/
```

The `build_mfixsolver` command creates a wrapper script `mfixsolver`, that runs
the case-specific MFIX solver (which is installed in the `lib` directory).

The `build` directory can be deleted, and the custom `mfixsolver` will still
run. However, leaving the `build` directory will greatly speed up rebuilds if
you want to edit the UDFs and run `build_mfixsolver` again.

From the GUI, in the [Run Dialog](USER_GUIDE.html#run-dialog), select the
mfixsolver file you have just built.

From the command line, run the custom solver with:

```shell
> ./mfixsolver -f MY_EXAMPLE.mfx
```

### Running mfixsolver from the command line




# Appendix

## Building MFIX (for developers)

This appendix describes how developers build MFIX packages for distribution. Most users will not need to do this.

- Get a source distribution of MFIX

<!-- TODO: replace with link to website for actual release -->
Building from source requires the MFIX source tarball: [mfix-17.1.tar.gz](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build%3Asrc). Extract the tarball and go into the top level source directory:

```shell
> tar xzf mfix-17.1.tar.gz
> cd mfix-17.1
```

- Install the conda build dependencies listed in `build-aux/meta.yaml`.

- Install the [build dependencies](Install build dependencies) for your platform.

- Install `conda-build`:
```shell
> conda install conda-build
```


You can now build the conda package with:

```shell
> conda build build-aux/meta.yaml
```

If you are a Python developer familiar with distutils/setuptools, you can also run `setup.py` directly:

```shell
> python setup.py --help-commands
> python setup.py install
```

This can be used to install MFiX in non-Anaconda Python environments, if you know what you are doing.

## Running non-interactive Fortran-only solver

You can still build and run a solver from source without Python from the command line (without installing the GUI).

- Get a source distribution of MFIX

<!-- TODO: replace with link to website for actual release -->
Building from source requires the MFIX source tarball: [mfix-17.1.tar.gz](https://mfix.netl.doe.gov/gitlab/develop/mfix/builds/artifacts/develop/download?job=build%3Asrc). Extract the tarball and go into the top level source directory:

```shell
> tar xzf mfix-17.1.tar.gz
> cd mfix-17.1
```

- Install a Fortran compiler for your platform

- Configure and build MFIX solver (as in MFIX 2016-1):

```shell
> cd tutorials/FluidBed_DES
> ../../configure_mfix
> make
```

Note that the command-line solver executable is now called `mfixsolver`, not
`mfix` as it was in MFiX 2016.1 and earlier. Also note that `mfix.dat` is now `<RUN_NAME>.mfx`.

```shell
> ./mfixsolver -f DES_FB1.mfx
```
