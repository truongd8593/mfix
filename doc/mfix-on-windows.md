# Building, Installing, and Running MFIX for Windows

## Building MFIX package

The binary Python extension (mfixsolver.so on Linux, mfixsolver.pyd on Windows)
can be built using a virtual machine running a free version of Windows 10.

If you have access to a Windows build host, you can skip ahead to "Install Anaconda".

### Set up the Windows 10 VM instance in VirtualBox

   Install [VirtualBox](virtualbox.org)

   Download 64-bit Windows 10 VM image, and import it into VirtualBox:
   
   https://www.microsoft.com/en-us/evalcenter/evaluate-windows-10-enterprise

   It is convenient to enable 'bidirectional copy and paste'
   (Settings/General/Advanced)

   It is also very helpful to enable a shared folder, which requires
   installation of VirtualBox "Guest additions".  This way you can use
   your exisiting Unix tools (editor, git) to work with files on the
   Windows filesystem.

   These instructions assume that you have created a shared folder as E:
   which points to a directory containing a clean git checkout of 'mfix'
   (gui branch)

   On the VM host:
```shell
   $ mkdir /somewhere/VBOX_SHARE
   $ cd /somewhere/VBOX_SHARE
   $ git clone -b gui <working_dir>/mfix  # can clone local mfix repo
```

   Then setup a shared folder pointing to VBOX_SHARE, and you should see
   E:\mfix in the Windows file browser


###  Disable Windows auto-updates

   Windows updates (BITS = "Background Intelligent Transfer Service")
   will use up all available bandwidth if allowed to.  This makes it
   difficult to download other needed packages.  To disable this, open
   the "Edit group policy" control panel, navigate to Computer
   Configuration/Administrative Templates/Network/BITS, select "Limit the
   maximim network bandwidth" and set the limits to 0.

   As the help text states: "If you disable or do not configure this policy
   setting, BITS uses all available unused bandwidth."

   Microsoft Search Service also uses a lot of bandwidth. To disable that,
   run "System Configuration" from the Start menu, select "Services"
   and disable "Windows Search"

   If you're going to do any serious work it is also helpful to disable Cortana:

   http://www.pcworld.com/article/2949759/windows/killing-cortana-how-to-disable-windows-10s-info-hungry-digital-assistant.html


### Install Anaconda
   As of this writing the current version is Anaconda2 4.2.0,
   Make sure to install the 64-bit Python 2.7 version, and
   say "No" to "Anaconda cloud" (not needed)

   https://www.continuum.io/downloads


### Install Cygwin

   Cygwin is needed for MinGW, GNU make, and autotools in order to build MFIX.

   https://cygwin.com/setup-x86_64.exe

   You need the following packages, in addition to the minimal 'base'
   packages.  Versions used are listed here for reference, but you
   should use the latest available packages.

   Cygwin package                  |    (version as of this writing)
-----------------------------------|----------------------------------
   autoconf                        |    13.1
   automake                        |    9.1
   make                            |    4.2.1-1
   binutils                        |    2.25-4
   mingw64-x86_64-gcc-fortran      |    5.4.0-3

   Be sure to install 86_64 packages not i686 (which are 32-bit).

   Do *not* install Cygwin Python.  We are using, and linking against,
   Anaconda.



### Install Anaconda gcc compatibility lib (needed to build Python
   extensions, since we're not using MSVC):

   In a Cygwin terminal (bash):
```shell
   $ conda install libpython
```

   This will probably trigger some other updates, say "yes" to everything.


### Fix up compiler names for f2py

   'make' and 'autoconf' can deal with the compiler name being prefixed
   with a host triple, but f2py expects the compilers to be called
   'gcc' and 'gfortran'

```shell
   $ pushd /usr/bin
   $ cp x86_64-w64-mingw32-gcc.exe gcc.exe
   $ cp x86_64-w64-mingw32-gfortran.exe gfortran.exe
   $ popd
```

   (Since f2py executes in Anaconda Python, which is *not* part of cygwin, it
   cannot handle cygwin-style symlinks, so you must copy, not link these files)

### (Optional) Add --verbose when calling F2PY in Makefile.usr.am command (helpful for debugging errors)

The line in Makefile.usr.am should read:
```shell
     python -m numpy.f2py --verbose -c --f77exec=...
```

### Configure mfix

```shell
   $ ./configure_mfix --enable-python --host=x86_64-w64-mingw32
```

### Build mfix

There is no "-static-libquadmath" linker option in GCC yet, so we use the workaround mentioned here:
https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539

The LDFLAGS arguments (and linking with gcc instead of gfortran) statically
compile libgcc, libgfortran, and libquadmath so those DLLs do not need to be
redistributed with the mfixsolver.pyd extension.

```shell
   $ make LDFLAGS='-static-libgcc -Wl,-Bstatic -lgfortran -lquadmath -Wl,-Bdynamic -lm -shared' LD=gcc
```

### Build mfixgui package as an exe installer (tested on Windows 10 VM, works)
```shell
   $ python setup.py bdist_wininst
```

### Build mfixgui package as an MSI installer (tested on Windows 10 VM, did not work)
```shell
   $ python setup.py bdist_msi
```

### Build mfixgui as a Python wheel package (tested on Windows 10 VM, works)
```shell
   $ python setup.py bdist_wheel
```


## Installing MFIX

If installing on the same system that was used for building, the runtime
dependency requirements are already satisfied. Cygwin is not a runtime
dependency, so the commands are shown for the Window Command Prompt
(CMD.exe).

### Install dependency (Anaconda)

   As of this writing the current version is Anaconda2 4.2.0,
   Make sure to install the 64-bit Python 2.7 version, and
   say "No" to "Anaconda cloud" (not needed)

   https://www.continuum.io/downloads

   TODO: document dependencies when installing with miniconda

### Install

   To install the mfixgui package from the exe file, double click on it in Windows Explorer.
          mfixgui-0.1.0.win-amd64-py2.7.exe

   To install the mfixgui package from the msi file, double click on it in Windows Explorer. (NOTE: MSI installer did not work for me)
          mfixgui-0.1.0.win-amd64-py2.7.msi

   To install the mfixgui package from the wheel file:
```shell
   C:\> pip install mfixgui-0.1.0-cp27-cp27m-win_amd64.whl
```

### Test install

   Note that if you want to start an interactive Anaconda python on
   Cygwin, you need to specify '-i'

```shell
   C:\> python -c 'import mfixsolver; print(mfixsolver.__file__)'
```

   If this fails, and you are have Cygwin installed, it is helpful to use Cygwin 'strace':

```shell
   $ strace python -c 'import mfixsolver'
```
   will show all file accesses and helps resolve dependencies

### Uninstall
   If you want to uninstall the mfixgui package:
```shell
   C:\> pip uninstall mfixgui
```

Or, from the Start menu go to "Add or remove program" and to uninstall mfixgui.

## Running MFIX

### Run and test
   To run the GUI:
```shell
   C:\> mfixgui
```

   Once the GUI is running, you can create a case and run it from the GUI.

   To run pymfix directly from the command line:
```shell
   C:\> pymfix -f mfix.dat -s
```
