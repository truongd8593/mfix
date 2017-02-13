# MFIX User Guide

This document explains how to run MFIX 17.1, using either the GUI or the command line.

This document assumes MFIX is already installed. For information on building or
installing MFIX, please see the setup guide: [INSTALL.md](INSTALL.md)

Everything in this document applies to each platform (Linux, Mac, Windows)
unless otherwise noted.

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

# Tutorial on Running MFIX with the GUI

The MFIX install should put an `mfixgui` binary in your PATH. To start the MFIX GUI, run:

```shell
> mfixgui
```

<img alt="command line" src="doc/tutorial_0.png" style="width:400;height:300" />

 - Create a new project by double-clicking on 3dFluidbed.
 - Enter a new filename "my_3d_fluidbed" and click Ok.
 - Select a directory. "my_3d_fluidbed/my_3d_fluidbed.mfx" will be created in that directory.

<img alt="create project" src="doc/tutorial_1.png" style="width:800;height:600" />

 - Click the Start button ![Start button](mfixgui/icons/play.png) to start the MFIX simulation.

<img alt="command line" src="doc/tutorial_2.png" style="width:800;height:600" />

 - Click Ok to use the default mfixsolver installed with MFIX.

<img alt="command line" src="doc/tutorial_3.png" style="width:800;height:600" />

 - The simulation runs with output in the Terminal Window.

<img alt="command line" src="doc/tutorial_4.png" style="width:800;height:600" />

# User Interface Reference

## Toolbar

|                                                                                                                              |
|------------------------------------------------------------------------------------------------------------------------------|
| ![File menu](mfixgui/icons/menu.png)               Shows the  file menu for creating, opening, and saving MFIX project files.|
| ![Save button](mfixgui/icons/save.png)             The save button saves the MFIX project file.                              |
| ![Start button](mfixgui/icons/play.png)            The Start button displays the Run dialog, or unpauses an MFIX simulation. |
| ![Pause button](mfixgui/icons/pause.png)           The Pause button pauses an MFIX simulation.                               |
| ![Stop button](mfixgui/icons/stop.png)             The Stop button stops an MFIX simulation.                                 |
| ![Clear button](mfixgui/icons/restore_delete.png)  The Clear button deletes MFIX simulation data.                            |
| ![Parameters button](mfixgui/icons/functions.png)  The Parameters menu adjusts MFIX parameters.                              |

If no MFIX job is currently running, the Start button shows the run dialog.
If an MFIX job is running, you can pause it with the pause button and unpause it with the start button.
You can stop a running MFIX job with the stop button. A stopped job leaves restart (\*.RES) files to resume the simulation.
Starting a job with \*.RES files present will resume the job where it stopped.
The Clear button will delete the \*.RES files, and the next time the job is run it starts from the beginning.

### Run dialog

The run dialog

 - Threads (number of threads used for OpenMP)
 - NODESI (number divisions in X direction)
 - NODESJ (number divisions in Y direction)
 - NODESK (number divisions in Z direction)

 $NODESI \times NODESJ \times NODESK = n$ where $n$ is the number of MPI processes running. If not using MPI, $NODESI=NODESJ=NODESK=1$.

 - Submit to Queue (submit to queueing system, e.g. Grid Engine, PBS, SLURM. Currently only supported on Joule)
 - Run local MFIX executable (Runs MFIX solver on same system running the GUI)

 To run the MFIX solver, select an executable. Usually the default
 `pymfixsolver` command in PATH should be sufficient. If running a case with
 UDFs, you need to first build a case-specific MFIX as described in
 the [setup guide](INSTALL.md#building-for-udfs). You may want to build your own
 solver for other reasons.

 Click "Run" in the Run dialog starts the simulation.

## File menu

### Project Info

Displays metadata about the current project file.

 - Project Version
 - MFIX Release version that created this project file
 - Author
 - Modified By
 - Last Modified
 - Created
 - Notes

### New

Create a new project file from a list of templates. The list of templates can be filtered by:

|                                           |                                 |
|-------------------------------------------|---------------------------------|
| ![single](mfixgui/icons/single.png)       | single phase                    |
| ![tfm](mfixgui/icons/tfm.png)             | two fluid model                 |
| ![pic](mfixgui/icons/pic.png)             | particle in cell                |
| ![dem](mfixgui/icons/dem.png)             | discrete element model          |
| ![hybrid](mfixgui/icons/hybrid.png)       | hybrid two fluid and DEM model  |
| ![geometry](mfixgui/icons/geometry.png)   | template has cartesian geometry |
| ![chemistry](mfixgui/icons/chemistry.png) | template has chemistry          |


### Open

Open a existing project. You can import mfix.dat files from previous releases of MFIX, but the GUI will save them as a new filename with a \*.mfx extension.

### Save

Saves the current project.

### Save As

Saves the current project as a new filename.

### Export Project
### Settings

 - Style
 - Enable animations
 - Animation Speed
 - Enable Developer Tools

### Help
### About

Displays the current MFIX version.

### Quit

Exits MFIX. Will as for confirmation if project is unsaved or if a job is running.

## Model panes

Each pane in the main window allow editing of different options for an MFIX simulation.

### Model Setup

The Model pane is used to specify overall options for the simulation. Depending
on what is selected, other panes may be enabled or disabled.

 - Solver (Single Phase, Two-Fluid Model, DEM, PIC, Hybrid)
 - Gravity
 - Drag Model

### Geometry
### Mesh

 - Background
 - Mesher

### Regions

The Regions pane is used to define spatial regions of the simulation space that
are referred to in the Initial Conditions and Boundary Conditions panes.

### Fluid

The fluid pane is used to defined the physical properties of each fluid phase.

### Solids

The solids phase is used to defined the physical properties of each solids phase.

#### TFM
#### DEM
#### PIC

### Initial Conditions

The initial conditions pane is used to define initial conditions for each Region
(defined on each Region pane) for each phase (defined on Fluid or Solids panes).

 - Volume fraction
 - Temperature
 - Pressure
 - Velocity

### Boundary Conditions

The boundary conditions pane is used to define boundary conditions for each Region
(defined on each Region pane) for each phase (defined on Fluid or Solids panes).

 - Volume fraction
 - Pressure
 - Velocity

### Point Sources
### Internal Surfaces
### Chemistry
### Numerics

 - Residuals
 - Discretization
 - Linear Solver
 - Preconditioner
 - Advanced

### Outputs

 - Basic
 - VTK
 - SPx
 - NetCDF

### Monitors
### Run

The Run pane is used to define parameters related to how long the simulation runs.

 - Stop time
 - Time step
 - Maximum time step
 - Minimum time step
 - Time step factor

### Post-processing

## Visualization window

The visualization window provides a 3D image of the simulation boundary
condition. The visualization window can also be used to graph live statistics
about the simulation as it is running.

 - ![overscan](mfixgui/icons/overscan.png)
 - ![xy](mfixgui/icons/xy.png)
 - ![xz](mfixgui/icons/xz.png)
 - ![yz](mfixgui/icons/yz.png)
 - ![perspective](mfixgui/icons/perspective.png)
 - ![c](mfixgui/icons/camera.png)
 - ![v](mfixgui/icons/visibility.png)

## Terminal window

The terminal window displays the output of the mfixsolver job that would be displayed when running the solver on the command line.

Error messages and warnings are colored in red.

Informational messages from the GUI unrelated to the solver are colored in blue.

## Status bar

 - Modeler
 - Workflow


# Running MFIX with the command line

The command line version of MFIX works the same as in previous MFIX releases.
The main difference is that it is now called `mfixsolver` to distinguish the
command from the `mfixgui`

You may still want to run the GUI if you do not have Python on your platform, or
if you want to use features not yet supported by the GUI.