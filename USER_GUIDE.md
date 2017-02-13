# MFIX User Guide

This document explains how to run MFIX 17.1, using either the GUI or the command line.

This document assumes MFIX is already installed. For information on building or
installing MFIX, please see the setup guide: [INSTALL.html](INSTALL.html)

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

Potential users may find reviewing the Frequently Asked Questions section of the
MFIX website useful before downloading the code.

## About the GUI
The GUI is written in pure Python, leveraging the strengths of Python for quick
code development, extensive existing libraries, and flexibility. The GUI will
run on any operating system that Python can be installed in including Linux,
Windows, and Mac. The GUI uses Qt either through PyQt4, PyQt5 or PySide
wrappers as the GUI library. The
[Visualization Toolkit (VTK)](http://www.vtk.org/) is used to visualize and
manipulate the input geometry.


## Development state of MFIX models

MFIX provides a suite of models that treat the carrier phase (typically
the gas phase) and disperse phase (typically the solids phase)
differently. Their current state of development is summarized in the
tables below.

**MFIX-TFM (Two-Fluid Model)** is an Eulerian-Eulerian model, which
supports a broad range of capabilities for dense, reacting, multiphase
flows by representing the fluid and solids as interpenetrating continua.
This is the most mature MFIX model and is capable of modeling multiphase
reactors ranging in size from benchtop to industry-scale. Approximation
of the solid phase as a continuum typically allows for faster simulation
time than Lagrangian techniques; however, it also introduces the need
for accurate mathematical models to capture realistic solids phase
behavior. This includes transport properties, heterogeneous reaction
kinetics, and constitutive relations for interaction between fluid and
solid phases, e.g., solids phase drag and interphase heat transfer.

<div>

<div style="float:left">

|                      | Serial   | ^†^DMP   | ^‡^SMP   |
| -------------------- | -------- | -------- | -------- |
| Momentum Equations   | ●        | ●        | ●        |
| Energy Equations     | ●        | ●        | ●        |
| Species Equations    | ●        | ●        | ●        |
| Chemical Reactions   | ●        | ●        |          |
| Cartesian cut-cell   | ●        | ●        | **□**    |

</div>

<div>

[]{#TFM_pic .anchor}![](doc/media/devstate_tfm.png){width="3.022222222222222in"
height="1.8520833333333333in"}

</div>

</div>

**MFIX-DEM (Discrete Element Model)** is an Eulerian-Lagrangian model
that treats the fluid phase as a continuum and models the individual
particles of the solid phase. This is a relatively new variation on
MFIX. While the treatment of individual particles can provide higher
fidelity over a broad range of flow regimes (from dilute to packed), it
also very challenging when dealing with very large numbers of particles
for large-scale simulations. These large-scale applications will require
high performance computing (HPC) resources and large amounts of computer
time. Code optimization and speed up are critical research fronts to
support industrial scale applications.

<div>

<div style="float:left">

|                      | Serial   | ^†^DMP   | ^‡^SMP   |
| -------------------- | -------- | -------- | -------- |
| Momentum Equations   | ●        | ●        | ●        |
| Energy Equations     | ●        | ●        |          |
| Species Equations    | ●        | ●        |          |
| Chemical Reactions   | ●        | ●        |          |
| Cartesian cut-cell   | ○        | ○        |          |

</div>

<div>

[]{#DEM_pic .anchor}![](doc/media/devstate_dem.png){width="3.0125in"
height="1.5895833333333333in"}

</div>

</div>

**MFIX-PIC (Multiphase Particle in Cell)** is another
Eulerian-Lagrangian model that represents the fluid as a continuum while
using "parcels" to represent groups of real particles with similar
physical characteristics. The MFIX-PIC approach offers reduced
computational cost over MFIX-DEM as there are typically few parcels to
track and parcel collisions are not resolved. However, the added
modeling approximations influence the overall accuracy of the method.
Development, validation, and optimization of modeling approximations are
critical research fronts.

<div>

<div style="float:left">

|                      | Serial   | ^†^DMP   | ^‡^SMP   |
| -------------------- | -------- | -------- | -------- |
| Momentum Equations   | ●        |          | ○        |
| Energy Equations     |          |          |          |
| Species Equations    |          |          |          |
| Chemical Reactions   |          |          |          |
| Cartesian cut-cell   | ○        | □        |          |

</div>

<div>

[]{#MPPIC_pic
.anchor}![](doc/media/devstate_pic.png){width="3.0104166666666665in"
height="1.7055555555555555in"}

</div>

</div>

**MFIX-Hybrid (Eulerian-Lagrangian-Eulerian)** is a blend of MFIX-TFM
and MFIX-DEM that represents the fluid as a continuum and models solids
as either a continuous phase (TFM) or discrete particles (DEM). This
technique is presently restricted to solving only the momentum equations
to yield hydrodynamic predictions. This model is still in its infancy
and has seen only limited testing.

<div>

<div style="float:left">

|                      | Serial   | ^†^DMP   | ^‡^SMP|
| -------------------- | -------- | -------- | ------|
| Momentum Equations   | ○        | ○        | ○     |
| Energy Equations     |          |          |       |
| Species Equations    |          |          |       |
| Chemical Reactions   |          |          |       |
| Cartesian cut-cell   | ○        | ○        | ○     |

</div>

[]{#Hybrid_pic .anchor}![](doc/media/devstate_hybrid.png){width="3.017361111111111in" height="1.601388888888889in"}

<div>

</div>

</div>

● – implemented and fully tested

○ – implemented with limited testing

**□** – not tested or status unknown

† Models not extended to DMP-parallel are only available for serial
runs.

‡ Models not extended to SMP-parallel are available for SMP runs but do
not scale with thread count.

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

## Main Toolbar

| Icon                                              | Description                                                               |
|---------------------------------------------------|---------------------------------------------------------------------------|
| ![File menu](mfixgui/icons/menu.png)              | Shows the  file menu for creating, opening, and saving project files.     |
| ![Save button](mfixgui/icons/save.png)            | The save button saves the project file.                                   |
| ![Start button](mfixgui/icons/play.png)           | The Start button displays the Run dialog, or un-pauses a simulation.       |
| ![Pause button](mfixgui/icons/pause.png)          | The Pause button pauses a simulation.                                     |
| ![Stop button](mfixgui/icons/stop.png)            | The Stop button stops a simulation.                                       |
| ![Rest button](mfixgui/icons/restore_delete.png)  | The Reset button deletes output data.                                     |
| ![Parameters button](mfixgui/icons/functions.png) | The Parameters menu allows changing of parameters.                        |

If no MFIX job is currently running, the Start button shows the run dialog.
If an MFIX job is running, you can pause it with the pause button and un-pause it with the start button.
You can stop a running MFIX job with the stop button. A stopped job leaves restart (\*.RES) files to resume the simulation.
Starting a job with \*.RES files present will resume the job where it stopped.
The Clear button will delete the \*.RES files, and the next time the job is run it starts from the beginning.

### Run dialog
The run dialog allows for running the simulation from the GUI. The following options can be changed in this
dialog:

 - Threads (number of threads used for OpenMP)
 - NODESI (number divisions in X direction)
 - NODESJ (number divisions in Y direction)
 - NODESK (number divisions in Z direction)

 $NODESI \times NODESJ \times NODESK = n$ where $n$ is the number of MPI processes running. If not using MPI, $NODESI=NODESJ=NODESK=1$.

The GUI supports running both locally as well as submitting to a queue.

#### Run Local
To run locally, select the "Run local MFiX executable" tab. Select an executable from the dropdown list or
click the browse button to specify an executable that is not in the list. Usually the default
`pymfixsolver` command in PATH should be sufficient. If running a case with UDFs, you need to first build
a case-specific MFIX as described in the [setup guide](INSTALL.md#building-for-udfs). You may want to build
your own solver for other reasons, such as specifying various compiler flags to optimize the executable
for your specific hardware.

Click "Run" in the Run dialog to start the simulation.

#### Submit to Queue
To submit a job to a queue (submit to queueing system, e.g. Grid Engine, PBS, SLURM), select the "Submit to Queue"
tab. Select an executable from the dropdown list or click the browse button to specify an executable that is not
in the list. Next, select a template from the dropdown list or click the browse button to choose a template that is
not in the dropdown list. Based on the selected template, the widgets in the "queue options" section will update
automatically. Once the options are specified, click "submit" in the run dialog to submit the job to the queue.

Custom queue scripts are supported. The format for this script is described in the
[Queue Templates](#queue-templates) section.

## File menu
The file menu allows for opening, creating, copying, and exporting projects as well
as viewing project meta data, changing settings, and viewing available help
documentation, including this document.

### Project Info
Displays metadata about the current project file.

| Label                     | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| Project Version           | a version that is increased by 1 everytime the project is saved |
| Created with MFiX Version | the version of MFiX that the project was created with           |
| Author                    | the user that created the project                               |
| Modified By               | list of users that edited the project                           |
| Last Modified             | the date that the project was last modified                     |
| Created                   | the date the project was created                                |
| Notes                     | area where the user can record notes about the project          |

### New
Create a new project file from a template file. The list of templates can be filtered
by selecting one or more of the following model types:

| Icon                                      | Description                            |
|-------------------------------------------|----------------------------------------|
| ![single](mfixgui/icons/single.png)       | Single phase Model                     |
| ![tfm](mfixgui/icons/tfm.png)             | Two Fuild Model (TFM)                  |
| ![pic](mfixgui/icons/pic.png)             | Particle in Cell Model (PIC)           |
| ![dem](mfixgui/icons/dem.png)             | Discrete Element Model (DEM)           |
| ![hybrid](mfixgui/icons/hybrid.png)       | Hybrid Model (TFM + DEM)               |
| ![geometry](mfixgui/icons/geometry.png)   | Cartesian cut-cell (complex) geometry  |
| ![chemistry](mfixgui/icons/chemistry.png) | Chemistry                              |


### Open

Open a existing project. You can import mfix.dat files from previous releases of
MFiX, but the GUI will save them as a new filename with a \*.mfx extension. The
GUI also performs a number of conversions, including converting old keywords to
new keywords and conversion from CGS units to SI units.

### Save

Saves the current project.

### Save As

Saves the current project to a new directory and/or as a new filename. The
project is then opened in that location and with that filename.

### Export Project

Export the current project to a new directory and/or as a new filename, bit keep
the original project opened.

### Settings

Change settings that affect how the GUI acts.

| Option | Description |
|--------|-------------|
| Style | Change the application style |
| Enable animations | enable or disable animated widgets |
| Animation Speed | set the speed at which animations occur |
| Enable Developer Tools | hide/show widgets that are mainly for developers of the GUI |

### Help

Shows available documentation (including this document) and tutorial videos
(if connected to the internet) to help users utilize the features of the GUI.

### About

Displays the current MFIX version as well as the version of the various
libraries that are currently being used.

### Quit

Exits MFIX. Will ask for confirmation if project is unsaved or if a job is
running.

## Model panes

Each pane in the main window allows editing of different options for an MFIX
simulation.

### Model Setup

The Model pane is used to specify overall options for the simulation. Depending
on what is selected, other panes may be enabled or disabled. Options that are
specified on this pane include:

 - Solver (Single Phase, Two-Fluid Model, DEM, PIC, Hybrid)
 - Option to disable the fluid phase
 - Option to enable thermal energy equations
 - Option to enable turbulence, if the fluid phase is enabled
 - Gravity in teh x, y, and z directions
 - Drag Model including parameters for the selected drag model

Other advanced options that can be selected include:

 - Momentum formulation (Model a, Model B, Jackson, or Ishii)
 - Subgrid model (only avaliable with TFM, Wen-Yu drag model, etc...)
 - Subgird filter size
 - Subgrid wall correction

### Geometry

The Geometry pane allows the specification of the model geometry. This includes
specifiying the domain extents (xmin, xmax, ymin, ymax, zmin, zmax) and 2D/3D
selection. If there is complex geometry, the "Autosize" button can automatically
set the extents to encompass the geometry.

The geometry section provides tools for adding, applying filters, using
automated wizards to create and copy geometry, remove, copy, and perform boolean
opartions on the geometry. All the geometry operations and visualizations are
perfromed using the [Visualization Toolkit (VTK)](http://www.vtk.org/)'s methods
and functions.

Geometry toolbar icons:

| Icon                                        | Description                                       |
|---------------------------------------------|---------------------------------------------------|
| ![geometry](mfixgui/icons/geometry.png)     | add geometry to model                             |
| ![filter](mfixgui/icons/filter.png)         | apply a filter to the currently selected geometry |
| ![wizard](mfixgui/icons/wand.png)           | use a wizard to create or copy geometry           |
| ![remove](mfixgui/icons/remove.png)         | remove the selected geometry                      |
| ![copy](mfixgui/icons/copy.png)             | copy the selected geometry                        |
| ![union](mfixgui/icons/union.png)           | perform a union of the selected geometry          |
| ![intersect](mfixgui/icons/intersect.png)   | perform a intersection of the selected geometry   |
| ![difference](mfixgui/icons/difference.png) | perform a difference of the selected geometry     |


#### Adding Geometry
Geometry can be added by pressing the add geometry icon and selecting the
geometry to add. There are two types of geometric objects supported in the GUI,
polydata (triangles, like stl files) and implicit functions (quadrics). Boolean
operations can not be perfromed between polydata and implicit geometry types.
The implicit function needs to be converted to polydata by using the `sample
implicit` filter. Converting the implicit function also needs to be done in
order for the GUI to export a STL file that the mifxsolver can use.

In the geomtry tree, the geometry object is described with an icon:

| Icon                                        | Geometry Type |
|---------------------------------------------|---------------|
| ![geometry](mfixgui/icons/geometry.png)     | polydata      |
| ![function](mfixgui/icons/function.png)     | implicit      |
| ![filter](mfixgui/icons/filter.png)         | filter        |
| ![union](mfixgui/icons/union.png)           | union         |
| ![intersect](mfixgui/icons/intersect.png)   | intersect     |
| ![difference](mfixgui/icons/difference.png) | difference    |


The following geometric objects can be added:
- STL file(s)
- Implicit (Quadric) functions
    - sphere
    - box
    - cylinder
    - cone
    - quadric
    - superquadric
- Primitives
    - sphere
    - box
    - cylinder
    - cone
- Parametrics (torus, boy, conic spiral, etc.)

#### Applying Filters
Filters can be applied to selected geometry by first, selecting the geometry the
filter will be applied to, next pressing the filter icon
and finally, selecting a filter from the filter menu. The filter options can be
edited in the parameter section. The following filters are included:

| Filter                | Description                                                          | vtk class                     |
|-----------------------|----------------------------------------------------------------------|-------------------------------|
| sample implicit       | converts an implicit function to polydata                            | vtkSampleFunction             |
| transform             | rotate, scale, translate polydata                                    | vtkTransformPolyDataFilter    |
| clean                 | merge duplicate points and remove unused points and degenerate cells | vtkCleanPolyData              |
| fill holes            | fill holes                                                           | vtkFillHolesFilter            |
| triangle              | make sure all polys are triangles                                    | vtkTriangleFilter             |
| decimate              | reduce the number of triangles                                       | vtkDecimatePro                |
| quadric decimation    | reduce the number of triangles                                       | vtkQuadricDecimation          |
| quadric clustering    | reduce the number of triangles                                       | vtkQuadricClustering          |
| linear subdivision    | subdivide based on a linear scheme                                   | vtkLinearSubdivisionFilter    |
| loop subdivision      | subdivide based on the Loop scheme                                   | vtkLoopSubdivisionFilter      |
| butterfly subdivision | subdivide based on 8-point butterfly scheme                          | vtkButterflySubdivisionFilter |
| smooth                | move points based on Laplacian smoothing                             | vtkSmoothPolyDataFilter       |
| windowed sinc         | move points based on a windowed sinc function interpolation kernel   | vtkWindowedSincPolyDataFilter |
| reverse sense         | reverse order and/or normals of triangles                            | vtkReverseSense               |


#### Wizards
Four wizards have been included to perform routine tasks when setting up
multiphase flow simulations. These wizards allow for creating cyclones,
reactors, and hoppers. The distributed wizard can also be used to distribute one
geometry inside another geometry with random, cubic, or body centered cubic
positions. Random rotations can also be applied with the wizard.

#### Boolean Operations
Boolean operations can be performed with geomtry objects of the same type
(implicit, polydata). Boolean operations can not be perfromed between polydata
and implicit geometry types. The implicit object needs to be first converted to
a polydata object using the sample implicit filter.

> Note: boolean operation between two polydata objects can crash the GUI due to
> segfaults in vtk.

### Mesh
The Mesh pane allows for the specification of the background mesh, add/remove
control points in the x, y, and z directions as well as change various
cut-cell tolerances.

#### Background Mesh
On the Background sub-pane, a uniform mesh can be
specified by entering the numner of cells in the x, y, and z directions. Control
points can be added by pressing the ![add](mfixgui/icons/add.png) button. Once a
control point has been added, the position, numner of cells, stretch, first and
last parameters can be changed. A control point can be split evenly by
`right-click` on the control point to be split and selecting split. This
operation will create a new control point at the midpoint between the previous
control point and the selected control point, dividing the cells evenly between
the two. Control points can be removed by pressing the
![remove](mfixgui/icons/remove.png) button.

The stretch parameter is a value that will apply a non-uniform grid spacing to
the cells. The value is defined as ${Last Width} \over {First Width}$. a value
larger than 1 stretches the grid as x increases, while a value smaller than one
compresses the grid as x increases. A value of 1 will keep the spacing uniform.

The first and last values allow the specification of the first and last widths,
the other cell widths adjust accordingly. If a negative value is specified for
first, the last width from the previous grid segment is copied. If a negative
value is specified, the first width from the next segment is copied.

#### Mesher
The Mesher sub-pane exposes options to adjust the cut-cell mesher. These options
include:

| Option                 | Description                                                                                                                                                                          |
|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| External flow          | select internal or external flow. Note: this depends on which way the normals are pointing on the stl file. If they are pointing out of the geometry, then the text will be correct. |
| Small cell tolerance   | tolerance to detect, and remove small cells                                                                                                                                          |
| Small area tolerance   | tolerance to detect, and remove cells with small faces                                                                                                                               |
| Merge tolerance        | tolerance used to merge duplicate nodes                                                                                                                                              |
| Snap tolerance         | tolerance to move an intersection point to an existing cell corner                                                                                                                   |
| Allocation Factor      | factor used in allocation of cut-cell arrays                                                                                                                                         |
| Maximum iterations     | maximum number of iterations used to find intersection points                                                                                                                        |
| Intersection tolerance | tolerance used to find intersection of background mesh and stl triangles                                                                                                             |
| Facet angle tolerance  | ignore stl facets that have an angle less than this tolerance                                                                                                                        |
| Dot product tolerance  | tolerance used to determine if a point lies in a facet                                                                                                                               |
| Max facets per cell    | maximum number of facets allowed in a cell                                                                                                                                           |


### Regions

The Regions pane is used to define spatial regions (points, lines, planes,
boxes, or stls) of the simulation space that are refered to later when
creating:
- [Initial Conditions](#initial-conditions)
- [Boundary Conditions](#boundary-conditions)
- [Point Sources](#point-sources)
- [Outputs](#outputs)

Tool bar: 

| Icon                                        | Description|
|---------------------------------------------|---------------|
| ![add](mfixgui/icons/add.png)     | create a new region     |
| ![remove](mfixgui/icons/remove.png) | delete the selected region |
| ![copy](mfixgui/icons/copy.png) | duplicate the selected region |
| ![all region](mfixgui/icons/all_region.png) | create a region the encompasses the entire domain |
| ![left region](mfixgui/icons/left_region.png) | create a region on the left side of the domain |
| ![right region](mfixgui/icons/right_region.png) | create a region on the right side of the domain |
| ![top region](mfixgui/icons/top_region.png) | create a region on the top side of the domain |
| ![bottom region](mfixgui/icons/bottom_region.png) | create a region on the bottom side of the domain |
| ![front region](mfixgui/icons/front_region.png) | create a region on the front side of the domain |
| ![back region](mfixgui/icons/back_region.png) | create a region on the back side of the domain |


•	Specify an alias for easy referencing (e.g., outlet, solids-bed).
•	Specify region extents (xmin, xmax, ymin, ymax, zmin, zmax )


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

# Appendix

## Queue Templates
Custom queue templates allow users to customize the functionality of the GUI for their specific system. Queue templates
included with the source can be found in the `MFIX_HOME\queue_templates` directory. Queue templates are files that contain two sections. The first section is the configuration section that tells the GUI what widgets to display as well as
various options for those widgets. The second section is the actual script that is executed.

Variables can be used throughout the template, including with the widgets, and are reference by `${my_variable}`. There are a couple of built in variables including:

| Variable     | Description                                   |
|--------------|-----------------------------------------------|
| SCRIPT       | the path of this script, in the run directory |
| PROJECT_NAME | the name of the project                       |
| JOB_ID       | the job id extracted from the submit command  |
| COMMAND      | the command that executes mfix                |

The configuration section starts with `## CONFIG` and ends with `## END CONFIG`. The format follows the `Microsoft Windows INI` format. Sections are defined with a `[section]` header, followed by a collection of `key=value` or `key: value` entries for defining parameters. For example:

```
## CONFIG
[My Section]
foodir: %(dir)s/whatever
dir=frob
long: this value continues
   in the next line
## END CONFIG
```

The configuration section has a special section called `[options]` where the following options can be specified:

| Key          | Description                                                                      |
|--------------|----------------------------------------------------------------------------------|
| name         | name of the template, this is displayed in the template drop-down box in the gui |
| job_id_regex | regular expression to extract the job id from the output of the submit command   |
| status_regex | regular expression to extract the job status from the status command             |
| submit       | the submission command                                                           |
| delete       | the job cancel or delete command                                                 |
| status       | the job status command                                                           |

An example of values that work with Grid Engine:

```
[options]
name: Joule
job_id_regex: (\d+)
status_regex: ([rqw])
submit: qsub ${SCRIPT}
delete: qdel ${JOB_ID}
status: qstat -j ${JOB_ID}
```
To define a new variable and widget to edit that variable in the GUI, create a new section:

```
[my_value]
```

The widget and options for that widget can then be selected by specifying various parameters including:

| Parameter | Description                                         | Values                                                         |
|-----------|-----------------------------------------------------|----------------------------------------------------------------|
| widget    | the widget to be used                               | `lineedit`, `combobox`, `checkbox`, `spinbox`, `doublespinbox` |
| label     | text to be placed beside the widget                 | `any string`                                                   |
| value     | default value                                       | a value such as `1`, `10.3`, `True`, `some text`               |
| items     | list of items for the combobox                      | items delimited by &#124; character                            |
| help      | text to be displayed in the tooltip for that widget | `this widget does this`                                        |
| true      | value to be returned if a checkbox is checked       | a value such as `1`, `10.3`, `True`, `some text`               |
| false     | value to be returned if a checkbox is un-checked    | a value such as `1`, `10.3`, `True`, `some text`               |

An example defining a combo box:

```
[my_email]
widget: combobox
label: email
value: you@mail.com
items: you@mail.com|
       me@mail.net|
       hi@mail.org
help: email to send notifications to
```

The value of this widget can now be referenced throughout the template with `${my_email}`

The rest of the configuration file, outside of the  `## CONFIG` to `## END CONFIG` section is the script that needs to
be executed to submit a job to your specific queue system. For example, with Grid Engine on Joule, the following script
specifies the job name, job type, cores, and queue as well as loads the required modules and finally runs mfix with
`${COMMAND}`.

```
## Change into the current working directory
#$ -cwd
##
## The name for the job. It will be displayed this way on qstat
#$ -N ${JOB_NAME}
##
## Number of cores to request
#$ -pe ${JOB_TYPE} ${CORES}
##
## Queue Name
#$ -q ${QUEUE}
##

##Load Modules
module load ${MODULES}
##Run the job
${COMMAND}
```
