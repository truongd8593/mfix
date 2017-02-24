# Tutorials

## Two-Fluid Models (MFiX-TFM)

### 2D Fluid Bed

This tutorial shows how to create a two dimensional fluidized bed simulation
using the two fluid model. The model setup is:

| Property       | Value                  |
|----------------|------------------------|
| geometry       | 10 cm x 30 cm          |
| mesh           | 20 x 60                |
| solid diameter | 200 microns (200e-6 m) |
| solid density  | 2500 kg/m2             |
| gas velocity   | 1 m/s                  |
| temperature    | 298 K                  |
| pressure       | 101325 Pa              |


#### Step 1. Create a new project
- On the file menu click on the ![new](../mfixgui/icons/newfolder.png) button
- Create a new project by double-clicking on "Blank" template.
- Enter a project name and browse to a location for the new project.

> Note: A new project directory will be created in the location directory, with
> the name being the project name.

<img alt="create project" src="media/gui_new_project.png" style="width:800;height:600" />

#### Step 2. Select model parameters
- On the `Model` pane, enter a descriptive text in the `Description` field
- Select "Two-Fluid Model (MFiX-TFM)" in the `Solver` combo-box.

<img alt="create project" src="media/gui_tfm_2d_model.png" style="width:800;height:600" />

#### Step 3. Enter the geometry
- On the `Geometry` pane select the `2 Dimensional` checkbox
- Enter "10/100" meters for the maximum x value
- Enter "30/100" meters for the maximum y value

<img alt="create project" src="media/gui_tfm_2d_geometry.png" style="width:800;height:600" />

#### Step 4. Enter the mesh
- On the `Mesh` pane, `Background` sub-pane
  - Enter "20" for the x cell value
  - Enter "60" for the y cell value

<img alt="create project" src="media/gui_tfm_2d_mesh.png" style="width:800;height:600" />

#### Step 5. Create regions for initial and boundary condition specification
- click the ![new](../mfixgui/icons/add.png) button to create a new region to be used for the bed initial condition.
  - Enter a name for the region in the `Name` field
  - Change the color by pressing the `Color` button
  - Enter "xmin" or "min" in the `From X` field
  - Enter "xmax" or "max" in the `To X` field
  - Enter "ymin" or "min" in the `From Y` field
  - Enter "ymax/2" or "max" in the `To Y` feild
  - Enter "zmin" or "min" in the `From Z` field
  - Enter "zmax" or "max" in the `To Z` feild

<img alt="create project" src="media/gui_tfm_2d_region1.png" style="width:800;height:600" />

- Click the ![new](../mfixgui/icons/bottom_region.png) button to create a new region with the `From` and `To` fields already filled out for a region at the bottom of the domain, to be used by the gas inlet boundary condition. `From Y` should equal `To Y`, defining an XZ-plane.
  - Enter a name for the region in the `Name` field

<img alt="create project" src="media/gui_tfm_2d_region2.png" style="width:800;height:600" />

- Click the ![new](../mfixgui/icons/top_region.png) button to create a new region with the `From` and `To` fields already filled out for a region at the top of the domain, to be used by the pressure outlet boundary condition. `From Y` should equal `To Y`, defining an XZ-plane.
  - Enter a name for the region in the `Name` field

<img alt="create project" src="media/gui_tfm_2d_region3.png" style="width:800;height:600" />

#### Step 6. Create a solid

- Click the ![new](../mfixgui/icons/add.png) button to create a new solid
- Enter a descriptive name in the `Name` field
- Enter the particle diameter of "200e-6" in the `Diameter` field
- ENter the particle density of "2500" in the `Density` field

<img alt="create project" src="media/gui_tfm_2d_solids.png" style="width:800;height:600" />

## Multiphase Particle in Cell Models (MFiX-PIC)

## Discrete Element Models (MFiX-DEM)

## Eulerian-Lagrangian-Eulerian (MFiX-Hybrid)
