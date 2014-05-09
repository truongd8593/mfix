!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: CARTESIAN_GRID_INIT_NAMELIST                           C
!     Purpose: initialize the cartesian_grid-namelist                     C
!                                                                         C
!                                                                         C
!     Author: Jeff Dietiker                              Date: 26-Aug-08  C
!     Reviewer:                                          Date:            C
!     Comments:                                                           C
!                                                                         C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST

      USE param1
      USE quadric
      USE cutcell
      USE polygon
      USE vtk
      USE progress_bar
      USE dashboard
      Use stl

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!-----------------------------------------------
!
!
      INCLUDE 'cartesian_grid/cartesian_grid_namelist.inc'

!<keyword category="cartesian grid" required="false">
!  <description>activate cartesian grid cut cell technique.</description>
!  <dependent keyword="COORDINATES" value="CARTESIAN"/>
!  <conflict keyword="COORDINATES" value="CYLINDRICAL"/>
!  <valid value=".false." note="do not use cartesian grid cut cell technique."/>
!  <valid value=".true." note="use cartesian grid cut cell technique. one of the following methods must be used to define the geometry:"/>
      CARTESIAN_GRID = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of quadric surfaces defining the boundaries (<=100).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <range min="0" max="100" />
      N_QUADRIC = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>use stl file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <valid value=".false." note="do not use stl file."/>
!  <valid value=".true." note="read triangulated geometry (for 3d geometry only) from geometry.stl."/>
      USE_STL = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>use .msh file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <valid value=".false." note="do not use .msh file."/>
!  <valid value=".true." note="read geometry (for 3d geometry only) from geometry.msh."/>
      USE_MSH = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>use polygons to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <valid value=".false." note="do not use polygons."/>
!  <valid value=".true." note="read polygon data (for 2d geometry only) from poly.dat."/>
      USE_POLYGON = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of user-defined functions (currently limited to 0 or 1). if set to 1, the geometry is defined in the user subroutine eval_usr_fct.f.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <valid value="0" note="Do not use user-defined function" alias=""/>                                     
!  <valid value="1" note="Use one user-defined function" alias=""/>                                     
!  <range min="0" max="1" />                                                
      N_USR_DEF = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>form of the quadric surface equation.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <valid value="normal" note="use normal form, as defined in equation (1). The lamdba's and d must be defined"/>
!  <valid value="plane" note="Plane. Needs to define n_x,n_y,n_z (unit normal vector pointing away from fluid cells)."/>
!  <valid value="x_cyl_int" note="Cylinder aligned with x-axis, internal flow. Needs to define radius(QID)."/>
!  <valid value="x_cyl_ext" note="Cylinder aligned with x-axis, external flow. Needs to define radius(QID)."/>
!  <valid value="y_cyl_int" note="Cylinder aligned with y-axis, internal flow. Needs to define radius(QID)."/>
!  <valid value="y_cyl_ext" note="Cylinder aligned with y-axis, external flow. Needs to define radius(QID)."/>
!  <valid value="z_cyl_int" note="Cylinder aligned with z-axis, internal flow. Needs to define radius(QID)."/>
!  <valid value="z_cyl_ext" note="Cylinder aligned with z-axis, external flow. Needs to define radius(QID)."/>
!  <valid value="x_cone" note="Cone aligned with x-axis, internal flow. Needs to define half_angle(QID)."/>
!  <valid value="y_cone" note="Cone aligned with y-axis, internal flow. Needs to define half_angle(QID)."/>
!  <valid value="z_cone" note="Cone aligned with z-axis, internal flow. Needs to define half_angle(QID)."/>
!  <valid value="sphere_int" note="Sphere, internal flow. Needs to define radius(QID)."/>
!  <valid value="sphere_ext" note="Sphere, external flow. Needs to define radius(QID)."/>
!  <valid value="C2C" note="Cylinder-to-cylinder conical junction, internal flow. Needs to be defined between two cylinders."/>
      quadric_form = 'NORMAL'
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>scaling factor, applied to all quadric geometry parameters. Must be a positive number</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>                                
!  <range min="0.0" max="" />  
      quadric_scale = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>coefficient lambda_x in equation (1) ('normal' form) or x-component of normal vector defining plane in equation (5) ('degenerate' form).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      lambda_x = ZERO
!</keyword>


!<keyword category="category name" required="true/false"                  
!                                    legacy="true/false">                 
!  <description></description>                                            
!  <arg index="" id="" min="" max=""/>                                    
!  <dependent keyword="" value="DEFINED"/>                                
!  <conflict keyword="" value="DEFINED"/>                                 
!  <valid value="" note="" alias=""/>                                     
!  <range min="" max="" />                                                
!  MFIX_KEYWORD=INIT_VALUE                                                
!</keyword> 


!<keyword category="cartesian grid" required="false">
!  <description>coefficient lambda_y in equation (1) ('normal' form) or y-component of normal vector defining plane in equation (5) ('degenerate' form).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      lambda_y = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>coefficient lambda_z in equation (1) ('normal' form) or z-component of normal vector defining plane in equation (5) ('degenerate' form).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>       
      lambda_z = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>coefficient d in equation (1).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      dquadric = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>rotation angle with respect to x-axis (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      theta_x = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>rotation angle with respect to y-axis (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      theta_y = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>rotation angle with respect to z-axis (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      theta_z = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>cylinder radius (used when quadric_form = *_cyl_***)</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      Radius = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>cone half angle, expressed in degrees (used when quadric_form = *_cone)</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      Half_angle = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>x-component of normal vector defining the plane (used when quadric_form = plane)</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      n_x = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>y-component of normal vector defining the plane (used when quadric_form = plane)</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      n_y = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>z-component of normal vector defining the plane (used when quadric_form = plane)</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      n_z = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in x-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      t_x = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in y-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      t_y = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in z-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      t_z = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>lower x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_xmin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>upper x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_xmax =   LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>lower y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_ymin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>upper y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_ymax =   LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>lower z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_zmin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>upper z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      clip_zmax =   LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>lower x-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_xmin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>upper z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_xmax =   LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Lower y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_ymin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Upper y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_ymax =   LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Lower z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_zmin = - LARGE_NUMBER
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Upper z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      piece_zmax =   LARGE_NUMBER
!</keyword>


!<keyword category="cartesian grid" required="false">
!  <description>flag defining the type of cells that are outside of the zone defined by [clip_xmin;clip_xmax], [clip_ymin;clip_ymax],[clip_zmin;clip_zmax].</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
!  <valid value=".false." note="remove cells from computational domain."/>
!  <valid value=".true." note="treat cells as fluid cells."/>
      FLUID_IN_CLIPPED_REGION = .TRUE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>boundary condition flag</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>                                    
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/> 
      BC_ID_Q = UNDEFINED_I
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of group(s) of quadrics (<=50).</description>
      N_GROUP = 1
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Number of quadrics in the group.</description>
      GROUP_SIZE(1) = 1
!</keyword>
      GROUP_SIZE(2:DIM_GROUP) = 0

!<keyword category="cartesian grid" required="false">
!  <description>NA</description>
      GROUP_Q = 0
!</keyword>

      GROUP_Q(1,1) = 1
!<keyword category="cartesian grid" required="false">
!  <description>relation among quadrics of a same group.</description>
!  <valid value="or" note="a point belongs to the computational domain if at least one of f(x,y,z) among all quadrics is negative"/>
!  <valid value="and" note="a point belongs to the computational domain if all of f(x,y,z) among all quadrics are negative"/>
      GROUP_RELATION = 'OR'
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>relation between current group and combination of all previous groups.</description>
!  <valid value="or" note="a point belongs to the computational domain if f-value for the current group or f-value for the combination of previous groups is negative"/>
!  <valid value="and" note="a point belongs to the computational domain if f-value for the current group and f-value for the combination of previous groups is negative"/>
      RELATION_WITH_PREVIOUS = 'OR'
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Tolerance used to snap an intersection point onto an existing cell corner (expressed as a fraction of edge length, between 0.0 and 0.5). For stretched grids, three values can be entered in the x, y and z directions.</description>
!  <range min="0.0" max="0.5" />   
      TOL_SNAP(1)    = 0.00D0  ! 0% of original edge length
!</keyword>
      TOL_SNAP(2)    = UNDEFINED
      TOL_SNAP(3)    = UNDEFINED

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to limit acceptable values of normal distance to the wall (expressed as a fraction of cell diagonal, between 0.0 and 1.0).</description>
      TOL_DELH       = 0.00D0  ! 0% of original Diagonal
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to detect small cells (expressed as a fraction of cell volume, between 0.0 and 1.0).</description>
      TOL_SMALL_CELL = 0.01D0  ! 1% of original cell volume
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to remove duplicate nodes (expressed as a fraction of cell diagonal, between 0.0 and 1.0).</description>
      TOL_MERGE      = 1.0D-12 ! fraction of original cell diagonal
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to detect small faces (expressed as a fraction of original face area, between 0.0 and 1.0).</description>
      TOL_SMALL_AREA = 0.01D0  ! 1% of original face area
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>maximum acceptable value of interpolation correction factor.</description>
      ALPHA_MAX      = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to find intersection of quadric surfaces or user-defined function with background grid.</description>
      TOL_F     = 1.0D-9
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to find intersection of polygon with background grid.</description>
      TOL_POLY  = 1.0D-9
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>maximum number of iterations used to find intersection points.</description>
      ITERMAX_INT = 10000
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to find intersection of stl triangles with background grid.</description>
      TOL_STL = 1.0D-6        ! Settings for STL file
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>smallest angle accepted for valid stl triangles (in degrees). triangles having an angle smaller that this value will be ignored.</description>
      STL_SMALL_ANGLE = 5.0   ! Degrees
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Dot product tolerance when determining if a point lies in a facet.</description>
      TOL_STL_DP = 1.0D-3        ! Settings for STL file
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Maximum number of STL facets per cell.</description>
      DIM_FACETS_PER_CELL  = 10        !
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Maximum number of STL facets per cell for des data arrays.</description>
      MAX_FACETS_PER_CELL_DES  = 24
!</keyword>


!<keyword category="cartesian grid" required="false">
!  <description>defines value of f_stl outside of the stl geometry. a value of 1.0 means the domain outside of the stl geometry is excluded from computation, i.e., an internal flow is computed.</description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
      OUT_STL_VALUE = 1.0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>boundary condition flag for the stl geometry</description>
      STL_BC_ID = UNDEFINED_I
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in x-direction, applied to the stl geometry.</description>
      TX_STL = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in y-direction, applied to the stl geometry.</description>
      TY_STL = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in z-direction, applied to the stl geometry.</description>
      TZ_STL = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>scaling factor, applied to the stl geometry. note that translation occurs after scaling.</description>
      SCALE_STL = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>tolerance used to find intersection of .msh file with background grid.</description>
      TOL_MSH = 1.0D-6        ! Settings for MSH file
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>defines value of f outside of the .msh geometry. a value of 1.0 means the domain outside of the .msh geometry is excluded from computation, i.e., an internal flow is computed.</description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
      OUT_MSH_VALUE = 1.0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in x-direction, applied to the .msh geometry.</description>
      TX_MSH = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in y-direction, applied to the .msh geometry.</description>
      TY_MSH = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>translation in z-direction, applied to the .msh geometry.</description>
      TZ_MSH = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>scaling factor, applied to the .msh geometry. note that translation occurs after scaling.</description>
      SCALE_MSH = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>ray propagation order used to determine whether any point is located inside or outside of the stl surface. a value of ijk means the propagation occurs in the i, followed by j, and k directions. other available orders are jki and kij</description>
!  <valid value="ijk" note="propagation occurs in the i, followed by j, and k directions"/>
!  <valid value="jki" note="propagation occurs in the j, followed by k, and i directions"/>
!  <valid value="kij" note="propagation occurs in the k, followed by i, and j directions"/>
      CAD_PROPAGATE_ORDER = '   '
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>NA</description>
      RAY_DIR = 'X-'
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>flag to detect and treat corner cells the same way as in the original mfix version (i.e. without cut cells). if set to .true., some cut cells may be treated as corner cells.</description>
!  <valid value=".true." note="some cut cells may be treated as corner cells."/>
!  <valid value=".false." note=""/>
      SET_CORNER_CELLS = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>factor used to allocate some cut cell arrays (expressed as a fraction of dimension_3g)</description>
      FAC_DIM_MAX_CUT_CELL = 0.25
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>write vtk files at regular intervals.</description>
!  <valid value=".false." note="do not write vtk files. if there are cut cells, they will not be displayed from the usual .res file"/>
!  <valid value=".true." note="valid only if cartesian_grid = .true."/>
!  <dependent keyword="cartesian_grid" value=".true."/>
      WRITE_VTK_FILES = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>use time-dependent vtk file names</description>
!  <valid value=".false." note="the vtk file overwrites the previous file (recommended for steady-state computation)."/>
!  <valid value=".true." note="a sequential integer is appended to the vtk filenames as they are written to create a series of files (recommended for transient computation)."/>
      TIME_DEPENDENT_FILENAME = .TRUE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>interval (expressed in seconds of simulation time) at which vtk files are written.</description>
      VTK_DT     = UNDEFINED
!</keyword>

      VTK_VAR(1) = 1
      VTK_VAR(2) = 2
!      VTK_VAR(3) = 3
!      VTK_VAR(4) = 4
      VTK_VAR(3:20) = UNDEFINED_I

!<keyword category="cartesian grid" required="false">
!  <description>NA</description>
      FRAME = -1
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>directory where vtk files are stored (default is run directory)</description>
      VTU_DIR = '.'
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>option for pressure gradient computation in cut cells.</description>
!  <valid value="1" note="use maximum of (east/west), (north/south), and (top/bottom) pairs of velocity cells."/>
!  <valid value="2" note="use both (east/west), (north/south), and (top/bottom) areas of velocity cells."/>
!  <valid value="0" note="use east, north and top areas of pressure cell (same as standard cells)."/>
      PG_OPTION = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>run code in safe mode.</description>
!  <valid value="1" note="performs initial preprocessing but use all original mfix subroutines during flow solution (using only cell volumes and areas of cut cells)."/>
!  <valid value="0" note="runs the code with modified subroutines for cut cell treatment."/>
      CG_SAFE_MODE = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>prints any warning message encountered during pre-processing on the screen.</description>
      PRINT_WARNINGS = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>under-relaxation factor used in cut cells (only cg_ur_fac(2) is used).</description>
      CG_UR_FAC = 1.0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>print a progress bar during each major step of pre-processing stage.</description>
      PRINT_PROGRESS_BAR = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>width of the progress bar (complete status), expressed in number of characters (between 10 and 80).</description>
!  <range min="10" max="80" />    
      BAR_WIDTH = 50
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>character used to create the progress bar.</description>
      BAR_CHAR  = '='
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>update frequency of progress bar, expressed in percent of total length (between 1.0 and 100.0).</description>
!  <range min="1.0" max="100.0" />  
      BAR_RESOLUTION = 5.0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>writes the file dashboard.txt at regular intervals. the file shows a summary of the simulation progress.</description>
      WRITE_DASHBOARD = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>frequency, expressed in terms of iterations, at which the dashboard is updated.</description>
      F_DASHBOARD = 1
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>location of control points in x-direction.</description>
      CPX      = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of cells within a segment (x-direction).</description>
      NCX      = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>expansion ratio (last dx/first dx) in a segment (x-direction).</description>
      ERX      = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of first dx in a segment (x-direction). a negative value will copy dx from previous segment (if available).</description>
      FIRST_DX = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of last dx in a segment (x-direction). a negative value will copy dx from next segment (if available).</description>
      LAST_DX  = ZERO
!</keyword>


!<keyword category="cartesian grid" required="false">
!  <description>location of control points in y-direction.</description>
      CPY      = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of cells within a segment (y-direction).</description>
      NCY      = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>expansion ratio (last dy/first dy) in a segment (y-direction).</description>
      ERY      = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of first dy in a segment (y-direction). a negative value will copy dy from previous segment (if available).</description>
      FIRST_DY = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of last dy in a segment (y-direction). a negative value will copy dy from next segment (if available).</description>
      LAST_DY  = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>location of control points in z-direction.</description>
      CPZ      = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>number of cells within a segment (z-direction).</description>
      NCZ      = 0
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>expansion ratio (last dz/first dz) in a segment (z-direction).</description>
      ERZ      = ONE
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of first dz in a segment (z-direction). a negative value will copy dz from previous segment (if available).</description>
      FIRST_DZ = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>value of last dz in a segment (z-direction). a negative value will copy dz from next segment (if available).</description>
      LAST_DZ  = ZERO
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Turns on the re-indexing of cells. When true, inactive (dead) cells are removed from computational domain. </description>
      RE_INDEXING = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Attemps to adjust grid partition. Each processor will be assigned its own size to minimize load imbalance  </description>
      ADJUST_PROC_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Attemps to adjust grid partition. Each processor will be assigned its own size to minimize load imbalance  </description>
      REPORT_BEST_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description> Temporary setting used in serial run to report best domain size for parallel run </description>
      NODESI_REPORT = 1 
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description> Temporary setting used in serial run to report best domain size for parallel run </description>
      NODESJ_REPORT = 1 
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description> Temporary setting used in serial run to report best domain size for parallel run </description>
      NODESK_REPORT = 1 
!</keyword>

!<keyword category="cartesian grid" required="false">
!  <description>Attemps to minimize the size of the send/receive layers </description>
      MINIMIZE_SEND_RECV = .TRUE.
!</keyword>

      RETURN
      END SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST
