!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: CARTESIAN_GRID_INIT_NAMELIST                           C
!     Purpose: initialize the cartesian_grid-namelist                     C
!                                                                         C
!                                                                         C
!     Author: Jeff Dietiker                              Date: 26-Aug-08  C
!     Reviewer:                                          Date:            C
!     Comments:                                                           C
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

      CARTESIAN_GRID = .FALSE.

      N_QUADRIC = 0

      USE_POLYGON = .FALSE.

      N_USR_DEF = 0

      quadric_form = 'NORMAL'

      lambda_x = ZERO
      lambda_y = ZERO
      lambda_z = ZERO
      dquadric = ZERO

      theta_x = ZERO
      theta_y = ZERO
      theta_z = ZERO

      Radius = ZERO

      Half_angle = ZERO

      n_x = ZERO
      n_y = ZERO
      n_z = ZERO

      t_x = ZERO
      t_y = ZERO
      t_z = ZERO

      clip_xmin = - LARGE_NUMBER
      clip_xmax =   LARGE_NUMBER
      clip_ymin = - LARGE_NUMBER
      clip_ymax =   LARGE_NUMBER
      clip_zmin = - LARGE_NUMBER
      clip_zmax =   LARGE_NUMBER

      FLUID_IN_CLIPPED_REGION = .TRUE.

      BC_ID_Q = UNDEFINED_I

      N_GROUP = 1

      GROUP_SIZE(1) = 1
      GROUP_SIZE(2:DIM_GROUP) = 0
      GROUP_Q = 0
      GROUP_Q(1,1) = 1
      GROUP_RELATION = 'OR'
      RELATION_WITH_PREVIOUS = 'OR'
    
      TOL_SNAP       = 0.00D0  ! 0% of original edge length
      TOL_DELH       = 0.00D0  ! 0% of original Diagonal
      TOL_SMALL_CELL = 0.01D0  ! 1% of original cell volume

      TOL_SMALL_AREA = 0.01D0  ! 1% of original face area
      ALPHA_MAX      = ONE

      TOL_F     = 1.0D-9
      TOL_POLY  = 1.0D-9

      ITERMAX_INT = 10000

      SET_CORNER_CELLS = .FALSE.

      FAC_DIM_MAX_CUT_CELL = 0.25

      WRITE_VTK_FILES = .FALSE.
      TIME_DEPENDENT_FILENAME = .TRUE.
      VTK_DT     = UNDEFINED
      VTK_VAR(1) = 1
      VTK_VAR(2) = 2
      VTK_VAR(3) = 3
      VTK_VAR(4) = 4
      VTK_VAR(5:20) = UNDEFINED_I

      FRAME = -1

      PG_OPTION = 0

      CG_SAFE_MODE = .FALSE.
      PRINT_WARNINGS = .FALSE.

      CG_UR_FAC = 1.0

      PRINT_PROGRESS_BAR = .FALSE.
      BAR_WIDTH = 50
      BAR_CHAR  = '=' 
      BAR_RESOLUTION = 5.0

      WRITE_DASHBOARD = .FALSE.
      F_DASHBOARD = 1
      

      RETURN
      END SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST
