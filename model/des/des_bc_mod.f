!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_INLET                                              !
!                                                                      !
!  Purpose: Common elements needed for the des mass inflow boundary    !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_BC

      USE param

! Logicals
      LOGICAL DES_MIO  ! either inlet or outlet exists
      LOGICAL DES_MI   ! inlet exists
      LOGICAL DES_MI_X, DES_MI_Y, DES_MI_Z  ! inlet exists on IJK face
      LOGICAL DES_MO_X, DES_MO_Y, DES_MO_Z  ! outlet exists on IJK face

! Total DES mass inlets and outlets
      INTEGER DES_BCMI   ! inlet count
      INTEGER DES_BCMO   ! outlet count

! Physical injection location
      DOUBLE PRECISION DES_BC_X_w(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_X_e(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Y_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Y_n(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Z_b(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Z_t(DIMENSION_BC)

! Logical that indicates whether a boundary has been defined by the user
      LOGICAL DES_BC_DEFINED(DIMENSION_BC)

! Specification for inflow into the system 
      DOUBLE PRECISION DES_BC_VOLFLOW_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION DES_BC_MASSFLOW_s(DIMENSION_BC, DIM_M)
      CHARACTER*16 DES_BC_TYPE(DIMENSION_BC)

! DES mass inlet velocity at the BC plane
      DOUBLE PRECISION DES_BC_U_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_V_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_W_s(DIMENSION_BC)

! Macroscopic density of solids phases in a specified boundary region
      DOUBLE PRECISION DES_BC_ROP_s (DIMENSION_BC, DIM_M)

! DES specification for solids phase velocity for WALL boundary
! conditions. The current setup is fairly limited. The specified
! boundary velocities are assigned to the indicated wall where a wall
! corresponds to one of the six planes in a cubic domain. Each wall
! corresponds to a number as follows west=1, east=2, bottom=3, top=4,
! south=5, north=6. See cfwallposvel for details. To specify a y or z
! velocity to the west wall set des_bc_vw_s(1,M) or des_bc_ww_s(1,M), 
! respectively (note an x velocity is not valid for a west or east wall).
! Since these are user input, they are allocated here with a constant
! preset size, but their actual size is represented by (nwalls, mmax)
      DOUBLE PRECISION DES_BC_Uw_s(DIMENSION_BC, DIM_M) 
      DOUBLE PRECISION DES_BC_Vw_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION DES_BC_Ww_s(DIMENSION_BC, DIM_M)


! Limit on the total number of divisions (fineness) used to represent
! the particle number distribution at an inlet.
      INTEGER, PARAMETER :: NUMFRAC_LIMIT = 10000

! Indicates if the boundary condition is mono or ploydisperse
! .F. = monodisperse and .T. = polydisperse
      LOGICAL, DIMENSION(:), ALLOCATABLE :: DES_BC_POLY    !(DES_BCMI)

! This array contains integers representing the mass/solid phase indices
! present at a specific boundary condtion in proportion to their
! respective number fraction at the inlet (i.e., it represents the
! particle number distribution of incoming solids at the inlet).  The 
! array is scaled in size according to the parameter NUMFRAC_LIMIT.  
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: DES_BC_POLY_LAYOUT   !(DES_BCMI,NUMFRAC_LIMIT)

! DES boundary condition ID number (boundary number used in mfix.dat
! when the boundary is originally defined)
      INTEGER, DIMENSION(:), ALLOCATABLE :: DES_BC_MI_ID   !(DES_BCMI)
      INTEGER, DIMENSION(:), ALLOCATABLE :: DES_BC_MO_ID   !(DES_BCMO)

! Boundary classification
! Assigned a value of the face/edge associated with inlet in 
! des_init_bc (see des_mi_classify)
!     Possible values in 2D: 'YN', 'YS', 'XE', 'XW'
!     Possible values in 3D: 'XZs', 'XZn', 'XYb', 'XYt', 'YZw', 'YZe'
      CHARACTER*3, DIMENSION(:), ALLOCATABLE :: DES_MI_CLASS   !(DES_BCMI)
! Assigned a value of the face/edge associated with outlet in 
! des_init_bc (see des_mo_classify)
!     Possible values in 2D: 'XW', 'XE', 'YN', 'YS'
!     Possible values in 3D: 'XW', 'XE', 'YN', 'YS', 'ZB', 'ZT'
      CHARACTER*2, DIMENSION(:), ALLOCATABLE :: DES_MO_CLASS   !(DES_BCMO)
! Assigned a value in des_init_bc (see des_mi_layout):
!     When the inlet bc velocity is sufficiently high it is assigned
!     'RAND', otherwise it is assigned 'ORDR'      
      CHARACTER*4, DIMENSION(:), ALLOCATABLE :: PARTICLE_PLCMNT
! Logical that can be flagged in the mfix.dat file to force the inlet
! to operate with an ordered boundary condition.  This may be useful
! during long simulations or if the inlet appears to be taking a long
! time to randomly place particles.
      LOGICAL FORCE_ORD_BC

! Particle injection factor; how many solid time steps (dtsolid) pass 
! before the next injection of a particle. if pi_count is greater than
! 1, then pi_factor is set to 1 (i.e. multiple particles enter every
! solids time step). 
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_FACTOR   !(DES_BCMI)

! Particle injection count (injection number); how many particles are
! injected in one solids time step. pi_count is set to one if
! less than 1 particle enters per solids time step.     
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_COUNT   !(DES_BCMI)


! Particle injection time scale; used when pi_factor > 1 to keep track
! of time needed for next injection     
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_MI_TIME
      

! Order inlet condition variables: used when the variable 
! particle_plcmnt is assigned 'ORDR'; in this case the inlet boundary is
! divided into grids and the particles are initially injected into an 
! available grid cell
! ----------------------------------------      
! keeps track of the number of grid cells that have been injected
! with a particle; is reset to 1 when all cells have been used.
      INTEGER, DIMENSION(:), ALLOCATABLE :: MI_FACTOR   !(DES_BCMI)
! the inlet edge/face is gridded and this variable stores the length
! dimension of the grid (no smaller than a particle diameter)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MI_WINDOW   !(DES_BCMI)
! Distance to offset incoming particles in ghost cells
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_BC_OFFSET   !(DES_BCMI)      

! the dimension of this variable is equal to the number of grid
! cells in the inlet edge/face
      TYPE dmi
         INTEGER, DIMENSION(:), POINTER :: VALUE   
      END TYPE dmi
! an array that stores integer values from 0 to a calculated factor in 
! sequential order.  
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: I_OF_MI
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: J_OF_MI   !(DES_BCMI)
! Construct an array of integers in values from 1 to a calculated factor
! in a random order, which is used when placing new particles.        
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: MI_ORDER   !(DES_BCMI)
! ----------------------------------------           

! Loop counter arrays: identifies the bounding i, j, k indices on the 
! Eulerian/Continuum mesh (defined by imax, jmax, kmax) that encompasses 
! an inlet boundary
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GS_ARRAY   !(DES_BCMI, 6)


     
      END MODULE DES_BC

