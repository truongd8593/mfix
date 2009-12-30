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
      LOGICAL DES_MI
      LOGICAL DES_MO_X, DES_MO_Y, DES_MO_Z

! Physical injection location
      DOUBLE PRECISION DES_BC_X_w(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_X_e(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Y_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Y_n(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Z_b(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_Z_t(DIMENSION_BC)

! Specification for inflow into the system 
      DOUBLE PRECISION DES_BC_VOLFLOW_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_MASSFLOW_s(DIMENSION_BC)
      CHARACTER*16 DES_BC_TYPE(DIMENSION_BC)

! DES mass INLET solids-phase velocity at the BC plane
      DOUBLE PRECISION DES_BC_U_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_V_s(DIMENSION_BC)
      DOUBLE PRECISION DES_BC_W_s(DIMENSION_BC)

! DES specification for solids phase velocity for WALL boundary conditions    
! Currently limited setup as des_bc_uw_s(nwalls,1), etc...
      DOUBLE PRECISION DES_BC_Uw_s(DIMENSION_BC, DIM_M)   !(number of des boundaries, number of solids phases)
      DOUBLE PRECISION DES_BC_Vw_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION DES_BC_Ww_s(DIMENSION_BC, DIM_M)

! DES boundary condition ID number (boundary number used in mfix.dat
! when the boundary is originally defined)
      INTEGER, DIMENSION(:), ALLOCATABLE :: DES_BC_MI_ID   !(number of MI boundaries)
      INTEGER, DIMENSION(:), ALLOCATABLE :: DES_BC_MO_ID   !(number of MO boundaries)

! Boundary classification
! Assigned a value of the face/edge associated with inlet in check_des_bc
! (see des_mi_classify)
!     Possible values in 2D: 'YN', 'YS', 'XE', 'XW'
!     Possible values in 3D: 'XZs', 'XZn', 'XYb', 'XYt', 'YZw', 'YZe'
      CHARACTER*4, DIMENSION(:), ALLOCATABLE :: DES_MI_CLASS   !(number of MI boundaries)
! Assigned a value of the face/edge associated with outlet in check_des_bc
! (see des_mo_classify)
!     Possible values in 2D: 'Ysn', 'Xwe'
!     Possible values in 3D: 'Xwe', 'Ysn', 'Zbt'
      CHARACTER*3, DIMENSION(:), ALLOCATABLE :: DES_MO_CLASS   !(number of MO boundaries)
! Assigned a value in check_des_bc (see des_mi_layout):
!     When the inlet bc velocity is sufficiently high it is assigned
!     'RAND', otherwise it is assigned 'ORDR'      
      CHARACTER*4, DIMENSION(:), ALLOCATABLE :: PARTICLE_PLCMNT

! Particle injection factor; how many solid time steps (dtsolid) pass 
! before the next injection of a particle. if pi_count is greater than
! 1, then pi_factor is set to 1 (i.e. multiple particles enter every
! solids time step). 
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_FACTOR   !(number of MI boundaries)

! Particle injection count (injection number); how many particles are
! injected in one solids time step. pi_count is set to one if
! less than 1 particle enters per solids time step.     
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_COUNT   !(number of MI boundaries)


! Particle injection time scale; used when pi_factor > 1 to keep track of
! time needed for next injection     
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_MI_TIME
      

! Order inlet condition variables: used when the variable 
! particle_plcmnt is assigned 'ORDR'; in this case the inlet boundary is
! divided into grids and the particles are initially injected into an 
! available grid cell
! ----------------------------------------      
! keeps track of the number of grid cells that have been injected with a
! particle; is reset to 1 when all cells have been used.
      INTEGER, DIMENSION(:), ALLOCATABLE :: MI_FACTOR   !(number of MI boundaries)
! the inlet edge/face is gridded and this variable stores the length
! dimension of the grid (no smaller than a particle diameter)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MI_WINDOW   !(number of MI boundaries)
! the dimension of this variable is equal to the number of grid cells
! in the inlet edge/face
      TYPE dmi
         INTEGER, DIMENSION(:), POINTER :: VALUE   
      END TYPE dmi
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: I_OF_MI
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: J_OF_MI   !(number of MI boundaries)
! Construct an array of integers from 0 to a calculated factor in a random
! order that is used when placing new particles.        
      TYPE(dmi), DIMENSION(:), ALLOCATABLE :: MI_ORDER   !(number of MI boundaries)
! ----------------------------------------           

! Grid search loop counter arrays
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: GS_ARRAY   !(number of MI boundaries, 6)


     
      END MODULE DES_BC

