!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GET_DATA                                                C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_DATA 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE funits 
      USE compar      
      USE gridmap
      USE discretelement
      USE des_thermo
      USE des_rxns
      USE leqsol
      USE parallel
      USE qmom_kinetic_equation
      USE mfix_pic
      USE cutcell
      USE dashboard
      USE des_stl_functions, only: des_stl_preprocessing, allocate_des_stl_arrays
      USE visc_g, only: L_SCALE
      USE constant, only: L_SCALE0

      USE error_manager


      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! shift DX, DY and DZ values 
      LOGICAL, PARAMETER :: SHIFT = .TRUE.


! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST 
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0) 

! Initialize the error manager. This call occurs after the mfix.dat
! is read so that message verbosity can be set and the .LOG file 
! can be opened.
      CALL INIT_ERROR_MANAGER

! Write header in the .LOG file and to screen.
! Not sure if the values are correct or useful
      CALL WRITE_HEADER

! Open files
      CALL OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX)

! These checks verify that sufficient information was provided
! to setup the domain indices and DMP gridmap.
      CALL CHECK_GEOMETRY_PREREQS
      CALL CHECK_DMP_PREREQS

! Set up the physical domain indicies (cell index max/min values).
      CALL SET_MAX2

! Set constants
      CALL SET_CONSTANTS 
      
! Adjust partition for better load balance (done when RE_INDEXING is .TRUE.)
      CALL ADJUST_IJK_SIZE

! Partition the domain and set indices
      CALL GRIDMAP_INIT

! Basic geometry checks.
      CALL CHECK_GEOMETRY(SHIFT)
! Set grid spacing variables.
      CALL SET_GEOMETRY 

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_INTERNAL_SURFACES
      CALL CHECK_POINT_SOURCES

      CALL CHECK_CHEMICAL_RXNS

!--------------------------  ARRAY ALLOCATION -----------------------!

! Allocate array storage.
      CALL ALLOCATE_ARRAYS
      IF (QMOMK) CALL QMOMK_ALLOCATE_ARRAYS


!--------------------------  GEOMETRY CONTROLS -----------------------!





!----------------------  DOMAIN SPECIFIC CHECKS  --------------------!

      CALL CHECK_DATA_ODEPACK


! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG

! Compute area of boundary surfaces.
      CALL GET_BC_AREA 

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS 

! JFD: cartesian grid implementation
      CALL CHECK_DATA_CARTESIAN
      IF(CARTESIAN_GRID) THEN
         CALL CUT_CELL_PREPROCESSING
      ELSE
         CALL ALLOCATE_DUMMY_CUT_CELL_ARRAYS
      ENDIF

      IF(DISCRETE_ELEMENT) then 
         !RG. Allocate the DES arrays dimensioned by Eulerian grid parameters 
         !Right now the below routine allocates and defines xe, yn, and zt that are used 
         !for des_stl_preprocessing. 
         CALL DES_ALLOCATE_ARRAYS_EULERIAN_GEOM
         
         ! if using stl representation for particle-wall interactions in discrete model
         ! then allocate the required arrays. This was earlier done in CG routines, 
         ! but to extend this capabilty to non-CG setups, allocate these arrays independently.
         IF(USE_STL_DES) then 
            CALL ALLOCATE_DES_STL_ARRAYS 
         
            ! do the stl preprocessing for discrete model if it is specified as true
            ! or forced to true through earlier consistency checks 
            CALL DES_STL_PREPROCESSING
         
         end IF
      end IF
      
! Initialize all field variables as undefined
      CALL INIT_FVARS 

! This is all that happens in SET_L_SCALE so it needs moved, maybe
! this should go in int_fluid_var.?
!     CALL SET_L_SCALE 
      L_SCALE(:) = L_SCALE0 


!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER
            
      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER
             
      NIT_MIN = MAX_NIT
      NIT_MAX = 0
             
      N_DASHBOARD = 0


      RETURN  

      END SUBROUTINE GET_DATA
