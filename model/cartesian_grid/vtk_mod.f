      MODULE vtk


      Use param
      Use param1
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of gas phase species
      use param, only: DIM_N_g
! Maximum number of solids phase species
      use param, only: DIM_N_s
! Maximum number of scalar equations
      use param, only: DIM_Scalar


      INTEGER NUMBER_OF_CELLS
      INTEGER NUMBER_OF_CUT_CELLS
      INTEGER NUMBER_OF_BLOCKED_CELLS
      INTEGER NUMBER_OF_STANDARD_CELLS
      INTEGER NUMBER_OF_VTK_CELLS

      LOGICAL             :: WRITE_VTK_FILES
      LOGICAL             :: TIME_DEPENDENT_FILENAME
      LOGICAL             :: RESET_FRAME_AT_TIME_ZERO=.TRUE.
      CHARACTER (LEN=64)  :: VTU_DIR
      CHARACTER (LEN=64)  :: VTK_FILENAME,FRAME_CHAR,VTU_FILENAME,PVD_FILENAME,PVTU_FILENAME
      CHARACTER (LEN=64)  :: VTU_FRAME_FILENAME='VTU_FRAME_INDEX.TXT'
      CHARACTER (LEN=512) :: BUFFER

      CHARACTER (LEN=1), PARAMETER  :: END_REC = CHAR(10)

      INTEGER :: BOUNDARY_UNIT=122
      INTEGER :: VTK_UNIT=123,VTU_UNIT=124,PVD_UNIT=125
      INTEGER :: PVTU_UNIT=126,VTU_FRAME_UNIT=127
      INTEGER, PARAMETER :: DIM_VTK_VAR = 20
      INTEGER, DIMENSION(DIM_VTK_VAR) :: VTK_VAR


      INTEGER :: POLY_COUNTER,NUMBER_OF_POINTS


      integer, allocatable :: GLOBAL_I_OF(:)
      integer, allocatable :: GLOBAL_J_OF(:)
      integer, allocatable :: GLOBAL_K_OF(:)
      integer, allocatable :: GLOBAL_CONNECTIVITY(:,:)
      integer, allocatable :: GLOBAL_CLEANED_CONNECTIVITY(:,:)
      integer, allocatable :: GLOBAL_NUMBER_OF_NODES(:)

      REAL, allocatable :: GLOBAL_COORDS_OF_POINTS(:,:)

      LOGICAL, allocatable :: GLOBAL_INTERIOR_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_BLOCKED_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_STANDARD_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_CUT_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_SNAP(:)
      DOUBLE PRECISION, allocatable :: GLOBAL_F_AT(:)

      double precision, allocatable :: GLOBAL_X_NEW_POINT(:)
      double precision, allocatable :: GLOBAL_Y_NEW_POINT(:)
      double precision, allocatable :: GLOBAL_Z_NEW_POINT(:)

      INTEGER :: GLOBAL_NUMBER_OF_NEW_POINTS


      LOGICAL :: GLOBAL_VAR_ALLOCATED

      LOGICAL :: GRID_INFO_PRINTED_ON_SCREEN

      LOGICAL :: WRITE_ANI_CUTCELL


      INTEGER :: VTU_offset

      LOGICAL, allocatable :: BELONGS_TO_VTK_SUBDOMAIN(:)

      INTEGER, PARAMETER :: DIMENSION_VTK = 100

! Current VTK region
      INTEGER :: VTK_REGION

! Time interval at which vtk files are saved    
      DOUBLE PRECISION :: VTK_DT(DIMENSION_VTK)

! Current vtk time      
      DOUBLE PRECISION :: VTK_TIME(DIMENSION_VTK)

! FRAME index of vtk file      
      INTEGER :: FRAME(DIMENSION_VTK)

! PVD file initialization flag
      LOGICAL :: PVD_FILE_INITIALIZED(DIMENSION_VTK)=.FALSE.

! Logical variable to determine whether an vtk region is defined
      LOGICAL :: VTK_DEFINED (DIMENSION_VTK)

! VTK region West face, X-coordinate
      DOUBLE PRECISION :: VTK_X_w (DIMENSION_VTK)

! VTK region East face, X-coordinate
      DOUBLE PRECISION :: VTK_X_e (DIMENSION_VTK)

! VTK region South face, Y-coordinate
      DOUBLE PRECISION :: VTK_Y_s (DIMENSION_VTK)

! VTK region North face, Y-coordinate
      DOUBLE PRECISION :: VTK_Y_n (DIMENSION_VTK)

! VTK region Bottom face, Z-coordinate
      DOUBLE PRECISION :: VTK_Z_b (DIMENSION_VTK)

! VTK region Top face, Z-coordinate
      DOUBLE PRECISION :: VTK_Z_t (DIMENSION_VTK)

! VTK filename base 
      CHARACTER(LEN=64) :: VTK_FILEBASE(DIMENSION_VTK)

! Gas phase volume fraction
      LOGICAL :: VTK_EP_g (DIMENSION_VTK)

! Gas pressure
      LOGICAL :: VTK_P_g (DIMENSION_VTK)

! Solids pressure
      LOGICAL :: VTK_P_star(DIMENSION_VTK)

! Macroscopic density of solids phases
      LOGICAL :: VTK_ROP_s(DIMENSION_VTK, DIM_M)

! Solids phase volume fraction
      LOGICAL :: VTK_EP_s (DIMENSION_VTK, DIM_M)

! Gas phase temperature
      LOGICAL :: VTK_T_g(DIMENSION_VTK)

! Solids phase temperature
      LOGICAL :: VTK_T_s(DIMENSION_VTK, DIM_M)

! Granular temperature
      LOGICAL :: VTK_Theta_m(DIMENSION_VTK, DIM_M)

! X-component of gas velocity
      LOGICAL :: VTK_U_g(DIMENSION_VTK)

! X-component of solids phase velocity
      LOGICAL :: VTK_U_s(DIMENSION_VTK, DIM_M)

! Y-component of gas velocity
      LOGICAL :: VTK_V_g(DIMENSION_VTK)

! Y-component of solids phase velocity
      LOGICAL :: VTK_V_s(DIMENSION_VTK, DIM_M)

! Z-component of gas velocity
      LOGICAL :: VTK_W_g(DIMENSION_VTK)

! Z-component of solids phase velocity
      LOGICAL :: VTK_W_s(DIMENSION_VTK, DIM_M)

! Gas species mass fractions
      LOGICAL :: VTK_X_g(DIMENSION_VTK, DIM_N_g)

! Solids species mass fractions
      LOGICAL :: VTK_X_s(DIMENSION_VTK, DIM_M, DIM_N_s)

! Scalar value
      LOGICAL :: VTK_Scalar(DIMENSION_VTK, DIM_scalar)

! K & Epsilon values
      LOGICAL :: VTK_K_Turb_G(DIMENSION_VTK)
      LOGICAL :: VTK_E_Turb_G(DIMENSION_VTK)


      END MODULE vtk
