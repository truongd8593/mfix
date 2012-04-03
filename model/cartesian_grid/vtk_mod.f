      MODULE vtk
 

      Use param 
      Use param1

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
      INTEGER :: FRAME,VTK_UNIT=123,VTU_UNIT=124,PVD_UNIT=125
      INTEGER :: PVTU_UNIT=126,VTU_FRAME_UNIT=127
      INTEGER, PARAMETER :: DIM_VTK_VAR = 20
      INTEGER, DIMENSION(DIM_VTK_VAR) :: VTK_VAR


      INTEGER :: POLY_COUNTER,NUMBER_OF_POINTS

      DOUBLE PRECISION :: VTK_DT, VTK_TIME

      integer, allocatable :: GLOBAL_I_OF(:)    
      integer, allocatable :: GLOBAL_J_OF(:)    
      integer, allocatable :: GLOBAL_K_OF(:)    
      integer, allocatable :: GLOBAL_CONNECTIVITY(:,:)    
      integer, allocatable :: GLOBAL_NUMBER_OF_NODES(:)   


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

      LOGICAL :: PVD_FILE_INITIALIZED=.FALSE.

      INTEGER :: VTU_offset

      END MODULE vtk
