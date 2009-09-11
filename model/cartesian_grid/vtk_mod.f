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
      CHARACTER (LEN=32)  :: VTK_FILENAME,FRAME_CHAR
      CHARACTER (LEN=512) :: BUFFER

      CHARACTER (LEN=1), PARAMETER  :: END_REC = CHAR(10)

      INTEGER :: FRAME,VTK_UNIT
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

      double precision, allocatable :: GLOBAL_X_NEW_POINT(:)
      double precision, allocatable :: GLOBAL_Y_NEW_POINT(:) 
      double precision, allocatable :: GLOBAL_Z_NEW_POINT(:)  

      INTEGER :: GLOBAL_NUMBER_OF_NEW_POINTS

 
      LOGICAL :: GLOBAL_VAR_ALLOCATED

      LOGICAL :: GRID_INFO_PRINTED_ON_SCREEN

      LOGICAL :: WRITE_ANI_CUTCELL

      END MODULE vtk
