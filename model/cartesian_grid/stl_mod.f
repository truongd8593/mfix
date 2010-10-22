      MODULE stl
 

      Use param 
      Use param1

!     Maximum of the number of polygons that can be read 
      INTEGER, PARAMETER          :: DIM_STL = 100000
!     Number of facets
      INTEGER                     :: N_FACETS
!     Vertex Coordinates X ,Y and Z
      DOUBLE PRECISION, DIMENSION(DIM_STL,3,3) :: VERTEX
!     Face normal vector (normalized)
      DOUBLE PRECISION, DIMENSION(DIM_STL,3) :: NORM_FACE
!     TRANSLATION COMPONENTS
      DOUBLE PRECISION :: TX_STL,TY_STL,TZ_STL  
!     SCALING FACTOR
      DOUBLE PRECISION :: SCALE_STL
!     RANGE OF STL FILE:
      DOUBLE PRECISION :: XMIN_STL,XMAX_STL
      DOUBLE PRECISION :: YMIN_STL,YMAX_STL
      DOUBLE PRECISION :: ZMIN_STL,ZMAX_STL
!     VALUE OF F_STL OUTSIDE OF STL BOUNDING BOX
      DOUBLE PRECISION :: OUT_STL_VALUE
!     SMALLEST ANGLE FOR DETECTION OF SMALL TRIANGLES
      DOUBLE PRECISION :: STL_SMALL_ANGLE
!     DIRECTION OF RAY TRACED TO DETERMINE WHETHER A POINT IS INSIDE/OUTSIDE
      CHARACTER(LEN=3) :: STL_RAY_DIR
!     Tolerance for polygone edge detection
      DOUBLE PRECISION :: TOL_STL
!     Boundary condition ID
      INTEGER :: STL_BC_ID


      END MODULE stl
