      MODULE quadric
 

      Use param 
      Use param1

!     Maximum of the number of quadrics that can be read 
      INTEGER, PARAMETER          :: DIM_QUADRIC = 100
!     Nnumber of quadrics
      INTEGER                     :: N_QUADRIC 
!     Current Quadric
      INTEGER :: QUADRIC_ID
!     form of quadric : 'normal' or one of the pre-defined quadrics
      CHARACTER (LEN=10), DIMENSION(DIM_QUADRIC) :: quadric_form
!     Scale factor for quadrics
      DOUBLE PRECISION :: quadric_scale
!     Characteristic values of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Lambda_x,Lambda_y,Lambda_z
!     d - coefficient of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: dquadric
!     Translation components of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: t_x,t_y,t_z
!     Rotation angles (Deg) of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: theta_x,theta_y,theta_z
!     Radius for either Spere or Cylinder (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Radius
!     Radii for Torus 
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Torus_R1, Torus_R2
!     Half-angle for cone (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Half_angle
!     Normal vector components for plane (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: n_x,n_y,n_z
!     A-matrices of the quadrics
      DOUBLE PRECISION, DIMENSION(3,3,DIM_QUADRIC) :: A_QUADRIC
!     Translation-matrices of the quadrics
      DOUBLE PRECISION, DIMENSION(1,3,DIM_QUADRIC) :: T_QUADRIC
!     Clipping range  of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: clip_xmin,clip_xmax,clip_ymin,clip_ymax,clip_zmin,clip_zmax
!     Piecewise range  of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: piece_xmin,piece_xmax,piece_ymin,piece_ymax,piece_zmin,piece_zmax
!     Clip flag
      LOGICAL, DIMENSION(DIM_QUADRIC) :: FLUID_IN_CLIPPED_REGION
!     Boundary condition ID
      INTEGER :: BC_ID_Q(DIM_QUADRIC)
!     Maximum number of groups
      INTEGER, PARAMETER :: DIM_GROUP = 50
!     Number of groups
      INTEGER :: N_GROUP
!     Number of quadric in each group
      INTEGER,DIMENSION(DIM_GROUP) :: GROUP_SIZE
!     Quadric ID list in each group
      INTEGER,DIMENSION(DIM_GROUP,DIM_QUADRIC) :: GROUP_Q
!     Quadric relation in each group
      CHARACTER(LEN=9),DIMENSION(DIM_GROUP) :: GROUP_RELATION
!     Relation between groups
      CHARACTER(LEN=9),DIMENSION(DIM_GROUP) :: RELATION_WITH_PREVIOUS
!     Tolerance for intersection between quadrics and planes
      DOUBLE PRECISION :: TOL_F
!     Maximum number of iterations while finding intersection between geometry and grid
      INTEGER :: ITERMAX_INT

      END MODULE quadric
