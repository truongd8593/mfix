!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: geometry.inc                                           C
!  Purpose: Common block containing geometry and discretization data   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add variables DTODX_E, DToXDX_E, DToDY_N, DToDZ_T for      C
!           Variable Grid Size Capability                              C
!  Author: W. Rogers                                  Date: 16-APR-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      MODULE geometry


      Use param
      Use param1


!
!                      Coordinates: CARTESIAN, CYLINDRICAL
      CHARACTER(LEN=16)     COORDINATES
!
!                      Indicates whether x or r direction is not
!                      considered
      LOGICAL          NO_I
!
!                      Indicates whether x or r direction is
!                      considered
      LOGICAL          DO_I
!
!                      Starting index in the x or r direction
      INTEGER          IMIN1
!
!                      Number of cells in the x or r direction
      INTEGER          IMAX
!
!                      Number of cells in the x or r direction + 1
      INTEGER          IMAX1
!
!                      Number of cells in the x or r direction + 2
      INTEGER          IMAX2
!
!                      Cell sizes in the x or r direction
      DOUBLE PRECISION DX (0:DIM_I)
!
!                      Starting value of X.  This quantity is useful for
!                      simulating an annular cylindrical region.
      DOUBLE PRECISION XMIN
!
!                      Reactor length in the x or r direction
      DOUBLE PRECISION XLENGTH
!
!                      Indicates whether y direction is not considered
      LOGICAL          NO_J
!
!                      Indicates whether y direction is considered
      LOGICAL          DO_J
!
!                      Starting index in the y direction
      INTEGER          JMIN1
!
!                      Number of cells in the y direction
      INTEGER          JMAX
!
!                      Number of cells in the y direction + 1
      INTEGER          JMAX1
!
!                      Number of cells in the y direction + 2
      INTEGER          JMAX2
!
!                      Cell sizes in the y direction
      DOUBLE PRECISION DY (0:DIM_J)
!
!                      Reactor length in the y direction
      DOUBLE PRECISION YLENGTH
!
!                      Indicates whether z or theta direction is not
!                      considered
      LOGICAL          NO_K
!
!                      Indicates whether z or theta direction is
!                      considered
      LOGICAL          DO_K
!
!                      Starting index in the z or theta direction
      INTEGER          KMIN1
!
!                      Number of cells in the z or theta direction
      INTEGER          KMAX
!
!                      Number of cells in the z or theta direction + 1
      INTEGER          KMAX1
!
!                      Number of cells in the z or theta direction + 2
      INTEGER          KMAX2
!
!                      Cell sizes in the z or theta direction
      DOUBLE PRECISION DZ (0:DIM_K)
!
!                      Reactor length in the z or theta direction
      DOUBLE PRECISION ZLENGTH
!
!                      IMAX2 * JMAX2
      INTEGER          IJMAX2
!
!                      IMAX2 * JMAX2 * KMAX2
      INTEGER          IJKMAX2
!
!                      IMAX2 * JMAX2 * KMAX2
      INTEGER          IJKMAX3
!
!                      IJMAX2 + 1
      INTEGER          IJKMIN1
!
!                      IJKMAX2 - IJMAX2
      INTEGER          IJKMAX1
!
!                      Cell flags.
      INTEGER, DIMENSION(:), ALLOCATABLE ::           FLAG
!
!                      Cell flags with 3rd layer.
      INTEGER, DIMENSION(:), ALLOCATABLE ::           FLAG3
!
!                      Flag for the East surface
      INTEGER, DIMENSION(:), ALLOCATABLE ::           FLAG_E
!
!                      Flag for North surface
      INTEGER, DIMENSION(:), ALLOCATABLE ::           FLAG_N
!
!                      Flag for Top surface
      INTEGER, DIMENSION(:), ALLOCATABLE ::           FLAG_T
!
!                      Cell flags (bc/ic conditions)
!//PG allocatable type causes PG internal error, Ed's soln: pointers
!      CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: ICBC_FLAG
      character(LEN=3),  dimension(:), pointer :: icbc_flag
!
!                      1 / dx_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDX
!
!                      1 / dy_j
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDY
!
!                      1 / dz_k
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDZ
!
!                      1 / dx_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDX_E
!
!                      1 / dy_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDY_N
!
!                      1 / dz_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oDZ_T
!
!                      Radial location at cell center (x_i).
!                      X = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  X
      
!                      For cylindrical_2d simulation      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  cyl_X       
!
!                      Radial location at East face (x_i+1/2).
!                      X_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  X_E

!                      For cylindrical_2d simulation      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  cyl_X_E      
!
!                      Reciprocal of radial location at cell center (1/x_i).
!                      oX = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oX
!
!                      Reciprocal of radial location at East face (1/x_i+1/2).
!                      oX_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oX_E
!
!                      Azimuthal location at cell center (z_k).
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Z
!
!                      Azimuthal location at top face (z_k+1/2).
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Z_T
!
!                      one or more periodic boundary condition is used
      LOGICAL          CYCLIC
!
!                      Variable to flag periodic boundary condition in X
      LOGICAL          CYCLIC_X
!
!                      Variable to flag periodic boundary condition in Y
      LOGICAL          CYCLIC_Y
!
!                      Variable to flag periodic boundary condition in Z
      LOGICAL          CYCLIC_Z
!
!                      Variable to flag periodic bc with pressure drop in X
      LOGICAL          CYCLIC_X_PD
!
!                      Variable to flag periodic bc with pressure drop in Y
      LOGICAL          CYCLIC_Y_PD
!
!                      Variable to flag periodic bc with pressure drop in Z
      LOGICAL          CYCLIC_Z_PD
!
!                      Variable to flag periodic bc with mass flux in X
      LOGICAL          CYCLIC_X_MF
!
!                      Variable to flag periodic bc with mass flux in Y
      LOGICAL          CYCLIC_Y_MF
!
!                      Variable to flag periodic bc with mass flux in Z
      LOGICAL          CYCLIC_Z_MF
!
!                      Variable to flag cylindrical coordinates
      LOGICAL          CYLINDRICAL
      
!                      Variables for cylindrical_2d simulation
!		       Turn on the cylindrical_2d simulation 		
      logical          CYLINDRICAL_2D
!		       Variables for cylindrical_2d simulation      
!		       Half width of the plate in term of cell count 	      
      integer          I_CYL_NUM
!		       Variables for cylindrical_2d simulation      
!		       Cell number used to smooth the transition from plate to wedge 	      
      integer          I_CYL_TRANSITION      
      
!
!                      Factor for x direction averaging of U velocity: FX_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FX
!
!                      1 - FX_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FX_bar
!
!                      Factor for x direction averaging of scalars: FX_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FX_E
!
!                      1 - FX_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FX_E_bar
!
!                      Factor for y direction averaging of scalars: FY_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FY_N
!
!                      1 - FY_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FY_N_bar
!
!                      Factor for z direction averaging of scalars: FZ_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FZ_T
!
!                      1 - FZ_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  FZ_T_bar
!
!                      East face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AYZ
!
!                      North face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXZ
!
!                      Top face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXY
!
!                      Cell volume - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  VOL

!                      Total volume of cell's DES stencil neighbors
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  VOL_SURR
!
!                      East face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AYZ_U
!
!                      North face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXZ_U
!
!                      Top face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXY_U
!
!                      Cell volume - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  VOL_U
!
!                      East face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AYZ_V
!
!                      North face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXZ_V
!
!                      Top face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXY_V
!
!                      Cell volume - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  VOL_V
!
!                      East face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AYZ_W
!
!                      North face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXZ_W
!
!                      Top face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AXY_W
!
!                      Cell volume - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  VOL_W
!
!
          common /geometry_i/imin1, imax1, jmin1, jmax1 !for Linux

!     ADDED FOLLOWING VARIABLES TO TAKE CARE OF THE NEW CONVENTION - Pannala - 08/11/99

      INTEGER IMIN3,JMIN3,KMIN3,IMAX3,JMAX3,KMAX3, IMIN2, JMIN2, KMIN2

!     ADDED FOLLOWING VARIABLES TO TAKE CARE OF 4th order discretization in parallel

      INTEGER IMIN4,JMIN4,KMIN4,IMAX4,JMAX4,KMAX4, IJKMAX4, IJKMIN4

!!!HPF$ align FLAG(:) with TT(:)
!!!HPF$ align FLAG_E(:) with TT(:)
!!!HPF$ align FLAG_N(:) with TT(:)
!!!HPF$ align FLAG_T(:) with TT(:)
!!!!HPF$ align ICBC_FLAG(:) with TT(:)
!!!HPF$ align AYZ(:) with TT(:)
!!!HPF$ align AXZ(:) with TT(:)
!!!HPF$ align AXY(:) with TT(:)
!!!HPF$ align VOL(:) with TT(:)
!!!HPF$ align AYZ_U(:) with TT(:)
!!!HPF$ align AXZ_U(:) with TT(:)
!!!HPF$ align AXY_U(:) with TT(:)
!!!HPF$ align VOL_U(:) with TT(:)
!!!HPF$ align AYZ_V(:) with TT(:)
!!!HPF$ align AXZ_V(:) with TT(:)
!!!HPF$ align AXY_V(:) with TT(:)
!!!HPF$ align VOL_V(:) with TT(:)
!!!HPF$ align AYZ_W(:) with TT(:)
!!!HPF$ align AXZ_W(:) with TT(:)
!!!HPF$ align AXY_W(:) with TT(:)
!!!HPF$ align VOL_W(:) with TT(:)

      END MODULE geometry

!// Comments on the modifications for DMP version implementation
!//PG allocatable type causes PG internal error, Ed's soln: pointers
