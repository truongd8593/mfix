CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: geometry.inc                                           C
C  Purpose: Common block containing geometry and discretization data   C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Add variables DTODX_E, DToXDX_E, DToDY_N, DToDZ_T for      C
C           Variable Grid Size Capability                              C
C  Author: W. Rogers                                  Date: 16-APR-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References: None                                C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      Coordinates: CARTESIAN, CYLINDRICAL
      CHARACTER*16     COORDINATES
C
C                      Indicates whether x or r direction is not
C                      considered
      LOGICAL          NO_I
C
C                      Indicates whether x or r direction is
C                      considered
      LOGICAL          DO_I
C
C                      Starting index in the x or r direction
      INTEGER          IMIN1
C
C                      Number of cells in the x or r direction
      INTEGER          IMAX
C
C                      Number of cells in the x or r direction + 1
      INTEGER          IMAX1
C
C                      Number of cells in the x or r direction + 2
      INTEGER          IMAX2
C
C                      Cell sizes in the x or r direction
      DOUBLE PRECISION DX (DIMENSION_I)
C
C                      Starting value of X.  This quantity is useful for
C                      simulating an annular cylindrical region.
      DOUBLE PRECISION XMIN
C
C                      Reactor length in the x or r direction
      DOUBLE PRECISION XLENGTH
C
C                      Indicates whether y direction is not considered
      LOGICAL          NO_J
C
C                      Indicates whether y direction is considered
      LOGICAL          DO_J
C
C                      Starting index in the y direction
      INTEGER          JMIN1
C
C                      Number of cells in the y direction
      INTEGER          JMAX
C
C                      Number of cells in the y direction + 1
      INTEGER          JMAX1
C
C                      Number of cells in the y direction + 2
      INTEGER          JMAX2
C
C                      Cell sizes in the y direction
      DOUBLE PRECISION DY (DIMENSION_J)
C
C                      Reactor length in the y direction
      DOUBLE PRECISION YLENGTH
C
C                      Indicates whether z or theta direction is not
C                      considered
      LOGICAL          NO_K
C
C                      Indicates whether z or theta direction is
C                      considered
      LOGICAL          DO_K
C
C                      Starting index in the z or theta direction
      INTEGER          KMIN1
C
C                      Number of cells in the z or theta direction
      INTEGER          KMAX
C
C                      Number of cells in the z or theta direction + 1
      INTEGER          KMAX1
C
C                      Number of cells in the z or theta direction + 2
      INTEGER          KMAX2
C
C                      Cell sizes in the z or theta direction
      DOUBLE PRECISION DZ (DIMENSION_K)
C
C                      Reactor length in the z or theta direction
      DOUBLE PRECISION ZLENGTH
C
C                      IMAX2 * JMAX2
      INTEGER          IJMAX2
C
C                      IMAX2 * JMAX2 * KMAX2
      INTEGER          IJKMAX2
C
C                      IJMAX2 + 1
      INTEGER          IJKMIN1
C
C                      IJKMAX2 - IJMAX2
      INTEGER          IJKMAX1
C
C                      Cell flags.
      INTEGER          FLAG (DIMENSION_3)
C
C                      Flag for the East surface
      INTEGER          FLAG_E (DIMENSION_3)
C
C                      Flag for North surface
      INTEGER          FLAG_N (DIMENSION_3)
C
C                      Flag for Top surface
      INTEGER          FLAG_T (DIMENSION_3)
C
C                      Cell flags (bc/ic conditions)
      CHARACTER*3      ICBC_FLAG (DIMENSION_3)
C
C                      1 / dx_i
      DOUBLE PRECISION oDX (DIMENSION_I)
C
C                      1 / dy_j
      DOUBLE PRECISION oDY (DIMENSION_J)
C
C                      1 / dz_k
      DOUBLE PRECISION oDZ (DIMENSION_K)
C
C                      1 / dx_i+1/2
      DOUBLE PRECISION oDX_E (DIMENSION_I)
C
C                      1 / dy_j+1/2
      DOUBLE PRECISION oDY_N (DIMENSION_J)
C
C                      1 / dz_k+1/2
      DOUBLE PRECISION oDZ_T (DIMENSION_K)
C
C                      Radial location at cell center (x_i).
C                      X = 1 in Cartesian coordinates.
      DOUBLE PRECISION X (DIMENSION_I)
C
C                      Radial location at East face (x_i+1/2).
C                      X_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION X_E (DIMENSION_I)
C
C                      Reciprocal of radial location at cell center (1/x_i).
C                      oX = 1 in Cartesian coordinates.
      DOUBLE PRECISION oX (DIMENSION_I)
C
C                      Reciprocal of radial location at East face (1/x_i+1/2).
C                      oX_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION oX_E (DIMENSION_I)
C
C                      one or more periodic boundary condition is used
      LOGICAL          CYCLIC
C
C                      Variable to flag periodic boundary condition in X
      LOGICAL          CYCLIC_X
C
C                      Variable to flag periodic boundary condition in Y
      LOGICAL          CYCLIC_Y
C
C                      Variable to flag periodic boundary condition in Z
      LOGICAL          CYCLIC_Z
C
C                      Variable to flag periodic bc with pressure drop in X
      LOGICAL          CYCLIC_X_PD
C
C                      Variable to flag periodic bc with pressure drop in Y
      LOGICAL          CYCLIC_Y_PD
C
C                      Variable to flag periodic bc with pressure drop in Z
      LOGICAL          CYCLIC_Z_PD
C
C                      Variable to flag cylindrical coordinates
      LOGICAL          CYLINDRICAL
C
C                      Factor for x direction averaging of U velocity: FX_i
      DOUBLE PRECISION FX (DIMENSION_I)
C
C                      1 - FX_i
      DOUBLE PRECISION FX_bar (DIMENSION_I)
C
C                      Factor for x direction averaging of scalars: FX_i+1/2
      DOUBLE PRECISION FX_E (DIMENSION_I)
C
C                      1 - FX_i+1/2
      DOUBLE PRECISION FX_E_bar (DIMENSION_I)
C
C                      Factor for y direction averaging of scalars: FY_j+1/2
      DOUBLE PRECISION FY_N (DIMENSION_J)
C
C                      1 - FY_j+1/2
      DOUBLE PRECISION FY_N_bar (DIMENSION_J)
C
C                      Factor for z direction averaging of scalars: FZ_k+1/2
      DOUBLE PRECISION FZ_T (DIMENSION_K)
C
C                      1 - FZ_k+1/2
      DOUBLE PRECISION FZ_T_bar (DIMENSION_K)
C
C                      East face area - scalar cell
      DOUBLE PRECISION AYZ (DIMENSION_3)
C
C                      North face area - scalar cell
      DOUBLE PRECISION AXZ (DIMENSION_I, DIMENSION_K)
C
C                      Top face area - scalar cell
      DOUBLE PRECISION AXY (DIMENSION_I, DIMENSION_J)
C
C                      Cell volume - scalar cell
      DOUBLE PRECISION VOL (DIMENSION_3)
C
C                      East face area - U cell
      DOUBLE PRECISION AYZ_E (DIMENSION_3)
C
C                      North face area - U cell
      DOUBLE PRECISION AXZ_E (DIMENSION_I, DIMENSION_K)
C
C                      Top face area - U cell
      DOUBLE PRECISION AXY_E (DIMENSION_I, DIMENSION_J)
C
C                      Cell volume - U cell
      DOUBLE PRECISION VOL_E (DIMENSION_3)
C
C                      East face area - V cell
      DOUBLE PRECISION AYZ_N (DIMENSION_3)
C
C                      North face area - V cell
      DOUBLE PRECISION AXZ_N (DIMENSION_I, DIMENSION_K)
C
C                      Top face area - V cell
      DOUBLE PRECISION AXY_N (DIMENSION_I, DIMENSION_J)
C
C                      Cell volume - V cell
      DOUBLE PRECISION VOL_N (DIMENSION_3)
C
C                      East face area - W cell
      DOUBLE PRECISION AYZ_T (DIMENSION_3)
C
C                      North face area - W cell
      DOUBLE PRECISION AXZ_T (DIMENSION_I, DIMENSION_K)
C
C                      Top face area - W cell
      DOUBLE PRECISION AXY_T (DIMENSION_I, DIMENSION_J)
C
C                      Cell volume - W cell
      DOUBLE PRECISION VOL_T (DIMENSION_3)
CC
C
      COMMON / GEOMETRY_DP /
     &       DX, XMIN, XLENGTH, DY, YLENGTH, DZ, ZLENGTH, 
     &       X, X_E, oX, oX_E,
     &       FX, FX_bar, FX_E, FX_E_bar, FY_N, FY_N_bar, FZ_T, FZ_T_bar,
     &       oDX, oDY, oDZ, oDX_E, oDY_N, oDZ_T,
     &       AYZ  , AXZ  , AXY  , VOL,
     &       AYZ_E, AXZ_E, AXY_E, VOL_E,
     &       AYZ_N, AXZ_N, AXY_N, VOL_N,
     &       AYZ_T, AXZ_T, AXY_T, VOL_T
C
      COMMON / GEOMETRY_I /
     &       IMIN1, IMAX, IMAX1, IMAX2, JMIN1, JMAX, JMAX1, JMAX2,
     &       KMIN1, KMAX, KMAX1, KMAX2, IJMAX2, IJKMAX2,
     &       IJKMIN1, IJKMAX1,
     &       FLAG, FLAG_E, FLAG_N, FLAG_T
C
      COMMON / GEOMETRY_L /
     &       NO_I, NO_J, NO_K, CYLINDRICAL, 
     &       DO_I, DO_J, DO_K, CYCLIC_X, CYCLIC_Y, CYCLIC_Z,
     &       CYCLIC_X_PD, CYCLIC_Y_PD, CYCLIC_Z_PD, CYCLIC
C
      COMMON / GEOMETRY_C /
     &       COORDINATES, ICBC_FLAG
C
