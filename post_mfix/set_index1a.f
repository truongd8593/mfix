CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK,      C
C                            IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,       C
C                            IJKB, IJKT)                               C
C  Purpose: Set the indices of the neighbors of cell ijk (brute force) C
C                                                                      C
C  Author: M. Syamlal                                 Date: 21-JAN-92  C
C  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
C           W. Rogers                                                  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Modify index computations for K for setting periodic       C
C           boundary conditions in a cylindrical geometry where z goes C
C           from 0 to 2 pi                                             C
C  Author: M. Syamlal                                 Date: 10-MAR-92  C
C  Revision Number: 2                                                  C
C  Purpose:  Calculate only the nearest neighbor indices.( for code    C
C            optimization)                                             C
C  Author: M. Syamlal                                 Date: 23-SEP-92  C
C  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: I, J, K, IJK                                  C
C                                                                      C
C  Variables modified: IJKM, IJMK, IMJK, IPJK, IJPK, IJKP, IJKW, IJKE, C
C                      IJKS, IJKN, IJKB, IJKT                          C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SET_INDEX1A(I, J, K, IJK, IMJK, IPJK, IJMK, IJPK,
     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
     &                       IJKB, IJKT)
C
      IMPLICIT NONE
C
C  Include param.inc file to specify parameter values
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C     Physical and Numerical Parameters Section
C
      INCLUDE 'physprop.inc'
C
C     Geometry and Discretization Section
C
      INCLUDE 'geometry.inc'
C
      INCLUDE 'constant.inc'
C
C     Field Variables
C
      INCLUDE 'fldvar.inc'
C
C     Indices
C
      INCLUDE 'indices.inc'
C
C  Function subroutines
C
      LOGICAL COMPARE
C
C  Local variables
C
C                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,
     &                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT
CC
      INCLUDE 'function.inc'
C
C
      IMJK    = BOUND_FUNIJK(I-1, J, K)
      IPJK    = BOUND_FUNIJK(I+1, J, K)
      IJMK    = BOUND_FUNIJK(I, J-1, K)
      IJPK    = BOUND_FUNIJK(I, J+1, K)
      IJKM    = BOUND_FUNIJK(I, J, K-1)
      IJKP    = BOUND_FUNIJK(I, J, K+1)
C
C Modify indices as needed to allow cyclic boundary conditions
C
      IF(CYCLIC_X) THEN
        IF(I .EQ. IMAX1)
     &    IPJK    = BOUND_FUNIJK(IMIN1, J, K)
        IF(I .EQ. IMIN1)
     &    IMJK    = BOUND_FUNIJK(IMAX1, J, K)
      ENDIF
      IF(CYCLIC_Y) THEN
        IF(J .EQ. JMAX1)
     &    IJPK    = BOUND_FUNIJK(I, JMIN1, K)
        IF(J .EQ. JMIN1)
     &    IJMK    = BOUND_FUNIJK(I, JMAX1, K)
      ENDIF
      IF(CYCLIC_Z) THEN
        IF(K .EQ. KMAX1)
     &    IJKP    = BOUND_FUNIJK(I, J, KMIN1)
        IF(K .EQ. KMIN1)
     &    IJKM    = BOUND_FUNIJK(I, J, KMAX1)
      ENDIF
C
C  IJKW
C
      IF ( WALL_AT( IMJK ) ) THEN
        IJKW = IJK
      ELSE
        IJKW = IMJK
      ENDIF
C
C  IJKE
C
      IF ( WALL_AT( IPJK ) ) THEN
        IJKE  = IJK
      ELSE
        IJKE = IPJK
      ENDIF
C
C  IJKS
C
      IF ( WALL_AT( IJMK ) ) THEN
        IJKS = IJK
      ELSE
        IJKS = IJMK
      ENDIF
C
C  IJKN
C
      IF ( WALL_AT( IJPK ) ) THEN
        IJKN  = IJK
      ELSE
        IJKN = IJPK
      ENDIF
C
C  IJKB
C
      IF ( WALL_AT( IJKM ) ) THEN
        IJKB = IJK
      ELSE
        IJKB = IJKM
      ENDIF
C
C  IJKT
C
      IF ( WALL_AT( IJKP ) ) THEN
        IJKT  = IJK
      ELSE
        IJKT = IJKP
      ENDIF
C
      RETURN
      END
