CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_RO_g(L)                                           C
C  Purpose: Calculate gas density                                      C
C                                                                      C
C  Author: M. Syamlal                                 Date: 22-NOV-93  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      REAL FUNCTION CALC_RO_g(L)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
C              Passed value of IJK index
      INTEGER  L, IJK
C
      DOUBLE PRECISION MW
C
C     Function subroutines
C
      DOUBLE PRECISION CALC_MW, EOSG
C
      INCLUDE 'function.inc'
C
        IF(RO_g0 .EQ. UNDEFINED .AND. .NOT.WALL_AT(L)) THEN
          IF(MW_AVG .EQ. UNDEFINED) THEN
            MW = CALC_MW(X_g, DIMENSION_3, L, NMAX(0), MW_g)
          ELSE
            MW = MW_AVG
          ENDIF
          CALC_RO_g = EOSG(MW, P_g(L), T_g(L))
        ELSE
          CALC_RO_g = RO_g0
        ENDIF
C
      RETURN
      END
