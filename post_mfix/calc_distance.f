CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_DISTANCE (DIST_MIN,D,DIM_MAX2,DIST_SC,DIST_VEC)   C
C  Purpose: Calculate the scalar and vector distance for the passed    C
C           direction.                                                 C
C                                                                      C
C  Author: P. Nicoletti                               Date: 07-JUL-92  C
C  Reviewer:                                                           C
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
C  Local variables: I                                                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DISTANCE (DIST_MIN, D,DIM_MAX2,DIST_SC,DIST_VEC)
C
      IMPLICIT NONE
C
      REAL               DIST_SC(*) , DIST_VEC(*)
      DOUBLE PRECISION   DIST_MIN, D(*)
      INTEGER            DIM_MAX2 , I
C
      DIST_SC(1)  = DIST_MIN - 0.5 * D(1)
      DIST_VEC(1) = DIST_MIN
C
      DO I = 2,DIM_MAX2
         DIST_SC(I)  = DIST_VEC(I-1) + 0.5 * D(I)
         DIST_VEC(I) = DIST_VEC(I-1) + D(I)
      END DO
C
      RETURN
      END
