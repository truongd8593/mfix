CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: IK_AVG (VAR,VOLUME,AVG,JUSE)                           C
C  Purpose: Calculate the value of a variable at a specified height -  C
C           averaged over X & Z (weighted by volume of the cell)       C
C                                                                      C
C  Author: P. Nicoletti                               Date: 20-JUL-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMIN1, IMAX1, KMIN1, KMAX1, FLAG              C
C  Variables modified: I,K,IJK                                         C
C                                                                      C
C  Local variables: TOTVOL                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE IK_AVG (VAR,VOLUME,AVG,JUSE)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'geometry.inc'
C
C     passed arguments
C
      DOUBLE PRECISION  VAR(*)
      REAL              VOLUME(*) , AVG
      INTEGER           JUSE
      INTEGER           I, K, IJK
C
C     local variables
C
      REAL              TOTVOL
C
      INCLUDE 'function.inc'
C
      TOTVOL = 0.0
      AVG    = 0.0
      DO K = KMIN1,KMAX1
         DO I = IMIN1,IMAX1
            IJK = FUNIJK(I,JUSE,K)
            IF (FLAG(IJK).EQ.1) THEN
               AVG = AVG + VAR(IJK) * VOLUME(IJK)
               TOTVOL = TOTVOL + VOLUME(IJK)
            END IF
         END DO
      END DO
C
      AVG = AVG / TOTVOL
C
      RETURN
      END
