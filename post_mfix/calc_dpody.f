CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_DPoDY (TAVG,NTAVG,NPOINTS)                        C
C  Purpose: Update the time-averaged DPoDY calculation                 C
C           Should FLAG be referenced ?????                            C
C                                                                      C
C  Author: P. Nicoletti, M. Syamlal                   Date: 16-FEB-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1      C
C                        P_g, DY                                       C
C  Variables modified: I,J,K,IJK,IJPK                                  C
C                                                                      C
C  Local variables: DPoDY, NCOUNT                                      C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DPoDY(TAVG,NTAVG,NPOINTS)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
C     Passed arguments
C
      DOUBLE PRECISION TAVG(DIMENSION_3,*)
      INTEGER          NPOINTS , NTAVG
      INTEGER        I, J, K, IJK, IJPK
C
C     local variables
C
      INTEGER          NCOUNT
      DOUBLE PRECISION DPoDY
C
      INCLUDE 'function.inc'
C
      DO J = JMIN1,JMAX1
         DPoDY = 0.0
         NCOUNT = 0
         DO K = KMIN1,KMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               IJPK = FUNIJK(I,J+1,K)
               NCOUNT = NCOUNT + 1
               DPoDY = DPoDY 
     &            -( P_g(IJPK) - P_g(IJK) ) * 2.0 / (DY(J+1)+DY(J))
            END DO
         END DO
         TAVG(J,NTAVG) = TAVG(J,NTAVG) + DPoDY / REAL(NCOUNT)
      END DO
      NPOINTS = NPOINTS + 1
C
      RETURN
      END
