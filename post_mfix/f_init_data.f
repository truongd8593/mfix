CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: F_INIT_DATA                                            C
C  Purpose:                                                            C
C                                                                      C
C  Author: P.Nicoletti                                Date: 05-JUN-95  C
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
      SUBROUTINE F_INIT_DATA
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'xforms.inc'
C
      INTEGER I
C
      I_MAX = IMAX2
      J_MAX = JMAX2
      K_MAX = KMAX2
      M_MAX = MMAX
C
      DO I = 0,MMAX
         N_MAX(I) = NMAX(I)
      END DO
C
      RN_SPX = REAL(N_SPX)
C
      IF (C_E.EQ.UNDEFINED) THEN
         E_PASS = -1.0
         HAVE_E = .FALSE.
      ELSE
         E_PASS = C_E
         HAVE_E = .TRUE.
      END IF
C
      IF (RO_g0.EQ.UNDEFINED) THEN
         HAVE_RO_g0 = .FALSE.
      ELSE
         HAVE_RO_g0 = .TRUE.
      END IF
C
      IF (MW_avg.EQ.UNDEFINED) THEN
         HAVE_MW_avg = .FALSE.
      ELSE
         HAVE_MW_avg = .TRUE.
      END IF
C
      RETURN
      END
      subroutine do_pan
      write (*,*) ' hello'
      return
      end
