CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_VOL                                               C
C  Purpose: Calculate the volume of the cells                          C
C                                                                      C
C  Author: P. Nicoletti                               Date: 01-AUG-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: DX,DY,DZ,IMIN1,IMAX1,IMAX2,JMIN1,JMAX1,JMAX2  C
C                        KMIN1,KMAX1,KMAX2,XDIST_SC,XDIST_VEC          C
C                        YDIST_SC,YDIST_VEC,ZDIST_SC,ZDIST_VEC         C
C  Variables modified: I,J,K,VOL_SC,VOL_U,VOL_V,VOL_W                  C
C                                                                      C
C  Local variables: LEN_K,LEN_K_KP,KEN_J,LEN_J_JP,LEN_I,LEN_I_IP       C
C                   XI,XIH                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_VOL
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'post3d.inc'
C
      REAL     DX_E, DY_N, DZ_T
      INTEGER  I, J, K, IJK
C
      INCLUDE 'function.inc'


      DZ_T(K) = HALF * (DZ(K) + DZ(Kp1(K)))
      DY_N(J) = HALF * (DY(J) + DY(Jp1(J)))
      DX_E(I) = HALF * (DX(I) + DX(Ip1(I)))
C
      DO K = 1,KMAX2
      DO J = 1,JMAX2
      DO I = 1,IMAX2
        IJK = FUNIJK(I,J,K)
        VOL_SC(IJK) = DX(I)   * DY(J)   * X(I)*DZ(K)
        VOL_U(IJK)  = DX_E(I) * DY(J)   * X_E(I)*DZ(K)
        VOL_V(IJK)  = DX(I)   * DY_N(J) * X(I)*DZ(K)
        VOL_W(IJK)  = DX(I)   * DY(J)   * X(I)*DZ_T(K)
      END DO
      END DO
      END DO
C
      RETURN
      END
