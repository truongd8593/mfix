!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VOL                                               C
!  Purpose: Calculate the volume of the cells                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 01-AUG-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: DX,DY,DZ,IMIN1,IMAX1,IMAX2,JMIN1,JMAX1,JMAX2  C
!                        KMIN1,KMAX1,KMAX2,XDIST_SC,XDIST_VEC          C
!                        YDIST_SC,YDIST_VEC,ZDIST_SC,ZDIST_VEC         C
!  Variables modified: I,J,K,VOL,VOL_U,VOL_V,VOL_W                  C
!                                                                      C
!  Local variables: LEN_K,LEN_K_KP,KEN_J,LEN_J_JP,LEN_I,LEN_I_IP       C
!                   XI,XIH                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_VOL
!
!
      Use param
      Use param1
      Use geometry
      Use indices
      Use fldvar
      Use physprop
      Use post3d
      Use compar

      IMPLICIT NONE
!
      REAL     DX_E, DY_N, DZ_T
      INTEGER  I, J, K, IJK
!
      INCLUDE 'function.inc'


      DZ_T(K) = HALF * (DZ(K) + DZ(Kp1(K)))
      DY_N(J) = HALF * (DY(J) + DY(Jp1(J)))
      DX_E(I) = HALF * (DX(I) + DX(Ip1(I)))
!
      DO K = 1,KMAX2
      DO J = 1,JMAX2
      DO I = 1,IMAX2
        IJK = FUNIJK(I,J,K)
        VOL(IJK)    = DX(I)   * DY(J)   * X(I)*DZ(K)
        VOL_U(IJK)  = DX_E(I) * DY(J)   * X_E(I)*DZ(K)
        VOL_V(IJK)  = DX(I)   * DY_N(J) * X(I)*DZ(K)
        VOL_W(IJK)  = DX(I)   * DY(J)   * X(I)*DZ_T(K)
      END DO
      END DO
      END DO
!
      RETURN
      END
