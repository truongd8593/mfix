!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL(L, II, VRl, TANGNT)                           C
!  Purpose: DES - Calculate relative velocity between a particle pair  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFRELVEL(L, II, VRl, TANGNT)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER L, KK, II
      DOUBLE PRECISION VRl(NDIM), TANGNT(NDIM), LVEL, RVEL

!-----------------------------------------------------------------------

      LVEL = 0.0
      RVEL = 0.0

      DO KK = 1, DIMN
         LVEL = (DES_VEL_NEW(KK,L) - DES_VEL_NEW(KK,II))
         RVEL = 0.0
         RVEL =(OMEGA_NEW(KK,L)*DES_RADIUS(L)+&
                OMEGA_NEW(KK,II)*DES_RADIUS(II))*TANGNT(KK)
         VRl(KK) = LVEL - RVEL
      END DO

!     PRINT *,'relative velocity', VRl(1), VRl(2)

      RETURN
      END SUBROUTINE CFRELVEL


