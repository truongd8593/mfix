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
      USE param1
      IMPLICIT NONE

      INTEGER L, KK, II
      DOUBLE PRECISION VRl(DIMN), TANGNT(DIMN), LVEL(DIMN), RVEL(DIMN)

!-----------------------------------------------------------------------

      LVEL(:) = ZERO
      RVEL(:) = ZERO

         LVEL(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))
         IF(DIMN.EQ.3) THEN
            RVEL(:) =(OMEGA_NEW(L,:)*DES_RADIUS(L)+&
                      OMEGA_NEW(II,:)*DES_RADIUS(II))*TANGNT(:)
         ELSE 
            RVEL(:) =(OMEGA_NEW(L,1)*DES_RADIUS(L)+&
                      OMEGA_NEW(II,1)*DES_RADIUS(II))*TANGNT(:)
         END IF        
         VRl(:) = LVEL(:) - RVEL(:)
         
      RETURN
      END SUBROUTINE CFRELVEL


