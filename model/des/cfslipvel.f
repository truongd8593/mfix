!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIPVEL(SL, SII, VRl, SVEL, NORML)                   C
!  Purpose: DES - Calculate relative velocity between a particle pair  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Implements Eqn 8 from the following paper                 C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFSLIPVEL(L, II, SVL, VRl, Vno, NORML)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      INTEGER L, KK, II
      DOUBLE PRECISION LVEL(DIMN), RVEL(DIMN), SVL(DIMN), Vno
      DOUBLE PRECISION VRl(DIMN), NORML(DIMN), OMEGA_SUM(DIMN)

!-----------------------------------------------------------------------

      RVEL(:) = ZERO

      OMEGA_SUM(:) = OMEGA_NEW(L,:) + OMEGA_NEW(II,:)
      CALL DES_CROSSPRDCT(RVEL, OMEGA_SUM, NORML)
      RVEL(:) =DES_RADIUS(L)*RVEL(:)
 
      SVL(:) = VRl(:) - Vno*NORML(:) + RVEL(:)      
  
      RETURN
      END SUBROUTINE CFSLIPVEL


