!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFSLIPVEL(L, II, VSLIP, VRELTRANS, VRN, NORM)          C
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

      SUBROUTINE CFSLIPVEL(L, II, VSLIP, VRELTRANS, VRN, NORM)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      INTEGER L, KK, II
      DOUBLE PRECISION LVEL(DIMN), V_ROT(DIMN), VSLIP(DIMN), VRN
      DOUBLE PRECISION VRELTRANS(DIMN), NORM(DIMN), OMEGA_SUM(DIMN)

!-----------------------------------------------------------------------

      V_ROT(:) = ZERO

      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DES_RADIUS(L)+ OMEGA_NEW(II,:)*DES_RADIUS(II)
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DES_RADIUS(L)+ OMEGA_NEW(II,1)*DES_RADIUS(II)
         OMEGA_SUM(2) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)
 
      VSLIP(:) = VRELTRANS(:) - VRN*NORM(:) + V_ROT(:)      
  
      RETURN
      END SUBROUTINE CFSLIPVEL


