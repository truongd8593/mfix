!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL(L, II, VRELTRANS)                             C
!  Purpose: DES - Calculate relative velocity between a particle pair  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: The relative velocity definition is now corrected         C
!            The tangent, tangential and normal velocity components    C
!            are now calculated in this routine only                   C
!                                                                      C
!  Comments: Relative (translational) velocity required                C
!  for Eqn 6  from the following paper                                 C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFRELVEL(L, II, VRELTRANS,VRN, VRT,  TANGNT, NORM)
      
      USE discretelement
      USE param1
      IMPLICIT NONE
      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      
      INTEGER L, KK, II
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION TANMOD, TANGNT(DIMN), NORM(DIMN), VSLIP(DIMN), V_ROT(DIMN), OMEGA_SUM(DIMN), VRN, VRT
!
!-----------------------------------------------------------------------

      
      V_ROT(:) = ZERO
      
      
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DES_RADIUS(L)+ OMEGA_NEW(II,:)*DES_RADIUS(II)
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DES_RADIUS(L)+ OMEGA_NEW(II,1)*DES_RADIUS(II)
         OMEGA_SUM(2) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)
      
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

      VRN = DES_DOTPRDCT(VRELTRANS,NORM)
      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
      TANMOD = ZERO
      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.0.) THEN
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      END IF
      
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)

      RETURN
      END SUBROUTINE CFRELVEL


