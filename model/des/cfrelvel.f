!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL(L, II, VRN, VRT, TANGNT, NORM, WALLCONTACT)   C
!
!  Purpose: DES - Calculate relative velocity between a particle pair  
!
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: The relative velocity definition is now corrected         C
!            The tangent, tangential and normal velocity components    C
!            are now calculated in this routine only                   C
!                                                                      C
!>
!!  Comments: Relative (translational) velocity required                
!!  for Eqn 6  from the following paper                                
!!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical       
!!  simulation of plug glow of cohesionless particles in a           
!!  horizontal pipe", Powder technology, 71, 239-250, 1992          
!<
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFRELVEL(L, II, VRN, VRT, TANGNT, NORM, WALLCONTACT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!-----------------------------------------------         
      INTEGER L, II

! Marker for particle/wall contact (1=p/w)
      INTEGER WALLCONTACT

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

!-----------------------------------------------      
! Functions
!-----------------------------------------------   
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT         
!-----------------------------------------------   


! translational relative velocity 
      IF (WALLCONTACT.EQ.1) THEN
         VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_WALL_VEL(II,:))
      ELSE
         VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))
      ENDIF

! rotational contribution  : v_rot
      IF (WALLCONTACT.EQ.1) THEN
         IF(DIMN.EQ.3) THEN
            OMEGA_SUM(:) = OMEGA_NEW(L,:)*DES_RADIUS(L)
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(L,1)*DES_RADIUS(L)
            OMEGA_SUM(2) = ZERO
         ENDIF
      ELSE              
         IF(DIMN.EQ.3) THEN
            OMEGA_SUM(:) = OMEGA_NEW(L,:)*DES_RADIUS(L)+ OMEGA_NEW(II,:)*DES_RADIUS(II)
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(L,1)*DES_RADIUS(L)+ OMEGA_NEW(II,1)*DES_RADIUS(II)
            OMEGA_SUM(2) = ZERO
         ENDIF
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! relative surface velocity in tangential direction 
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
! the magnitude of the tangential vector      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction  
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)

      RETURN
      END SUBROUTINE CFRELVEL

