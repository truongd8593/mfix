!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL(I, J, VRELTRANS)                             C
!>
!!  Purpose: DES - Calculate relative velocity between a particle pair  
!<
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
      SUBROUTINE CFRELVEL(I, J, VRN, VRT, TANGNT, NORM)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      
      INTEGER I, J
      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)
            
!-----------------------------------------------------------------------

! translational relative velocity 
      VRELTRANS(:) = (DES_VEL_NEW(I,:) - DES_VEL_NEW(J,:))

! rotational contribution  : v_rot
      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(I,:)*DES_RADIUS(I)+ OMEGA_NEW(J,:)*DES_RADIUS(J)
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(I,1)*DES_RADIUS(I)+ OMEGA_NEW(J,1)*DES_RADIUS(J)
         OMEGA_SUM(2) = ZERO
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
      END IF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)

      RETURN
      END SUBROUTINE CFRELVEL
