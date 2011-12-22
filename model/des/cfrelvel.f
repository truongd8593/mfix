!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL                                               C
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
!
!  Comments: Relative (translational) velocity required for eqn 6 
!  from the following paper:           
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical       
!    simulation of plug glow of cohesionless particles in a           
!    horizontal pipe", Powder technology, 71, 239-250, 1992          
!
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFRELVEL(L, II, VRN, VRT, TANGNT, NORM, DIST_LI)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L, II

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 

! translational relative velocity 
         VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2       
         DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
            (2.d0*DIST_LI)
         DIST_CI = DIST_LI - DIST_CL
         IF(DIMN.EQ.3) THEN
            OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
               OMEGA_NEW(II,:)*DIST_CI
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL + &
               OMEGA_NEW(II,1)*DIST_CI
            OMEGA_SUM(2) = ZERO
         ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
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


!==================================================================
      SUBROUTINE CFRELVEL_WALL(l,w_vel_l,VRN,VRT,TANGNT,NORM,DIST_LI)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L

      DOUBLE PRECISION w_vel_l(dimn)
      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL 

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 

! translational relative velocity 
         VRELTRANS(:) = DES_VEL_NEW(L,:) - w_vel_l(:)

! calculate the distance from the particle center to the wall
         DIST_CL = DIST_LI - DES_RADIUS(L)
         IF(DIMN.EQ.3) THEN
            OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL
            OMEGA_SUM(2) = ZERO
         ENDIF


      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
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
      END SUBROUTINE CFRELVEL_WALL
