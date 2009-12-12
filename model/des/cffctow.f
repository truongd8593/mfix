!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOW                                                C
!
!  Purpose: DES - Calculate the total force and torque on a particle
!
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: 2-D case torque calculation corrected                     C
!
!  Comments: Implement eqns 13 & 14 from the following paper:       
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical       
!    simulation of plug glow of cohesionless particles in a             
!    horizontal pipe", Powder technology, 71, 239-250, 1992              
!
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOW(L, II,  NORM, DIST_LI)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER L, II
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN)

! temporary variable for particle L tangential force
      DOUBLE PRECISION FT_TMP(DIMN) 
! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!-----------------------------------------------

      FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 
      FC(II,:) = FC(II,:) - FN(L,:) - FT(L,:)

! temporary holder of tangential force         
      FT_TMP(:) = FT(L,:)

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_li^2+a^2=ri^2;  dist_lj^2+a^2=rj^2       
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
         (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FT_TMP)
         TOW(L,:)  = TOW(L,:)  + DIST_CL*CROSSP(:)
         TOW(II,:) = TOW(II,:) + DIST_CI*CROSSP(:)
! remember torque is R cross FT, which, compared to I particle, are
! both negative for the J particle.  Therefore, the toqrue, unlike tangential
! and normal contact forces, will be in the same direction for both the 
! particles making the pair 
      ELSE 
         CROSSP(1) = NORM(1)*FT_TMP(2) - NORM(2)*FT_TMP(1)
         TOW(L,1)  = TOW(L,1)  + DIST_CL*CROSSP(1)
         TOW(II,1) = TOW(II,1) + DIST_CI*CROSSP(1)
      ENDIF 


      RETURN
      END SUBROUTINE CFFCTOW


