!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFFCTOWALL                                             C
!
!  Purpose: DES - Calculate the total force and torque on a particle 
!
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: 2-D case torque calculation corrected                     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFFCTOWALL(L, NORM, DIST_LI)
      
      USE param1 
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L
      DOUBLE PRECISION NORM(DIMN), CROSSP(DIMN)
! temporary variable for particle L tangential force
      DOUBLE PRECISION FT_TMP(DIMN) 
! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL      
!---------------------------------------------------------------------

      FC(L,:) = FC(L,:) + FN(L,:) + FT(L,:) 

! temporary holder of tangential force         
      FT_TMP(:) = FT(L,:)
      
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FT_TMP)
         TOW(L,:) = TOW(L,:) + DIST_CL*CROSSP(:)
      ELSE 
         CROSSP(1) = NORM(1)*FT_TMP(2) - NORM(2)*FT_TMP(1)
         TOW(L,1) =  TOW(L,1) + DIST_CL*CROSSP(1)
      ENDIF 

      RETURN
      END SUBROUTINE CFFCTOWALL


