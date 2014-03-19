!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFFCTOWALL
!  Purpose: Calculate the total force and torque on a particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 02-Aug-07!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE CFFCTOWALL(L, NORM, DIST_LI)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle index
      INTEGER, INTENT(IN) :: L
! distance between particle center and wall
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to wall
      DOUBLE PRECISION, INTENT(IN) :: NORM(DIMN)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable for calculating torque on particle
      DOUBLE PRECISION :: CROSSP(DIMN)
! temporary variable for particle L tangential force
      DOUBLE PRECISION FT_TMP(DIMN)
! distance from the contact point to the particle center
      DOUBLE PRECISION DIST_CL
!------------------------------------------------

! total contact force
      FC(:,L) = FC(:,L) + FN(:,L) + FT(:,L)

! temporary holder of tangential force
      FT_TMP(:) = FT(:,L)

! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)

! total torque
      IF(DIMN.EQ.3) THEN
         CALL DES_CROSSPRDCT(CROSSP, NORM, FT_TMP)
         TOW(:,L) = TOW(:,L) + DIST_CL*CROSSP(:)
      ELSE
         CROSSP(1) = NORM(1)*FT_TMP(2) - NORM(2)*FT_TMP(1)
         TOW(1,L) =  TOW(1,L) + DIST_CL*CROSSP(1)
      ENDIF

      RETURN
      END SUBROUTINE CFFCTOWALL


