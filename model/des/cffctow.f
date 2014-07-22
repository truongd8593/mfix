!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFFCTOW
!  Purpose: Calculate the total force and torque on a particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 02-Aug-07
!
!  Comments: Implement eqns 13 & 14 from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFFCTOW(L, II,  NORM, DIST_LI, FC_tmp, FN_tmp, FT_tmp, TOW_tmp)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: DES_RADIUS
      use geometry, only: DO_K
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! indices of particle-particle contact pair
      INTEGER, INTENT(IN) :: L, II
! distance between particle centers
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
      DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT) :: FC_tmp, FN_tmp, FT_tmp, TOW_tmp
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable for calculating torque on particle
      DOUBLE PRECISION :: CROSSP(3)
! distance from the contact point to the particle center
      DOUBLE PRECISION :: DIST_CL
!------------------------------------------------

! total contact force
      FC_tmp(:) = FC_tmp(:) + FN_tmp(:) + FT_tmp(:)

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
         (2.d0*DIST_LI)

! total torque
      IF(DO_K) THEN
         CALL DES_CROSSPRDCT(CROSSP, NORM, FT_TMP)
         TOW_tmp(:)  = TOW_tmp(:)  + DIST_CL*CROSSP(:)
      ELSE
         CROSSP(1) = NORM(1)*FT_TMP(2) - NORM(2)*FT_TMP(1)
         TOW_tmp(1)  = TOW_tmp(1)  + DIST_CL*CROSSP(1)
      ENDIF

      RETURN
      END SUBROUTINE CFFCTOW
