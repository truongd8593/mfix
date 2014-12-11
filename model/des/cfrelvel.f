!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFRELVEL
!  Purpose: Calculate the normal and tangential components of the
!           relative velocity between contacting particles
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 01-Aug-07
!
!  Comments: Relative (translational) velocity required for eqn 6
!  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFRELVEL(L, II, VRN, VSLIP, NORM, DIST_LI)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use geometry, only: DO_K
      USE discretelement, only: DES_VEL_NEW, DES_RADIUS, OMEGA_NEW
      USE param1, only: ZERO
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
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(OUT) :: VSLIP(3)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(OUT) :: VRN
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! translational relative velocity
      DOUBLE PRECISION :: VRELTRANS(3)
! rotational velocity at point of contact
      DOUBLE PRECISION :: V_ROT(3), OMEGA_SUM(3)
! distance from the contact point to the particle centers
      DOUBLE PRECISION :: DIST_CL, DIST_CI
!-----------------------------------------------

! translational relative velocity
         VRELTRANS(:) = (DES_VEL_NEW(:,L) - DES_VEL_NEW(:,II))

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
         DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
            (2.d0*DIST_LI)
         DIST_CI = DIST_LI - DIST_CL

         IF(DO_K) THEN
            OMEGA_SUM(:) = OMEGA_NEW(:,L)*DIST_CL + &
               OMEGA_NEW(:,II)*DIST_CI
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(1,L)*DIST_CL + &
               OMEGA_NEW(1,II)*DIST_CI
            OMEGA_SUM(2) = ZERO
            OMEGA_SUM(3) = ZERO
         ENDIF

! calculate the rotational relative velocity
      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFRELVEL_WALL
!  Purpose: Calculate the normal and tangential components of the
!           relative velocity between a particle and wall contact
!
!  Comments: this subroutine is the same as above but it used for
!            particle-wall contact rather than particle-particle
!            contact. so the wall velocity is passed to this routine
!            rather than index of a contacting particle
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFRELVEL_WALL(L,WALL_VEL,VRN,VSLIP,NORM,DIST_LI)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only: DES_VEL_NEW, DES_RADIUS, OMEGA_NEW
      use geometry, only: DO_K
      USE param1, only: ZERO
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! indices of particle
      INTEGER, INTENT(IN) :: L
! wall velocity
      DOUBLE PRECISION, INTENT(IN) :: WALL_VEL(3)
! distance between particle centers
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to wall
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(INOUT) :: VSLIP(3)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(INOUT) :: VRN
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! translational relative velocity
      DOUBLE PRECISION :: VRELTRANS(3)
! rotational velocity at point of contact
      DOUBLE PRECISION :: V_ROT(3), OMEGA_SUM(3)
! distance from the contact point to the particle centers
      DOUBLE PRECISION :: DIST_CL
!-----------------------------------------------

! translational relative velocity
         VRELTRANS(:) = DES_VEL_NEW(:,L) - WALL_VEL(:)

! calculate the distance from the particle center to the wall
         DIST_CL = DIST_LI - DES_RADIUS(L)
         IF(DO_K) THEN
            OMEGA_SUM(:) = OMEGA_NEW(:,L)*DIST_CL
         ELSE
            OMEGA_SUM(1) = OMEGA_NEW(1,L)*DIST_CL
            OMEGA_SUM(2) = ZERO
            OMEGA_SUM(3) = ZERO
         ENDIF

! calculate the rotational relative velocity
      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! slip velocity of the contact point
! Equation (8) in Tsuji et al. 1992
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL
