!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Calculate the particle stress.                             !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC


! Global Variables:
!---------------------------------------------------------------------//
! Flag to use Snider's particle stress model
      use mfix_pic, only: MPPIC_SOLID_STRESS_SNIDER
! Particle stress
      use mfix_pic, only: PIC_P_S

! Module procedures:
!---------------------------------------------------------------------//
      use sendrecv, only: SEND_RECV

       use rxns

      IMPLICIT NONE

!......................................................................!


      IF(MPPIC_SOLID_STRESS_SNIDER) THEN
         CALL CALC_PS_PIC_SNIDER
      ELSE
         CALL CALC_PS_PIC_GARG
      ENDIF

      CALL SEND_RECV(PIC_P_S,1)

      ReactionRates(:,1) = PIC_P_S(:,1)

      RETURN
      END SUBROUTINE CALC_PS_PIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Evaluate the particle stress model of Snider.              !
!                                                                      !
!  REF: D.M. Snider, "Three-Dimensional Multiphase Particle-in-Cell    !
!     Model for Dense Particle Flows," Journal of Computational        !
!     Physics, Vol. 170, No. 2, pp. 523-549, 2001.                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC_SNIDER

! Global Variables:
!---------------------------------------------------------------------//
! Model parameters for Snider particle stress model
      use mfix_pic, only: FRIC_EXP_PIC
      use mfix_pic, only: PSFAC_FRIC_PIC
      use mfix_pic, only: FRIC_NON_SING_FAC
! Calculated particle stress
      use mfix_pic, only: PIC_P_S
! Fluid phase volume fraction
      use fldvar, only: EP_G
! Fluid volume fraction at close-pack
      use constant, only: EP_STAR
! Domain bounds
      use compar, only: IJKSTART3, IJKEND3
! Double precision parameters
      use param1, only: ONE


! Module procedures:
!---------------------------------------------------------------------//
      use functions, only: FLUID_AT

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: IJK
! Volume fraction of cell, modified for wall cells.
      DOUBLE PRECISION :: lEPg
! Volume fraction assigned to wall cells to prevent parcels from leaving
! the domain.
      DOUBLE PRECISION, PARAMETER :: WALL_EPg = 0.1d0

      double precision :: epg_min, ps_max
!......................................................................!

      epg_min = 1.0
      ps_max = -1.0

      DO IJK = IJKSTART3, IJKEND3
! Set a high value for the volume fraction in wall cells.
         lEPG = merge(EP_G(IJK), WALL_EPg, FLUID_AT(IJK))
! Particle stress :: Snider (Eq 33)
         PIC_P_S(IJK,1) = PSFAC_FRIC_PIC *((ONE - lEPg)**FRIC_EXP_PIC)/&
            MAX(lEPg - EP_STAR, FRIC_NON_SING_FAC*lEPg)

         if(fluid_at(ijk)) then
            epg_min = min(epg_min, ep_g(ijk))
            ps_max = max(ps_max, PIC_P_S(IJK,1))
         endif

      ENDDO

      write(*,"(/3x,'Epg Min: ',f15.4,/3x,'Ps  Max: ',f15.4)") &
         epg_min, ps_max

      RETURN
      END SUBROUTINE CALC_PS_PIC_SNIDER


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_PS_PIC                                             !
!  Author: R. Garg                                                     !
!                                                                      !
!  Purpose: Evaluate the particle stress as a coloring function:       !
!     X=0.0 :: cells below packing limit                               !
!     X=EPg :: cells above packing limit                               !
!     X=1.0 :: wall cells (far avobe max packing)                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_PS_PIC_GARG


! Global Variables:
!---------------------------------------------------------------------//
! Calculated particle stress
      use mfix_pic, only: PIC_P_S
! Resulting particle stress force
      use mfix_pic, only: PS_FORCE_PIC
! Fluid phase volume fraction
      use fldvar, only: EP_G
! Fluid volume fraction at close-pack
      use constant, only: EP_STAR
! Domain bounds
      use compar, only: IJKSTART3, IJKEND3
! Double precision parameters
      use param1, only: ZERO, ONE

! Module procedures:
!---------------------------------------------------------------------//
      use functions, only: FLUID_AT


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: IJK
!......................................................................!


! The Garg model uses a coloring function approach. 
      DO IJK = IJKSTART3, IJKEND3
         PS_FORCE_PIC(IJK,:) = ZERO
         IF(FLUID_AT(IJK)) THEN
            IF(EP_G(IJK) < EP_STAR) THEN
               PIC_P_S(IJK,1) = (ONE - EP_G(IJK))
            ELSE
               PIC_P_S(IJK,1) = ZERO
            ENDIF
         ELSE
            PIC_P_S(IJK,1) = ONE
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_PS_PIC_GARG