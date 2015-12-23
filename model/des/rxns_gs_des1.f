#include "version.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CONV_GS_DES1                                           !
!  Author: J.Musser: 16-Jun-10                                         !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the gas-particle convective heat transfer.                          !
!                                                                      !
!  Comments: Explicitly coupled simulations use a stored convective    !
!  heat transfer coefficient. Otherwise, the convective heat transfer  !
!  coeff is calculated every time step and the total interphase energy !
!  transfered is 'stored' and used explictly in the gas phase. The     !
!  latter conserves all energy
!                                                                      !
!  REF: Zhou, Yu, and Zulli, "Particle scale study of heat transfer in !
!       packed and bubbling fluidized beds," AIChE Journal, Vol. 55,   !
!       no 4, pp 868-884, 2009.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE RXNS_GS_DES1

      use constant, only: Pi
! Flag: The fluid and discrete solids are explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DTSOLID
      use discretelement, only: PIJK
      use discretelement, only: MAX_PIP

      use particle_filter, only: DES_INTERP_ON

      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_SIZE

      use des_rxns, only: DES_R_gp, DES_R_gc
      use des_rxns, only: DES_R_PHASE, DES_SUM_R_g
      use des_rxns, only: DES_HOR_G

      use physprop, only: NMAX

      use functions, only: FLUID_AT
      use functions, only: IS_NORMAL

      use param1, only: ZERO
      use des_rxns, only: DES_R_s
      use des_thermo, only: RXNS_Qs

      use param1, only: DIMENSION_LM
      use param, only: DIMENSION_M

      IMPLICIT NONE

      INTEGER :: IJK, LC, NP


! Local gas phase values.
      DOUBLE PRECISION :: lRgp(NMAX(0)) ! Rate of species production
      DOUBLE PRECISION :: lRgc(NMAX(0)) ! Rate of species consumption
      DOUBLE PRECISION :: lHoRg         ! Heat of reaction
      DOUBLE PRECISION :: lSUMRg        ! lSUMRg

! Interphase mass transfer
      DOUBLE PRECISION :: lRPhase(DIMENSION_LM+DIMENSION_M-1)

      DOUBLE PRECISION :: WEIGHT

      DO NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP)) CYCLE

! Avoid convection calculations in cells without fluid (cut-cell)
         IF(.NOT.FLUID_AT(PIJK(NP,4))) THEN
            DES_R_s(NP,:) = ZERO
            RXNS_Qs(NP) = ZERO

! No additional calculations are needed for explictly coupled
         ELSEIF(.NOT.DES_EXPLICITLY_COUPLED) THEN

! Calculate the heat transfer coefficient.
            CALL CALC_RRATES_DES(NP, lRgp, lRgc, lRPhase, lHoRg, lSUMRg)

! Integrate over solids time step and store in global array. This
! needs updated when interpolation is reintroduced into thermo code.
!---------------------------------------------------------------------//
! Store the gas phase source terms.
            IF(DES_INTERP_ON) THEN
               DO LC=1,FILTER_SIZE
                  IJK = FILTER_CELL(LC,NP)
                  WEIGHT = DTSOLID*FILTER_WEIGHT(LC,NP)

                  DES_R_gp(IJK,:) = DES_R_gp(IJK,:)+lRgp*WEIGHT
                  DES_R_gc(IJK,:) = DES_R_gc(IJK,:)+lRgc*WEIGHT
                  DES_R_PHASE(IJK,:) = DES_R_PHASE(IJK,:)+lRPhase*WEIGHT
                  DES_HOR_G(IJK) = DES_HOR_G(IJK) + lHoRg*WEIGHT
                  DES_SUM_R_g(IJK) = DES_SUM_R_g(IJK) + lSUMRg*WEIGHT
               ENDDO
            ELSE
               IJK = PIJK(NP,4)
               WEIGHT = DTSOLID

               DES_R_gp(IJK,:) = DES_R_gp(IJK,:) + lRgp*WEIGHT
               DES_R_gc(IJK,:) = DES_R_gc(IJK,:) + lRgc*WEIGHT
               DES_R_PHASE(IJK,:) = DES_R_PHASE(IJK,:) + lRPhase*WEIGHT
               DES_HOR_G(IJK) = DES_HOR_G(IJK) + lHoRg*WEIGHT
               DES_SUM_R_g(IJK) = DES_SUM_R_g(IJK) + lSUMRg*WEIGHT
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE RXNS_GS_DES1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: RXNS_GS_GAS1                                            !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using         !
!  interpolated values for the gas velocity and volume fraction. The   !
!  The resulting sources are interpolated back to the fluid grid.      !
!                                                                      !
!  NOTE: The loop over particles includes ghost particles so that MPI  !
!  communications are needed to distribute overlapping force between   !
!  neighboring grid cells. This is possible because only cells "owned" !
!  by the current process will have non-zero weights.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE RXNS_GS_GAS1

! Size of particle array on this process.
      use discretelement, only: MAX_PIP
! Flag to use interpolation
      use particle_filter, only: DES_INTERP_ON
! Interpolation cells and weights
      use particle_filter, only: FILTER_CELL,FILTER_WEIGHT,FILTER_SIZE
! IJK of fluid cell containing particles center
      use discretelement, only: PIJK
! Gas phase mass, species, and energy equation sources
      use des_rxns, only: DES_R_gp, DES_R_gc, DES_SUM_R_g
      use des_rxns, only: DES_R_PHASE, DES_HOR_g
      use physprop, only: NMAX
! MPI wrapper for halo exchange.
      use sendrecv, only: SEND_RECV

      use functions, only: FLUID_AT
      use functions, only: IS_NONEXISTENT
      use functions, only: IS_ENTERING, IS_ENTERING_GHOST
      use functions, only: IS_EXITING, IS_EXITING_GHOST

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE
      use param1, only: DIMENSION_LM
      use param, only: DIMENSION_M

      IMPLICIT NONE

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Loop bound for filter
      DOUBLE PRECISION :: GAMMAxTp
! Local gas phase values.
      DOUBLE PRECISION :: lRgp(NMAX(0)) ! Rate of species production
      DOUBLE PRECISION :: lRgc(NMAX(0)) ! Rate of species consumption
      DOUBLE PRECISION :: lHoRg         ! Heat of reaction
      DOUBLE PRECISION :: lSUMRg        ! lSUMRg

! Interphase mass transfer
      DOUBLE PRECISION :: lRPhase(DIMENSION_LM+DIMENSION_M-1)

      DOUBLE PRECISION :: WEIGHT

! Initialize fluid cell values.
      DES_R_gp = ZERO
      DES_R_gc = ZERO
      DES_R_PHASE = ZERO
      DES_HOR_G = ZERO
      DES_SUM_R_g = ZERO

! Calculate the gas phase forces acting on each particle.

      DO NP=1,MAX_PIP

         IF(IS_NONEXISTENT(NP)) CYCLE
         IF(.NOT.FLUID_AT(PIJK(NP,4))) CYCLE

! The drag force is not calculated on entering or exiting particles
! as their velocities are fixed and may exist in 'non fluid' cells.
         IF(IS_ENTERING(NP) .OR. IS_ENTERING_GHOST(NP) .OR. &
            IS_EXITING(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE

! Calculate the rates of species formation/consumption.
         CALL CALC_RRATES_DES(NP, lRgp, lRgc, lRPhase, lHoRg, lSUMRg)

! Store the gas phase source terms.
         IF(DES_INTERP_ON) THEN
            DO LC=1,FILTER_SIZE
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)

               DES_R_gp(IJK,:) = DES_R_gp(IJK,:)+lRgp*WEIGHT
               DES_R_gc(IJK,:) = DES_R_gc(IJK,:)+lRgc*WEIGHT
               DES_R_PHASE(IJK,:) = DES_R_PHASE(IJK,:)+lRPhase*WEIGHT
               DES_HOR_G(IJK) = DES_HOR_G(IJK) + lHoRg*WEIGHT
               DES_SUM_R_g(IJK) = DES_SUM_R_g(IJK) + lSUMRg*WEIGHT
            ENDDO
         ELSE
            IJK = PIJK(NP,4)

            DES_R_gp(IJK,:) = DES_R_gp(IJK,:) + lRgp
            DES_R_gc(IJK,:) = DES_R_gc(IJK,:) + lRgc
            DES_R_PHASE(IJK,:) = DES_R_PHASE(IJK,:) + lRPhase
            DES_HOR_G(IJK) = DES_HOR_G(IJK) + lHoRg
            DES_SUM_R_g(IJK) = DES_SUM_R_g(IJK) + lSUMRg
         ENDIF
      ENDDO

! Update the species mass sources in ghost layers.
      CALL SEND_RECV(DES_R_gp, 2)
      CALL SEND_RECV(DES_R_gc, 2)
      CALL SEND_RECV(DES_R_PHASE, 2)
      CALL SEND_RECV(DES_HOR_g, 2)
      CALL SEND_RECV(DES_SUM_R_g, 2)

      RETURN
      END SUBROUTINE RXNS_GS_GAS1
