
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: CALC_RRATE_DES                                     !
!                                                                      !
!  Purpose: This routine manages gas-solid reactions for the continuum !
!  phase.                                                              !
!                                                                      !
!  Author: J.Musser                                   Date: 16-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_RRATE_DES

      USE discretelement
      USE interpolation
      USE param1
      USE rxns

      USE energy

      INTEGER PC, NP
      INTEGER INTERP_IJK(2**DIMN)

      DOUBLE PRECISION INTERP_WEIGHTS(2**DIMN)

! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS


!------------------------------------------------------------------------->>>>
! Force the gas field to be uniform.
      R_gp(:,:) = ZERO
      RoX_gc(:,:) = ZERO
      R_PHASE(IJK,LM) = ZERO
      HOR_G(:) = ZERO
      SUM_R_G(:) = ZERO
      RETURN
!-------------------------------------------------------------------------<<<<

! Initialize the particle counter
      PC = 1
! Loop over the particles in the system
      DO NP = 1, MAX_PIP
! Exit the loop if all the particles in the system have been looped over
         IF(PC .GT. PIP) EXIT
! Cycle the loop if there is no particle associated with this index
         IF(.NOT.PEA(NP,1)) CYCLE
! Set debug flag
         FOCUS = .FALSE.
         IF(NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.

! Calculate time dependent physical properties
         CALL DES_PHYSICAL_PROP(NP, FOCUS)

! Determine the IJK values for the surrounding fluid cells and calculate
! the cell-centered interpolation weights. This is only needed if 
! the DES interpolation modules have been invoked.
         IF(DES_INTERP_ON) THEN
            INTERP_IJK(:) = -1
            INTERP_WEIGHTS(:) = ZERO
            CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS)
         ENDIF

! Obtain the convective heat transfer coefficient of the particle
         CALL DES_RRATES(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS, 'GAS')

      ENDDO

      RETURN

      END SUBROUTINE CALC_RRATE_DES
