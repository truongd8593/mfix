!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: REINITIALIZE                                            !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!  Author: P. Nicoletti                               Date: 04-DEC-91  !
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REINITIALIZE

      use run, only: REINITIALIZING

      use error_manager

      IMPLICIT NONE

      INTEGER :: IER

      WRITE(ERR_MSG,"('Reinitializing.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      IER = 0
      REINITIALIZING = .TRUE.

! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(2)

      CALL REINITIALIZE0(IER)

      REINITIALIZING = .FALSE.

      IF(IER /=0) THEN
         WRITE(ERR_MSG,"('Reinitialization failed!')")
      ELSE
         WRITE(ERR_MSG,"('Done.')")
      ENDIF
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN
      END SUBROUTINE REINITIALIZE




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: REINITIALIZE0                                           !
!  Purpose: read and verify input data, open files                     !
!                                                                      !
!  Author: P. Nicoletti                               Date: 04-DEC-91  !
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REINITIALIZE0(pIER)

      USE cutcell, only: CARTESIAN_GRID
      use coeff, only: INIT_COEFF

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: pIER
      INTEGER :: IER

! Set the default error flag to ERROR state.
      pIER = 1

      CALL CHECK_RUN_CONTROL(IER)
      IF(REINIT_ERROR()) RETURN

      CALL CHECK_NUMERICS(IER)
      IF(REINIT_ERROR()) RETURN

      CALL CHECK_OUTPUT_CONTROL(IER)
      IF(REINIT_ERROR()) RETURN

      CALL CHECK_GAS_PHASE
      IF(REINIT_ERROR()) RETURN

      CALL CHECK_SOLIDS_PHASES
      IF(REINIT_ERROR()) RETURN

      CALL CHECK_INITIAL_CONDITIONS
      IF(REINIT_ERROR()) RETURN
      CALL CHECK_BOUNDARY_CONDITIONS
      IF(REINIT_ERROR()) RETURN
      CALL CHECK_INTERNAL_SURFACES
      IF(REINIT_ERROR()) RETURN
      CALL CHECK_POINT_SOURCES

      CALL CHECK_CHEMICAL_RXNS
      IF(REINIT_ERROR()) RETURN
      CALL CHECK_ODEPACK_STIFF_CHEM
      IF(REINIT_ERROR()) RETURN

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW
      IF(REINIT_ERROR()) RETURN

! This is all that happens in SET_L_SCALE so it needs moved, maybe
! this should go in int_fluid_var.?
!     CALL SET_L_SCALE
!      L_SCALE(:) = L_SCALE0

! Set constant physical properties
      CALL SET_CONSTPROP
      IF(REINIT_ERROR()) RETURN

! Set initial conditions
      CALL SET_IC
      IF(REINIT_ERROR()) RETURN

! Set point sources.
      CALL SET_PS
      IF(REINIT_ERROR()) RETURN

! Set boundary conditions
      CALL ZERO_NORM_VEL
      IF(REINIT_ERROR()) RETURN
      CALL SET_BC0
      IF(REINIT_ERROR()) RETURN
      IF(CARTESIAN_GRID) CALL CG_SET_BC0
      IF(REINIT_ERROR()) RETURN

! Set gas mixture molecular weight
      CALL SET_MW_MIX_G
      IF(REINIT_ERROR()) RETURN

! Initialize densities.
      CALL SET_RO_G
      IF(REINIT_ERROR()) RETURN
      CALL SET_RO_S
      IF(REINIT_ERROR()) RETURN

! Initialize time dependent boundary conditions
      CALL SET_BC1
      IF(REINIT_ERROR()) RETURN

! Check the field variable data and report errors.
      IF(.NOT.CARTESIAN_GRID) CALL CHECK_DATA_20
      IF(REINIT_ERROR()) RETURN

! Parse residual strings
      CALL PARSE_RESID_STRING (IER)
      IF(REINIT_ERROR()) RETURN

      CALL RRATES_INIT(IER)
      IF(REINIT_ERROR()) RETURN

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)
      IF(REINIT_ERROR()) RETURN

! Made it here without error.
      pIER = 0

      RETURN
      END SUBROUTINE REINITIALIZE0
