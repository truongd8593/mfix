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
      USE cutcell, only: CARTESIAN_GRID
      use coeff, only: INIT_COEFF

      use error_manager

      IMPLICIT NONE

      INTEGER :: IER

      WRITE(ERR_MSG,"('Reinitializing.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      REINITIALIZING = .TRUE.

! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0)

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_INTERNAL_SURFACES
      CALL CHECK_POINT_SOURCES

      CALL CHECK_CHEMICAL_RXNS
      CALL CHECK_ODEPACK_STIFF_CHEM

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! This is all that happens in SET_L_SCALE so it needs moved, maybe
! this should go in int_fluid_var.?
!     CALL SET_L_SCALE
!      L_SCALE(:) = L_SCALE0

! Set constant physical properties
      CALL SET_CONSTPROP

! Set initial conditions
      CALL SET_IC

! Set point sources.
      CALL SET_PS

! Set boundary conditions
      CALL ZERO_NORM_VEL
      CALL SET_BC0
      IF(CARTESIAN_GRID) CALL CG_SET_BC0

! Set gas mixture molecular weight
      CALL SET_MW_MIX_G

! Initialize densities.
      CALL SET_RO_G
      CALL SET_RO_S

! Initialize time dependent boundary conditions
      CALL SET_BC1

! Check the field variable data and report errors.
      IF(.NOT.CARTESIAN_GRID) CALL CHECK_DATA_20

! Parse residual strings
      CALL PARSE_RESID_STRING (IER)

      CALL RRATES_INIT(IER)

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)

      CALL WRITE_OUT0

      REINITIALIZING = .FALSE.

      WRITE(ERR_MSG,"('Done.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)


      RETURN
      END SUBROUTINE REINITIALIZE
