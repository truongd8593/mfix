!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_NUMERICS                                          !
!  Purpose: Check the numerics control namelist section                !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_NUMERICS


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Use the Fourth-order scheme
      use run, only: FPFOI
! Flag: Use the Chi-scheme
      use run, only: CHI_SCHEME
! Flag: Set optimal LEQ solver parameters for parallel runs.
      use leqsol, only: OPT_PARALLEL
! Discretization scheme for various equations
      USE run, only: DISCRETIZE
! Solve system transpose
      use leqsol, only: DO_TRANSPOSE
! Minimize dot porducts in BiCGSTAB
      use leqsol, only: MINIMIZE_DOTPRODUCTS
! Report solver stats.
      use leqsol, only: SOLVER_STATISTICS
! Controls reduction of global sums for residuals.
      use run, only: DEBUG_RESID
! Linear equation, preconditioner sweep method.
      use leqsol, only: LEQ_SWEEP
! Linear equation solution method.
      use leqsol, only: LEQ_METHOD
! Calculate dot-products more efficiently for serial runs.
      use parallel, only: IS_SERIAL


! Global Parameters:
!---------------------------------------------------------------------//
! NONE

! Global Module proceedures:
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: L


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_NUMERICS")


! Check  Fourth-order scheme requirements.
      IF (FPFOI) THEN
         DO L = 1,8
            IF(DISCRETIZE(L) <= 1) THEN
               WRITE(ERR_MSG,2000)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
 2000 FORMAT('Error 2000: Fourth-order scheme (FPFOI) requireds ',     &
         'DISCRETIZE >= 2',/'for all equations. Please correct the ',  &
         'mfix.dat file.')
      ENDIF


! Check chi scheme requirements.
      IF(CHI_SCHEME)THEN
         IF(DISCRETIZE(7) .NE. 3 .AND. DISCRETIZE(7).NE.6) THEN
            WRITE(ERR_MSG,2001)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 2001 FORMAT('Error 2001: CHI_SCHEME for species equations is only ',  &
         'implemented',/'for SMART and MUSCL discretization schemes ', &
         '[DISCRTIZE(7)].',/'Please correct the mfix.dat file.')
        ENDIF
      ENDIF


! Set the optimizations for DMP runs.
      IF (OPT_PARALLEL) THEN
         IS_SERIAL = .FALSE.
         DO_TRANSPOSE = .FALSE.
         MINIMIZE_DOTPRODUCTS = .TRUE.
         SOLVER_STATISTICS = .TRUE.
         DEBUG_RESID = .FALSE.
         LEQ_SWEEP(1:2) = 'ASAS'
         LEQ_METHOD(1:2) = 2
         LEQ_METHOD(3:9) = 1
      ENDIF

! Finalize the error msg.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_NUMERICS
