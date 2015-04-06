!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ADJUST_DT(IER, NIT)                                    !
!  Author: M. Syamlal                                 Date: FEB-10-97  !
!                                                                      !
!  Purpose: Automatically adjust time step.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION ADJUST_DT (IER, NIT)

! Global Variables:
!---------------------------------------------------------------------//
! User defined type of run: new or restart
      use run, only: RUN_TYPE
! User defined aximum number of iterations
      use leqsol, only: MAX_NIT
! User defined: min, max DT and adjustment factor
      use run, only: DT_MIN, DT_MAX, DT_FAC
! Flag: Use stored DT value for advancing TIME
      use run, only: USE_DT_PREV
! Flag: 2nd order time implementation
      use run, only: CN_ON
! Flag: Continue to run at DT_MIN
      use run, only: PERSISTENT_MODE
! The current number of time steps (value at restart).
      use run, only: NSTEP, NSTEPRST
! Current DT (1/DT) and direction of last change (+/-)
      use run, only: DT, oDT, DT_DIR


! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, ONE, UNDEFINED


! Module proceedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Integer flag: 0=Good, 100=initialize, otherwise bad.
      INTEGER, INTENT(IN) :: IER
! Number of iterations for current time step
      INTEGER, INTENT(IN) :: NIT


! Local Variables:
!---------------------------------------------------------------------//
! Number of steps in between DT adjustments.
      INTEGER, PARAMETER :: STEPS_MIN = 5
! Number of time steps since last DT adjustment
      INTEGER, SAVE :: STEPS_TOT=0
! number of iterations since last DT adjustment
      INTEGER, SAVE :: NIT_TOT=0
! Iterations per second for last dt
      DOUBLE PRECISION, SAVE :: NIToS=0.0
! Current number of iterations per second
      DOUBLE PRECISION :: NITOS_NEW
! Flag to half/double the current time step
      LOGICAL :: CN_ADJUST_DT
!......................................................................!


! Initialize the function result.
      ADJUST_DT = .FALSE.
      USE_DT_PREV = .FALSE.

! Steady-state simulation.
      IF (DT==UNDEFINED .OR. DT<ZERO) RETURN

! Local flag for adjusting the time step for CN implementation.
      CN_ADJUST_DT = CN_ON .AND. ((RUN_TYPE=='NEW' .AND. NSTEP>1) .OR. &
         (RUN_TYPE/='NEW' .AND. NSTEP >= (NSTEPRST+1)))


! Initialize first call from time march.
      IF (IER == 100) THEN
         DT_DIR = -1
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Iterate successfully converged.
!---------------------------------------------------------------------//
      ELSE IF(IER == 0) THEN

! Set back the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0D0*DT

! Calculate a new DT every STEPS_MIN time steps.
         IF(STEPS_TOT >= STEPS_MIN) THEN
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT)
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1)
            STEPS_TOT = 0
            NITOS = NITOS_NEW
            NIT_TOT = 0
            IF (DT_DIR > 0) THEN
               IF(NIT < MAX_NIT) DT = MIN(DT_MAX,DT/DT_FAC)
            ELSE
               DT = DT*DT_FAC
               IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ENDIF

! DT was modified. Use the stored DT should be used to update TIME.
            USE_DT_PREV = .TRUE.

! Write the convergence stats to the screen/log file.
            WRITE(ERR_MSG,"('DT=',A,3x,'NIT/s=',A)")  &
               trim(iVal(DT)), trim(iVal(nint(NITOS)))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

         ELSE
            STEPS_TOT = STEPS_TOT + 1
            NIT_TOT = NIT_TOT + NIT
         ENDIF
! No need to iterate again
         ADJUST_DT = .FALSE.
! Cut the timestep into half for 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Iterate failed to converge.
!---------------------------------------------------------------------//
      ELSE

! Reset the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0d0*DT

! Reset counters.
         DT_DIR = -1
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Reduce the step size.
         DT = DT*DT_FAC

! The step size has decreased to the minimum.
         IF (DT > DT_MIN) THEN

            WRITE(ERR_MSG,"(3X,'Recovered: Dt=',G11.5,' :-)')") DT
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            CALL RESET_NEW

! Iterate again with new dt
            ADJUST_DT = .TRUE.

! Cut the timestep for 2nd order accurate time implementation.
            IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Set the return flag stop iterating.
         ELSE

! Prevent DT from dropping below DT_MIN.
            IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ADJUST_DT = .FALSE.
         ENDIF
      ENDIF

! Update ONE/DT variable.
      ODT = ONE/DT

      RETURN
      END FUNCTION ADJUST_DT
