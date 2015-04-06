      MODULE INTERACTIVE

      use param1, only: UNDEFINED_I, UNDEFINED

      IMPLICIT NONE

      PUBLIC :: INTERACT
      PRIVATE

!----------------------------------------------------------------------!

      INTEGER :: ACTION
! Do nothing (release interactive control)
      INTEGER, PARAMETER :: NULL_ENUM=0
! Pause the simulation to await further commands.
      INTEGER, PARAMETER :: WAIT_ENUM=1
! Exit the simulation cleanly.
      INTEGER, PARAMETER :: EXIT_ENUM=2
! Take N time-steps then wait for next command.
      INTEGER, PARAMETER :: NSTEPS_ENUM=3
! Take N time-steps then wait for next command.
      INTEGER, PARAMETER :: WAIT_AT_ENUM=4


      INTEGER, PARAMETER :: WRITE_RES_ENUM=10

! Refresh the simulation by processing the mfix.dat file.
      INTEGER, PARAMETER :: REINITIALIZE_ENUM=100

! Integer input for action items
      INTEGER :: ACTION_INT
! Float input for action items
      DOUBLE PRECISION :: ACTION_DP


      INTEGER :: USER_DEFINED_STEPS = UNDEFINED_I
      INTEGER :: STEP_CNT

      DOUBLE PRECISION :: USER_WAIT_AT_TIME = UNDEFINED

      CONTAINS
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACT                                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACT(pEXIT_SIGNAL)

      use compar, only: myPE, PE_IO

      use mpi_utility
      use error_manager

      LOGICAL :: INTERACTING

      LOGICAL, INTENT(INOUT) :: pEXIT_SIGNAL

      CALL CHECK_INTERACT_SIGNAL(INTERACTING)

      DO WHILE(INTERACTING)

         CALL GET_INTERACTIVE_DATA

         SELECT CASE(ACTION)
         CASE(NULL_ENUM)
            INTERACTING = .FALSE.

         CASE(WAIT_ENUM)
            INTERACTING = .TRUE.

         CASE(EXIT_ENUM)
            pEXIT_SIGNAL = .TRUE.
            INTERACTING = .FALSE.

         CASE(NSTEPS_ENUM)
            CALL CHECK_INTERACTIVE_NSTEPS(INTERACTING)

         CASE(WAIT_AT_ENUM)
            CALL CHECK_INTERACTIVE_WAIT_AT(INTERACTING)

         CASE(REINITIALIZE_ENUM)
            CALL REINITIALIZE
            INTERACTING = .TRUE.

         CASE DEFAULT
            CALL UNKNOWN_INTERACTIVE_CMD
            INTERACTING = .TRUE.

         END SELECT

         IF(INTERACTING) CALL INTERACTIVE_WAIT
      ENDDO

      RETURN
      END SUBROUTINE INTERACT


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACT                                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACT_SIGNAL(pINTERACT)

      use compar, only: myPE, PE_IO
      use param1, only: UNDEFINED_I, UNDEFINED

      use run, only: TIME, DT

      use mpi_utility

      LOGICAL, INTENT(OUT) :: pINTERACT

! Root checks if there is there is information from the user.
      IF(myPE == PE_IO) THEN

         INQUIRE(file="interact.dat", exist=pINTERACT)

         IF(USER_DEFINED_STEPS /= UNDEFINED_I) THEN
            STEP_CNT = STEP_CNT + 1
            IF(STEP_CNT >= USER_DEFINED_STEPS) THEN
               STEP_CNT = UNDEFINED_I
               USER_DEFINED_STEPS = UNDEFINED_I
               pINTERACT = .TRUE.
               CALL INTERACTIVE_WAIT
            ENDIF
         ENDIF

         IF(USER_WAIT_AT_TIME /= UNDEFINED) THEN
            IF(TIME + 0.1d0*DT >= USER_WAIT_AT_TIME) THEN
               USER_WAIT_AT_TIME = UNDEFINED
               pINTERACT = .TRUE.
               CALL INTERACTIVE_WAIT
            ENDIF
         ENDIF

      ENDIF

      CALL BCAST(pINTERACT, PE_IO)

      RETURN
      END SUBROUTINE CHECK_INTERACT_SIGNAL


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: UNKNOWN_INTERACTIVE_CMD                                 !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE UNKNOWN_INTERACTIVE_CMD

      use compar, only: myPE, PE_IO

      IF(myPE == PE_IO) THEN
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' Unknown interactive command.'
         WRITE(*,*) ' Enter a different command.'
         WRITE(*,*) ' '
      ENDIF

      RETURN
      END SUBROUTINE UNKNOWN_INTERACTIVE_CMD

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_INTERACTIVE_NSTEPS                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACTIVE_NSTEPS(pINTERACT)

      use compar, only: myPE, PE_IO

      LOGICAL, INTENT(OUT) :: pINTERACT

      IF(ACTION_INT == UNDEFINED_I) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' ACTION: NSTEPS is not fully specified.'
            WRITE(*,*) ' I will wait until you fix it...'
            WRITE(*,*) ' '
         ENDIF
         pINTERACT = .TRUE.
      ELSE
         USER_DEFINED_STEPS = ACTION_INT
         STEP_CNT = 0
         pINTERACT = .FALSE.
      ENDIF

      RETURN
      END SUBROUTINE CHECK_INTERACTIVE_NSTEPS

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_INTERACTIVE_WAIT_AT                               !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACTIVE_WAIT_AT(pINTERACT)

      use compar, only: myPE, PE_IO

      LOGICAL, INTENT(OUT) :: pINTERACT

      IF(ACTION_DP == UNDEFINED_I) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' ACTION: WAIT_AT is not fully specified.'
            WRITE(*,*) ' I will wait until you fix it...'
            WRITE(*,*) ' '
         ENDIF
         pINTERACT = .TRUE.
      ELSE
         USER_WAIT_AT_TIME = ACTION_DP
         pINTERACT = .FALSE.
      ENDIF

      RETURN
      END SUBROUTINE CHECK_INTERACTIVE_WAIT_AT

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: GET_INTERACTIVE_DATA                                    !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GET_INTERACTIVE_DATA

      CALL INTERACTIVE_DATA_INIT
      CALL INTERACTIVE_DATA_READ
      CALL INTERACTIVE_DATA_FINL

      RETURN
      END SUBROUTINE GET_INTERACTIVE_DATA


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_DATA_INIT                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_INIT

      use error_manager

      ACTION = NULL_ENUM

      ACTION_INT = UNDEFINED_I
      ACTION_DP = UNDEFINED

      RETURN
      END SUBROUTINE INTERACTIVE_DATA_INIT


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_DATA_READ                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_READ

      use funits, only: UNIT_DAT
      use error_manager

      INTEGER :: IOS
      CHARACTER(LEN=512) :: LINE_STRING, STRING

      NAMELIST / INTERACT_INPUT / ACTION, ACTION_INT, ACTION_DP

      OPEN(UNIT=UNIT_DAT, FILE='interact.dat', STATUS='OLD', IOSTAT=IOS)

      READ_LP: DO
         READ(UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP
! Setup the namelist format as a string.
         STRING=''; STRING = '&INTERACT_INPUT '//&
            trim(adjustl(LINE_STRING))//'/'
! Process the input line.
         READ(STRING, NML=INTERACT_INPUT, IOSTAT=IOS)
      ENDDO READ_LP

      CLOSE(UNIT_DAT)

      RETURN
      END SUBROUTINE INTERACTIVE_DATA_READ


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INIT_INTERACTIVE_DATA                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_FINL

      use compar, only: myPE, PE_IO
      use mpi_utility

      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF(myPE == PE_IO) CALL SYSTEM('rm interact.dat')


      RETURN
      END SUBROUTINE INTERACTIVE_DATA_FINL

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_WAIT                                        !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_WAIT

      use compar, only: myPE, PE_IO
      use mpi_utility

      LOGICAL :: FOUND

      IF(myPE == PE_IO) THEN
         DO; INQUIRE(file="interact.dat", exist=FOUND)
            IF(FOUND) EXIT
            CALL SLEEP(1)
         ENDDO
      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)

      RETURN
      END SUBROUTINE INTERACTIVE_WAIT


      END MODULE INTERACTIVE
