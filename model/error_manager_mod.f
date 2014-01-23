!----------------------------------------------------------------------!
! Module: ERROR_MANAGER                                                !
!                                                                      !
! Purpose: Unify error message handeling.                              !
!                                                                      !
!----------------------------------------------------------------------!
      MODULE ERROR_MANAGER

      implicit none

! Maximum number of lines a message can have before a flush is needed.
      INTEGER, PARAMETER :: LINE_COUNT  = 32
! Maximum number of characters per line.
      INTEGER, PARAMETER :: LINE_LENGTH = 256

! Character string for storing the error message.
      CHARACTER(LEN=LINE_LENGTH), DIMENSION(LINE_COUNT) :: ERR_MSG

! The name of the calling routine. Set by calling: INIT_ERR_MSG
      CHARACTER(LEN=128), PRIVATE :: PRIVATE_CALLER

! Lock. The lock must be set/unset for each calling routine. This
! is done by calling the start and end routines.
      LOGICAL, PRIVATE :: LOCKED

! Flag for writing messages to the screen.
      LOGICAL, PRIVATE :: SCR_LOG

      contains

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERROR_MANAGER                                       !
!                                                                      !
! Purpose: Initialize the error manager. This routine also opens the   !
! .LOG file(s) based on user input settings.                           !
!......................................................................!
      SUBROUTINE INIT_ERROR_MANAGER

! Global Variables:
!---------------------------------------------------------------------//
! Name given to current run.
      use run, only: RUN_NAME
! Flag: All ranks report errors.
      use output, only: ENABLE_DMP_LOG
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! Flag: Provide the full log.
      use output, only: FULL_LOG
! Rank ID of process
      use compar, only: myPE
! Rank ID for IO handeling
      use compar, only: PE_IO
! Number of ranks in parallel run.
      use compar, only: numPEs
! File unit for LOG messages.
      use funits, only: UNIT_LOG
! Undefined character string.
      use param1, only: UNDEFINED_C

! Global Routine Access:
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Log file name.
      CHARACTER(len=64) :: LOGFILE
      CHARACTER(len=64) :: FILE_NAME
! First non-blank character in run_name.
      INTEGER :: NB
! Integer error flag
      INTEGER :: IER(0:numPEs-1)

! Initizilae the error flag.
      IER = 0
! Initialize the lock.
      LOCKED = .FALSE.
! Clear the error message storage container.
      ERR_MSG = ''
! Clear the caller routine information.
      PRIVATE_CALLER = ''

! This turns on error messaging from all processes.
      DMP_LOG = (myPE == PE_IO) .OR. ENABLE_DMP_LOG
! Flag for printing screen messages.
      SCR_LOG = (myPE == PE_IO) .AND. FULL_LOG

! Verify the length of user-provided name.
      LOGFILE = ''
      NB = INDEX(RUN_NAME,' ')
! RUN_NAME length too short.
      IF(RUN_NAME == UNDEFINED_C .OR. NB <= 1) THEN
         IF(myPE  == PE_IO) WRITE (*, 1000) 'short'
         CALL MFIX_EXIT(myPE) 
! RUN_NAME length too long.
      ELSEIF(NB + 10 > LEN(LOGFILE)) THEN
         IF(myPE == PE_IO) WRITE (*, 1000) 'long'
         CALL MFIX_EXIT(myPE) 
! RUN_NAME legnth too just right.
      ELSE
! Specify the .LOG file name based on MPI Rank extenion.
         IF(numPEs == 1 .OR. .NOT.ENABLE_DMP_LOG) THEN
            WRITE(LOGFILE,"(A)")RUN_NAME(1:(NB-1))
         ELSEIF(numPEs <    10) THEN
            WRITE(LOGFILE,"(A,'_',I1.1)") RUN_NAME(1:(NB-1)), myPE
         ELSEIF(numPEs <   100) THEN
            WRITE(LOGFILE,"(A,'_',I2.2)") RUN_NAME(1:(NB-1)), myPE
         ELSEIF(numPEs <  1000) THEN
            WRITE(LOGFILE,"(A,'_',I3.3)") RUN_NAME(1:(NB-1)), myPE
         ELSEIF(numPEs < 10000) THEN
            WRITE(LOGFILE,"(A,'_',I4.4)") RUN_NAME(1:(NB-1)), myPE
         ELSE
            WRITE(LOGFILE,"(A,'_',I8.8)") RUN_NAME(1:(NB-1)), myPE
         ENDIF
      ENDIF

! Open the .LOG file. From here forward, all routines should store
! error messages (at a minimum) in the .LOG file.
      IF(DMP_LOG) THEN
         NB = len_trim(LOGFILE)+1
         CALL OPEN_FILE(LOGFILE, NB, UNIT_LOG, '.LOG', FILE_NAME,      &
            'APPEND', 'SEQUENTIAL', 'FORMATTED', 132,  IER(myPE))
      ENDIF

! Verify that the .LOG file was successfully opened. Otherwise, flag the
! error and abort.
      CALL GLOBAL_ALL_SUM(IER)
      IF(sum(IER) /= 0) THEN
         IF(myPE == PE_IO) WRITE(*,1001) trim(FILE_NAME)
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1000 FORMAT(2/,1X,70('*')/' From: INIT_ERROR_MANAGER',/               &
         ' Error 1000: RUN_NAME too ',A,'. Please correct the',        &
         ' mfix.dat file.',/1x,70('*'),2/)

 1001 FORMAT(2/,1X,70('*')/' From: INIT_ERROR_MANAGER',/               &
         ' Error 1001: Failed to open log file: ',A,/' Aborting run.'/,&
         1x,70('*'),2/)

      END SUBROUTINE INIT_ERROR_MANAGER


!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERR_MSG                                             !
!                                                                      !
! Purpose: Initialize the error manager for the local routine. This    !
! call is needed to set the caller routines name for error messages.   !
!......................................................................!
      SUBROUTINE INIT_ERR_MSG(CALLER)

! Rank ID of process
      use compar, only: myPE
! Rank ID for IO handeling
      use compar, only: PE_IO

      implicit none

      CHARACTER(LEN=*), intent(IN) :: CALLER

! Verify that the lock isn't set.
      IF(LOCKED)THEN
         IF(myPE == PE_IO) WRITE(*,1000) CALLER
         CALL MFIX_EXIT(myPE)
      ENDIF

! Set the error manager lock. This lock remains set until FNL_ERR_MSG
! is called. A hard stop is called if the lock is set.
      LOCKED = .TRUE.

! Clear out the error manager.
      ERR_MSG=''
! Store the caller routines name.
      PRIVATE_CALLER = trim(CALLER)

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> INIT_ERR_MSG',/     &
         ' Error 1000: LOCK set at call to initialize routine.',/      &
         ' Previous routine may not have called FINL_ERR_MSG.',/,      &
         ' CALLER: ',A,/1x,70('*'))

      END SUBROUTINE INIT_ERR_MSG


!``````````````````````````````````````````````````````````````````````!
! Subroutine: FINL_ERR_MSG                                             !
!                                                                      !
! Purpose: Finalize the error manager. The call is needed to clear out !
! old information and unset the lock.                                  !
!......................................................................!
      SUBROUTINE FINL_ERR_MSG

! Rank ID of process
      use compar, only: myPE
! Rank ID for IO handeling
      use compar, only: PE_IO

      implicit none

! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Line Counter
      INTEGER :: LC
! Number of non-empty lines.
      INTEGER :: COUNT

! Verify that the lock was set, otherwise flag the error.
      IF(.NOT.LOCKED)THEN
         IF(myPE == PE_IO) WRITE(*,1000) trim(PRIVATE_CALLER)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Verify that the error message container is empty.
      COUNT = 0
      DO LC = 1, LINE_COUNT
         LINE = ERR_MSG(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < 256 ) COUNT = COUNT + 1
      ENDDO

! If the error message container empty, finish the reset and exit.
      IF(COUNT /= 0) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,1001) trim(PRIVATE_CALLER)
! Write out the error message container contents.
            DO LC = 1, LINE_COUNT
               LINE = ERR_MSG(LC)
               LENGTH = len_trim(LINE)
               IF(0 < LENGTH .AND. LENGTH < 256 )                      &
                  WRITE(*,"(' LC ',I2.2,': LEN: ',I3.3,1x,A)")         &
                     LC, LENGTH, trim(LINE)
            ENDDO
            WRITE(*,"(1x,70('*'),2/)")
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

      PRIVATE_CALLER = ''
      ERR_MSG = ''
      LOCKED = .FALSE.

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1000: LOCK not set prior to finalize call.',/         &
         ' Previous routine may not have called INIT_ERR_MSG.',/,      &
         ' PRIVATE_CALLER: ',A,/1x,70('*'))

 1001 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1001: Error container ERR_MSG not empty.',/           &
         ' PRIVATE_CALLER: ',A,2/' Contents:')

      END SUBROUTINE FINL_ERR_MSG



!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE FLUSH_ERR_MSG(DEBUG, HEADER, FOOTER, ABORT)

! Rank ID of process
      use compar, only: myPE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Debug flag.
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG
! Flag to suppress the message header.
      LOGICAL, INTENT(IN), OPTIONAL :: HEADER
! Flag to suppress the message footer.
      LOGICAL, INTENT(IN), OPTIONAL :: FOOTER
! Flag to abort execution by invoking MFIX_EXIT.
      LOGICAL, INTENT(IN), OPTIONAL :: ABORT

! Local Variables:
!---------------------------------------------------------------------//
! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Line Counter
      INTEGER :: LC
! Local debug flag.
      LOGICAL :: D_FLAG
! Local flag to suppress writing the header.
      LOGICAL :: H_FLAG
! Local flag to suppress writing the footer.
      LOGICAL :: F_FLAG
! Local abort flag.
      LOGICAL :: A_FLAG

! Set the abort flag. Continue running by default.
      A_FLAG = merge(ABORT, .FALSE., PRESENT(ABORT))
! Set the local debug flag. Suppress debugging messages by default.
      D_FLAG = merge(DEBUG, .FALSE., PRESENT(DEBUG))
! Set the header flag. Write the header by default.
      H_FLAG = merge(HEADER, .TRUE., PRESENT(HEADER))
! Set the footer flag. Write the footer by default.
      F_FLAG = merge(FOOTER, .TRUE., PRESENT(FOOTER))

! Write out header infomration.
      IF(H_FLAG) THEN
         write(*,"(2/,'')")
         IF(D_FLAG) WRITE(*,"('--- HEADER --->')",advance="no")
         WRITE(*,1000)
         IF(D_FLAG) WRITE(*,"('--- HEADER --->')",advance="no")
         WRITE(*,1001) trim(PRIVATE_CALLER)
      ENDIF


! Write the message body.
      IF(D_FLAG)THEN
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(LENGTH == 0) THEN
               write(*,"('LC ',I2.2,': LEN: ',I3.3,1x,A)")             &
                  LC, LENGTH, "EMPTY."
            ELSEIF(LENGTH >=  LINE_LENGTH)THEN
               write(*,"('LC ',I2.2,': LEN: ',I3.3,1x,A)")             &
                  LC, LENGTH, "OVERFLOW."
            ELSE
               write(*,"('LC ',I2.2,': LEN: ',I3.3,1x,A)")             &
                  LC, LENGTH, trim(LINE)
            ENDIF
         ENDDO
      ELSE
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(0 < LENGTH .AND. LENGTH < 256 ) THEN
               write(*,"(1x,A)") trim(LINE)
            ENDIF
         ENDDO
      ENDIF

! Clear the message array.
      ERR_MSG=''

! Print footer.
      IF(F_FLAG) THEN
         IF(D_FLAG) WRITE(*,"('--- FOOTER --->')",advance="no")
         WRITE(*,1000)
      ENDIF

! Abort the run if specified.
      IF(A_FLAG) CALL MFIX_EXIT(myPE)

      RETURN

 1000 FORMAT(1x,70('*'))
 1001 FORMAT(1x,'From: ',A)

      END SUBROUTINE FLUSH_ERR_MSG

      END MODULE ERROR_MANAGER
