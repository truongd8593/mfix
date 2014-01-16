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
      LOGICAL, PRIVATE :: ER_LOCK

      contains

!``````````````````````````````````````````````````````````````````````!
! Module: INIT_ERR_MSG                                                 !
!                                                                      !
! Purpose: Initialize the error manager.                               !
!......................................................................!
      SUBROUTINE INIT_ERR_MGR

      ER_LOCK = .FALSE.

      END SUBROUTINE INIT_ERR_MGR


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE INIT_ERR_MSG(CALLER)

      CHARACTER(LEN=*), intent(IN) :: CALLER

      ERR_MSG=''

      PRIVATE_CALLER = trim(CALLER)

      RETURN
      END SUBROUTINE INIT_ERR_MSG



!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE FLUSH_ERR_MSG(DEBUG, HEADER, FOOTER)

! Debug flag.
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG
      LOGICAL, INTENT(IN), OPTIONAL :: HEADER
      LOGICAL, INTENT(IN), OPTIONAL :: FOOTER

! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Line Counter
      INTEGER :: LC
! Local debug flag.
      LOGICAL :: D_FLAG
      LOGICAL :: H_FLAG
      LOGICAL :: F_FLAG

! Set the local debug flag.
      D_FLAG = merge(DEBUG, .FALSE., PRESENT(DEBUG))
! Set the header flag.
      H_FLAG = merge(HEADER, .TRUE., PRESENT(HEADER))
! Set the footer flag.
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

      RETURN

 1000 FORMAT(1x,70('*'))
 1001 FORMAT(1x,'From: ',A)

      END SUBROUTINE FLUSH_ERR_MSG



      END MODULE ERROR_MANAGER
