!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MACHINE_CONS                                           C
!  Purpose: set the machine constants    ( SGI ONLY )                  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: OPEN_N1, NWORDS_DP, NWORDS_R, N_WORDS_I         C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MACHINE_CONS
!
!
      Use machine
      IMPLICIT NONE
!
      OPEN_N1   = 512
      NWORDS_DP =  64
      NWORDS_R  = 128
      NWORDS_I  = 128
      JUST_FLUSH = .false.
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_RUN_ID                                             C
!  Purpose: get the run id for this run                                C
!                                                                      C
!  Author: P. Nicoletti                               Date: 16-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: add ndoe name                                              C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
!                      ID_SECOND, ID_NODE                              C
!                                                                      C
!  Local variables: TIME_ARRAY                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_RUN_ID
!
      USE param
      USE run
      IMPLICIT NONE
!
!             temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER*10 DATE, TIM, ZONE

      CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)
      
!     For SGI only
      CALL GETHOSTNAME(ID_NODE,64)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CPU_TIME (CPU)                                         C
!  Purpose: get the CPU time for the run                               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: TA, XT                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CPU_TIME(CPU)
      USE compar      !// 001 Include MPI header file
!
      IMPLICIT NONE
!
! passed arguments
!
!                      cpu time since start of run
      DOUBLE PRECISION CPU
!
! local variables
!
!
!                      TA(1) = user cpu time   TA(2) = system cpu time
!     REAL             TA(2) 
!
!                      XT = TA(1) + TA(2)
!     REAL             XT
!
!                      ETIME is an SGI system function which returns
!                      the elasped CPU time
!     REAL             ETIME
!
!      XT = MPI_WTIME()
!      XT  = ETIME(TA)
!     CPU = XT
      
      
!-------------------------------------------F90
!                       clock cycle
       INTEGER           COUNT

!                       number of cycles per second
       INTEGER           COUNT_RATE
      
!                       max number of cycles, after which count is reset to 0
       INTEGER           COUNT_MAX

       CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
      
       CPU           = DBLE(COUNT)/DBLE(COUNT_RATE)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: START_LOG                                              C
!  Purpose: does nothing ... for VAX/VMS compatibility (SGI ONLY)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE START_LOG
      IMPLICIT NONE
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: END_LOG                                                C
!  Purpose: flushes the log file                                       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE END_LOG
      USE funits
      IMPLICIT NONE
      IF(DMP_LOG)CALL FLUSH (UNIT_LOG)
      RETURN
      END
!
      subroutine slumber
      return
      end
!
      subroutine pc_quickwin
      return
      end
!
!
      subroutine flush_bin(iunit)
!
      use funits
      use run
      use machine
!
      implicit none
!
      integer    :: iunit , lc  , nb , ier
      CHARACTER  :: EXT*4 , file_name*64 

!
! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
!
      DO LC = 1, LEN(RUN_NAME) 
         IF (RUN_NAME(LC:LC) == ' ') THEN 
            NB = LC 
            GO TO 125 
         ENDIF 
      END DO 
      WRITE (*, *) 'RUN_NAME TOOOOOOO LOOOONG' 
      STOP  
 125  continue
      lc  = iunit - unit_spx
      ext = '.SPx'
      WRITE (EXT(4:4), 1000) LC         
      close (unit=iunit)

      CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, &
              'old' , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER)
!
 1000 FORMAT(I1)
      return
      end
!
!
      subroutine flush_res(iunit)
!
      use funits
      use run
      use machine
!
      implicit none
!
      integer    :: iunit , lc , nb , ier
      CHARACTER  :: EXT*4 , file_name*64 
!
! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
!
      DO LC = 1, LEN(RUN_NAME) 
         IF (RUN_NAME(LC:LC) == ' ') THEN 
            NB = LC 
            GO TO 125 
         ENDIF 
      END DO 
      WRITE (*, *) 'RUN_NAME TOOOOOOO LOOOONG' 
      STOP  
 125  continue
      close (unit=iunit)
!
      ext = '.RES'
      CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, EXT, FILE_NAME, &
              'old' , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER)
!
 1000 FORMAT(I1)
      return
      end
