!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_FILES                                             C
!  Purpose: open all the files for this run                            C
!                                                                      C
!  Author: P. Nicoletti                               Date: 12-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: EXT, FILE_NAME, LC, NB                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE machine 
      USE funits 
      USE compar 
      USE cdist

      use error_manager

      
      IMPLICIT NONE

! Error index: 0 - no error, 1 could not open file
      INTEGER :: IER(0:numPEs-1)
! RUN_NAME (as specified in input file)
      CHARACTER*(*) :: RUN_NAME
! Run_type (as specified in input file)
      CHARACTER*(*) :: RUN_TYPE
! Number of single precision output files (param.inc)
      INTEGER :: N_SPX

! local variables
      CHARACTER     EXT*4
!
!                   run_name + extension
      CHARACTER     FILE_NAME*64
!
!
!                   Loop counter
      INTEGER       LC
!
!                   index to first blank character in run_name
      INTEGER       NB, NBL

      CHARACTER     EXT_END*35 , cstatus*10

! Character error code.
      CHARACTER(len=32) :: CER
!-----------------------------------------------


! Initialize the error manager.
      CALL INIT_ERR_MSG("OPEN_FILES")

! Initialize the error flag array.
      IER = 0

! Initialize the generic SPx extension.
      EXT = '.SPx' 

! Generic SPx end characters in order.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Get the length of RUN_NAME. Note that len_trim would allow the
! name to still contain spaces. The following approach truncates
! RUN_NAME at the first blank character.
      NB = INDEX(RUN_NAME,' ')

! Only PE_IO opens the RUN_NAME.OUT file.
      IF(myPE == PE_IO) CALL OPEN_FILE (RUN_NAME, NB, UNIT_OUT, '.OUT',&
         FILE_NAME, 'UNKNOWN', 'SEQUENTIAL','FORMATTED',132, IER(myPE))

! Check if there was an error opening the file.
      IF(ERROR_OPENING(IER)) THEN
         WRITE(ERR_MSG,3000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Open the RES and SPx files. By default, only PE_IO opens these files, 
! but all ranks open a rank-specific copy for distributed IO runs.
      SELECT CASE (TRIM(RUN_TYPE))  

! Open the RES and SPx files for a new run.
!......................................................................
      CASE ('NEW')  

         IF(myPE==PE_IO .OR.  bDist_IO) THEN

! Open the RES file.
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, &
               'NEW', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
            IF (IER(myPE) == 100) THEN
               WRITE(ERR_MSG, 1000)'RES', 'NEW', trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX 
               EXT(4:4) = ext_end(LC:LC)
               CALL OPEN_FILE(RUN_NAME, NB, UNIT_SPX+LC, EXT,FILE_NAME,&
                  'NEW', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE)) 
! Report errors.
               IF (IER(myPE) == 100) THEN
                  WRITE(ERR_MSG, 1000)EXT(2:4), 'NEW', trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            ENDDO 
         ENDIF


! Open the RES and SPx files for a typical restart run.
!......................................................................
      CASE ('RESTART_1')

! Open the RES file.
         IF(myPE == PE_IO .or. bDist_IO) THEN
            CALL OPEN_FILE(RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME,  &
               'OLD', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
            IF (IER(myPE) == 101) THEN
               WRITE(ERR_MSG, 1001)'RES', 'RESTART_1',trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX 
               EXT(4:4) = EXT_END(LC:LC)
               CALL OPEN_FILE (RUN_NAME,NB, UNIT_SPX+LC,EXT, FILE_NAME,&
                  'OLD', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE)) 
! Report errors.
               IF (IER(myPE) == 101) THEN
                  WRITE(ERR_MSG, 1001) EXT(2:4), 'RESTART_1',         &
                     trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            END DO 
         ENDIF


! Open the RES and SPx files for a typical restart run.
!......................................................................
      CASE ('RESTART_2')  
! Open the RES file.
         CSTATUS = 'OLD'
         IF(myPE == PE_IO .OR. bDist_IO) THEN
            IF(bStart_with_one_res) CSTATUS = 'UNKNOWN'
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, &
               CSTATUS,'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE)) 
! Report errors.
            IF (IER(myPE) == 101) THEN
               WRITE(ERR_MSG, 1001)'RES', 'RESTART_2',trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX 
               EXT(4:4) = EXT_END(LC:LC)
               CALL OPEN_FILE (RUN_NAME,NB,UNIT_SPX+LC, EXT, FILE_NAME,&
                  'NEW' , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE)) 
! Report errors.
               IF (IER(myPE) == 100) THEN
                  WRITE(ERR_MSG, 1000)EXT(2:4), 'RESTART_2',          &
                     trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            END DO
         ENDIF

      CASE DEFAULT 
         WRITE(ERR_MSG, 3000)
         CALL FLUSH_ERR_MSG
         GO TO 100

      END SELECT 

! If an error was detected, abort the run.
  100 IF(ERROR_OPENING(IER)) CALL MFIX_EXIT(myPE)

! Initialize the error manager.
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: ',A,' file detected but RUN_TYPE=',A/,       &
         'Cannot open file: ',A)

 1001 FORMAT('Error 1001: ',A,' file missing for RUN_TYPE=',A/,        &
         'Cannot open file: ',A)

 2000 FORMAT('Error 2000: Unknown error opening file ',A,/             &
         'Error code: ',A)

 3000 FORMAT('Error 3000: Unknown run type: ',A)


      CONTAINS


!``````````````````````````````````````````````````````````````````````!
! FUNCTION: ERROR_OPENING                                              !
! Purpose: Collect the error flags from all processes and sum them.    !
! RESULT: .TRUE.  :: Sum of IER over all processes is non-zero.        !
!         .FALSE. :: GLOBAL_ALL_SUM is zero.                           !
!                                                                      !
!......................................................................!
      LOGICAL FUNCTION ERROR_OPENING(IER_l)

! MPI Wrapper function.
      use mpi_utility, only: GLOBAL_ALL_SUM

! Array containing error flags from all ranks.
      INTEGER, INTENT(IN) :: IER_L(0:numPEs-1)
! Initialize error flags.
      ERROR_OPENING = .FALSE.
! Globally collect flags.
      CALL GLOBAL_ALL_SUM(IER)
! Report errors.
      IF(sum(IER_l) /= 0) ERROR_OPENING = .TRUE.

      RETURN
      END FUNCTION ERROR_OPENING

      END SUBROUTINE OPEN_FILES 


















      SUBROUTINE OPEN_PE_LOG (IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE machine 
      USE funits 
      USE run
      USE compar 
      
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                     Error index: 0 - no error, 1 could not open file
      INTEGER         IER
!
!                   run_name + extension
      CHARACTER     FILE_NAME*64
!
!
!                   Log file name: dmp mode adds processor no to file name
      CHARACTER     LOGFILE*60
!
!                   Loop counter
      INTEGER       LC
!
!                   index to first blank character in run_name
      INTEGER       NB, NBL

!-----------------------------------------------

!
! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
!

!//PAR_I/O all PEs must exec this check in order to avoid Bcast of NB
      DO LC = 1, LEN(RUN_NAME) 
         IF (RUN_NAME(LC:LC) == ' ') THEN 
            NB = LC 
            EXIT 
         ENDIF
         LOGFILE(LC:LC) = RUN_NAME(LC:LC) 
      END DO 
!
      NBL = NB
      write(LOGFILE(NB:NB+3),'(I3.3)') myPE
      NBL = NB + 3
!
      CALL OPEN_FILE (LOGFILE, NBL, UNIT_LOG, '.LOG', FILE_NAME, 'NEW', &
          'SEQUENTIAL', 'FORMATTED', 132, IER(myPE))
        
      RETURN
      END SUBROUTINE OPEN_PE_LOG 
