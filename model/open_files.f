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
!
      SUBROUTINE OPEN_FILES(RUN_NAME, RUN_TYPE, NN_SPX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE machine 
      USE funits 
      USE geometry           !//
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                     Error index: 0 - no error, 1 could not open file
      INTEGER         IER
!
!                   run_name (as specified in input file)
      CHARACTER*(*) RUN_NAME
!
!                   run_type (as specified in input file)
      CHARACTER*(*) RUN_TYPE
!
!                   number of single precision output files (param.inc)
      INTEGER       NN_SPX
!
! local variables
!
!                   Answer
      CHARACTER     ANS
!
!                   extension to filename
      CHARACTER     EXT*4
!
!                   run_name + extension
      CHARACTER     FILE_NAME*64
!
!                   Loop counter
      INTEGER       LC
!
!                   index to first blank character in run_name
      INTEGER       NB
!-----------------------------------------------
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
!
  125 CONTINUE 
      IF (NB + 4 > LEN(FILE_NAME)) THEN 
         WRITE (*, *) 'RUN_NAME TOOOOOOO LOOOONG' 
         STOP  
      ENDIF 
!
      CALL OPEN_FILE (RUN_NAME, NB, UNIT_LOG, '.LOG', FILE_NAME, 'NEW', &
         'SEQUENTIAL', 'FORMATTED', 132, IER) 
      IF (IER /= 0) THEN 
         CALL OPEN_FILE (RUN_NAME, NB, UNIT_LOG, '.LOG', FILE_NAME, 'OLD', &
            'SEQUENTIAL', 'FORMATTED', 132, IER) 
         IF (IER /= 0) GO TO 500
	 DO WHILE(IER ==0)
	   READ(UNIT_LOG,'(a)', IOSTAT = IER)ANS
	 ENDDO
	 BACKSPACE(UNIT_LOG) 
      ENDIF 
!
      CALL OPEN_FILE (RUN_NAME, NB, UNIT_OUT, '.OUT', FILE_NAME, 'UNKNOWN', &
         'SEQUENTIAL', 'FORMATTED', 132, IER) 
!
!
      EXT = '.SPx' 
      SELECT CASE (RUN_TYPE)  
      CASE ('NEW')  
         CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, 'NEW', &
            'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
         IF (IER /= 0) THEN 
            WRITE (*, 1001) FILE_NAME 
            GO TO 600 
         ENDIF 
!
         DO LC = 1, NN_SPX 
            WRITE (EXT(4:4), 1000) LC 
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, 'NEW'&
               , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
            IF (IER /= 0) GO TO 500 
         END DO 
      CASE ('RESTART_1')  
         CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, 'OLD', &
            'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
         IF (IER /= 0) THEN 
            WRITE (*, 1002) FILE_NAME 
            GO TO 600 
         ENDIF 
         DO LC = 1, NN_SPX 
            WRITE (EXT(4:4), 1000) LC 
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, 'OLD'&
               , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
            IF (IER /= 0) GO TO 500 
         END DO 
      CASE ('RESTART_2')  
         CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, 'OLD', &
            'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
         IF (IER /= 0) THEN 
            WRITE (*, 1002) FILE_NAME 
            GO TO 600 
         ENDIF 
         DO LC = 1, NN_SPX 
            WRITE (EXT(4:4), 1000) LC 
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, 'NEW'&
               , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
            IF (IER /= 0) GO TO 500 
         END DO 
      CASE DEFAULT 
         WRITE (*, *) ' OPEN_FILES: DO NOT KNOW HOW TO PROCESS' 
         WRITE (*, *) ' RUN_TYPE in the input file' 
         STOP  
      END SELECT 
      RETURN  
  500 CONTINUE 
      WRITE (*, 1100) FILE_NAME 
  600 CONTINUE 
      CALL SLUMBER 
      STOP  
!
 1000 FORMAT(I1) 
 1001 FORMAT(/70('*')//' From: OPEN_FILES',/&
         ' Error: NEW run -- .RES file should NOT be in the run directory'/&
         ' Cannot open new file -- ',A,/70('*')/) 
 1002 FORMAT(/70('*')//' From: OPEN_FILES',/&
         ' Error: RESTART run -- .RES file should be in the run directory'/&
         ' Cannot open existing file -- ',A,/70('*')/) 
 1100 FORMAT(/70('*')//' From: OPEN_FILES',/' Error: Cannot open file -- ',A,/&
         70('*')/) 
      END SUBROUTINE OPEN_FILES 
