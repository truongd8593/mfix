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
      SUBROUTINE OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE machine 
      USE funits 
      USE compar 
      USE cdist
      
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
      INTEGER       N_SPX
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
!
!                   Log file name: dmp mode adds processor no to file name
      CHARACTER     LOGFILE*60
!
!                   Loop counter
      INTEGER       LC
!
!                   index to first blank character in run_name
      INTEGER       NB, NBL

      CHARACTER     EXT_END*35 , cstatus*10
!-----------------------------------------------

      ext_end = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
! DETERMINE THE FIRST BLANK CHARCATER IN RUN_NAME
!

!//PAR_I/O all PEs must exec this check in order to avoid Bcast of NB
      DO LC = 1, LEN(RUN_NAME) 
         IF (RUN_NAME(LC:LC) == ' ') THEN 
            NB = LC 
            GO TO 125 
         ENDIF
         LOGFILE(LC:LC) = RUN_NAME(LC:LC) 
      END DO 
      WRITE (*, *) 'RUN_NAME TOOOOOOO LOOOONG' 
      call mfix_exit(myPE) 
!
  125 CONTINUE 
      IF (NB + 7 > LEN(FILE_NAME)) THEN 
         WRITE (*, *) 'RUN_NAME TOOOOOOO LOOOONG' 
         call mfix_exit(myPE) 
      ENDIF 
!
      NBL = NB
      if ( numPEs > 1 ) then
         write(LOGFILE(NB:NB+3),'(I3.3)') myPE
         NBL = NB + 3
      endif
!
      IF (DMP_LOG)Then
         CALL OPEN_FILE (LOGFILE, NBL, UNIT_LOG, '.LOG', FILE_NAME, 'NEW', &
          'SEQUENTIAL', 'FORMATTED', 132, IER) 

         IF (IER /= 0) THEN 
            CALL OPEN_FILE (LOGFILE, NBL, UNIT_LOG, '.LOG', FILE_NAME, 'OLD', &
                           'SEQUENTIAL', 'FORMATTED', 132, IER) 
            IF (IER /= 0) GO TO 500
            DO WHILE(IER ==0)
               READ(UNIT_LOG,'(a)', IOSTAT = IER)ANS
            ENDDO
            BACKSPACE(UNIT_LOG)
         ENDIF
      ENDIF

!//PAR_I/O only PE 0 opens the ASCI output (.out), restart (.res) and species (.spX) files 
      if ( myPE == PE_IO ) then  ! dem 25jun07
!
         CALL OPEN_FILE (RUN_NAME, NB, UNIT_OUT, '.OUT', FILE_NAME, 'UNKNOWN', &
                         'SEQUENTIAL', 'FORMATTED', 132, IER) 
      end if ! pan 27-jun-2007
!
!
      EXT = '.SPx' 
      SELECT CASE (TRIM(RUN_TYPE))  

      CASE ('NEW')  

        if (myPE==PE_IO .or.  bDist_IO) then 
           CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, 'NEW', &
                           'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
           IF (IER /= 0) THEN 
              WRITE (*, 1001) FILE_NAME 
              GO TO 600 
           ENDIF 
           DO LC = 1, N_SPX 
              ext(4:4) = ext_end(LC:LC)
              CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, 'NEW'&
                              , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
              IF (IER /= 0) GO TO 500 
           END DO 
        end if ! pan



      CASE ('RESTART_1') 

         if (myPE == PE_IO .or. bDist_IO) then ! pan  
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, 'OLD', &
                           'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
            IF (IER /= 0) THEN 
               WRITE (*, 1002) FILE_NAME 
               GO TO 600 
            ENDIF 

           DO LC = 1, N_SPX 
              ext(4:4) = ext_end(LC:LC)
              CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, 'OLD'&
                           , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
              IF (IER /= 0) GO TO 500 
           END DO 

        end if ! pan

      CASE ('RESTART_2')  

	   cstatus = 'old'

         if (myPE == PE_IO .or. bDist_IO) then ! pan  

	      if (bStart_with_one_res) cstatus = 'unknown'

            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, cstatus, &
                         'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 

            IF (IER /= 0) THEN 
               WRITE (*, 1002) FILE_NAME 
               GO TO 600 
            ENDIF 

            DO LC = 1, N_SPX 
               ext(4:4) = ext_end(LC:LC)
               CALL OPEN_FILE (RUN_NAME, NB, UNIT_SPX + LC, EXT, FILE_NAME, &
                              'NEW' , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER) 
               IF (IER /= 0) GO TO 500 
            END DO 

        end if ! pan



      CASE DEFAULT 
          WRITE (*, *) ' OPEN_FILES: DO NOT KNOW HOW TO PROCESS' 
          WRITE (*, *) ' RUN_TYPE in the input file' 
          call mfix_exit(myPE) 

      END SELECT 

 !     endif   ! end of myPE=PE_IO if block commented out pan ... 

      RETURN  
  500 CONTINUE 
      WRITE (*, 1100) myPE,FILE_NAME  !//PAR_I/O added myPE for output
  600 CONTINUE 
      CALL SLUMBER 
      call mfix_exit(myPE) 
!
 1000 FORMAT(I1) 
 1001 FORMAT(/70('*')//' From: OPEN_FILES',/&
         ' Error: NEW run -- .RES file should NOT be in the run directory'/&
         ' Cannot open new file -- ',A,/70('*')/) 
 1002 FORMAT(/70('*')//' From: OPEN_FILES',/&
         ' Error: RESTART run -- .RES file should be in the run directory'/&
         ' Cannot open existing file -- ',A,/70('*')/) 
 1100 FORMAT(/70('*')//'(PE ',I3,'): From: OPEN_FILES',/&
         ' Error: Cannot open file -- ',A,/70('*')/) 
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
          'SEQUENTIAL', 'FORMATTED', 132, IER)
        
      RETURN
      END SUBROUTINE OPEN_PE_LOG 
