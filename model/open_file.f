!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_FILE                                              C
!  Purpose: open a file                                                C
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
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OPEN_FILE(RUN_NAME, NB, IUNIT, EXT, FILE_NAME, OPEN_STAT, &
         OPEN_ACCESS, OPEN_FORM, IRECL, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                     Error index: 0 - no error, 1 could not open file
      INTEGER         IER
!
!                     run_name (without extension)
      CHARACTER*(*)   RUN_NAME
!
!                     extension
      CHARACTER*(*)   EXT
!
!                     run_name + extension
      CHARACTER*(*)   FILE_NAME
!
!                     open status ('NEW' or 'OLD')
      CHARACTER*(*)   OPEN_STAT
!
!                     open access ('SEQUENTIAL' or 'DIRECT')
      CHARACTER*(*)   OPEN_ACCESS
!
!                     open form ('FORMATTED' or 'UNFORMATTED')
      CHARACTER*(*)   OPEN_FORM
!
!                     index to first blank character in run_name
      INTEGER         NB
!
!                     unit number to open
      INTEGER         IUNIT
!
!                     record length
      INTEGER         IRECL
!-----------------------------------------------
!

      FILE_NAME = ' ' 
      FILE_NAME(1:NB-1) = RUN_NAME(1:NB-1) 
      FILE_NAME(NB:NB+3) = EXT(1:4) 
!
      IF (OPEN_ACCESS == 'DIRECT') THEN 
         OPEN(UNIT=IUNIT, FILE=FILE_NAME(1:NB+3), STATUS=OPEN_STAT, RECL=IRECL, ACCESS=&
            OPEN_ACCESS, FORM=OPEN_FORM, ERR=100) 
      ELSE 
         OPEN(UNIT=IUNIT, FILE=FILE_NAME(1:NB+3), STATUS=OPEN_STAT, ACCESS=OPEN_ACCESS&
            , FORM=OPEN_FORM, ERR=100) 
      ENDIF 
      IER = 0 
      RETURN  
!
  100 CONTINUE 
      IER = 1 
      RETURN  
      END SUBROUTINE OPEN_FILE 
