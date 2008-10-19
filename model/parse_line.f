!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG)            C
!  Purpose: Parse input line                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-JUN-97  C
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
      SUBROUTINE PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
      CHARACTER*(*)     LINE
!
!                      The part of LINE containing input 
      INTEGER          LMAX 
! 
!                      Cumulative value and sub value 
      DOUBLE PRECISION VALUE, SUB_VALUE 
! 
!                      Start and end locations for the arithmetic operation 
      INTEGER          LSTART, LEND 
! 
!                      Length of arithmetic operation string 
      INTEGER          LENGTH 
! 
!                      22 - LENGTH 
      INTEGER          LDIF 
! 
!                      Locations in SUB_STR, and LINE 
      INTEGER          LSUB, L 
! 
!                      Operator symbol (Legal values: *, /) 
      CHARACTER        OPERATION*1 
! 
!                      Substring taken from LINE 
      CHARACTER        SUB_STR*80 
! 
!                      Indicate whether currently reading rxns 
      LOGICAL          RXN_FLAG 
! 
!                      Indicate whether to do a namelist read on the line 
      LOGICAL          READ_FLAG 
  
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: START_RXN, END_RXN 
!-----------------------------------------------
!
!  Preliminary processing
!
      IF (LMAX == 0) THEN                        !Blank line -- no need to read 
         READ_FLAG = .FALSE. 
         RETURN  
      ENDIF 
!
      LSTART = INDEX(LINE,START_STR)             !Is there a string to parse? 
!
      IF (LSTART /= 0) THEN 
!
         LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),END_STR) 
         IF (LEND <= LSTART) THEN 
            WRITE (*, 1000) myPE,LINE(LSTART:LMAX) 
            CALL MFIX_EXIT(myPE)  
         ENDIF 
!
         IF (END_RXN(LINE(LSTART:LEND),LEND-LSTART)) THEN 
            IF (.NOT.RXN_FLAG) WRITE (*, 1010) myPE,LINE(1:LMAX) 
!
            IF (READING_RATE) CALL CLOSE_READING_RATE 
            IF (READING_RXN) CALL CLOSE_READING_RXN 
!
            RXN_FLAG = .FALSE. 
            READ_FLAG = .FALSE. 
            RETURN  
         ENDIF 
!
         IF (START_RXN(LINE(LSTART:LEND),LEND-LSTART)) THEN 
            RXN_FLAG = .TRUE. 
            READ_FLAG = .FALSE. 
!
            READING_RXN = .FALSE. 
            READING_RATE = .FALSE. 
            RETURN  
         ENDIF 
!
      ENDIF 
!
      IF (RXN_FLAG) THEN 
         CALL PARSE_RXN (LINE, LMAX) 
         READ_FLAG = .FALSE. 
         RETURN  
      ENDIF 
!
      LSTART = INDEX(LINE,START_STR)             !Arithmetic processing ? 
      IF (LSTART /= 0) CALL PARSE_ARITH (LINE, LMAX) 
      READ_FLAG = .TRUE. 
!
      RETURN  

 1000 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_LINE',/&
         ' Message: No ending ) found in the input line: ',/9X,A,/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_LINE',/&
         ' Message: END keyword before a start keyword in line: ',/9X,A,/1X,70(&
         '*')/) 
      END SUBROUTINE PARSE_LINE 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: START_RXN(LINE, LMAX)                                  C
!  Purpose: Check for the start of rxn block                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-JUN-97  C
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
      LOGICAL FUNCTION START_RXN (LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!C 
!                      The part of LINE containing input 
      INTEGER LMAX 
!                      Input line with arithmetic operations.  Out put
!                      line with completed arithmetic statements.
!
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------

      IF (INDEX(LINE(1:LMAX),RXN_BLK) == 0) THEN 
         START_RXN = .FALSE. 
      ELSE 
         START_RXN = .TRUE. 
      ENDIF 
!
      RETURN  
      END FUNCTION START_RXN 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: END_RXN(LINE, LMAX)                                    C
!  Purpose: Check for the end of rxn block                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-JUN-97  C
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
      LOGICAL FUNCTION END_RXN (LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!C 
!                      The part of LINE containing input 
      INTEGER LMAX 
!                      Input line with arithmetic operations.  Out put
!                      line with completed arithmetic statements.
!
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
      IF (INDEX(LINE(1:LMAX),END_BLK) == 0) THEN 
         END_RXN = .FALSE. 
      ELSE 
         END_RXN = .TRUE. 
      ENDIF 
!
      RETURN  
      END FUNCTION END_RXN 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_ARITH(LINE, LMAX)                                C
!  Purpose: Complete arithmetic operations and expand the line         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-AUG-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
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
      SUBROUTINE PARSE_ARITH(LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parse 
      USE compar     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!C 
!                      The part of LINE containing input 
      INTEGER LMAX 
!                      Input line with arithmetic operations.  Out put
!                      line with completed arithmetic statements.
!
      CHARACTER LINE*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      Value of pi 
      DOUBLE PRECISION PI 
! 
!                      Cumulative value and sub value 
      DOUBLE PRECISION VALUE, SUB_VALUE 
! 
!                      Start and end locations for the arithmetic operation 
      INTEGER          LSTART, LEND 
! 
!                      Length of arithmetic operation string 
      INTEGER          LENGTH 
! 
!                      22 - LENGTH 
      INTEGER          LDIF 
! 
!                      Locations in SUB_STR, and LINE 
      INTEGER          LSUB, L 
! 
!                      Operator symbol (Legal values: *, /) 
      CHARACTER        OPERATION*1 
! 
!                      Substring taken from LINE 
      CHARACTER        SUB_STR*80 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      INTEGER , EXTERNAL :: SEEK_END 
!-----------------------------------------------
!
!
      PI = 4.0D0*ATAN(ONE) 
!
!  Search for arithmetic operation
!
   10 CONTINUE 
      LMAX = SEEK_END(LINE,LEN(LINE)) 
!
      LSTART = INDEX(LINE,START_STR) 
!
      IF (LSTART == 0) RETURN  
!
      LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),END_STR) 
      IF (LEND <= LSTART) THEN 
         WRITE (*, 1000) myPE,LINE(LSTART:LMAX) 
         CALL MFIX_EXIT(myPE)  
      ENDIF 
!
!    Do the arithmetic
!
      VALUE = ONE 
      OPERATION = '*' 
      LSUB = 1 
      DO L = LSTART + 2, LEND 
         IF (LINE(L:L)=='*' .OR. LINE(L:L)=='/' .OR. LINE(L:L)==END_STR) THEN 
            IF (LSUB == 1) THEN 
               WRITE (*, 1015) myPE,LINE(LSTART:LEND) 
               CALL MFIX_EXIT(myPE)  
            ENDIF 
            IF (SUB_STR(1:LSUB-1) == 'PI') THEN 
               SUB_VALUE = PI 
            ELSE 
               READ (SUB_STR(1:LSUB-1), *, ERR=900) SUB_VALUE 
            ENDIF 
            IF (OPERATION == '*') THEN 
               VALUE = VALUE*SUB_VALUE 
            ELSE IF (OPERATION == '/') THEN 
               VALUE = VALUE/SUB_VALUE 
            ENDIF 
            LSUB = 1 
            OPERATION = LINE(L:L) 
         ELSE IF (LINE(L:L) == ' ') THEN 
         ELSE 
            SUB_STR(LSUB:LSUB) = LINE(L:L) 
            LSUB = LSUB + 1 
         ENDIF 
      END DO 
      LENGTH = LEND - LSTART + 1 
      IF (LENGTH > 22) THEN 
         DO L = LSTART + 22, LEND 
            LINE(L:L) = ' ' 
         END DO 
      ELSE IF (LENGTH < 22) THEN 
         LMAX = SEEK_END(LINE,LEN(LINE)) 
         LDIF = 22 - LENGTH 
         IF (LMAX + LDIF > LEN(LINE)) THEN 
            WRITE (*, 1020) myPE,LINE(1:80) 
            CALL MFIX_EXIT(myPE)  
         ENDIF 
         DO L = LMAX, LEND + 1, -1 
            LINE(L+LDIF:L+LDIF) = LINE(L:L) 
         END DO 
      ENDIF 
!
!  Transfer the value to LINE
!
      WRITE (SUB_STR, '(G22.15)') VALUE 
      L = LSTART 
      DO LSUB = 1, 22 
         LINE(L:L) = SUB_STR(LSUB:LSUB) 
         L = L + 1 
      END DO 
      GO TO 10 
!
  900 CONTINUE 
      WRITE (*, 1010) myPE, SUB_STR(1:LSUB-1) 
      CALL MFIX_EXIT(myPE)  
 1000 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: No ending ) found in the input line: ',/9X,A,/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Error reading the input string: ',/9X,A,/1X,70('*')/) 
 1015 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Invalid operator in the input string: ',/9X,A,/1X,70('*')/) 
 1020 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Too many arithmetic operations in the line: ',/1X,A,/1X,70(&
         '*')/) 
      END SUBROUTINE PARSE_ARITH 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!//PAR_I/O added myPE stamp in output
