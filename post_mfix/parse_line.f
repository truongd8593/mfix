CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG)            C
C  Purpose: Parse input line                                           C
C                                                                      C
C  Author: M. Syamlal                                 Date: 27-JUN-97  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'parse.inc'
C
C  Functions
C
      INTEGER SEEK_END
C
C  Passed variables
C
C                      Input line with arithmetic operations.  Out put
C                      line with completed arithmetic statements.
      CHARACTER*(*)     LINE
C
C  Local variables
C
C                      The part of LINE containing input
      INTEGER          LMAX
C
C                      Cumulative value and sub value
      DOUBLE PRECISION VALUE, SUB_VALUE
C
C                      Start and end locations for the arithmetic operation
      INTEGER          LSTART, LEND 
C
C                      Length of arithmetic operation string
      INTEGER          LENGTH
C
C                      22 - LENGTH
      INTEGER          LDIF
C
C                      Locations in SUB_STR, and LINE
      INTEGER          LSUB, L
C
C                      Operator symbol (Legal values: *, /)
      CHARACTER        OPERATION*1
C
C                      Substring taken from LINE
      CHARACTER        SUB_STR*80
C
C                      Indicate whether currently reading rxns
      LOGICAL          RXN_FLAG
C
C                      Indicate whether to do a namelist read on the line
      LOGICAL          READ_FLAG

C
C     Functions

      LOGICAL START_RXN, END_RXN
C
C  Preliminary processing
C
      IF(LMAX .EQ. 0) THEN                   !Blank line -- no need to read
        READ_FLAG = .FALSE.
        RETURN
      ENDIF

      LSTART = INDEX(LINE, START_STR)        !Is there a string to parse?      

      IF( LSTART .NE. 0)THEN

        LEND = LSTART - 1+ INDEX(LINE(LSTART:LMAX), END_STR)
        IF(LEND .LE. LSTART)THEN
          WRITE(*, 1000)LINE(LSTART:LMAX)
          STOP
        ENDIF

        IF( END_RXN( LINE(LSTART:LEND), (LEND-LSTART) ) )THEN
          IF(.NOT.RXN_FLAG)THEN
            WRITE(*,1010)LINE(1:LMAX)
          ENDIF

          IF(READING_RATE)CALL CLOSE_READING_RATE
          IF(READING_RXN)CALL CLOSE_READING_RXN

          RXN_FLAG = .FALSE.
          READ_FLAG = .FALSE.
          RETURN
        ENDIF
      
        IF( START_RXN( LINE(LSTART:LEND), (LEND-LSTART) ) )THEN
          RXN_FLAG = .TRUE.
          READ_FLAG = .FALSE.

          READING_RXN  = .FALSE.
          READING_RATE = .FALSE.
          RETURN
        ENDIF

      ENDIF

      IF(RXN_FLAG)THEN
        CALL PARSE_RXN(LINE, LMAX)
        READ_FLAG = .FALSE.
        RETURN
      ENDIF

      LSTART = INDEX(LINE, START_STR)         !Arithmetic processing ?
      IF( LSTART .NE. 0)CALL PARSE_ARITH(LINE, LMAX)
      READ_FLAG = .TRUE.
        
      RETURN
1000  FORMAT(/1X,70('*')//' From: PARSE_LINE',
     &/' Message: No ending ) found in the input line: ',/9X,A
     & ,/1X, 70('*')/)
1010  FORMAT(/1X,70('*')//' From: PARSE_LINE',
     &/' Message: END keyword before a start keyword in line: ',
     &/9X, A, /1X, 70('*')/)
      END

CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: START_RXN(LINE, LMAX)                                  C
C  Purpose: Check for the start of rxn block                           C
C                                                                      C
C  Author: M. Syamlal                                 Date: 27-JUN-97  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      LOGICAL FUNCTION START_RXN(LINE, LMAX)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'parse.inc'
C
C                      Input line with arithmetic operations.  Out put
C                      line with completed arithmetic statements.
      CHARACTER*(*)     LINE
C
C                      The part of LINE containing input
      INTEGER          LMAX

      IF(INDEX(LINE(1:LMAX), RXN_BLK) .EQ. 0)THEN
        START_RXN = .FALSE.
      ELSE
        START_RXN = .TRUE.
      ENDIF

      RETURN
      END


CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: END_RXN(LINE, LMAX)                                    C
C  Purpose: Check for the end of rxn block                             C
C                                                                      C
C  Author: M. Syamlal                                 Date: 27-JUN-97  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      LOGICAL FUNCTION END_RXN(LINE, LMAX)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'parse.inc'
C
C                      Input line with arithmetic operations.  Out put
C                      line with completed arithmetic statements.
      CHARACTER*(*)     LINE
C
C                      The part of LINE containing input
      INTEGER          LMAX

      IF(INDEX(LINE(1:LMAX), END_BLK) .EQ. 0)THEN
        END_RXN = .FALSE.
      ELSE
        END_RXN = .TRUE.
      ENDIF

      RETURN
      END

CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: PARSE_ARITH(LINE, LMAX)                                C
C  Purpose: Complete arithmetic operations and expand the line         C
C                                                                      C
C  Author: M. Syamlal                                 Date: 10-AUG-92  C
C  Reviewer: W. Rogers                                Date: 11-DEC-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE PARSE_ARITH(LINE, LMAX)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'parse.inc'
C
C  Functions
C
      INTEGER SEEK_END
C
C  Passed variables
C
C                      Input line with arithmetic operations.  Out put
C                      line with completed arithmetic statements.
      CHARACTER*(*)     LINE
C
C  Local variables
CC
C                      The part of LINE containing input
      INTEGER          LMAX
C
C                      Value of pi
      DOUBLE PRECISION PI
C
C                      Cumulative value and sub value
      DOUBLE PRECISION VALUE, SUB_VALUE
C
C                      Start and end locations for the arithmetic operation
      INTEGER          LSTART, LEND 
C
C                      Length of arithmetic operation string
      INTEGER          LENGTH
C
C                      22 - LENGTH
      INTEGER          LDIF
C
C                      Locations in SUB_STR, and LINE
      INTEGER          LSUB, L
C
C                      Operator symbol (Legal values: *, /)
      CHARACTER        OPERATION*1
C
C                      Substring taken from LINE
      CHARACTER        SUB_STR*80
C
      Pi = 4.0D0 * ATAN(ONE)
C
C  Search for arithmetic operation
C
10    LMAX = SEEK_END(LINE, LEN(LINE))
      
      LSTART = INDEX(LINE, START_STR)         

      IF( LSTART .EQ. 0)RETURN

      LEND = LSTART - 1+ INDEX(LINE(LSTART:LMAX), END_STR)
      IF(LEND .LE. LSTART)THEN
        WRITE(*, 1000)LINE(LSTART:LMAX)
        STOP
      ENDIF
C
C    Do the arithmetic
C
      VALUE = ONE
      OPERATION = '*'
      LSUB = 1
      DO 100 L = LSTART+2, LEND
        IF(LINE(L:L) .EQ. '*' .OR. LINE(L:L) .EQ. '/'
     &     .OR. LINE(L:L) .EQ. END_STR)THEN
          IF(LSUB .EQ. 1)THEN
            WRITE(*,1015)LINE(LSTART:LEND)
            STOP
          ENDIF
          IF(SUB_STR(1:LSUB-1) .EQ. 'PI')THEN
            SUB_VALUE = Pi
          ELSE
            READ(SUB_STR(1:LSUB-1),*, ERR=900)SUB_VALUE
          ENDIF
          IF(OPERATION .EQ. '*')THEN
            VALUE = VALUE * SUB_VALUE
          ELSEIF(OPERATION .EQ. '/')THEN
            VALUE = VALUE / SUB_VALUE
          ENDIF
          LSUB = 1
          OPERATION = LINE(L:L)
        ELSEIF(LINE(L:L) .EQ. ' ')THEN
        ELSE
          SUB_STR(LSUB:LSUB) = LINE(L:L)
          LSUB = LSUB + 1
        ENDIF
100   CONTINUE
C
C  Make space in LINE for the value
C
      LENGTH = LEND-LSTART+1
      IF(LENGTH .GT. 22)THEN
        DO 120 L = (LSTART + 22), LEND
          LINE(L:L) = ' '
120     CONTINUE
      ELSEIF(LENGTH .LT. 22)THEN
        LMAX = SEEK_END(LINE, LEN(LINE))
        LDIF = 22 - LENGTH
        IF( (LMAX+LDIF) .GT. LEN(LINE))THEN
          WRITE(*, 1020)LINE(1:80)
          STOP
        ENDIF
        DO 140 L = LMAX, (LEND+1), -1
          LINE(L+LDIF:L+LDIF) = LINE(L:L)
140     CONTINUE
      ENDIF
C
C  Transfer the value to LINE
C
      WRITE(SUB_STR,'(G22.15)') VALUE
      L = LSTART
      DO 200 LSUB = 1, 22
        LINE(L:L) = SUB_STR(LSUB:LSUB)
        L = L + 1
200   CONTINUE
      GOTO 10

900   WRITE(*, 1010)SUB_STR(1:LSUB-1)
      STOP
1000  FORMAT(/1X,70('*')//' From: PARSE_ARITH',
     &/' Message: No ending ) found in the input line: ',/9X,A
     & ,/1X, 70('*')/)
1010  FORMAT(/1X,70('*')//' From: PARSE_ARITH',
     &/' Message: Error reading the input string: ',/9X,A
     & ,/1X, 70('*')/)
1015  FORMAT(/1X,70('*')//' From: PARSE_ARITH',
     &/' Message: Invalid operator in the input string: ',/9X,A
     & ,/1X, 70('*')/)
1020  FORMAT(/1X,70('*')//' From: PARSE_ARITH',
     &/' Message: Too many arithmetic operations in the line: ',/1X,A
     & ,/1X, 70('*')/)
      END
