CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_FROM_SPX (AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)  C
C                                                                      C
C  Purpose: OUTARR type output for a specified time from the SPX files C
C           USES N_SPX                                                 C
C                                                                      C
C  Author: P. Nicoletti                               Date: 15-MAR-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: TIME_FOR_OUT, FILE_NAME, PRINTED_MESS, MORE_DATA   C
C                   ERROR, L                                           C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER,TIME_REAL)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'xforms.inc'
C
      REAL      TIME_FOUND
      LOGICAL   AT_EOF(*), READ_SPX(*)
      INTEGER   REC_POINTER(*)
      INTEGER   NSTEP_1
      REAL      TIME_REAL(*)
      LOGICAL   PRINTED_MESS(N_SPX)
      LOGICAL   ERROR
C
      INTEGER L
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A,$)') 'Enter time to retrieve from Spx files > '
         READ  (*,*) TIME_FOR_OUT
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF
C
      OPEN (UNIT=UNIT_OUT,FILE=TEMP_FILE,STATUS='UNKNOWN')
      CALL READ_RES0
C
      DO 10 L = 1,N_SPX
         READ_SPX(L) = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L) = .FALSE.
         PRINTED_MESS(L) = .FALSE.
10    CONTINUE
      DO 20 L = 1,N_SPX
         IF(.NOT.SPX_OPEN(L))GOTO 20
         READ_SPX(L) = .TRUE.
         IF(L .GT. 1)READ_SPX(L-1)=.FALSE.
         CALL SEEK_TIME(READ_SPX, TIME_FOR_OUT, REC_POINTER, TIME_FOUND)
         IF(TIME_FOUND .GE. ZERO) THEN
           WRITE (*,*) ' Found time in SPX file ', L
           PRINTED_MESS(L) = .TRUE.
         ENDIF
20    CONTINUE
      DO 30 L = 1,N_SPX
        IF(SPX_OPEN(L)) READ_SPX(L) = .TRUE.
30    CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL,NSTEP_1)
C
      ERROR = .FALSE.
      DO L = 1,N_SPX
         IF (.NOT.PRINTED_MESS(L)) THEN
            WRITE(UNIT_OUT,*) ' Did not find time in SPX file ' , L
            ERROR = .TRUE.
         END IF
      END DO
C
      TIME  = TIME_FOR_OUT
      CALL WRITE_OUT1
      CLOSE (UNIT=UNIT_OUT)
C
      IF (ERROR) THEN
        WRITE (*,*)' WARNING: Some variables were not found !'
      END IF
C
      RETURN  
      END
