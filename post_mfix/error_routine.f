CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: ERROR_ROUTINE (CALL_ROUTINE,MESSAGE,ACTION_CODE,       C
C                              MESSAGE_CODE)                           C
C  Purpose:  Assist in printing error messages during input processing C
C                                                                      C
C  Author: P.NICOLETTI                                Date: 25-NOV-91  C
C  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 29-JAN-92  C
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
C  Local variables: ABORT_CONT                                         C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ERROR_ROUTINE(CALL_ROUTINE,MESSAGE,ACTION_CODE,
     &                         MESSAGE_CODE)
C
      IMPLICIT NONE
      INCLUDE 'funits.inc'
C
      CHARACTER*(*)   CALL_ROUTINE , MESSAGE
      CHARACTER       ABORT_CONT*10
      INTEGER         ACTION_CODE , MESSAGE_CODE
C
C SET UP THE ABORT / CONTINUE MESSAGE
C
      IF (ACTION_CODE.EQ.0) THEN
         ABORT_CONT = 'continued'
      ELSE
         ABORT_CONT = 'aborted'
      END IF
C
C WRITE OUT HEADER INFO , UNLESS MESSAGE_CODE = 3
C
      IF (MESSAGE_CODE.NE.3) 
     &  WRITE(UNIT_LOG,1000) CALL_ROUTINE , MESSAGE
C
C WRITE OUT TRAILER INFO, UNLESS MESSAGE_CODE = 2
C
      IF (MESSAGE_CODE.NE.2) WRITE(UNIT_LOG,1100) ABORT_CONT
C
      IF (ACTION_CODE.EQ.0) THEN
         RETURN
      ELSE
         STOP
      END IF
C
1000  FORMAT(1X,70('*'),/,/,
     &       1X,'From : ' , A , / ,
     &       1X,'Message : ' , A)
1100  FORMAT(1X,'Program execution ' , A ,/,/,1x,70('*'))
C
      END
