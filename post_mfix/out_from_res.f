CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: OUT_FROM_RES                                           C
C  Purpose: OUTARR type output from the RES file                       C
C                                                                      C
C  Author: P. Nicoletti                               Date: 20-MAR-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: TIME                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: FILE_NAME                                          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE OUT_FROM_RES(FILE_NAME)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'xforms.inc'
C
      CHARACTER FILE_NAME*(*)
C
      IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(FILE_NAME)
      OPEN (UNIT=UNIT_OUT,FILE=FILE_NAME,STATUS='UNKNOWN')
      CALL READ_RES0
      CALL READ_RES1
      WRITE (*,*) ' time in RES file = ' , TIME
      CALL WRITE_OUT1
      CLOSE (UNIT=UNIT_OUT)
C
      RETURN
      END
