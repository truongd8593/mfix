CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: WRITE_SPX0(L)                                          C
C  Purpose: write out the initial restart records (REAL)               C
C                                                                      C
C  Author: P. Nicoletti                               Date: 13-DEC-91  C
C  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR  C
C                        ID_MINUTE, ID_SECOND                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: LC, VERSION                                        C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE WRITE_SPX0(L)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
C
C                SPx file number
      INTEGER    L
C
C                file version ID
      CHARACTER  VERSION*512
C
      VERSION = 'SPx = 02.00'
         WRITE (VERSION(3:3),1000) L
         WRITE (UNIT_SPX+L,REC=1) VERSION
         WRITE (UNIT_SPX+L,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,
     &          ID_HOUR,ID_MINUTE,ID_SECOND
C
C  The first field contains the pointer to the next record.
C  The second field contains the number of records written each time step
C  (The 4 and -1 will be overwritten in WRITE_SPX1)
C
         WRITE (UNIT_SPX+L,REC=3) 4, -1
         CALL FLUSH(UNIT_SPX+L)
1000  FORMAT(I1)
      RETURN
      END
