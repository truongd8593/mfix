CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: READ_SPX0(READ_SPX)                                    C
C  Purpose: Read the initial records (REAL)                            C
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
      SUBROUTINE READ_SPX0(READ_SPX)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'post3d.inc'
C
      logical read_spx(*)
C                loop counter
      INTEGER    LC,l1
      
C
C                file version ID
      CHARACTER  VERSION*512
C
      DO 100 LC = 1,N_SPX
        IF (READ_SPX(LC) .AND. SPX_OPEN(LC)) THEN
           READ (UNIT_SPX+LC,REC=1) VERSION
           READ(VERSION(6:512),*)VERSION_NUMBER
           READ (UNIT_SPX+LC,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,
     &          ID_HOUR,ID_MINUTE,ID_SECOND
C
C  The first field contains the pointer to the next record.
C  The second field contains the number of records written each time step
C
           READ (UNIT_SPX+LC,REC=3) LAST_REC(LC),NUM_REC(LC)
        ENDIF
100   CONTINUE
      RETURN
      END
