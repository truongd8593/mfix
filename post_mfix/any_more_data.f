CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: ANY_MORE_DATA (READ_SPX,AT_EOF)                        C
C  Purpose: Determine whether all data has been read from the          C
C           requested/needed SPX files ... uses N_SPX                  C
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
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: L                                                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      LOGICAL FUNCTION ANY_MORE_DATA(READ_SPX,AT_EOF)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C     passed arguments
C
      LOGICAL READ_SPX(*) , AT_EOF(*)
C
C     local variables
C
      INTEGER L
C
      ANY_MORE_DATA = .FALSE.
      DO 100 L = 1,N_SPX
         IF (READ_SPX(L) .AND. .NOT.AT_EOF(L)) ANY_MORE_DATA = .TRUE.
100   CONTINUE
C
      RETURN
      END
