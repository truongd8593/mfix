CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: PRINT_OUT                                              C
C  Purpose: Create ASCII outputs of all arrays                         C
C                                                                      C
C  Author: P. Nicoletti                               Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE PRINT_OUT
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'xforms.inc'
C
      REAL              TIME_REAL(N_SPX)
      INTEGER           REC_POINTER(N_SPX)
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
C
      IF (DO_XFORMS) THEN
         SELECTION = XCODE
         GOTO 20
      END IF
C
10    WRITE (*,*) 
     &   '  0   - Return to main menu'
      WRITE (*,*) 
     &   '  1   - Print out variables from RES file'
      WRITE (*,*) 
     &   '  2   - Print out variables from SPX files'
C
      CALL GET_SELECTION (SELECTION)
C
 20   CONTINUE
      IF(SELECTION .EQ. 0) THEN
        RETURN
      ELSEIF(SELECTION .EQ. 1) THEN
         CALL OUT_FROM_RES(TEMP_FILE)
      ELSEIF(SELECTION .EQ. 2)THEN
        CALL OUT_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
      ENDIF
      IF (DO_XFORMS) RETURN
      GOTO 10
C
      END
