CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: REMOVE_COMMENT (LINE, LSTART, MAXCOL)                  C
C  Purpose: Remove comments                                            C
C                                                                      C
C  Author: P.Nicoletti                                Date:            C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables: DIM_COMMENT, COMMENT_CHAR, L, COMMENT, L2          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE REMOVE_COMMENT (LINE, LSTART, MAXCOL)
C
      IMPLICIT NONE
C
C
C passed arguments
C
C                   input data line
      CHARACTER*(*) LINE
C
C                   start of comments
      INTEGER       LSTART
C
C                   maximum column of input data line to search
      INTEGER       MAXCOL
C
C local variables
C
C
C                   loop index
      INTEGER       L 
C
      DO 200 L = LSTART, MAXCOL
         LINE(L:L) = ' '
200   CONTINUE
C
C
      RETURN
      END
