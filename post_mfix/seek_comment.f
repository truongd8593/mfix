CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SEEK_COMMENT (LINE_MAXCOL)                             C
C  Purpose: determine if (and where) a comment character appears       C
C           in a data input line                                       C
C                                                                      C
C  Author: P.Nicoletti                                Date: 25-NOV-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: SEEK_COMMENT                                    C
C                                                                      C
C  Local variables: DIM_COMMENT, COMMENT_CHAR, L, COMMENT, L2          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      INTEGER FUNCTION SEEK_COMMENT (LINE,MAXCOL)
C
      IMPLICIT NONE
C
C     The function SEEK_COMMENT returns the index to where a comment 
C     character was found in the input data line.  Equals MAXCOL + 1 
C     if no-comment characters in the line
C
C passed arguments
C
C                   input data line
      CHARACTER*(*) LINE
C
C                   maximum column of input data line to search
      INTEGER       MAXCOL
C
C local variables
C
C                   the number of designated comment characters
      INTEGER       DIM_COMMENT
      PARAMETER     (DIM_COMMENT = 2)
C
C                   the comment characters
      CHARACTER     COMMENT_CHAR(DIM_COMMENT)*1
C
C                   loop indicies
      INTEGER       L , L2
C
      DATA          COMMENT_CHAR / '#' , '!' /
C
      DO 200 L = 1,MAXCOL
         DO 100 L2 = 1,DIM_COMMENT
            IF (LINE(L:L).EQ.COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            END IF
100      CONTINUE
200   CONTINUE
C
      SEEK_COMMENT = MAXCOL + 1
C
      RETURN
      END
