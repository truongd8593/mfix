CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: LINE_TOO_BIG (LINE,LINE_LEN,MAXCOL)                    C
C  Purpose: return an error condition if input data is located past    C
C           column MAXCOL in the data input file                       C
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
C  Variables modified: LINE_TOO_BIG                                    C
C                                                                      C
C  Local variables: L                                                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      INTEGER FUNCTION LINE_TOO_BIG (LINE,LINE_LEN,MAXCOL)
C
      IMPLICIT NONE
C
C     The function LINE_TOO_BIG returns a value greater than 0 to
C     indicate an error condition (data passed column MAXCOL in LINE)
C
C passed arguments
C
C                   input data line
      CHARACTER*(*) LINE
C
C                   length of input data line
      INTEGER       LINE_LEN
C
C                   maximum column that non-blank charcaters are
C                   are in the input data line
      INTEGER       MAXCOL
C
C local variables
C
C               loop index
      INTEGER   L
C
      DO 100 L = MAXCOL+1,LINE_LEN
         IF (LINE(L:L) .NE. ' ') THEN
            LINE_TOO_BIG = L
            RETURN
         END IF
100   CONTINUE
C
      LINE_TOO_BIG = 0
      RETURN
      END
