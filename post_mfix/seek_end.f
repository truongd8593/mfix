CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SEEK_END (LINE, MAXCOL)                                C
C  Purpose: determine where trailing blanks begin in a line            C
C                                                                      C
C  Author: P.Nicoletti, M. Syamlal                    Date: 7-AUG-92   C
C  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: SEEK_END                                        C
C                                                                      C
C  Local variables: L                                                  C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      INTEGER FUNCTION SEEK_END (LINE,MAXCOL)
C
      IMPLICIT NONE
C
C     The function SEEK_END returns the index to where the last 
C     character was found in the input data line.  Equals MAXCOL
C     if no trailing blank characters in the line
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
C                   loop indicies
      INTEGER       L
C
      SEEK_END = 0
      DO 200 L = 1, MAXCOL
        IF (LINE(L:L).NE.' ') SEEK_END = L
200   CONTINUE
C
      RETURN
      END
