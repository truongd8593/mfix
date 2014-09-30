!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SEEK_COMMENT (LINE_MAXCOL)                             C
!  Purpose: determine if (and where) a comment character appears       C
!           in a data input line                                       C
!                                                                      C
!  Author: P.Nicoletti                                Date: 25-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: SEEK_COMMENT                                    C
!                                                                      C
!  Local variables: DIM_COMMENT, COMMENT_CHAR, L, COMMENT, L2          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      INTEGER FUNCTION SEEK_COMMENT (LINE, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input data line
      CHARACTER*(*) LINE
!
!                   maximum column of input data line to search
      INTEGER       MAXCOL
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                   the number of designated comment characters
      INTEGER, PARAMETER :: DIM_COMMENT = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                   loop indicies
      INTEGER :: L, L2
!
!                   the comment characters
      CHARACTER, DIMENSION(DIM_COMMENT) :: COMMENT_CHAR
!-----------------------------------------------
!
!     The function SEEK_COMMENT returns the index to where a comment
!     character was found in the input data line.  Equals MAXCOL + 1
!     if no-comment characters in the line
!
!
      DATA COMMENT_CHAR/'#', '!'/
!
      DO L = 1, MAXCOL
         DO L2 = 1, DIM_COMMENT
            IF (LINE(L:L) == COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            ENDIF
         END DO
      END DO
      SEEK_COMMENT = MAXCOL + 1
!
      RETURN
      END FUNCTION SEEK_COMMENT
