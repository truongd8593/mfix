!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SEEK_END (LINE, MAXCOL)                                C
!  Purpose: determine where trailing blanks begin in a line            C
!                                                                      C
!  Author: P.Nicoletti, M. Syamlal                    Date: 7-AUG-92   C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: SEEK_END                                        C
!                                                                      C
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      INTEGER FUNCTION SEEK_END (LINE, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   maximum column of input data line to search
      INTEGER MAXCOL
!
!                   input data line
      CHARACTER LINE*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L
!-----------------------------------------------
!
!     The function SEEK_END returns the index to where the last
!     character was found in the input data line.  Equals MAXCOL
!     if no trailing blank characters in the line
!
!
      SEEK_END = 0
      DO L = 1, MAXCOL
         IF (LINE(L:L) /= ' ') SEEK_END = L
      END DO
      RETURN
      END FUNCTION SEEK_END
