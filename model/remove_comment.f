!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REMOVE_COMMENT (LINE, LSTART, MAXCOL)                  C
!  Purpose: Remove comments                                            C
!                                                                      C
!  Author: P.Nicoletti                                Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: DIM_COMMENT, COMMENT_CHAR, L, COMMENT, L2          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE REMOVE_COMMENT(LINE, LSTART, MAXCOL) 
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
!                   start of comments
      INTEGER       LSTART
!
!                   maximum column of input data line to search
      INTEGER       MAXCOL
!
! local variables
!
!
!                   loop index
      INTEGER       L 
!-----------------------------------------------
!
!
      DO L = LSTART, MAXCOL 
         LINE(L:L) = ' ' 
      END DO 
      RETURN  
      END SUBROUTINE REMOVE_COMMENT 
