!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_SELECTION (SELECTION)                              C
!  Purpose: GET A SELECTION INTEGER FROM THE USER                      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 05-FEB-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_SELECTION (SELECTION)
!
      IMPLICIT NONE
!
      INTEGER SELECTION
!
      WRITE (*,*) ' '
      WRITE (*,'(A,$)') ' Enter menu selection > '
!
      READ (*,'(I1)',ERR=10) SELECTION
      RETURN
10    SELECTION = -1000000
      RETURN
      END
