CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_SELECTION (SELECTION)                              C
C  Purpose: GET A SELECTION INTEGER FROM THE USER                      C
C                                                                      C
C  Author: P. Nicoletti                               Date: 05-FEB-92  C
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
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_SELECTION (SELECTION)
C
      IMPLICIT NONE
C
      INTEGER SELECTION
C
      WRITE (*,*) ' '
      WRITE (*,'(A,$)') ' Enter menu selection > '
C
      READ (*,'(I1)',ERR=10) SELECTION
      RETURN
10    SELECTION = -1000000
      RETURN
      END
