CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: GET_FILE_NAME(FILE_NAME)                               C
C  Purpose: GET A FILENAME FOR OUTPUT DATA                             C
C                                                                      C
C  Author: P. Nicoletti                               Date: 27-FEB-92  C
C  Reviewer:                                                           C
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
      SUBROUTINE GET_FILE_NAME (FILE_NAME)
C
      IMPLICIT NONE
C
      CHARACTER*(*) FILE_NAME
      LOGICAL       FILE_EXIST
      CHARACTER   ANSWER
C
100   WRITE (*,'(A,$)') 'Enter Filename for this data > '
      READ  (*,'(A)') FILE_NAME
      INQUIRE (FILE=FILE_NAME,EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
         WRITE (*,*) ' '
         WRITE (*,'(A,$)')
     &       'File already exists.  Overwrite ? (Y/N) > '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') RETURN
         GOTO 100
      END IF
C
      RETURN
      END
