CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: HEADER_MAIN                                            C
C  Purpose: Write out the main selection menu and get user selection   C
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
C  Variables referenced:                                               C
C  Variables modified: SELECTION                                       C
C                                                                      C
C  Local variables: NERROR                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE HEADER_MAIN
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'post3d.inc'
C
C             NUMBER OF INCORRECT INPUT SELECTIONS MADE
      INTEGER NERROR
C
      NERROR = -1
C
10    NERROR = NERROR + 1
      IF (NERROR.GT.10) THEN
         WRITE (*,*) ' HEADER_MAIN : TOO MANY INCORRECT INPUTS'
         STOP
      END IF
      WRITE (*,*)
     &  ' *************************************************'
      WRITE (*,*)
     &  '  0   - Exit POST_MFIX'     
      WRITE (*,*)
     &  '  1   - Examine/print data' 
      WRITE (*,*)
     &  '  2   - Write .RES from data in .SPx files'
      WRITE (*,*)
     &  '  3   - Write .RES for a new grid, using old data'
      WRITE (*,*)
     &  '  4   - Calculate miscellaneous quantities'
      WRITE (*,*)
     &  '  5   - Print out variables'
      WRITE (*,*)
     &  '  6   - Call user defined subroutine USR_POST'
      WRITE (*,*)
     &  '  7   - Write a new SPx file with selected records'
      WRITE (*,*)
     &  ' *************************************************'
C
      CALL GET_SELECTION (SELECTION)
C
      RETURN
      END
