!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEPRECIATED_OR_UNKNOWN(LINE_NO, INPUT)

      use compar, only: myPE, PE_IO
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LINE_NO
      CHARACTER(len=*), INTENT(IN) :: INPUT

      CHARACTER(len=256) :: STRING

! Old keyword for solids density :: Replaced by RO_s0
      DOUBLE PRECISION :: RO_s

      NAMELIST / DEPRECIATED / RO_s

      STRING = '&DEPRECIATED '//trim(adjustl(INPUT))//'/'
      READ(STRING,NML=DEPRECIATED,ERR=999, END=999)

      IF(myPE == 0) WRITE(*,2000) trim(iVAL(LINE_NO)), trim(INPUT)
      CALL MFIX_EXIT(myPE)

 2000 FORMAT(//1X,70('*')/' From DEPRECIATED_OR_UNKNOWN:',/1x,         &
         'Error 2000: A keyword pair on line ',A,' of the mfix.dat ',  &
         'file was',/' identified as being depreciated.',//3x,A,//1x,  &
         'Please see the user documentation and update the mfix.dat ', &
         'file.',/1X,70('*')//)

  999 if(myPE == 0)WRITE (*, 2001) trim(iVAL(LINE_NO)), trim(INPUT)
      CALL MFIX_EXIT(myPE)


 2001 FORMAT(//1X,70('*')/' From: DEPRECIATED_OR_UNKNOWN',/1x,         &
         'Error 2001: Unable to process line ',A,' of the mfix.dat ',  &
         'file.',2/3x,A,2/1x,'Possible causes are',/3x,'*  Incorrect', &
         ' or illegal keyword format',/3x'*  Unknown or mistyped name',&
         /3x,'*  The dimensioned item is too small (array overflow).', &
         2/1x,'Please see the user documentation and update the ',     &
         'mfix.dat file. ',/1X,70('*')//)

      END SUBROUTINE DEPRECIATED_OR_UNKNOWN
