!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_error(Name, Line, L)                             C                     C
!  Purpose: Write an error message                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_ERROR(NAME, LINE, LMAX) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE funits 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Subroutine name
      CHARACTER*(*)    Name
!
!                      Message
      CHARACTER*(*)    LINE(*)
!
!                      Dimension of message array
      INTEGER          LMAX
!
!                      Index
      INTEGER          L
!
!-----------------------------------------------
!

      CALL START_LOG 
      WRITE (UNIT_LOG, 1000) NAME 
      DO L = 1, LMAX 
         WRITE (UNIT_LOG, 1010) LINE(L) 
      END DO 
      WRITE (UNIT_LOG, 1020) 
      CALL END_LOG 
      RETURN  
 1000 FORMAT(1X,70('*'),/,/,1X,'From : ',A) 
 1010 FORMAT(1X,A) 
 1020 FORMAT(/,/,1X,70('*')) 
      END SUBROUTINE WRITE_ERROR 
