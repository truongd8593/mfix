!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ERROR_ROUTINE (CALL_ROUTINE,MESSAGE,ACTION_CODE,       C
!                              MESSAGE_CODE)                           C
!  Purpose:  Assist in printing error messages during input processing C
!                                                                      C
!  Author: P.NICOLETTI                                Date: 25-NOV-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 29-JAN-92  C
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
!  Local variables: ABORT_CONT                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ERROR_ROUTINE(CALL_ROUTINE,MESSAGE,ACTION_CODE,MESSAGE_CODE) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE funits 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER ACTION_CODE, MESSAGE_CODE 
      CHARACTER CALL_ROUTINE*(*), MESSAGE*(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER :: ABORT_CONT*10 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
! SET UP THE ABORT / CONTINUE MESSAGE
!
      IF (ACTION_CODE == 0) THEN 
         ABORT_CONT = 'continued' 
      ELSE 
         ABORT_CONT = 'aborted' 
      ENDIF 
!
! WRITE OUT HEADER INFO , UNLESS MESSAGE_CODE = 3
!
      IF (MESSAGE_CODE /= 3) WRITE (UNIT_LOG, 1000) CALL_ROUTINE, MESSAGE 
!
! WRITE OUT TRAILER INFO, UNLESS MESSAGE_CODE = 2
!
      IF (MESSAGE_CODE /= 2) WRITE (UNIT_LOG, 1100) ABORT_CONT 
!
      IF (ACTION_CODE == 0) THEN 
         RETURN  
      ELSE 
         STOP  
      ENDIF 
!
 1000 FORMAT(1X,70('*'),/,/,1X,'From : ',A,/,1X,'Message : ',A) 
 1100 FORMAT(1X,'Program execution ',A,/,/,1X,70('*')) 
!
      END SUBROUTINE ERROR_ROUTINE 
