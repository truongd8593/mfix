!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_OUT3                                             C
!  Purpose: To write cpu time used by the code                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-JAN-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_OUT3(CPU) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE funits 
      USE compar     !//d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION CPU 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!
      if (myPE.eq.PE_IO) then   !//d

         WRITE (*, 1000) 
         WRITE (*, *) ' ' 
         WRITE (*, *) ' Total CPU time used = ', CPU, 'seconds' 
         WRITE (*, *) ' ' 
         WRITE (*, 1000) 
!C
         WRITE (UNIT_OUT, 1000) 
         WRITE (UNIT_OUT, *) ' ' 
         WRITE (UNIT_OUT, *) ' Total CPU time used = ', CPU, 'seconds' 
         WRITE (UNIT_OUT, *) ' ' 
         WRITE (UNIT_OUT, 1000) 

      end if                      !//
!
      CALL START_LOG 
      IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' ' 
      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Total CPU time used = ', CPU, 'seconds' 
      IF(DMP_LOG)WRITE (UNIT_LOG, *) ' ' 
      IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
      CALL END_LOG 
!
      RETURN  
 1000 FORMAT(/,1X,70('*')) 
      END SUBROUTINE WRITE_OUT3 
