!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY (ARRAY, MESSAGE)                             C
!  Purpose: print out a 3D array to standard output                    C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: KMAX2                                         C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: IJK, K                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_ARRAY(ARRAY, MESSAGE) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE indices
      USE funits 
      USE tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       array to print out
      DOUBLE PRECISION  ARRAY(*)
!
!                       message to print out
      CHARACTER*(*)     MESSAGE
!
! local variables
!
!                       pointer into array (points to start of a k-plane)
      INTEGER           IJK
!
!                       loop counter
      INTEGER           K
!
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      call lock_tmp_array
      call convert_to_io_dp(array,array1,ijkmax2)

      DO K = 1, KMAX2 
         IJK = FUNIJK_IO(1,1,K) 
         WRITE (UNIT_OUT, 1100) MESSAGE, K 
         CALL OUT_ARRAY_K (ARRAY1(IJK)) 
      END DO 
 1100 FORMAT(/,1X,A,' at K = ',I4,/) 

      call unlock_tmp_array
!
      RETURN  
      END SUBROUTINE OUT_ARRAY 
