!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: bound_x(Array, IJKMAX2, IER)                           C
!
!  Purpose: bound the values of x array                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-SEP-98  C
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
      SUBROUTINE BOUND_X(ARRAY, IJKMAX2, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      cell index
      INTEGER          IJK
!
!                      Maximum dimension
      INTEGER          IJKMAX2
!
!                      Array
      DOUBLE PRECISION Array(DIMENSION_3)
!-----------------------------------------------
!
      IF (IJKMAX2 > 0) THEN 
!//? would it better to replace ARRAY(:) with ARRAY(ijkstart3:ijkend3)??      
!// 1207 rearrange array indices to remove IJKMAX2 dependency      
!         ARRAY(:IJKMAX2) = MIN(1.,MAX(0.,ARRAY(:IJKMAX2))) 
         ARRAY(:) = MIN(1.,MAX(0.,ARRAY(:))) 	 
      ENDIF 
      RETURN  
      END SUBROUTINE BOUND_X 
