!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Zero_array(Array, IJKMAX2, IER)                        C               !  Purpose: Zero out an array                                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-AUG-96  C
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
      SUBROUTINE ZERO_ARRAY(ARRAY, IJKMAX2, IER) 
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
!
!-----------------------------------------------
!
      IJK = 1 
      IF (IJKMAX2 > 0) THEN 
         ARRAY(:IJKMAX2) = ZERO 
         IJK = IJKMAX2 + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE ZERO_ARRAY 
