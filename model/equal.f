!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EQUAL (ARRAY1, IJK1, SIGN, ARRAY2, IJK2)               C
!  Purpose: Loop on the number of solids phases to set a variable      C
!           equal to the value or negative value of another variable   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX                                          C
!  Variables modified: M                                               C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE EQUAL(ARRAY1, IJK1, SIGN, ARRAY2, IJK2) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE indices
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      First array
      DOUBLE PRECISION ARRAY1 (DIMENSION_3, *)
!
!                      Second array
      DOUBLE PRECISION ARRAY2 (DIMENSION_3, *)
!
!                      IJK index for the first array
      INTEGER          IJK1
!
!                      IJK index for the second array
      INTEGER          IJK2
!
!                      Solids phase index
      INTEGER          M
!
!                      Sign to be used when setting ARRAY1.  Legal values
!                      are + or - 1.0.
      DOUBLE PRECISION SIGN
!
      M = 1 
      IF (MMAX > 0) THEN 
         ARRAY1(IJK1,:MMAX) = SIGN*ARRAY2(IJK2,:MMAX) 
         M = MMAX + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE EQUAL 
