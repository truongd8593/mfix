!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_ROP(ROP, IER)
!  Purpose: Remove small negative values of density.                   C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-AUG-96   C
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
      SUBROUTINE ADJUST_ROP(ROP, IER) 
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
      USE geometry
      USE indices
      USE compar       !//d
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
!                      Indices
      INTEGER          IJK
!
!                      density
      DOUBLE PRECISION ROP(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!HPF$ independent
      DO IJK = ijkstart3, ijkend3 
         IF (FLUID_AT(IJK)) ROP(IJK) = DMAX1(ZERO,ROP(IJK)) 
      END DO 
      RETURN  
      END SUBROUTINE ADJUST_ROP 
