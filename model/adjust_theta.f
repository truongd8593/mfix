!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_THETA (M, IER)                                  C
!  Purpose: Remove small negative values of theta caused by linear     C
!           solvers                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 02-APR-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ADJUST_THETA(M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE compar      !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IJK
!
!                      Solids phase
      INTEGER          M
      
!
!                      error indicator
      INTEGER          IER
!      CHARACTER*80     LINE
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IER = 0 
!
!HPF$ independent
      DO IJK = 1, IJKMAX2 
         IF (FLUID_AT(IJK)) THEN 
            IF (THETA_M(IJK,M) < ZERO_EP_S) THETA_M(IJK,M) = ZERO_EP_S 
!
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE ADJUST_THETA 
