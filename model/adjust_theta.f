!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_THETA (M, IER)                                  C
!  Purpose: Remove small negative values of theta caused by linear     C
!           solvers                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 02-APR-98  C
!  Reviewer:                                          Date:            C
!  Modified: S. Benyahia                              Date: 02-AUG-06  C
!  Purpose: check for small negative numbers at walls (not just fluid  C
!           cells)                                                     C
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
      USE compar
      
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
!                      Particle mass and diameter for use with IA theory only
!                      because theta definition includes mass of particle.
      DOUBLE PRECISION M_PM, D_PM
!
!                      small value of theta_m, 1 cm2/s2 = 1e-4 m2/s2
      DOUBLE PRECISION smallTheta
!
!                      error indicator
      INTEGER          IER
!      CHARACTER*80     LINE
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IER = 0 
      smallTheta = (to_SI)**4 * ZERO_EP_S
!
!!!HPF$ independent
      DO IJK = IJKSTART3, IJKEND3
         IF ( FLUID_AT(IJK) ) THEN 
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
	       D_PM = D_P(IJK,M)
	       M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
	       IF (THETA_M(IJK,M) < smallTheta*M_PM) THETA_M(IJK,M) = smallTheta*M_PM 
            ELSE
	       IF (THETA_M(IJK,M) < smallTheta) THETA_M(IJK,M) = smallTheta 
            ENDIF
!
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE ADJUST_THETA 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
