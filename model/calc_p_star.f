! TO DO:
! p_star calculation should be based on the sum of volume fractions of
! close-packed solids.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_P_star(EP_g, P_star, IER)                         C
!  Purpose: Calculate P_star in cells where solids continuity is solvedC
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-AUG-96  C
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
      SUBROUTINE CALC_P_STAR(EP_G, P_STAR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE pgcor
      USE pscor
      USE ur_facs 
      USE residual
      USE compar
!     USE fldvar
!     USE run
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Solids pressure
      DOUBLE PRECISION P_star(DIMENSION_3)
!
!                      Gas volume fraction
      DOUBLE PRECISION EP_g(DIMENSION_3)
      double precision calc_ep_star  !GERA

!!!HPF$ align P_star(:) with TT(:)
!!!HPF$ align EP_g(:) with TT(:)

!
!                      error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK

      DOUBLE PRECISION dPs
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
!     INCLUDE 'ep_s1.inc'
!     INCLUDE 'ep_s2.inc'
!
!
!
!!$omp parallel do private(ijk)
!!!HPF$ independent

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 
!GERA ******************
	    if (MMAX >= 2) EP_star = calc_ep_star(ijk, ier)
!END GERA***************
            IF (EP_G(IJK) < EP_STAR) THEN 
               P_STAR(IJK) = NEG_H(EP_G(IJK)) 
            ELSE 
               P_STAR(IJK) = ZERO 
            ENDIF 
         ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE CALC_P_STAR 

!TO DO:
! Extend the formulation for more than two sizes
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_ep_star(IJK,IER)                                  C
!  Purpose: calculte the local value of maximum packing                C
!                                                                      C
!  Author: D. Gera/M. Syamlal                         Date: 31-DEC-02  C
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
      Double Precision function CALC_ep_star(IJK,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE compar
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
      double precision xbar
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
 
 	if ((EP_s(IJK,1)+EP_s(IJK,2)) .NE. 0) THEN
	   xbar = EP_s(IJK,1)/(EP_s(IJK,1)+EP_s(IJK,2))

	   if (xbar .LE. ep_s_max_ratio(1,2)) THEN
	      CALC_EP_star =MAX(0.36d0, (1.-(((ep_s_max(1)-ep_s_max(2))+(1.-d_p_ratio(1,2))*(1-ep_s_max(1))*ep_s_max(2))*(ep_s_max(1)+(1-ep_s_max(1)) &
        	              *ep_s_max(2))*xbar/ep_s_max(1)+ep_s_max(2))))
     	   else
    	      CALC_EP_star =MAX(0.36d0, (1.-((1. -d_p_ratio(1,2))*(ep_s_max(1)+(1-ep_s_max(1))*ep_s_max(2))*(1. -xbar) +ep_s_max(1))))
	   end if
	else
	   CALC_EP_star = MIN(ep_s_max(1), ep_s_max(2))
	end if
 
	
      RETURN  
      END function CALC_ep_star

!
