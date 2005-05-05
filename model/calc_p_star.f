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
! sof (done)
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_ep_star(IJK,IER)                                  C
!  Purpose: calculte the local value of maximum packing                C
!                                                                      C
!  Author: D. Gera/M. Syamlal                         Date: 31-DEC-02  C
!  Reviewer:                                          Date:            C
!  Modified: Sof                                      Date: 02-May-05  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!    A.B. Yu and N. Standish. Powder Tech, 52 (1987) 233-241           C
! The commented version implemented by Dinesh has the following ref:   C
!     R.F. Fedors and R.F. Landel. Powder Tech, 23 (1979) 225-231      C
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
      USE toleranc
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
      INTEGER          IJK, I, J
      double precision xbar
!
! start sof modifications (02-May-05)
!
!                       maximum packing for the mixture
       DOUBLE PRECISION P_IT(MMAX)
!
!                       maximum packing fraction for a binary mixture
       DOUBLE PRECISION P_IJ(MMAX, MMAX)
!
!                       particle diameter ratio
       DOUBLE PRECISION R_IJ(MMAX, MMAX)
!      
!                       fractional solids volume corresponding to P_IJ
       DOUBLE PRECISION X_IJ(MMAX, MMAX)
!      
!                       fractional solids volume in a mixture
       DOUBLE PRECISION COMP_X_I(MMAX), SUM_LOCAL ! this is Xj in eq. 22 of Yu-Standish
!      
!                       tmp variables for rearanging solids phases from coarsest to finest
       DOUBLE PRECISION DP_TMP(MMAX), EPs_TMP(MMAX), EPs_max_TMP(MMAX), old_value
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
! rearanging phases to start from coarsest to finest particles
! this is the way the algorithm was written by Yu and Standish.
!
       DO I = 1, MMAX
         
	 DP_TMP(I) = D_P(I)
	 EPs_TMP(I) = EP_s(IJK,I)
	 EPs_max_TMP(I) = ep_s_max(I)
       END DO

       DO I = 1, MMAX	 
	 DO J = I, MMAX
	   
	   IF(DP_TMP(I) < DP_TMP(J)) THEN
	     
	     old_value = DP_TMP(I)
	     DP_TMP(I) = DP_TMP(J)
	     DP_TMP(J) = old_value
	     
	     old_value = EPs_TMP(I)
	     EPs_TMP(I) = EPs_TMP(J)
	     EPs_TMP(J) = old_value
	     
	     old_value = EPs_max_TMP(I)
	     EPs_max_TMP(I) = EPs_max_TMP(J)
	     EPs_max_TMP(J) = old_value

	   ENDIF
	   
	 ENDDO
       END DO
!
!end of solids phase rearangment
!
! compute equations 25 in Yu-Standish
!      
       DO I = 1, MMAX
         SUM_LOCAL = ZERO
	 DO J = 1, MMAX

	   IF( I .GE. J) THEN
	     R_IJ(I, J) = DP_TMP(I)/DP_TMP(J)
	   ELSE
	     R_IJ(I, J) = DP_TMP(J)/DP_TMP(I)
	   ENDIF
	 SUM_LOCAL = SUM_LOCAL + EPs_TMP(J)
	 END DO
	 
	 IF(SUM_LOCAL > DIL_EP_s) THEN
	   COMP_X_I(I) = EPs_TMP(I)/SUM_LOCAL ! fractional solids volume see eq. 20
	 ELSE
	   CALC_EP_star = EPs_max_TMP(1) !return first phase ep_s_max in case very dilute
	   RETURN
	 ENDIF
       
       END DO     
!
! compute equation 23-24 in Yu-Standish
!
       DO I = 1, MMAX
         DO J = 1, MMAX
	   
	   IF( R_IJ(I, J) .LE. 0.741) THEN

	     IF( J .LT. I ) THEN
	       X_IJ(I, J) = (1.0 - R_IJ(I, J)*R_IJ(I, J))/(2.0 -  EPs_max_TMP(I))
	     ELSE
	       X_IJ(I, J) = 1.0 - (1.0 - R_IJ(I, J)*R_IJ(I, J))/(2.0 -  EPs_max_TMP(I))
	     ENDIF

	     P_IJ(I, J) = EPs_max_TMP(I) + EPs_max_TMP(I)* (1.0-EPs_max_TMP(I)) *  &
	                  (1.0 - 2.35*R_IJ(I, J) + 1.35*R_IJ(I, J)*R_IJ(I, J))
	   ELSE
	   
	     P_IJ(I, J) = EPs_max_TMP(I)
	   ENDIF

	 END DO
       END DO
!
! Compute equation 22
!
	 CALC_EP_star = 1.0
	 DO I = 1, MMAX
           SUM_LOCAL = 0.0
         
	   IF( I .GE. 2) THEN
	     DO J = 1, (I-1)
	     
	       IF( P_IJ(I, J) == EPs_max_TMP(I) ) THEN
	         SUM_LOCAL = SUM_LOCAL
	       ELSE
	         SUM_LOCAL = SUM_LOCAL + (1.0 - EPs_max_TMP(I)/P_IJ(I, J))*COMP_X_I(J)/X_IJ(I, J)
	       ENDIF
	   
	     END DO
	   ENDIF
	   
	   IF( (I+1) .LE. MMAX) THEN  
	     DO J = (I+1), MMAX
	       
	       IF( P_IJ(I, J) == EPs_max_TMP(I) ) THEN
	         SUM_LOCAL = SUM_LOCAL
	       ELSE
	         SUM_LOCAL = SUM_LOCAL + (1.0 - EPs_max_TMP(I)/P_IJ(I, J))*COMP_X_I(J)/X_IJ(I, J)
	       ENDIF

	     END DO
	   ENDIF
	       
	   IF (SUM_LOCAL .NE. 0.0) THEN
	     P_IT(I) = EPs_max_TMP(I)/(1.0 - SUM_LOCAL)
	   ELSE
	     P_IT(I) = 1.0 ! do nothing if particles have same diameter
	   ENDIF
!   
	   CALC_EP_star = MIN(P_IT(I), CALC_EP_star)

         END DO !for I

! for the case of all phases having same diameter	 

	 IF (CALC_EP_star == 1.0) CALC_EP_star = EPs_max_TMP(1)
!
! end modifications by sof (May-02-2005)

!
! Part implemented by Dinesh for binary mixture, uncomment to use (Sof)
!
! 	if ((EP_s(IJK,1)+EP_s(IJK,2)) .NE. 0) THEN
!	   xbar = EP_s(IJK,1)/(EP_s(IJK,1)+EP_s(IJK,2))
!
!	   if (xbar .LE. ep_s_max_ratio(1,2)) THEN
!	      CALC_EP_star =MAX(0.36d0, (1.-(((ep_s_max(1)-ep_s_max(2))+&
!              (1.-d_p_ratio(1,2))*(1-ep_s_max(1))*ep_s_max(2))*(ep_s_max(1)+&
!              (1-ep_s_max(1)) *ep_s_max(2))*xbar/ep_s_max(1)+ep_s_max(2))))
!     	   else
!    	      CALC_EP_star =MAX(0.36d0, (1.-((1. -d_p_ratio(1,2))*(ep_s_max(1)&
!              +(1-ep_s_max(1))*ep_s_max(2))*(1. -xbar) +ep_s_max(1))))
!	   end if
!	else
!	   CALC_EP_star = MIN(ep_s_max(1), ep_s_max(2))
!	end if
! 
	
      RETURN  
      END function CALC_ep_star

!
