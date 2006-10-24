!     TO DO:
!     p_star calculation should be based on the sum of volume fractions of
!     close-packed solids.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: CALC_P_star(EP_g, P_star, IER)                         C
!     Purpose: Calculate P_star in cells where solids continuity is solvedC
!                                                                         C
!     Author: M. Syamlal                                 Date: 21-AUG-96  C
!     Reviewer:                                          Date:            C
!                                                                         C
!                                                                         C
!     Literature/Document References:                                     C
!                                                                         C
!     Variables referenced:                                               C
!     Variables modified:                                                 C
!                                                                         C
!     Local variables:                                                    C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE CALC_P_STAR(EP_G, P_STAR, IER) 
!...  Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...  Switches: -xf
!     
!-----------------------------------------------
!     M o d u l e s 
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
      USE run
      USE visc_s
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     Solids pressure
      DOUBLE PRECISION P_star(DIMENSION_3)
!     
!     Gas volume fraction
      DOUBLE PRECISION EP_g(DIMENSION_3)
      double precision calc_ep_star !GERA
      DOUBLE PRECISION , EXTERNAL :: BLEND_FUNCTION

!!!   HPF$ align P_star(:) with TT(:)
!!!   HPF$ align EP_g(:) with TT(:)

!     
!     error index
      INTEGER          IER
!     
!     Indices
      INTEGER          IJK

      DOUBLE PRECISION dPs
!     Blend Factor
      Double Precision blend
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
!     INCLUDE 'ep_s1.inc'
!     INCLUDE 'ep_s2.inc'
!     
!     
!     
!     !$omp parallel do private(ijk)
!!!   HPF$ independent

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 
!     GERA ******************
!     if Yu_Standish or Fedors_Landel correlations are not used, ep_star_array will not
!     be modified (sof Nov-16-2005)
!       changed blend_start to 0.99*ep_star from 0.95*ep_star
!       changed blend_end to 1.01*ep_star from 1.02*ep_star
!       [ceaf 2006-03-31]
!       added option for sigmoid function [sp 2006-10-24]
                 
	    if (YU_STANDISH .OR. FEDORS_LANDEL) THEN
              EP_star_array(ijk) = calc_ep_star(ijk, ier)
	      IF(BLENDING_STRESS.AND.TANH_BLEND) THEN
                ep_g_blend_start(ijk) = ep_star_array(ijk) * 0.99d0
                ep_g_blend_end(ijk)   = ep_star_array(ijk) * 1.01d0
              ELSE IF(BLENDING_STRESS.AND.SIGM_BLEND) THEN
                ep_g_blend_start(ijk) = ep_star * 0.97d0
                ep_g_blend_end(ijk) = ep_star * 1.01d0
              ELSE
                ep_g_blend_start(ijk) = ep_star_array(ijk)
                ep_g_blend_end(ijk)   = ep_star_array(ijk)
	      ENDIF
	    endif
!     END GERA***************
	    
            IF (EP_G(IJK) < EP_g_blend_end(ijk)) THEN 
               P_STAR(IJK) = NEG_H(EP_G(IJK),EP_g_blend_end(ijk))
               IF(BLENDING_STRESS) THEN
	         blend =  blend_function(IJK)
                 P_STAR (IJK) = (1.0d0-blend) * P_STAR (IJK)
               ENDIF
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
! The commented version has the following ref:                         C
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
      USE run
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
!                       true maximum packing for the mixture
       DOUBLE PRECISION EPs_max_local
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
! rearranging to start from coarsest to finest particles (see check_data_06.f)
! this is the way the algorithm was written by Yu and Standish (sof).
!
      IF (CALL_DQMOM) THEN
!
       DO I = 1, MMAX
	 DP_TMP(I) = D_P(IJK,I)
	 EPs_TMP(I) = EP_s(IJK,I)
	 EPs_max_TMP(I) = ep_s_max(I)
       END DO
!
       DO I = 1, MMAX	 
	 
	 DO J = I , MMAX
	   
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
      ELSE
       DO I = 1, MMAX 
	 DP_TMP(I) = D_P(IJK,M_MAX(I))
	 EPs_TMP(I) = EP_s(IJK,M_MAX(I))
	 EPs_max_TMP(I) = ep_s_max(M_MAX(I))
       END DO
!
      ENDIF !for dqmom
!
! compute equations 25 in Yu-Standish (this is also needed by Fedors_Landel)
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
	   CALC_EP_star = ONE - EPs_max_TMP(1) !return first phase ep_s_max in case very dilute
	   RETURN
	 ENDIF
       
       END DO 
!
       IF(YU_STANDISH) THEN  
!
! compute equation 23-24 in Yu-Standish
!
         DO I = 1, MMAX
           DO J = 1, MMAX
	   
	     IF( R_IJ(I, J) .LE. 0.741d0) THEN

	       IF( J .LT. I ) THEN
	         X_IJ(I, J) = (ONE - R_IJ(I, J)*R_IJ(I, J))/(2.0d0 -  EPs_max_TMP(I))
	       ELSE
	         X_IJ(I, J) = ONE - (ONE - R_IJ(I, J)*R_IJ(I, J))/(2.0d0 -  EPs_max_TMP(I))
	       ENDIF

	       P_IJ(I, J) = EPs_max_TMP(I) + EPs_max_TMP(I)* (ONE-EPs_max_TMP(I)) *  &
	                    (ONE - 2.35d0*R_IJ(I, J) + 1.35d0*R_IJ(I, J)*R_IJ(I, J))
	     ELSE
	   
	       P_IJ(I, J) = EPs_max_TMP(I)
	     ENDIF

	   END DO
         END DO
!
! Compute equation 22
!
	 EPs_max_local = ONE
	 DO I = 1, MMAX
           SUM_LOCAL = ZERO
         
	   IF( I .GE. 2) THEN
	     DO J = 1, (I-1)
	     
	       IF( P_IJ(I, J) == EPs_max_TMP(I) ) THEN
	         SUM_LOCAL = SUM_LOCAL
	       ELSE
	         SUM_LOCAL = SUM_LOCAL + (ONE - EPs_max_TMP(I)/P_IJ(I, J))*COMP_X_I(J)/X_IJ(I, J)
	       ENDIF
	   
	     END DO
	   ENDIF
	   
	   IF( (I+1) .LE. MMAX) THEN  
	     DO J = (I+1), MMAX
	       
	       IF( P_IJ(I, J) == EPs_max_TMP(I) ) THEN
	         SUM_LOCAL = SUM_LOCAL
	       ELSE
	         SUM_LOCAL = SUM_LOCAL + (ONE - EPs_max_TMP(I)/P_IJ(I, J))*COMP_X_I(J)/X_IJ(I, J)
	       ENDIF

	     END DO
	   ENDIF
	       
	   IF (SUM_LOCAL .NE. ZERO) THEN
	     P_IT(I) = EPs_max_TMP(I)/(ONE - SUM_LOCAL)
	   ELSE
	     P_IT(I) = ONE ! do nothing if particles have same diameter
	   ENDIF
!   
	   EPs_max_local = MIN(P_IT(I), EPs_max_local)

         END DO !for I

! for the case of all phases having same diameter	 

	 IF (EPs_max_local == ONE) EPs_max_local = EPs_max_TMP(1)
	 
	 CALC_EP_star = ONE - EPs_max_local
!
! end modifications by sof (May-02-2005)


!
! Part implemented by Dinesh for binary mixture, uncomment to use (Sof)
!
! 	if ((EP_s(IJK,1)+EP_s(IJK,2)) .NE. ZERO) THEN
!	   xbar = EP_s(IJK,1)/(EP_s(IJK,1)+EP_s(IJK,2))
!
!	   if (xbar .LE. ep_s_max_ratio(1,2)) THEN
!	      CALC_EP_star =MAX(0.36d0, (ONE-(((ep_s_max(1)-ep_s_max(2))+&
!              (ONE-d_p_ratio(1,2))*(ONE-ep_s_max(1))*ep_s_max(2))*(ep_s_max(1)+&
!              (ONE-ep_s_max(1)) *ep_s_max(2))*xbar/ep_s_max(1)+ep_s_max(2))))
!     	   else
!    	      CALC_EP_star =MAX(0.36d0, (ONE-((ONE -d_p_ratio(1,2))*(ep_s_max(1)&
!              +(ONE-ep_s_max(1))*ep_s_max(2))*(ONE -xbar) +ep_s_max(1))))
!	   end if
!	else
!	   CALC_EP_star = ONE - MIN(ep_s_max(1), ep_s_max(2)) !corrected by sof
!	end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I suggest to use the following code (below) instead of the above commented few lines
! of code because the phases were not rearranged and I didn't want to modify it (sof)
! If you don't understand what's going on, contact me: sof@fluent.com
!
! In the case of binary mixture (Fedors-Landel empirical correlation)
!
!
       ELSEIF(FEDORS_LANDEL) THEN
!       
         IF(COMP_X_I(1) .LE. (EPs_max_TMP(1)/(EPs_max_TMP(1)+ &
	   (ONE - EPs_max_TMP(1))*EPs_max_TMP(2)))) THEN
	 
	   CALC_EP_star = (EPs_max_TMP(1)-EPs_max_TMP(2)+(1-sqrt(R_IJ(2, 1)))* &
	               (ONE - EPs_max_TMP(1)) *EPs_max_TMP(2))*(EPs_max_TMP(1) &
		       +(ONE-EPs_max_TMP(1))*EPs_max_TMP(2))*  &
			  COMP_X_I(1)/EPs_max_TMP(1) + EPs_max_TMP(2)
	 ELSE
	   CALC_EP_star = (ONE-sqrt(R_IJ(2, 1)))*(EPs_max_TMP(1)+(ONE-EPs_max_TMP(1))* &
	               EPs_max_TMP(2))*(ONE - COMP_X_I(1)) + EPs_max_TMP(1)
         ENDIF
         CALC_EP_star = ONE - CALC_EP_star ! this is gas volume fraction at packing
       ENDIF ! for Yu_Standish and Fedors_Landel correlations
!	
      RETURN  
      END function CALC_ep_star

!


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: blend_function(IJK)                                    C
!  Purpose: To calculate blending function                             C
!                                                                      C
!  Author: S. Pannala                                 Date: 28-FEB-06  C
!  Reviewer:                                          Date:            C
!  Modified:                                          Date: 24-OCT-06  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      Double Precision function blend_function(IJK) 
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
      USE run
      USE visc_s
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
!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
!     Blend Factor
      Double Precision:: blend, blend_right
!     Scale Factor
      Double Precision, Save:: scale
!     Midpoint
      Double Precision, Save:: ep_mid_point
!
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

      IF(TANH_BLEND) then  ! Tan hyperbolic blending of stresses
!
         IF(EP_g(IJK) .LT. ep_g_blend_end(ijk).AND. EP_g(IJK) .GT. ep_g_blend_start(ijk)) THEN
            ep_mid_point = (ep_g_blend_end(IJK)+ep_g_blend_start(IJK))/2.0d0
            blend = tanh(2.0d0*pi*(ep_g(IJK)-ep_mid_point)/ &
            (ep_g_blend_end(IJK)-ep_g_blend_start(IJK)))
            blend = (blend+1.0d0)/2.0d0
         ELSE IF(EP_g(IJK) .GE. ep_g_blend_end(ijk)) THEN
            blend = 1.0d0
         ELSE IF(EP_g(IJK) .LE. ep_g_blend_start(ijk)) THEN
            blend = 0.0d0
         ENDIF
!
      ELSEIF(SIGM_BLEND) then !  Truncated and Scaled Sigmoidal blending of stresses

         IF(FIRST_PASS) THEN
            blend_right =  1.0d0/(1+0.01d0**((ep_g_blend_end(IJK)-ep_star_array(IJK))&
            /(ep_g_blend_end(IJK)-ep_g_blend_start(IJK))))
            blend_right = (blend_right+1.0d0)/2.0d0
            scale = 1.0d0/blend_right
            write(*,*) 'Blending value at end and scaling factor', blend_right, scale
            FIRST_PASS = .FALSE.
         ENDIF

         IF(EP_g(IJK) .LT. ep_g_blend_end(ijk)) THEN
            blend =  scale/(1+0.01d0**((ep_g(IJK)-ep_star_array(IJK))&
            /(ep_g_blend_end(IJK)-ep_g_blend_start(IJK))))
         ELSE IF(EP_g(IJK) .GE. ep_g_blend_end(ijk)) THEN
            blend = 1.0d0
         ENDIF
!
      ENDIF 
!     
      blend_function = blend
!     
      RETURN  
      END function blend_function

!
