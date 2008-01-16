!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Scalar_PROP(IER)                                       C
!  Purpose: Calculate diffusion coefficeint and sources for user-defined
!           scalars
!                                                                      C
!  Author:                                                    Date:    C
!  Reviewer:                                                  Date:    C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SCALAR_PROP( IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE scalars
      USE scales
      USE toleranc 
      USE trace
      USE constant
      USE compar
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                       Error index
      INTEGER          IER

      INTEGER          L, IJK, M
      
      DOUBLE PRECISION term0,term1,term2,term2p, term3, term4, &
                       term5, term6, term7, term8, &
		       Dterm1Dx, Dterm2Dx, Dterm3Dx, Dterm4Dx, Dterm5Dx, Dterm6Dx, &
		       Dterm7Dx, Dterm8Dx, trS, &
		       Fc,C_1,C_2, C3,C4,a0,a1,a2,sq2o3
      DOUBLE PRECISION voidRatio, voidRatioMin, voidRatioCrit
!                 
!-----------------------------------------------

      INCLUDE 'function.inc'

      IF(NScalar == 0) RETURN
      
      M = 1
      C_1 = -33.5d0
      C_2 = -341.4
      C3 = -339.7d0
      C4 = 446.5d0 !! Must be ZERO for true 1D cases (where S11=0) since tr(S*) must be zero
      a0 = 0.8d0 !0.8d0
      a1 = 0.31d0
      a2 = 0.02d0
      sq2o3 = (2d0/3d0)**2
      voidRatioMin = ep_star/(one-ep_star)
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do private(ijk, L)
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
	   
	   voidRatio = ep_g(ijk)/(one-ep_g(ijk))
	   voidRatioCrit = voidRatioMin + a1*exp(-a2*abs(scalar(ijk,1)+scalar(ijk,2))/(P_g(ijk)+P_ref))
	      
	   Fc = (one-a0)*(voidRatio - voidRatioMin)/(voidRatioCrit - voidRatioMin) + a0
	      
	   Fc = Fc * DSQRT(trD_s2(IJK,M)) !DSQRT(D_s11(IJK,M)**2 + D_s22(IJK,M)**2) ! should be
	   
	   trS = (scalar(ijk,1)+scalar(ijk,2))
	   if(trS < ZERO) THEN ! to make sure that trS is not within ]-zero_ep_s, +zero_ep_s[ (non-zero)
	     if(trS > -ZERO_EP_S) trS = -ZERO_EP_S
	   elseif(trS < ZERO_EP_S) THEN
	     trS = ZERO_EP_S
	   endif
           
	   DO L = 1, 3 
	   if(ep_g(ijk) < ONE) then
	      
	      Scalar_c(IJK, L) = ZERO
              Scalar_p(IJK, L) = ZERO

!            d (Scalar)/dt = S
!            S is linearized as S = Scalar_c - Scalar_p * Scalar
!            Scalar_c and Scalar_p must be >= 0  ! Scalar_c needs not be >0 for negative scalars (sof)
!            *** Uncomment next two lines ***
!
!            to compute the hypoplastic stresses, I used a formal linearization of the source term S
!            S = S(Scalar0) + d(S)/d(Scalar) (Scalarp-Scalar0)
!        ==> S = ( S(Scalar0) - d(S)/d(Scalar) Scalar0 + d(S)/d(Scalar) Scalarp )
!        
!           Following the MFIX convention that S = Scalar_c ""-"" Scalar_p * Scalar (note the minus sign)
!           so that: Scalar_c = S(Scalar0) - d(S)/d(Scalar) Scalar0  (this is from the "normal" linearization)
!           and      Scalar_p = - d(S)/d(Scalar) Scalarp      (this was factored with minus sign)
!           since Scalar_p must be >0, thus this linearization will be done only if d(S)/d(Scalar) < 0
!           otherwise I simply use: Scalar_p = zero and Scalar_c = S(Scalar0)
!
!           Note that I do this for all individual source terms instead of the sum of all terms to
!           increase stability as more terms contribute to the negative Am(0) center coefficient 
!          (diagonal term); hence all the if statements below.
!!!!!!!!!!
!!!!          Compute normal stress S11  in x-direction
!!!!!!!!!!              
!
              IF(L==1) THEN ! S11
	        term0 = 2d0*scalar(ijk,3)*vort_s12(IJK,M) ! vorticity is computed at end of calc_mu_s
		
		term1 = C_2*scalar(ijk,1)**2 * D_s11(IJK,M) / trS ! D_sij are computed at end of calc_mu_s
		Dterm1Dx = C_2*D_s11(IJK,M)*scalar(ijk,1)/trS * (2d0 - scalar(ijk,1)/trS)
		IF(Dterm1Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm1Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm1Dx * Scalar(IJK, L)
		ENDIF
		
		term2 = C_1*(scalar(ijk,1)+scalar(ijk,2)) * D_s11(IJK,M)
		Dterm2Dx = C_1*D_s11(IJK,M)
		IF(Dterm2Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm2Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm2Dx * Scalar(IJK, L)
		ENDIF
		
		term3 = C_2*scalar(ijk,1)*scalar(ijk,2) * D_s22(IJK,M) / trS
		Dterm3Dx = C_2*D_s22(IJK,M)*scalar(ijk,2)/trS* (ONE - scalar(ijk,1)/trS)
		IF(Dterm3Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm3Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm3Dx * Scalar(IJK, L)
		ENDIF
		
		term4 = 2d0*C_2*scalar(ijk,1)*scalar(ijk,3) * D_s12(IJK,M) / trS
		Dterm4Dx = 2d0*C_2*D_s12(IJK,M)*scalar(ijk,3)/trS* (ONE - scalar(ijk,1)/trS)
		IF(Dterm4Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm4Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm4Dx * Scalar(IJK, L)
		ENDIF
	        
		term5 = Fc*C3 * scalar(ijk,1)**2 / trS
		Dterm5Dx = Fc*C3 *scalar(ijk,1)/trS * (2d0 - scalar(ijk,1)/trS)
		IF(Dterm5Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm5Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm5Dx * Scalar(IJK, L)
		ENDIF
		
		term6 = Fc*C3 * scalar(ijk,3)**2 / trS
		Dterm6Dx = -Fc*C3 * scalar(ijk,3)**2/trS**2
		IF(Dterm6Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm6Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm6Dx * Scalar(IJK, L)
		ENDIF
	        
		term7 = Fc*C4 * (scalar(ijk,1)-scalar(ijk,2))**2/4d0 / trS
		Dterm7Dx = HALF*Fc*C4*(scalar(ijk,1)-scalar(ijk,2))/trS *(ONE-HALF* &
		          ((scalar(ijk,1)-scalar(ijk,2))/trS))
		IF(Dterm7Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm7Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm7Dx * Scalar(IJK, L)
		ENDIF
		
		term8 = Fc*C4 * scalar(ijk,3)**2 / trS
		Dterm8Dx = -Fc*C4*scalar(ijk,3)**2/trS**2
		IF(Dterm8Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm8Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm8Dx * Scalar(IJK, L)
		ENDIF
!!!!!!!!!!
!!!!          Compute normal stress S22  in y-direction
!!!!!!!!!!              
!
	      ELSEIF(L==2) THEN ! S22
	        term0 = -2d0*scalar(ijk,3)*vort_s12(IJK,M)
		
		term1 = C_2*scalar(ijk,2)**2 * D_s22(IJK,M) / trS
		Dterm1Dx = C_2*D_s22(IJK,M)*scalar(ijk,2)/trS * (2d0 - scalar(ijk,2)/trS)
		IF(Dterm1Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm1Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm1Dx * Scalar(IJK, L)
		ENDIF
		
		term2 = C_1*(scalar(ijk,1)+scalar(ijk,2)) * D_s22(IJK,M)
		Dterm2Dx = C_1*D_s22(IJK,M)
		IF(Dterm2Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm2Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm2Dx * Scalar(IJK, L)
		ENDIF
		
		term3 = C_2*scalar(ijk,1)*scalar(ijk,2) * D_s11(IJK,M) / trS
		Dterm3Dx = C_2*D_s11(IJK,M)*scalar(ijk,1)/trS* (ONE - scalar(ijk,2)/trS)
		IF(Dterm3Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm3Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm3Dx * Scalar(IJK, L)
		ENDIF
		
		term4 = 2d0*C_2*scalar(ijk,2)*scalar(ijk,3) * D_s12(IJK,M) / trS
		Dterm4Dx = 2d0*C_2*D_s12(IJK,M)*scalar(ijk,3)/trS* (ONE - scalar(ijk,2)/trS)
		IF(Dterm4Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm4Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm4Dx * Scalar(IJK, L)
		ENDIF
	        
		term5 = Fc*C3 * scalar(ijk,2)**2 / trS
		Dterm5Dx = Fc*C3 *scalar(ijk,2)/trS * (2d0 - scalar(ijk,2)/trS)
		IF(Dterm5Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm5Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm5Dx * Scalar(IJK, L)
		ENDIF
		
		term6 = Fc*C3 * scalar(ijk,3)**2 / trS
		Dterm6Dx = -Fc*C3 * scalar(ijk,3)**2/trS**2
		IF(Dterm6Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm6Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm6Dx * Scalar(IJK, L)
		ENDIF
	        
		term7 = Fc*C4 * (scalar(ijk,2)-scalar(ijk,1))**2/4d0 / trS
		Dterm7Dx = HALF*Fc*C4*(scalar(ijk,2)-scalar(ijk,1))/trS *(ONE-HALF* &
		          ((scalar(ijk,2)-scalar(ijk,1))/trS))
		IF(Dterm7Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm7Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm7Dx * Scalar(IJK, L)
		ENDIF
		
		term8 = Fc*C4 * scalar(ijk,3)**2 / trS
		Dterm8Dx = -Fc*C4*scalar(ijk,3)**2/trS**2
		IF(Dterm8Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm8Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm8Dx * Scalar(IJK, L)
		ENDIF
!!!!!!!!!!
!!!!          Compute Shear stress S12 = S21
!!!!!!!!!!              
              
	      ELSEIF(L==3) THEN ! S12
	        term0 = (scalar(ijk,2)-scalar(ijk,1))*vort_s12(IJK,M)
		
		term1 = 2d0*C_2*scalar(ijk,1)*scalar(ijk,3) * D_s11(IJK,M) / trS
		Dterm1Dx = 2d0*C_2*scalar(ijk,1)*D_s11(IJK,M) / trS
		IF(Dterm1Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm1Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm1Dx * Scalar(IJK, L)
		ENDIF
		
		term2 = 2d0*C_2*scalar(ijk,2)*scalar(ijk,3) * D_s22(IJK,M) / trS
		Dterm2Dx = 2d0*C_2*scalar(ijk,2)*D_s22(IJK,M) / trS
		IF(Dterm2Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm2Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm2Dx * Scalar(IJK, L)
		ENDIF
		
		term3= (2d0*C_2*scalar(ijk,3)**2) * D_s12(IJK,M) / trS
		Dterm3Dx = 4d0*C_2*scalar(ijk,3) * D_s12(IJK,M) / trS
		IF(Dterm3Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm3Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm3Dx * Scalar(IJK, L)
		ENDIF
		
		term4= (scalar(ijk,1)+scalar(ijk,2)) * D_s12(IJK,M)
	        
		term5 = Fc*C3 *scalar(ijk,3)  ! this simplification will not work in 3D
		Dterm5Dx = Fc*C3
		IF(Dterm5Dx < ZERO) THEN
		  Scalar_p(IJK, L) = Scalar_p(IJK, L) - Dterm5Dx
		  Scalar_c(IJK, L) = Scalar_c(IJK, L) - Dterm5Dx * Scalar(IJK, L)
		ENDIF
	        
		term6 = ZERO !S11* + S22* = ZERO for 2D cases
		term7 = ZERO
		term8 = ZERO
              
	      ENDIF
	      
	      
	      
	      Scalar_c(IJK, L) = Scalar_c(IJK, L) + term0+term1+term2+term3+term4+term5+term6+term7+term8
	      Scalar_p(IJK, L) = Scalar_p (IJK, L)
	      
	      Scalar_c(IJK, L) = Scalar_c (IJK, L) * ROP_S(IJK,M)
	      Scalar_p(IJK, L) = Scalar_p (IJK, L) * ROP_S(IJK,M)
	      
	! check for potential nan's in the computation (Not a Number, such as /zero etc.)
	      if(ISNAN(Scalar_c(IJK, L)) .OR. ISNAN(Scalar_p(IJK, L))) THEN
	      write(*,*) I_OF(IJK), J_OF(IJK), ep_g(ijk), voidRatioMin, voidRatio, voidRatioCrit, Fc,term1, scalar(ijk,1), &
	      (voidRatio - voidRatioMin)/(voidRatioCrit - voidRatioMin), DSQRT(trD_s2(IJK,M)), U_s(ijk,m), V_s(ijk,m)
	      stop
	      endif   
	! end of nan check
!
           else ! don't change source terms for empty cells
	      Scalar_c (IJK, L) = ZERO
              Scalar_p (IJK, L) = ZERO
	   
	   endif
!            Diffusion coefficient for User-defined Scalars
!            *** Uncomment next one line ***
              Dif_Scalar(IJK, L) = ZERO ! hypoplastic model has no diffusion terms.
	   END DO 
!
         ENDIF 
      END DO 
!\\Sendrecv operations - just to make sure all the variables computed are
!  are passed and updated locally - fool-proof approach - Sreekanth - 102199

!      call send_recv(Scalar_c,2)
!      call send_recv(Scalar_p,2)
!      call send_recv(Dif_Scalar,2)
      RETURN  
      END SUBROUTINE SCALAR_PROP 
