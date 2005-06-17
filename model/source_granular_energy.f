!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_GRANULAR_ENERGY(sourcelhs,sourcerhs,IJK,M,IER)  C
!  Purpose: Calculate the source terms in the granular energy equation C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 04-FEB-98  C
!  Reviewer: M. Syamlal                               Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!     Local variables: sourcelhs, sourcerhs                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE compar        !//d
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
      INTEGER          IJK, I, J, K, M, MM
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs
!
!                      Particle relaxation time
      DOUBLE PRECISION Tau_12_st
!
!                      Sum of eps*G_0
      DOUBLE PRECISION SUM_EpsGo
!
!                      Slip velocity
      DOUBLE PRECISION VSLIP
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
!
      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 
!
      SOURCERHS = (ZMAX(LAMBDA_S_C(IJK,M))*TRD_S_C(IJK,M)**2+2.*MU_S_C(IJK,M)*&
         TRD_S2(IJK,M)+P_S_C(IJK,M)*ZMAX((-TRD_S_C(IJK,M))))*VOL(IJK)  
!
      IF(SIMONIN) THEN
        SOURCERHS = SOURCERHS +  F_GS(IJK,M)*K_12(IJK)*VOL(IJK)
!
      ELSE IF(AHMADI) THEN
	IF(Ep_s(IJK,M) > DIL_EP_S .AND. F_GS(IJK,1) > small_number) THEN
          Tau_12_st = Ep_s(IJK,M)*RO_s(M)/F_GS(IJK,1)
	  SOURCERHS = SOURCERHS + 2.D+0*F_GS(IJK,M)* (ONE/(ONE+Tau_12_st/  &
                      (Tau_1(ijk)+small_number)))*K_Turb_G(IJK)*VOL(IJK)
	ELSE
	  SOURCERHS = SOURCERHS
	ENDIF
!
      ELSE ! no modifications done to the original KTGF
!
        VSLIP = (U_S(IJK,M)-U_G(IJK))**2 + (V_S(IJK,M)-V_G(IJK))**2 + (W_S(IJK,M)&
           -W_G(IJK))**2 
        VSLIP = DSQRT(VSLIP) 
!
        SOURCERHS = SOURCERHS + (SWITCH*81D0*EP_S(IJK,M)*(MU_G(IJK)*VSLIP)**2D0/(&
           G_0(IJK,M,M)*D_P(IJK,M)**3D0*RO_S(M)*(PI*THETA_M(IJK,M)+SMALL_NUMBER)**&
           0.5D0))*VOL(IJK) 
      ENDIF
!
!---------------------------------------------------------------------
!  The following lines are commented out since Kphi_s has been set to zero
!  in subroutine calc_mu_s.  To activate the feature uncomment the following
!  lines and the lines in calc_mu_s.
!
!
!      sourcerhs=sourcerhs+
!     &    AVG_X_h(Kphi_s(IJK,M),Kphi_s(EAST_OF(IJK),M),I)
!     &    *(EP_s(EAST_OF(IJK),M)-EP_s(IJK,M))*oDX_E(I)
!     &    *AYZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(WEST_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,I_OF(WEST_OF(IJK)))*(EP_s(IJK,M)-EP_s(WEST_OF(IJK),M))
!     &    *oDX_E(I_OF(WEST_OF(IJK)))*AYZ(WEST_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Y_h(Kphi_s(IJK,M),Kphi_s(NORTH_OF(IJK),M),J)
!     &    *(EP_s(NORTH_OF(IJK),M)-EP_s(IJK,M))*oDY_N(J)
!     &    *AXZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(SOUTH_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,J_OF(SOUTH_OF(IJK)))*(EP_s(IJK,M)-EP_s(SOUTH_OF(IJK),M))
!     &    *oDY_N(J_OF(SOUTH_OF(IJK)))
!     &    *AXZ(JM_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Z_h(Kphi_s(IJK,M),Kphi_s(TOP_OF(IJK),M),K)
!     &    *(EP_s(TOP_OF(IJK),M)-EP_s(IJK,M))
!     &    *oX(I)*oDZ_T(K)
!     &    *AXY(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_Z_h(Kphi_s(BOTTOM_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,K_OF(BOTTOM_OF(IJK)))*(EP_s(IJK,M)-
!     &    EP_s(BOTTOM_OF(IJK),M))
!     &    *oX(I)*oDZ_T(K_OF(BOTTOM_OF(IJK)))
!     &    *AXY(KM_OF(IJK))
!---------------------------------------------------------------------
!
! Changes needed for multitype particles, sof June 16 2005
! Sum of eps*G_0 is used instead of Eps*G_0
!
      SUM_EpsGo = ZERO
      DO MM = 1, MMAX
        SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,MM,MM)
      ENDDO
      SOURCELHS = ((48./DSQRT(PI))*ETA*(1.-ETA)*ROP_S(IJK&
         ,M)*SUM_EpsGo*DSQRT(THETA_M(IJK,M))/D_P(IJK,M)+ &
	 P_S_C(IJK,M)*ZMAX((TRD_S_C(IJK,M)))/(THETA_M(IJK,M)+SMALL_NUMBER) &
	 +ZMAX((-LAMBDA_S_C(IJK,M)))*TRD_S_C(IJK,M)**2/(THETA_M(IJK,M)+&
         SMALL_NUMBER))*VOL(IJK)  
!
      IF(SIMONIN .OR. AHMADI) THEN
        SOURCELHS = SOURCELHS +  3.0 *F_GS(IJK,M)*VOL(IJK)
!
      ELSE ! no modifications done to the original KTGF
!
        SOURCELHS = SOURCELHS + SWITCH *3.0 *F_GS(IJK,M)*VOL(IJK)
      ENDIF
!
      RETURN  
      END SUBROUTINE SOURCE_GRANULAR_ENERGY 
