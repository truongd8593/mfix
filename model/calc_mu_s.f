!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!     C
!     Module name: CALC_MU_s(M, IER)                                      C
!     Purpose: Calculate granular stress terms: THETA, P_s, LAMBDA_s, MU_sC
!     C
!     Author: W. Rogers                                  Date: 04-mar-92  C
!     Reviewer: M. Syamlal                               Date: 16-MAR-92  C
!     C
!     Revision Number: 1                                                  C
!     Purpose: Modifications for cylindrical geometry                     C
!     Author: M. Syamlal                                 Date: 15-MAY-92  C
!     Revision Number: 2                                                  C
!     Purpose: Add volume-weighted averaging statement functions for      C
!     variable grid capability                                   C
!     Author:  W. Rogers                                 Date: 21-JUL-92  C
!     Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!     Revision Number: 3                                                  C
!     Purpose: Add frictional-flow stress terms                           C
!     Author: M. Syamlal                                 Date: 10-FEB-93  C
!     Revision Number: 4                                                  C
!     Purpose: Add Boyle-Massoudi stress terms                            C
!     Author: M. Syamlal                                 Date: 2-NOV-95   C
!     Revision Number: 5                                                  C
!     Purpose: MFIX 2.0 mods  (old name CALC_THETA)                       C
!     Author: M. Syamlal                                 Date: 24-APR_96  C
!     Author: Kapil Agrawal, Princeton University        Date: 6-FEB-98   C
!     Revision Number: 6                                                  C
!     Purpose: Add calculation of viscosities and conductivities for use  C
!     with granular temperature PDE. New common block contained  C
!     in 'trace.inc' contains trD_s_C(DIMENSION_3, DIMENSION_M)  C
!     and trD_s2(DIMENSION_3, DIMENSION_M)                       C
!     Author: Anuj Srivastava, Princeton University      Date: 20-APR-98  C
!     Revision Number:7                                                   C
!     Purpose: Add calculation of frictional stress terms                 C
!     C
!     Author: Sofiane Benyahia, Fluent Inc.      Date: 02-01-05           C
!     Revision Number:8                                                   C
!     Purpose: Add Simonin and Ahmadi models                              C
!     C
!     Literature/Document References:                                     C
!     1- Simonin, O., 1996. Combustion and turbulence in two-phase flows. C
!     Von Karman institute for fluid dynamics, lecture series 1996-02  C
!     2- Balzer, G., Simonin, O., Boelle, A., and Lavieville, J., 1996.   C
!     A unifying modelling approach for the numerical prediction of    C
!     dilute and dense gas-solid two phase flow. CFB5, 5th int. conf.  C
!     on circulating fluidized beds, Beijing, China.                   C
!     3- Cao, J. and Ahmadi, G., 1995. Gas-particle two-phase turbulent   C
!     flow in a vertical duct. Int. J. Multiphase Flow, vol. 21 No. 6  C
!     pp. 1203-1228.                                                   C
!     C
!     Author: Sreekanth Pannala, ORNL            Date: 10-08-05           C
!     Revision Number:9                                                   C
!     Purpose: Rewrite different modules to increase modularity           C
!     C
!     Variables referenced: U_s, V_s, W_s, IMAX2, JMAX2, KMAX2, DX, DY,   C
!     DZ, IMJPK, IMJK, IPJMK, IPJK, IJMK, IJKP,     C
!     IMJKP, IPJKM, IJKM, IJMKP, IJPK, IJPKM, IJMK, C
!     M,  RO_s, C_e, D_p, Pi, G_0, X                C
!     C
!     Variables modified: I, J, K, IJK, MU_s, LAMBDA_s, P_s               C
!     C
!     Local variables: K_1m, K_2m, K_3m, K_4m, D_s, U_s_N, U_s_S, V_s_E,  C
!     V_s_W, U_s_T, U_s_B, W_s_E, W_s_W, V_s_T, V_s_B,   C
!     W_s_N, W_s_S, trD_s_C, W_s_C                       C
!     trD_s2, EP_s2xTHETA, EP_sxSQRTHETA, I1, I2, U_s_C, C
!     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE CALC_MU_s(M, IER)

!     
      USE run
      USE vshear
      USE visc_s
      USE physprop
      USE constant
      USE fldvar
      USE compar
      USE indices
      USE geometry
      Implicit NONE

!     Local Variables
!     Error index
      INTEGER          IER
!     Solids Phase
      INTEGER          M, IJK
!     Blend Factor
      Double Precision blend

      DOUBLE PRECISION , EXTERNAL :: BLEND_FUNCTION
     
      Include 'function.inc'


! GHD Theory is called only for the mixture granular energy, i.e. for m == mmax
      IF (TRIM(KT_TYPE) == 'GHD' .AND. M /= MMAX) RETURN
! end of GHD theory

      IF (SHEAR) CALL add_shear(M)

      CALL init_mu_s(M, IER)    ! initializing/calculating all the quantities needed for various options
     
      IF (SHEAR) call remove_shear(M)

      IF(MU_s0 /= UNDEFINED) RETURN ! constant solids viscosity case 

!     GRANULAR_ENERGY
!     .FALSE.
!     EP_g < EP_star   -->    friction_schaeffer
!     EP_g >= EP_star  -->    viscous (algebraic)
      
!     GRANULAR_ENERGY
!     .TRUE.
!     
!     FRICTION
!     .TRUE.
!     EP_s(IJK,M) > EPS_f_min  -->  friction + viscous(pde)
!     EP_s(IJK,M) < EP_f_min   -->  viscous (pde)
      
!     FRICTION
!     .FALSE.
!     EP_g < EP_star  -->  friction_schaeffer + viscous(pde)
!     EP_g >= EP_star -->  viscous (pde)
     
!     
!     Viscous-flow stress tensor
!     
      IF(.NOT.GRANULAR_ENERGY) then
         call gt_algebraic(M,IER) !algebraic granular energy equation
      ELSE                        !granular energy transport equation
          IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               CALL gt_pde_ia_nonep(M,IER)   ! complete polydisperse IA theory
          ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN
               CALL gt_pde_gd_99(M,IER)      ! monodisperse GD theory
          ELSEIF (TRIM(KT_TYPE) == 'GHD') THEN
               CALL TRANSPORT_COEFF_GHD(M,IER) ! GHD theory for mixture temperature               
          ELSE
               CALL gt_pde(M,IER)   ! This is also used whith Simonin or Ahmadi models
          END IF
      ENDIF
!     
!     Frictional stress tensors
!     
      IF (SCHAEFFER .AND. CLOSE_PACKED(M)) call friction_schaeffer(M,IER) ! Schaeffer's Frictional Formulation
     
      IF(FRICTION .AND. CLOSE_PACKED(M)) call friction_princeton(M,IER) ! Princeton's frictional implementation
      
      IF(BLENDING_STRESS) THEN

         DO 200 IJK = ijkstart3, ijkend3
    
            blend =  blend_function(IJK)
           
            Mu_s_c(IJK,M) = Mu_s_v(IJK)
            Mu_s(IJK,M) = (1.0d0-blend)*Mu_s_p(IJK) &
            + blend*Mu_s_v(IJK) + Mu_s_f(IJK)
            
!     bulk viscosity in Mth solids phase   (add to plastic part)
            
            LAMBDA_s_c(IJK,M)= Lambda_s_v(IJK)
            LAMBDA_s(IJK,M) = (1.0d0-blend)*LAMBDA_s_p(IJK) &
            + blend*Lambda_s_v(IJK) + Lambda_s_f(IJK)
            
            P_s_c(IJK,M) = P_s_v(IJK)
            P_s(IJK,M) = (1.0d0-blend)*P_s_p(IJK) + blend*P_s_v(IJK) &
            + P_s_f(IJK)        !add to P_s
     
            ALPHA_s(IJK,M) = (1.0d0-blend)*ALPHA_s_p(IJK) &
            + blend*ALPHA_s_v(IJK) + ALPHA_s_f(IJK)
     
 200     ENDDO

      ELSE                      ! Blending Stress

         Mu_s_c(:,M) = Mu_s_v(:)
         Mu_s(:,M) = Mu_s_p(:) + Mu_s_v(:) + Mu_s_f(:)
     
!     bulk viscosity in Mth solids phase   (add to plastic part)
         
         LAMBDA_s_c(:,M)= Lambda_s_v(:)
         LAMBDA_s(:,M) = LAMBDA_s_p(:) + Lambda_s_v(:) + Lambda_s_f(:)
         
         P_s_c(:,M) = P_s_v(:)
         P_s(:,M) = P_s_p(:) + P_s_v(:) + P_s_f(:) !add to P_s
     
      ENDIF                     ! Blending Stress
      
      RETURN
      END





      Subroutine friction_schaeffer(M, IER)
     
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      IMPLICIT NONE

!     Local Variables
!     Index
      INTEGER          IJK
!     Solids phase
      INTEGER          M, MM
!     Maximum value of solids viscosity in poise
      DOUBLE PRECISION MAX_MU_s
      PARAMETER (MAX_MU_s = 1000.D0)
     
!     Sum of all solids volume fractions
      DOUBLE PRECISION   SUM_EPS_CP
     
!     Factor in frictional-flow stress terms
      DOUBLE PRECISION qxP_s
!     Error index
      INTEGER          IER     
!     Include statement functions
     
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 's_pr2.inc'
     
      DO 200 IJK = ijkstart3, ijkend3       
     
         IF ( FLUID_AT(IJK) ) THEN
            
!     added closed pack, this has to be consistent with the normal frictional force
!     see for example source_v_s.f. Tardos Powder Tech. 92 (1997) 61-74 explains in
!     his equation (3) that solids normal and shear frictional stresses have to be 
!     treated consistently. --> sof May 24 2005. 
!     
            IF(EP_g(IJK) .LT. EP_g_blend_end(IJK)) THEN
!     part copied from source_v_s.f (sof)
               SUM_EPS_CP=0.0
               DO MM=1,MMAX
                  IF (CLOSE_PACKED(MM)) SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
               END DO
!     end of part copied
!     
!      P_star(IJK) = Neg_H(EP_g(IJK),EP_star_array(IJK))
!     
!     Frictional-flow stress tensor
!     
!-----------------------------------------------------------------------
!     Gray and Stiles (1988)
!     IF(Sin2_Phi .GT. SMALL_NUMBER) THEN
!     qxP_s = SQRT( (4. * Sin2_Phi) * I2_devD_s(IJK)
!     &                       + trD_s_C(IJK,M) * trD_s_C(IJK,M))
!     MU_s(IJK, M)     = P_star(IJK) * Sin2_Phi
!     &                         / (qxP_s + SMALL_NUMBER)
!     MU_s(IJK, M)     = MIN(MU_s(IJK, M), MAX_MU_s)
!     LAMBDA_s(IJK, M) = P_star(IJK) * F_Phi
!     &                         / (qxP_s + SMALL_NUMBER)
!     LAMBDA_s(IJK, M) = MIN(LAMBDA_s(IJK, M), MAX_MU_s)
!     ELSE
!     MU_s(IJK, M)     = ZERO
!     LAMBDA_s(IJK, M) = ZERO
!     ENDIF
!-----------------------------------------------------------------------
!     Schaeffer (1987)
!     
               qxP_s            = SQRT( (4.D0 * Sin2_Phi) * I2_devD_s(IJK))
               MU_s_p(IJK)     = P_star(IJK) * Sin2_Phi&
               / (qxP_s + SMALL_NUMBER) &
               *(EP_S(IJK,M)/SUM_EPS_CP) ! added by sof for consistency
                                ! with solids pressure treatment
               MU_s_p(IJK)     = MIN(MU_s_p(IJK), to_SI*MAX_MU_s)
               
               LAMBDA_s_p(IJK) = ZERO
               ALPHA_s_p(IJK)  = ZERO
               P_s_p(IJK)  = ZERO
     
!     when solving for the granular energy equation (PDE) setting theta = 0 is done 
!     in solve_granular_energy.f to avoid convergence problems. (sof)
               IF(.NOT.GRANULAR_ENERGY) THETA_m(IJK, M) = ZERO
            ENDIF
         ENDIF ! Fluid_at
!     end Schaeffer

 200  Continue                  ! outer IJK loop
      
      Return
      End




      Subroutine Gt_algebraic (M, IER)
!     
!     
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      USE trace
      IMPLICIT NONE

!     Local Variables
!     Index
      INTEGER          IJK
!     Solids phase
      INTEGER          M, MM
!     
!     Sum of all solids volume fractions
      DOUBLE PRECISION   SUM_EPS_CP
!     Error index
      INTEGER          IER     
!     
!     Coefficients of quadratic equation
      DOUBLE PRECISION aq, bq, cq
!     
!     Constant in equation for mth solids phase pressure
      DOUBLE PRECISION K_1m
!     
!     Constant in equation for mth solids phase bulk viscosity
      DOUBLE PRECISION K_2m
!     
!     Constant in equation for mth solids phase viscosity
      DOUBLE PRECISION K_3m
!     
!     Constant in equation for mth solids phase dissipation
      DOUBLE PRECISION K_4m
!     
!     Constant in Boyle-Massoudi stress term
      DOUBLE PRECISION K_5m
!     
!     Factor in frictional-flow stress terms
      DOUBLE PRECISION qxP_s
!     
!     Value of EP_s * SQRT( THETA )for Mth solids phase
!     continuum
      DOUBLE PRECISION EP_sxSQRTHETA
!     
!     Value of EP_s * EP_s * THETA for Mth solids phase
!     continuum
      DOUBLE PRECISION EP_s2xTHETA
!     
!
!     Function subroutines
!     
      DOUBLE PRECISION G_0
!
!     Include statement functions
!     
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'  


      DO 200 IJK = ijkstart3, ijkend3   
!     
         IF ( FLUID_AT(IJK) ) THEN

            IF(EP_g(IJK) .GE. EP_g_blend_start(IJK)) THEN
!     
!     
!     Calculate K_1m, K_2m, K_3m, K_4m
               K_1m = 2.D0 * (ONE + C_e) * RO_s(M) * G_0(IJK, M,M)
               K_3m = HALF * D_p(IJK,M) * RO_s(M) * (&
               ( (SQRT_PI / (3.D0*(3.D0 - C_e))) *&
               (HALF*(3d0*C_e+ONE) + 0.4D0*(ONE + C_e)*(3.D0*C_e - ONE)*&
               EP_s(IJK,M)*G_0(IJK, M,M)) ) + 8.D0*EP_s(IJK,M)&
               *G_0(IJK, M,M)*(ONE + C_e)/ (5.D0*SQRT_PI) )
               K_2m = 4.D0 * D_p(IJK,M) * RO_s(M) * (ONE + C_e) *&
               EP_s(IJK,M) * G_0(IJK, M,M) / (3.D0 * SQRT_PI) - 2.D0/3.D0 * K_3m
               K_4m = 12.D0 * (ONE - C_e*C_e) *&
               RO_s(M) * G_0(IJK, M,M) / (D_p(IJK,M) * SQRT_PI)
               aq   = K_4m*EP_s(IJK,M)
               bq   = K_1m*EP_s(IJK,M)*trD_s_C(IJK,M)
               cq   = -(K_2m*trD_s_C(IJK,M)*trD_s_C(IJK,M)&
               + 2.D0*K_3m*trD_s2(IJK,M))
!     
!     Boyle-Massoudi Stress term
!     
               IF(V_ex .NE. ZERO) THEN
                  K_5m = 0.4 * (ONE + C_e) * G_0(IJK, M,M) * RO_s(M) *&
                  ( (V_ex * D_p(IJK,M)) / (ONE - EP_s(IJK,M) * V_ex) )**2
                  bq   = bq + EP_s(IJK,M) * K_5m * (trM_s(IJK) + 2.D0 * trDM_s(IJK))
               ELSE
                  K_5m = ZERO
               ENDIF
!     
!     Calculate EP_sxSQRTHETA and EP_s2xTHETA
               EP_sxSQRTHETA = (-bq + SQRT(bq**2 - 4.D0 * aq * cq ))&
               / ( 2.D0 * K_4m )
               EP_s2xTHETA = EP_sxSQRTHETA * EP_sxSQRTHETA
               
               IF(EP_s(IJK,M) > SMALL_NUMBER)THEN
!     start      kapil&anuj 01/19/98
!     Find pseudo-thermal temperature in the Mth solids phase
                  THETA_m(IJK,M) = EP_s2xTHETA/(EP_s(IJK,M)*EP_s(IJK,M))
!     end      kapil&anuj 01/19/98
               ELSE
                  THETA_m(IJK,M) = ZERO
               ENDIF
     
!     Find pressure in the Mth solids phase
               P_s_v(IJK) = K_1m * EP_s2xTHETA
     
!     bulk viscosity in Mth solids phase
               LAMBDA_s_v(IJK) = K_2m * EP_sxSQRTHETA
     
!     shear viscosity in Mth solids phase
               MU_s_v(IJK) = K_3m * EP_sxSQRTHETA
     
!     Boyle-Massoudi stress coefficient
               ALPHA_s(IJK, M) = -K_5m * EP_s2xTHETA
            ENDIF

         Endif                  ! Fluid_at
 200  Continue                  ! outer IJK loop
      
      Return
      End



      Subroutine gt_pde (M, IER)

      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      USE toleranc
      USE turb
      USE drag
      IMPLICIT NONE

!     Local Variables
!     
!     Index
      INTEGER          IJK, I, J, K
!     
!     Solids phase
      INTEGER          M, MM
!     
!     Error index
      INTEGER          IER     
!     
!     Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION Mu_star, Kth, Kth_star
!     
!     defining parametrs for Simonin and Ahmadi models
      DOUBLE PRECISION Tau_12_st, Tau_2_c, Tau_2, Zeta_r, C_Beta
      DOUBLE PRECISION Sigma_c, Zeta_c, Omega_c, Zeta_c_2, C_mu, X_21, Nu_t
      DOUBLE PRECISION MU_2_T_Kin, Mu_2_Col, Kappa_kin, Kappa_Col
      DOUBLE PRECISION Tmp_Ahmadi_Const
      DOUBLE PRECISION DGA, C_d, Re
!     
!     Sum of ep_s * g_0
      DOUBLE PRECISION   SUM_EpsGo
      
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
!     (Same for conductivity)
      
!     dg0/dep
      DOUBLE PRECISION DG_0DNU,SRT

!
!     Function subroutines
!     
      DOUBLE PRECISION G_0
!     
!     Include statement functions
!     
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!     
!     
      DO 200 IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)       
     
         IF ( FLUID_AT(IJK) ) THEN
     
!     Defining a single particle drag coefficient (similar to one defined in drag_gs)
     
            Re = D_p(IJK,M)*VREL_array(IJK)*ROP_G(IJK)/(MU_G(IJK) + small_number)
            IF(Re .LE. 1000D0)THEN
               C_d = (24.D0/(Re+SMALL_NUMBER)) * (ONE + 0.15D0 * Re**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF
!     This is from Wen-Yu correlation, you can put here your own single particle drag
     
            DgA = 0.75D0 * C_d * VREL_array(IJK) * ROP_g(IJK) / D_p(IJK,M)
            IF(VREL_array(IJK) == ZERO) DgA = LARGE_NUMBER !for 1st iteration and 1st time step
     
!     Define some time scales and constants related to Simonin and Ahmadi models
     
            IF(SIMONIN .OR. AHMADI) THEN
               C_mu = 9.0D-02
!     particle relaxation time. For very dilute flows avoid singularity by
!     redefining the drag as single partilce drag
     
               IF(Ep_s(IJK,M) > DIL_EP_S .AND. F_GS(IJK,1) > small_number) THEN
                  Tau_12_st = Ep_s(IJK,M)*RO_s(M)/F_GS(IJK,1)
               ELSE             !for dilute flows, drag equals single particle drag law
                  Tau_12_st = RO_s(M)/DgA
               ENDIF            !for dilute flows
               
               
!     time scale of turbulent eddies
               Tau_1(ijk) = 3.d0/2.d0*C_MU*K_Turb_G(IJK)/(E_Turb_G(IJK)+small_number)
            ENDIF
     
!     Define some time scales and constants and K_12 related to Simonin model only	  
     
            IF(SIMONIN) THEN
!     This is Zeta_r**2 as defined by Simonin
               Zeta_r = 3.0d0 * VREL_array(IJK)**2 / (2.0d0*K_Turb_G(IJK)+small_number)
     
!     parameters for defining Tau_12: time-scale of the fluid turbulent motion
!     viewed by the particles (crossing trajectory effect)
     
               C_Beta = 1.8d0 - 1.35d0*Cos_Theta(IJK)**2
     
!     Lagrangian Integral time scale: Tau_12	    
               Tau_12(ijk) = Tau_1(ijk)/sqrt(ONE+C_Beta*Zeta_r)
     
!     Defining the inter-particle collision time
     
               IF(Ep_s(IJK,M) > DIL_EP_S) THEN
                  Tau_2_c = D_p(IJK,M)/(6.d0*Ep_s(IJK,M)*G_0(IJK,M,M) &
                  *DSQRT(16.d0*(Theta_m(ijk,m)+Small_number)/PI))
               ELSE             ! assign it a large number
                  Tau_2_c = LARGE_NUMBER
               ENDIF
     
               Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
     
!     Zeta_c: const. to be used in the K_2 Diffusion coefficient.
               Zeta_c  = (ONE+ C_e)*(49.d0-33.d0*C_e)/100.d0

               Omega_c = 3.d0*(ONE+ C_e)**2 *(2.d0*C_e-ONE)/5.d0

               Zeta_c_2= 2./5.*(ONE+ C_e)*(3.d0*C_e-ONE)

!     mixed time scale in the generalized Simonin theory (switch between dilute
!     and kinetic theory formulation of the stresses)
               Tau_2 = ONE/(2./Tau_12_st+Sigma_c/Tau_2_c)
     
!     The ratio of densities
               X_21 = Ep_s(IJK,M)*RO_s(M)/(EP_g(IJK)*RO_g(IJK))
     
!     The ratio of these two time scales.
               Nu_t =  Tau_12(ijk)/Tau_12_st
     
!     Definition of an "algebraic" form of of Simonin K_12 PDE. This is obtained
!     by equating the dissipation term to the exchange terms in the PDE and 
!     neglecting all other terms, i.e. production, convection and diffusion.
!     This works because Tau_12 is very small for heavy particles

               K_12(ijk) = Nu_t / (ONE+Nu_t*(ONE+X_21)) * &
               (2.d+0 *K_Turb_G(IJK) + 3.d+0 *X_21*theta_m(ijk,m))
!     Realizability Criteria         
               IF(K_12(ijk) > DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))) THEN
                  K_12(ijk) = DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))
               ENDIF
     
            ENDIF               ! for Simonin

!     This is added for consistency of multi-particles kinetic theory. Solids pressure,
!     viscosity and conductivity must be additive. Thus non-linear terms (eps^2) are 
!     corrected so the stresses of two identical solids phases are equal to those
!     of a single solids phase. sof June 15 2005.
     
            SUM_EpsGo = ZERO
            DO MM = 1, MMAX
               SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,M,MM)
            ENDDO 
     
!     Find pressure in the Mth solids phase
            P_s_v(IJK) = ROP_s(IJK,M)*(1d0+ 4.D0 * Eta *&
            SUM_EpsGo)*Theta_m(IJK,M)
     
!     implement Simonin (same as granular) and Ahmadi solids pressures
            IF(SIMONIN) THEN
               P_s_v(IJK) = P_s_v(IJK) ! no changes to solids pressure
            ELSE IF(AHMADI) THEN
               P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M) * ( (ONE + 4.0D0* &
               SUM_EpsGo ) + HALF*(ONE - C_e*C_e) )
            ENDIF

!     find bulk and shear viscosity
     
            Mu_s_v(IJK) = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(IJK,M)*RO_s(M))/96d0
            
            Mu_b_v(IJK) = (256d0*Mu_s_v(IJK)*EP_s(IJK,M)*SUM_EpsGo)&
            /(5d0*Pi)

!     added Ro_g = 0 for granular flows (no gas). sof Aug-02-2005 
            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN !sof modifications (May 20 2005)
               Mu_star = Mu_s_v(IJK)
               
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Mu_star = ZERO
               
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               
               
               Mu_star = RO_S(M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)* Mu_s_v(IJK)/ &
               (RO_S(M)*SUM_EpsGo*Theta_m(IJK,M) &
               + 2.0d0*SWITCH*DgA/RO_S(M)* Mu_s_v(IJK))
               
            ELSE
               Mu_star = RO_S(M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)*Mu_s_v(IJK)/ &
               (RO_S(M)*SUM_EpsGo*Theta_m(IJK,M)+ &
               (2d0*SWITCH*F_gs(IJK,M)*Mu_s_v(IJK)/(RO_S(M)*EP_s(IJK,M))) )
            ENDIF
     
!     shear viscosity in Mth solids phase  (add to frictional part)
            Mu_s_v(IJK) =&
            ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
            G_0(IJK,M,M)))*(ONE+1.6d0*Eta*SUM_EpsGo)&
            *(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
            SUM_EpsGo)+(0.6d0*Mu_b_v(IJK)*Eta))
     
!     implement Simonin and Ahmadi solids viscosities
            IF(SIMONIN) THEN
     
!     Defining Simonin solids turbulent Kinetic (MU_2_T_Kin) and collisional (Mu_2_Col)
!     viscosities
               MU_2_T_Kin = (2.0d0/3.0d0*K_12(ijk)*Nu_t + Theta_m(IJK,M) * &
               (ONE+ zeta_c_2*EP_s(IJK,M)*G_0(IJK,M,M)))*Tau_2
     
               Mu_2_Col = 8.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* (MU_2_T_Kin+ &
               D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
     
               Mu_b_v(IJK) = 5.d0/3.d0*EP_s(IJK,M)*RO_s(M)*Mu_2_Col
     
               Mu_s_v(IJK) = EP_s(IJK,M)*RO_s(M)*(MU_2_T_Kin + Mu_2_Col)
               
            ELSE IF(AHMADI) THEN

              IF(EP_s(IJK,M) < (ONE-EP_star_array(ijk))) THEN
                Tmp_Ahmadi_Const = &
                  ONE/(ONE+ Tau_1(ijk)/Tau_12_st * (ONE-EP_s(IJK,M)/(ONE-EP_star_array(ijk)))**3)
              ELSE
                Tmp_Ahmadi_Const = ONE
              ENDIF
     
!     Defining Ahmadi shear and bulk viscosities. Ahmadi coefficient 0.0853 in C_mu
!     was replaced by 0.1567 to include 3/2*sqrt(3/2) because K = 3/2 Theta_m
               Mu_s_v(IJK) = Tmp_Ahmadi_Const &
               *0.1045d0*(ONE/G_0(IJK,M,M)+3.2d0*EP_s(IJK,M)+12.1824d0*   &
               G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M))*D_p(IJK,M)*RO_s(M)*  &
               DSQRT(Theta_m(IJK,M))
     
!     This is a guess of what Mu_b might be by taking 5/3 of the collisional viscosity
!     contribution. In this case col. visc. is the eps^2 contribution to Mu_s_v(IJK). This
!     might be changed later if communications with Ahmadi reveals a diffrent appoach
               Mu_b_v(IJK) = 5.d0/3.d0* Tmp_Ahmadi_Const                  &
               *0.1045d0*(12.1824d0*G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M)) &
               *D_p(IJK,M)*RO_s(M)* DSQRT(Theta_m(IJK,M))
               
            ENDIF               !for simonin or ahmadi viscosity
            
            
            Kth=75d0*RO_s(M)*D_p(IJK,M)*DSQRT(Pi*Theta_m(IJK,M))/&
               (48d0*Eta*(41d0-33d0*Eta))
            
            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN ! sof modifications (May 20 2005)
               Kth_star=Kth
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Kth_star = ZERO
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               Kth_star = RO_S(M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)* Kth/ &
               (RO_S(M)*SUM_EpsGo*Theta_m(IJK,M) &
               + 1.2d0*SWITCH*DgA/RO_S(M)* Kth)
            ELSE
               Kth_star = RO_S(M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)*Kth/ &
               (RO_S(M)*SUM_EpsGo*Theta_m(IJK,M)+ &
               (1.2d0*SWITCH*F_gs(IJK,M)*Kth/(RO_S(M)*EP_s(IJK,M))) )
            ENDIF
     
!     granular conductivity in Mth solids phase
            Kth_s(IJK,M) = Kth_star/G_0(IJK,M,M)*(&
            ( ONE + (12d0/5.d0)*Eta*SUM_EpsGo )&
            * ( ONE + (12d0/5.d0)*Eta*Eta*(4d0*Eta-3d0)* SUM_EpsGo )&
            + (64d0/(25d0*Pi)) * (41d0-33d0*Eta) * (Eta*SUM_EpsGo)**2 )
     
!     implement Simonin and Ahmadi solids conductivities
            IF(SIMONIN) THEN

!     Defining Simonin's Solids Turbulent Kinetic diffusivity: Kappa
               Kappa_kin = (9.d0/10.d0*K_12(ijk)*Nu_t + 3.0D0/2.0D0 * &
               Theta_m(IJK,M)*(ONE+ Omega_c*EP_s(IJK,M)*G_0(IJK,M,M)))/&
               (9.d0/(5.d0*Tau_12_st) + zeta_c/Tau_2_c)
               
               Kappa_Col = 18.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* & 
               (Kappa_kin+ 5.d0/9.d0*D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
               
               Kth_s(IJK,M) =  EP_s(IJK,M)*RO_s(M)*(Kappa_kin + Kappa_Col)
     
            ELSE IF(AHMADI) THEN

!     Defining Ahmadi conductivity from his equation 42 in Cao and Ahmadi 1995 paper
!     note the constant 0.0711 is now 0.1306 because K = 3/2 theta_m
               Kth_s(IJK,M) = 0.1306D0*RO_s(M)*D_p(IJK,M)*(ONE+C_e**2)* &
               (ONE/G_0(IJK,M,M)+4.8D0*EP_s(IJK,M)+12.1184D0 &
               *EP_s(IJK,M)*EP_s(IJK,M)*G_0(IJK,M,M) )  &
               *DSQRT(Theta_m(IJK,M))
     
            ENDIF

            LAMBDA_S_V(IJK) = Eta*Mu_b_v(IJK) - (2d0*Mu_s_v(IJK))/3d0
     
!     granular 'conductivity' in the Mth solids phase associated
!     with gradient in volume fraction
            
!--------------------------------------------------------------------
!     Kphi_s has been set to zero.  To activate the feature uncomment the
!     following lines and also the lines in source_granular_energy.
            Kphi_s(IJK,M) = ZERO
!     &            (Kth_star/(G_0(IJK,M,M)))*(12d0/5.)*Eta*(Eta-1.)*
!     &            (2.*Eta-1.)*(1.+(12d0/5.)*Eta*EP_s(IJK,M)*
!     &            G_0(IJK,M,M))*(EP_s(IJK,M)*
!     &            DG_0DNU(EP_s(IJK,M))
!     &            + 2*G_0(IJK,M,M))*Theta_m(IJK,M)
!--------------------------------------------------------------------
     
!     Boyle-Massoudi stress coefficient
            ALPHA_s_v(IJK) = ZERO

         Endif                  ! Fluid_at
 200  Continue                  ! outer IJK loop
      

      Return
      End
!----------------------------------------------- 



!     JEG Added 
!     Implement kinetic theory of Garzo and Dufty (1999) 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      Subroutine gt_pde_gd_99 (M, IER)
!
!-----------------------------------------------
!     Modules 
!-----------------------------------------------  
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      USE toleranc
      USE turb
      USE drag
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Index
      INTEGER          IJK, I, J, K
!     
!                      Solids phase
      INTEGER          M, L 
!     
!     Error index
      INTEGER          IER     
!     
!     Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION Mu_star, Kth_star
!
!
      DOUBLE PRECISION DGA, C_d, Re

      DOUBLE PRECISION D_PM, M_PM, NU_PM
      DOUBLE PRECISION c_star, zeta0_star, nu_eta_star, &
                       gamma_star, eta_k_star, eta_star, eta0, &
                       kappa0, nu_kappa_star, kappa_k_star, &
                       qmu_k_star, qmu_star, kappa_star, press_star

!
!     SWITCH enables us to turn on/off the modification to the
!     particulate phase viscosity. If we want to simulate gas-particle
!     flow then SWITCH=1 to incorporate the effect of drag on the
!     particle viscosity. If we want to simulate granular flow
!     without the effects of an interstitial gas, SWITCH=0.
!     (Same for conductivity)
      
!     dg0/dep
      DOUBLE PRECISION DG_0DNU,SRT
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION G_0
!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------   
     
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

          IF ( FLUID_AT(IJK) ) THEN
    
!              Defining a single particle drag coefficient (similar to one defined in drag_gs)

               Re = D_p(IJK,M)*VREL_array(IJK)*ROP_G(IJK)/&
                    (MU_G(IJK) + small_number)
               IF(Re .LE. 1000D0)THEN
                    C_d = (24.D0/(Re+SMALL_NUMBER)) * &
                         (ONE + 0.15D0 * Re**0.687D0)
               ELSE
                    C_d = 0.44D0
               ENDIF

!              This is from Wen-Yu correlation, you can put here your own single particle drag
    
               DgA = 0.75D0 * C_d * VREL_array(IJK) * ROP_g(IJK) / D_p(IJK,M)
               IF (VREL_array(IJK) == ZERO) THEN
                    DgA = LARGE_NUMBER     ! for 1st iteration and 1st time step
               ENDIF

    
!              Pressure/Viscosity/Bulk Viscosity
!              Note: k_boltz = M_PM
!-----------------------------------
               D_PM = D_P(IJK,M)
               M_PM = (PI/6.d0)*D_PM**3 * RO_S(M)
               NU_PM = ROP_S(IJK,M)/M_PM

     
!              Find pressure in the Mth solids phase
               press_star = 1.d0 + 2.d0*(1.d0+C_E)*EP_s(IJK,M)*G_0(IJK,M,M)

!              n*k_boltz = n*m = ep_s*ro_s
               P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M)*press_star
    

!              find bulk and shear viscosity
               c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
                    / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

               zeta0_star = (5.d0/12.d0)*G_0(IJK,M,M)*(1.d0 - C_E*C_E) &
                    * (1.d0 + (3.d0/32.d0)*c_star)

               nu_eta_star = G_0(IJK,M,M)*(1.d0 - 0.25d0*(1.d0-C_E)*(1.d0-C_E)) &
                    * (1.d0-(c_star/64.d0))

               gamma_star = (4.d0/5.d0)*(32.d0/PI)*EP_s(IJK,M)*EP_S(IJK,M) &
                    * G_0(IJK,M,M)*(1.d0+C_E) * (1.d0 - (c_star/32.d0))

               eta_k_star = (1.d0 - (2.d0/5.d0)*(1.d0+C_E)*(1.d0-3.d0*C_E) &
                    * EP_s(IJK,M)*G_0(IJK,M,M) ) / (nu_eta_star - 0.5d0*zeta0_star)

               eta_star = eta_k_star*(1.d0 + (4.d0/5.d0)*EP_s(IJK,M)*G_0(IJK,M) &
                    * (1.d0+C_E) ) + (3.d0/5.d0)*gamma_star

               eta0 = 5.0d0*M_PM*DSQRT(Theta_m(IJK,M)/PI) / (16.d0*D_PM*D_PM)
    
!              added Ro_g = 0 for granular flows (no gas). 
               IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN 
                    Mu_star = eta0
               
               ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                    Mu_star = ZERO
               
               ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
                    Mu_star = RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M)*eta0 / &
                         ( RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M) + &
                         2.d0*DgA*eta0/RO_S(M) )
               
               ELSE
                    Mu_star = RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M)*eta0 / &
                         ( RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M) + &
                         (2.d0*F_gs(IJK,M)*eta0/(RO_S(M)*EP_s(IJK,M))) )
               ENDIF

!              shear viscosity in Mth solids phase  (add to frictional part)
               Mu_s_v(IJK) = Mu_star * eta_star

               Mu_b_v(IJK) = Mu_star * gamma_star 

!              second viscosity
               LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)

 
!              Granular Conductivity
!-----------------------------------
               kappa0 = (15.d0/4.d0)*eta0

               nu_kappa_star = (G_0(IJK,M,M)/3.d0)*(1.d0+C_E) * ( 1.d0 + &
                    (33.d0/16.d0)*(1.d0-C_E) + ((19.d0-3.d0*C_E)/1024.d0)*c_star)
!              nu_mu_star = nu_kappa_star

               kappa_k_star = (2.d0/3.d0)*(1.d0 + 0.5d0*(1.d0+press_star)*c_star + &
                    (3.d0/5.d0)*EP_s(IJK,M)*G_0(IJK,M,M)*(1.d0+C_E)*(1.d0+C_E) * &
                    (2.d0*C_E - 1.d0 + ( 0.5d0*(1.d0+C_E) - 5.d0/(3*(1.d0+C_E))) * &
                     c_star ) ) / (nu_kappa_star - 2.d0*zeta0_star)

               kappa_star = kappa_k_star * (1.d0 + (6.d0/5.d0)*EP_s(IJK,M)* &
                    G_0(IJK,M,M)*(1.d0+C_E) ) + (256.d0/25.d0)*(EP_s(IJK,M)* &
                    EP_s(IJK,M)/PI)*G_0(IJK,M,M)*(1.d0+C_E)*(1.d0+(7.d0/32.d0)* &
                    c_star)

               IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN ! sof modifications (May 20 2005)
                    Kth_star= kappa0
               
               ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                    Kth_star = ZERO
               
               ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
                    Kth_star = RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M)*kappa0/ &
                         (RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M) + &
                         1.2d0*DgA*kappa0/RO_S(M) )

               ELSE
                    Kth_star = RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M)*kappa0/ &
                         (RO_S(M)*EP_s(IJK,M)*G_0(IJK,M,M)*Theta_m(IJK,M)+ &
                         (1.2d0*F_gs(IJK,M)*kappa0/(RO_S(M)*EP_s(IJK,M))) )
               ENDIF

!              granular conductivity in Mth solids phase
               Kth_s(IJK,M) = Kth_star * kappa_star
   
!              transport coefficient of the Mth solids phase associated
!              with gradient in volume fraction in heat flux
               qmu_k_star = 2.d0*( (1.d0+EP_s(IJK,M)*DG_0DNU(EP_s(IJK,M)))* &
                    zeta0_star*kappa_k_star + ( (press_star/3.d0) + (2.d0/3.d0)* &
                    EP_s(IJK,M)*(1.d0+C_E) * (G_0(IJK,M,M)+EP_s(IJK,M)* &
                    DG_0DNU(EP_s(IJK,M))) )*c_star - (4.d0/5.d0)*EP_s(IJK,M)* &
                    G_0(IJK,M,M)* (1.d0+(EP_s(IJK,M)/2.d0)*DG_0DNU(EP_s(IJK,M)))* &
                    (1.d0+C_E) * ( C_E*(1.d0-C_E)+0.25d0*((4.d0/3.d0)+C_E* &
                    (1.d0-C_E))*c_star ) ) / (2.d0*nu_kappa_star-3.d0*zeta0_star)

               qmu_star = qmu_k_star*(1.d0+(6.d0/5.d0)*EP_s(IJK,M)*G_0(IJK,M,M)*&
                    (1.d0+C_E) )


               IF (EP_S(IJK,M) .LT. SMALL_NUMBER) THEN
                    Kphi_s(IJK,M) = ZERO
               ELSE   
                    Kphi_s(IJK,M) = (Theta_m(IJK,M)*Kth_star/NU_PM)*qmu_star
               ENDIF

!              Boyle-Massoudi stress coefficient
               ALPHA_s_v(IJK) = ZERO

          ENDIF     ! Fluid_at
 200  Continue     ! outer IJK loop

      Return
      End
!----------------------------------------------- 



!     JEG Added 
!     Implement complete kinetic theory of Iddir & Arastoopour (2005) 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      Subroutine gt_pde_ia_nonep (M, IER)
!
!-----------------------------------------------
!     Modules 
!-----------------------------------------------  
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      USE toleranc
      USE turb
      USE drag
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Index
      INTEGER          IJK, I, J, K
!     
!                      Solids phase
      INTEGER          M, L 
!     
!                      Error index
      INTEGER          IER     
!     
!                      Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION Mu_star, Mu_s_dil, Kth_star, K_s_dil, XI_star
!
!                      variables for Iddir equipartition model
      DOUBLE PRECISION P_s_sum, P_s_MM, P_s_LM
      DOUBLE PRECISION MU_common_term, K_common_term
      DOUBLE PRECISION Mu_sM_sum, MU_s_MM, MU_s_LM, MU_sM_LM, MU_sL_LM
      DOUBLE PRECISION XI_sM_sum, XI_s_v
      DOUBLE PRECISION M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION Ap_lm, Dp_lm, R0p_lm, R1p_lm, R8p_lm, R9p_lm, Bp_lm,&
                       R5p_lm, R6p_lm, R7p_lm
      DOUBLE PRECISION K_s_sum, K_s_MM, K_s_LM
      DOUBLE PRECISION Re, C_d, DgA
!     
!     Sum of ep_s * g_0
      DOUBLE PRECISION   SUM_EpsGo
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION G_0
!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------     
     
      DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)       
     
          IF ( FLUID_AT(IJK) ) THEN
     
! Defining a single particle drag coefficient
! (similar to one defined in drag_gs)

               Re = D_p(IJK,M)*VREL_array(IJK)*ROP_G(IJK)/&
                    (MU_G(IJK) + small_number)
               IF(Re .LE. 1000D0)THEN
                    C_d = (24.D0/(Re+SMALL_NUMBER)) * &
                         (ONE + 0.15D0 * Re**0.687D0)
               ELSE
                    C_d = 0.44D0
               ENDIF

! This is from Wen-Yu correlation, you can put here your
! own single particle drag
     
               DgA = 0.75D0 * C_d * VREL_array(IJK) * ROP_g(IJK) / D_p(IJK,M)
               IF (VREL_array(IJK) == ZERO) THEN
                    DgA = LARGE_NUMBER     ! for 1st iteration and 1st time step
               ENDIF
  
! Added for concistancy of IA KT: 2 identical solids phases must yield same
! results as one solids phase. Both Mus and Kths are modified. 
! This is an ad-hoc modification as there are other possible ways of doing this.
! SWITCH_IA can be set to false in constant_mod to use standard IA theory.
! sof Dec 05 2006. 
   
               IF(SWITCH_IA) THEN
                  SUM_EpsGo = ZERO
                  DO L = 1, MMAX
                     SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,L)*G_0(IJK,M,L)
                  ENDDO 
               ELSE
                  SUM_EpsGo =  EP_s(IJK,M)*G_0(IJK,M,M) 
               ENDIF
     
               P_s_sum = ZERO
               Mu_sM_sum = ZERO
               XI_sM_sum = ZERO

               D_PM = D_P(IJK,M)
               M_PM = (PI/6.d0)*D_PM**3 * RO_S(M)
               NU_PM = ROP_S(IJK,M)/M_PM

               P_s_MM = NU_PM*Theta_m(IJK,M)

               MU_s_dil = (5.d0/96.d0)*D_PM* RO_S(M)*&
                          DSQRT(PI*Theta_m(IJK,M)/M_PM)

               IF(.NOT.SWITCH_IA .OR. RO_G(IJK) == ZERO) THEN 
                    Mu_star = MU_s_dil ! do nothing... granular flow
               ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
                    Mu_star = ZERO

               ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
                    Mu_star = MU_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                             (SUM_EpsGo + 2.0d0*DgA*MU_s_dil &
                             / (RO_S(M)**2 *(Theta_m(IJK,M)/M_PM)))
               ELSE
                    Mu_star = MU_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                      (SUM_EpsGo + 2.0d0*F_gs(IJK,M)*MU_s_dil &
                      / (RO_S(M)**2 *EP_s(IJK,M)*(Theta_m(IJK,M)/M_PM)))
               ENDIF
   
               MU_s_MM = (Mu_star/G_0(IJK,M,M))*&
                    (1.d0+(4.d0/5.d0)*(1.d0+C_E)*SUM_EpsGo)**2

               DO L = 1, MMAX

                    D_PL = D_P(IJK,L)
                    M_PL = (PI/6.d0)*D_PL**3 * RO_S(L)
                    MPSUM = M_PM + M_PL

                    DPSUMo2 = (D_PM+D_PL)/2.d0
                    NU_PL = ROP_S(IJK,L)/M_PL

                    IF ( L .eq. M) THEN
                         Ap_lm = MPSUM/(2.d0)
                         Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
                         R0p_lm = ONE/( Ap_lm**1.5 * Dp_lm**2.5 )
                         R1p_lm = ONE/( Ap_lm**1.5 * Dp_lm**3 )
 
                         P_s_LM = PI*(DPSUMo2**3 / 48.d0)*G_0(IJK,M,L)*&
                            (M_PM*M_PL/MPSUM)* (M_PM*M_PL)**1.5 *&
                            NU_PM*NU_PL*(1.d0+C_E)*R0p_lm*Theta_m(IJK,M)


                         MU_s_LM = DSQRT(PI)*( DPSUMo2**4 / 240d0 )*&
                             G_0(IJK,M,L)*(M_PL*M_PM/MPSUM)**2 *&
                             (M_PL*M_PM)**1.5 * NU_PM*NU_PL*&
                             (1.d0+C_E) * R1p_lm * DSQRT(Theta_m(IJK,M))

! This is Mu_i_1 as defined in eq (16) of Galvin document
                         MU_sM_ip(IJK,M,L) = (MU_s_MM + MU_s_LM)

! This is Mu_ii_2 as defined in eq (17) of Galvin document
                         MU_sL_ip(IJK,M,L) = MU_s_LM

! solids phase viscosity associated with the trace of
! solids phase M (eq. 18 from Galvin theory document)
                         XI_sM_ip(IJK,M,L) = (5.d0/3.d0)*MU_s_LM

! solids phase viscosity associated with the trace 
! of (sum of) all solids phases (eq. 19)
                         XI_sL_ip(IJK,M,L) = (5.d0/3.d0)*MU_s_LM
                    ELSE
 
                         Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*&
                                 Theta_m(IJK,M))/2.d0
                         Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                                 Theta_m(IJK,M) ))/(2.d0*MPSUM)
                         Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+&
                                 M_PL*Theta_m(IJK,L) ))/&
                                 (2.d0*MPSUM*MPSUM)
 
                         R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+ &
                              ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 *&
                              Dp_lm**3.5 ) )+&
                              ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 * &
                              Dp_lm**4.5 ) )
      
                         R1p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) ) )+ &
                              ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * &
                              Dp_lm**4 ) )+&
                              ( (30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 * &
                              Dp_lm**5 ) )
  
                         P_s_LM = PI*(DPSUMo2**3 / 48.d0)*G_0(IJK,M,L)*&
                              (M_PM*M_PL/MPSUM)* (M_PM*M_PL)**1.5 *&
                              NU_PM*NU_PL*(1.d0+C_E)*R0p_lm* &
                              (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

                         MU_common_term = DSQRT(PI)*( DPSUMo2**4 / 240d0 )*&
                              G_0(IJK,M,L)*(M_PL*M_PM/MPSUM)**2 *&
                              (M_PL*M_PM)**1.5 * NU_PM*NU_PL*&
                              (1.d0+C_E) * R1p_lm

                         MU_sM_LM = MU_common_term * Theta_m(IJK,M)**2 *&
                              Theta_m(IJK,L)**3

                         MU_sL_LM = MU_common_term * Theta_m(IJK,L)**2 *&
                              Theta_m(IJK,M)**3

! solids phase 'viscosity' associated with the divergence
! of solids phase M. defined in eq (16) of Galvin document   
                         MU_sM_ip(IJK,M,L) = MU_sM_LM

! solids phase 'viscosity' associated with the divergence
! of all solids phases. defined in eq (17) of Galvin document
                         MU_sL_ip(IJK,M,L) = MU_sL_LM

! solids phase viscosity associated with the trace of
! solids phase M
                         XI_sM_ip(IJK,M,L) = (5.d0/3.d0)*MU_sM_LM

! solids phase viscosity associated with the trace 
! of all solids phases 
                         XI_sL_ip(IJK,M,L) = (5.d0/3.d0)*MU_sL_LM
                    ENDIF

                    P_s_sum = P_s_sum + P_s_LM

                    MU_sM_sum = MU_sM_sum + MU_sM_ip(IJK,M,L)

                    XI_sM_sum = XI_sM_sum + XI_sM_ip(IJK,M,L)

               END DO

! Find the term proportional to the identity matrix
! (pressure in the Mth solids phase)
               P_s_v(IJK) = P_s_sum + P_S_MM

! Find the term proportional to the gradient in velocity
! of phase M  (shear viscosity in the Mth solids phase)
               MU_s_v(IJK) = MU_sM_sum + MU_sL_ip(IJK,M,M)

               XI_s_v = XI_sM_sum + XI_sL_ip(IJK,M,M)
     
! bulk viscosity in the Mth solids phase
               LAMBDA_s_v(IJK) = -(2.d0/3.d0)*Mu_s_v(IJK) + XI_s_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! find the granular conductivity in Mth solids phase      
               K_s_sum = ZERO

               K_s_dil = (75.d0/384.d0)*D_PM* RO_S(M)*&
                    DSQRT(PI*Theta_m(IJK,M)/M_PM)
               
               IF(.NOT.SWITCH_IA .OR. RO_G(IJK) == ZERO) THEN 
                    Kth_star = K_s_dil ! do nothing... granular flow
               ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
                    Kth_star = ZERO
        
               ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
                    Kth_star = K_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                         (SUM_EpsGo+ 1.2d0*DgA*K_s_dil &
                         / (RO_S(M)**2 *(Theta_m(IJK,M)/M_PM)))
               ELSE
                    Kth_star = K_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                         (SUM_EpsGo+ 1.2d0*F_gs(IJK,M)*K_s_dil &
                         / (RO_S(M)**2 *EP_s(IJK,M)*(Theta_m(IJK,M)/M_PM)))
               ENDIF

! Kth doesn't include the mass.      
              K_s_MM = (Kth_star/(M_PM*G_0(IJK,M,M)))*&  
                   (1.d0+(3.d0/5.d0)*(1.d0+C_E)*(1.d0+C_E)*SUM_EpsGo)**2

               DO L = 1, MMAX

                    D_PL = D_P(IJK,L)
                    M_PL = (PI/6.d0)*D_PL**3 *RO_S(L)
                    MPSUM = M_PM + M_PL

                    DPSUMo2 = (D_PM+D_PL)/2.d0
                    NU_PL = ROP_S(IJK,L)/M_PL

                    IF ( L .eq. M) THEN 

! solids phase 'conductivity' associated with the
! difference in velocity. again these terms cancel when
! added together so do not explicity include them
! in calculations
                         Kvel_s_ip(IJK,M,L) = ZERO

                         !     K_common_term*NU_PM*NU_PL*&
                         !     (3.d0*PI/10.d0)*R0p_lm*Theta_m(IJK,M)

! solids phase 'conductivity' associated with the 
! difference in the gradient in number densities.
! again these terms cancel so do not explicity include
! them in calculations
                         Knu_sL_ip(IJK,M,L) = ZERO

                         !     K_common_term*NU_PM*&
                         !     (PI*DPSUMo2/6.d0)*R1p_lm*(Theta_m(IJK,M)**(3./2.))

                         Knu_sM_ip(IJK,M,L) = ZERO

                         !     K_common_term*NU_PL*&
                         !     (PI*DPSUMo2/6.d0)*R1p_lm*(Theta_m(IJK,M)**(3./2.))

                         K_s_sum = K_s_sum + K_s_MM

                    ELSE

                         Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*&
                              Theta_m(IJK,M))/2.d0
                         Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                              Theta_m(IJK,M) ))/(2.d0*MPSUM)
                         Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+&
                              M_PL*Theta_m(IJK,L) ))/&
                              (2.d0*MPSUM*MPSUM)

                         R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+&
                              ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 *&
                              Dp_lm**3.5 ) )+&
                              ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 *&
                              Dp_lm**4.5 ) )
                             
                         R1p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) ) )+ &
                              ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 *&
                              Dp_lm**4 ) )+&
                              ( (30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 *&
                              Dp_lm**5 ) )
                         
                         R5p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**3 ) )+ &
                              ( (5.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * &
                              Dp_lm**4 ) )+&
                              ( (14.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5 ) )
                         
                         R6p_lm = ( 1.d0/( Ap_lm**3.5 * Dp_lm**3 ) )+ &
                              ( (7.d0*Bp_lm*Bp_lm)/( Ap_lm**4.5 * Dp_lm**4 ) )+&
                              ( (126.d0*Bp_lm**4)/( 5.d0*Ap_lm**5.5 * Dp_lm**5 ) )
                         
                         R7p_lm = ( 3.d0/( 2.d0*Ap_lm**2.5 * Dp_lm**4 ) )+ &
                              ( (10.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**5 ) )+&
                              ( (35.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**6 ) )
                         
                         R8p_lm = ( 1.d0/( 2.d0*Ap_lm**1.5 * Dp_lm**4 ) )+ &
                              ( (6.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**5 ) )+&
                              ( (25.d0*Bp_lm**4)/( Ap_lm**3.5 * Dp_lm**6 ) )
                         
                         R9p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**3 ) )+ &
                              ( (15.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**4 ) )+&
                              ( (70.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5 ) )

                         K_common_term = DPSUMo2**3 * M_PL*M_PM/(2.d0*MPSUM)*&
                              (1.d0+C_E)*G_0(IJK,M,L) * (M_PM*M_PL)**1.5

! solids phase 'conductivity' associated with the 
! gradient in granular temperature of species M
                    K_s_LM = - K_common_term*NU_PM*NU_PL*(&
                        ((DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)+&
                        ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                        (3.d0/2.d0)*R1p_lm)-(&
                        ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM/8.d0)*Bp_lm*R6p_lm)+&
                        ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                        8.d0)*M_PM*R9p_lm)+&
                        ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
                        M_PL*Bp_lm*R7p_lm)+&
                        ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                        2.d0)*(M_PM/MPSUM)**2 * M_PL*R8p_lm)+&
                        ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
                        R9p_lm)-&
                        ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
                        (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*Theta_M(IJK,L) )*&
                        (Theta_M(IJK,M)**2 * Theta_M(IJK,L)**3)

! solids phase 'conductivity' associated with the 
! gradient in granular temperature of species L
                    Kth_sL_ip(IJK,M,L) = K_common_term*NU_PM*NU_PL*(&
                         (-(DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)-&
                         ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                         (3.d0/2.d0)*R1p_lm)+(&
                         ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL/8.d0)*Bp_lm*R6p_lm)+&
                         ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                         8.d0)*M_PL*R9p_lm)+&
                         ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
                         M_PM*Bp_lm*R7p_lm)+&
                         ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                         2.d0)*(M_PM/MPSUM)**2 * M_PM*R8p_lm)+&
                         ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
                         R9p_lm)+&
                         ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
                         (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*Theta_M(IJK,M) )*&
                         (Theta_M(IJK,L)**2 * Theta_M(IJK,M)**3)

! solids phase 'conductivity' associated with the
! difference in velocity
                    Kvel_s_ip(IJK,M,L) = K_common_term*NU_PM*NU_PL*&
                         (M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(3.d0*PI/10.d0)*&
                         R0p_lm * (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

! solids phase 'conductivity' associated with the 
! difference in the gradient in number densities
                         
                    Knu_sL_ip(IJK,M,L) = K_common_term*(&
                         ((DPSUMo2*DSQRT(PI)/16.d0)*Bp_lm*R5p_lm)+&
                         ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                         R1p_lm) )* (Theta_m(IJK,M)*Theta_m(IJK,L))**3

! to avoid recomputing Knu_sL_ip, sof.
                    Knu_sM_ip(IJK,M,L) = NU_PL * Knu_sL_ip(IJK,M,L)

                    Knu_sL_ip(IJK,M,L) = NU_PM * Knu_sL_ip(IJK,M,L)

                    K_s_sum = K_s_sum + K_s_LM
                    ENDIF

               ENDDO

!              granular conductivity in Mth solids phase
               Kth_s(IJK,M) = K_s_sum

!              Boyle-Massoudi stress coefficient
               ALPHA_s_v(IJK) = ZERO

          ENDIF     ! Fluid_at
 200  Continue     ! outer IJK loop

      Return
      End
!----------------------------------------------- 




      Subroutine Friction_princeton(M, IER)
!     
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      USE visc_s
      USE physprop
      USE run
      USE constant
      USE trace
      IMPLICIT NONE

!     Local Variables
!     Index
      INTEGER          IJK
!
!     Solids phase
      INTEGER          M, MM
!
!     Used to compute frictional terms
      DOUBLE PRECISION Chi, Pc, Mu_zeta,Phin,PfoPc, N_Pff
!
      DOUBLE PRECISION ZETA
!     
!     Sum of all solids volume fractions
      DOUBLE PRECISION   SUM_EPS_CP
!
!     Error index
      INTEGER          IER     

!     Function subroutines
      DOUBLE PRECISION G_0

!     Include statement functions
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      DO 200 IJK = ijkstart3, ijkend3       
     
         IF ( FLUID_AT(IJK) ) THEN

!     close_packed was added for concistency with the Schaeffer model
!     I'm also extending this model in the case where more than 1 solids
!     phase are used (not completed yet). sof May24 2005.

               IF (EP_g(IJK) .LT. (ONE-eps_f_min)) THEN
     
!     part copied from source_v_s.f (sof)
                  SUM_EPS_CP=0.0
                  DO MM=1,SMAX
                     IF (CLOSE_PACKED(MM)) SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
                  END DO
!     end of part copied
     
                  IF (SAVAGE.EQ.1) THEN !form of Savage
                     Mu_zeta =&
                     ((2d0+ALPHA)/3d0)*((Mu_s_v(IJK)/(Eta*(2d0-Eta)*&
                     G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                     G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                     EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b_v(IJK)*Eta))
                     
                     ZETA =&
                     ((48d0*Eta*(1d0-Eta)*RO_s(M)*EP_s(IJK,M)*&
                     EP_s(IJK,M)*G_0(IJK,M,M)*&
                     (Theta_m(IJK,M)**1.5d0))/&
                     (SQRT_Pi*D_p(IJK,M)*2d0*Mu_zeta))**0.5d0
                     
                  ELSEIF (SAVAGE.EQ.0) THEN !S:S form
                     ZETA = (SMALL_NUMBER +&
                     trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                     trD_s_C(IJK,M))/3.d0))**0.5d0
                     
                  ELSE          !combined form
                     ZETA = ((Theta_m(IJK,M)/(D_p(IJK,M)*D_p(IJK,M))) +&
                     (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                     trD_s_C(IJK,M))/3.d0)))**0.5d0
                     
                  ENDIF
                  
                  IF ((ONE-EP_G(IJK)) .GT. (ONE-ep_star_array(ijk))) THEN
                     Pc = 1d25*(((ONE-EP_G(IJK))- (ONE-ep_star_array(ijk)))**10d0)
                  ELSE
                     Pc = Fr*(((ONE-EP_G(IJK)) - EPS_f_min)**N_Pc)/&
                     (((ONE-ep_star_array(ijk)) - (ONE-EP_G(IJK)) +&
                     SMALL_NUMBER)**D_Pc)
                  ENDIF
 
                  IF (trD_s_C(IJK,M) .GE. ZERO) THEN
                     N_Pff = DSQRT(3d0)/(2d0*Sin_Phi) !dilatation
                  ELSE
                     N_Pff = N_Pf !compaction
                  ENDIF
                  
                  IF ((trD_s_C(IJK,M)/(ZETA*N_Pff*DSQRT(2d0)&
                  *Sin_Phi)) .GT. 1d0) THEN
                    P_s_f(IJK) =ZERO
                    PfoPc = ZERO
                  
                  ELSEIF(trD_s_C(IJK,M) == ZERO) THEN
                    
                    P_s_f(IJK) = Pc
                    PfoPc = ONE
                  
                  ELSE
                  
                    P_s_f(IJK) = Pc*(1d0 - (trD_s_C(IJK,M)/(ZETA&
                    *N_Pff*DSQRT(2d0)*Sin_Phi)))**(N_Pff-1d0)
                  
                    PfoPc = (1d0 - (trD_s_C(IJK,M)/(ZETA&
                    *N_Pff*DSQRT(2d0)*Sin_Phi)))**(N_Pff-1d0)
                 ENDIF
               
                 Chi = DSQRT(2d0)*P_s_f(IJK)*Sin_Phi*(N_Pff - (N_Pff-1d0)*&
                 (PfoPc)**(1d0/(N_Pff-1d0)))
               
                 IF (Chi < ZERO) THEN
                    P_s_f(IJK) = Pc*((N_Pff/(N_Pff-1d0))**(N_Pff-1d0))
                    Chi = ZERO
                 ENDIF
               
                 Mu_s_f(IJK) = Chi/(2d0*ZETA)
                 Lambda_s_f(IJK) = - 2d0*Mu_s_f(IJK)/3d0
     
!     modification of the stresses in case of more than one solids phase are used (sof)
                 P_s_f(IJK) = P_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                 Mu_s_f(IJK) = Mu_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                 Lambda_s_f(IJK) = Lambda_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)

               
               ENDIF

         Endif                  ! Fluid_at
 200  Continue                  ! outer IJK loop
      

      Return
      End




      Subroutine add_shear(M) 
!
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      IMPLICIT NONE

!     Local Variables
!     Index
      INTEGER          IJK
!     Solids phase
      INTEGER          M
     
!     Include statement functions
      INCLUDE 'function.inc'

!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3         
         IF (FLUID_AT(IJK)) THEN  
            V_s(ijk,m)=V_s(IJK,m)+VSH(IJK)
         END IF
      ENDDO

      Return
      End                       ! add_shear



      Subroutine remove_shear(M)
!
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      IMPLICIT NONE

!     Local Variables
!     Index
      INTEGER          IJK
!     Solids phase
      INTEGER          M
     
!     Include statement functions
      INCLUDE 'function.inc'

!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
            V_s(IJK,m)=V_s(IJK,m)-VSH(IJK)
         ENDIF
      ENDDO

      Return
      End                       ! remove_shear




      Subroutine init_mu_s (M,IER)
!     
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE drag
      USE run
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE is
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      Use vshear
      USE compar
      USE sendrecv
      IMPLICIT NONE
      
!     Function subroutines
      DOUBLE PRECISION G_0
      double precision calc_ep_star

!     Local Variables
!     Error index
      INTEGER          IER
!     
!     Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_s(3,3), D_sl(3,3)
!     
!     U_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_s_N, Usl_N
!     
!     U_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_s_S, Usl_S
!     
!     U_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_s_T, Usl_T
!     
!     U_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_s_B, Usl_B
!     
!     U_s at the center of the THETA cell-(i, j, k)
!     Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_s_C, Usl_C
!     
!     V_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_s_E, Vsl_E
!     
!     V_s at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_s_W, Vsl_W
!     
!     V_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_s_T, Vsl_T
!     
!     V_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_s_B, Vsl_B
!     
!     W_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_s_E, Wsl_E
!     
!     W_s at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_s_W, Wsl_W
!     
!     W_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_s_N, Wsl_N
!     
!     W_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_s_S, Wsl_S
!     
!     W_s at the center of the THETA cell-(i, j, k).
!     Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_s_C, Wsl_C
!     
!     Cell center value of solids and gas velocities 
      DOUBLE PRECISION USCM, UGC, VSCM, VGC, WSCM, WGC 
!     
!     Local DO-LOOP counters and phase index
      INTEGER          I1, I2, MM
!     
!     d(EP_sm)/dX
      DOUBLE PRECISION DEP_soDX
!     
!     d(EP_sm)/dY
      DOUBLE PRECISION DEP_soDY
!     
!     d(EP_sm)/XdZ
      DOUBLE PRECISION DEP_soXDZ
!     
!     Solids volume fraction gradient tensor
      DOUBLE PRECISION M_s(3,3)
!     
!     Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,&
      IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,&
      IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP,&
      IJMKM, IJPKM
!     
!     Solids phase
      INTEGER          M, L
!     
!     Shear related reciprocal time scale
      DOUBLE PRECISION SRT
     
!     Include statement functions
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'

!     !$omp  parallel do &
!     !$omp& private(IMJPK, I, J, K, IJK,  IMJK, IPJK, IJMK, IJPK, IJKM, &
!     !$omp&  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM, &
!     !$omp&  U_s_N, U_s_S, U_s_T, U_s_B, V_s_E, V_s_W, V_s_T, V_s_B, W_s_N, &
!     !$omp&  W_s_S, W_s_E, W_s_W, U_s_C, W_s_C, D_s, I2_devD_s, trD_s_C, &
!     !$omp&  qxP_s, trD_s2, K_1m, K_2m, K_3m, K_4m, K_5m, aq, bq, cq, &
!     !$omp&  DEP_soDX, DEP_soDY, DEP_soXDZ, M_s, I1, I2, &
!     !$omp&  KTH_STAR,KTH,CHI,PFOPC,PC,ZETA,MU_ZETA,PF,&
!     !$omp&  MU_STAR,MU_B,MU,M,IJPKM,IJMKM,IJMKP,IPJMK,IPJKM,IMJKM,IMJKP,IMJMK, &
!     !$omp&  EP_sxSQRTHETA, EP_s2xTHETA )  
      
      Mu_s(:,M)     = ZERO
      LAMBDA_s(:,M) = ZERO
      ALPHA_s(:,M)  = ZERO
      Mu_s_v(:)     = ZERO
      Mu_b_v(:)     = ZERO
      LAMBDA_s_v(:) = ZERO
      ALPHA_s_v(M)  = ZERO
      Mu_s_p(:)     = ZERO
      LAMBDA_s_p(:) = ZERO
      ALPHA_s_p(M)  = ZERO
      Mu_s_f(:)     = ZERO
      LAMBDA_s_f(:) = ZERO
      ALPHA_s_f(M)  = ZERO

      P_s(:,M) = ZERO
      P_s_v(:) = ZERO
      P_s_f(:) = ZERO
      P_s_p(:) = ZERO

      IF (SHEAR) SRT=(2d0*V_sh/XLENGTH)

      DO 200 IJK = ijkstart3, ijkend3       
     
         IF ( FLUID_AT(IJK) ) THEN
     
!------------------------------------------------------------------------
!     CALL SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,
!     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
!     &                       IJKB, IJKT, IM, JM, KM)
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IM = Im1(I)
            JM = Jm1(J)
            KM = Km1(K)
            IJKW  = WEST_OF(IJK)
            IJKE  = EAST_OF(IJK)
            IJKS  = SOUTH_OF(IJK)
            IJKN  = NORTH_OF(IJK)
            IJKB  = BOTTOM_OF(IJK)
            IJKT  = TOP_OF(IJK)
            IMJK  = IM_OF(IJK)
            IPJK  = IP_OF(IJK)
            IJMK  = JM_OF(IJK)
            IJPK  = JP_OF(IJK)
            IJKM  = KM_OF(IJK)
            IJKP  = KP_OF(IJK)
            IMJPK = IM_OF(IJPK)
            IMJMK = IM_OF(IJMK)
            IMJKP = IM_OF(IJKP)
            IMJKM = IM_OF(IJKM)
            IPJKM = IP_OF(IJKM)
            IPJMK = IP_OF(IJMK)
            IJMKP = JM_OF(IJKP)
            IJMKM = JM_OF(IJKM)
            IJPKM = JP_OF(IJKM)

            U_s_N = AVG_Y(                                   & !i, j+1/2, k
            AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
            AVG_X_E(U_s(IMJPK, M), U_s(IJPK, M), I), J)

            U_s_S = AVG_Y(                                   & !i, j-1/2, k
            AVG_X_E(U_s(IMJMK, M), U_s(IJMK, M), I),&
            AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), JM)

            U_s_T = AVG_Z(                                 & !i, j, k+1/2
            AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
            AVG_X_E(U_s(IMJKP, M), U_s(IJKP, M), I), K)

            U_s_B = AVG_Z(                                 & !i, j, k-1/2
            AVG_X_E(U_s(IMJKM, M), U_s(IJKM, M), I),&
            AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), KM)

            IF (SHEAR)  THEN
               
               V_s_E = AVG_X(                                 & !i+1/2, j, k
               AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
               AVG_Y_N((V_s(IPJMK, M)-VSH(IPJMK)+VSH(IJMK)&
               +SRT*1d0/oDX_E(I)),&
               (V_s(IPJK, M)-VSH(IPJK)+VSH(IJK)&
               +SRT*1d0/oDX_E(I))), I)

               V_s_W = AVG_X(                                 & !i-1/2, j, k
               AVG_Y_N((V_s(IMJMK, M)-VSH(IMJMK)+VSH(IJMK)&
               -SRT*1d0/oDX_E(IM1(I))),&
               (V_s(IMJK, M)-VSH(IMJK)+VSH(IJK)&
               -SRT*1d0/oDX_E(IM1(I)))),&
               AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM)

            ELSE
               V_s_E = AVG_X(                                 & !i+1/2, j, k
               AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
               AVG_Y_N(V_s(IPJMK, M), V_s(IPJK, M)), I )

               V_s_W = AVG_X(                                 & !i-1/2, j, k
               AVG_Y_N(V_s(IMJMK, M), V_s(IMJK, M)),&
               AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM )
            END IF

            V_s_T = AVG_Z(                                 & !i, j, k+1/2
            AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
            AVG_Y_N(V_s(IJMKP, M), V_s(IJKP, M)), K )

            V_s_B = AVG_Z(                                 & !i, j, k-1/2
            AVG_Y_N(V_s(IJMKM, M), V_s(IJKM, M)),&
            AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), KM )

            W_s_N = AVG_Y(                                 & !i, j+1/2, k
            AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
            AVG_Z_T(W_s(IJPKM, M), W_s(IJPK, M)), J )

            W_s_S = AVG_Y(                                 & !i, j-1/2, k
            AVG_Z_T(W_s(IJMKM, M), W_s(IJMK, M)),&
            AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), JM )

            W_s_E = AVG_X(                                 & !i+1/2, j, k
            AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
            AVG_Z_T(W_s(IPJKM, M), W_s(IPJK, M)), I)

            W_s_W = AVG_X(                                 & !i-1/2, j, k
            AVG_Z_T(W_s(IMJKM, M), W_s(IMJK, M)),&
            AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), IM )
     
            IF(CYLINDRICAL) THEN
               U_s_C = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I) !i, j, k
               W_s_C = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)) !i, j, k
            ELSE
               U_s_C = ZERO
               W_s_C = ZERO
            ENDIF
     
!     Check for IS surfaces and modify solids velocity-comp accordingly
            IF(ANY_IS_DEFINED) THEN
              IF(IS_AT_N(IJK) .AND. .NOT.WALL_AT(IJPK)) U_s_N = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
              IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) U_s_S = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
              IF(IS_AT_T(IJK) .AND. .NOT.WALL_AT(IJKP)) U_s_T = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
              IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) U_s_B = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
              IF(IS_AT_E(IJK) .AND. .NOT.WALL_AT(IPJK)) V_s_E = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
              IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) V_s_W = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
              IF(IS_AT_T(IJK) .AND. .NOT.WALL_AT(IJKP)) V_s_T = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
              IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) V_s_B = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
              IF(IS_AT_N(IJK) .AND. .NOT.WALL_AT(IJPK)) W_s_N = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
              IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) W_s_S = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
              IF(IS_AT_E(IJK) .AND. .NOT.WALL_AT(IPJK)) W_s_E = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
              IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) W_s_W = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
            ENDIF
     
!     Find components of Mth solids phase continuum strain rate
!     tensor, D_s, at center of THETA cell-(i, j, k)
            D_s(1,1) = ( U_s(IJK,M) - U_s(IMJK,M) ) * oDX(I)

            D_s(1,2) = HALF * ( (U_s_N - U_s_S) * oDY(J) +&
            (V_s_E - V_s_W) * oDX(I) )

            D_s(1,3) = HALF * ( (W_s_E - W_s_W) * oDX(I) +&
            (U_s_T - U_s_B) * (oX(I)*oDZ(K)) - W_s_C * oX(I) )

            D_s(2,1) = D_s(1,2)

            D_s(2,2) = ( V_s(IJK,M) - V_s(IJMK,M) ) * oDY(J)

            D_s(2,3) = HALF * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +&
            (W_s_N - W_s_S) * oDY(J) )

            D_s(3,1) = D_s(1,3)

            D_s(3,2) = D_s(2,3)

            D_s(3,3) = ( W_s(IJK,M) - W_s(IJKM,M) ) * (oX(I)*oDZ(K)) +&
            U_s_C * oX(I)

!     Calculate the trace of D_s
            trD_s_C(IJK,M) = D_s(1,1) + D_s(2,2) + D_s(3,3)

!     Calculate trace of the square of D_s
            trD_s2(IJK,M) = 0.D0 !Initialize the totalizer
            DO 20 I1 = 1,3
               DO 10 I2 = 1,3
                  trD_s2(IJK,M) = trD_s2(IJK,M) + D_s(I1,I2)*D_s(I1,I2)
 10            CONTINUE
 20         CONTINUE
!     use this fact to prevent underflow during theta calculation
             if(trD_s2(IJK,M) == zero)trD_s_C(IJK,M) = zero 


!     The trace of D_sm dot D_sl is required in the implementation of
!     Iddir's (2004) kinetic theory 
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               DO L = 1,MMAX
                  IF (L .ne. M) THEN
                     IF (L > M) THEN !done because trD_s2_ip(IJK,M,L) is symmetric, sof.
                         Usl_N = AVG_Y(&                    !i, j+1/2, k
                              AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I),&
                              AVG_X_E(U_s(IMJPK, L), U_s(IJPK, L), I), J)
      
                         Usl_S = AVG_Y(&                    !i, j-1/2, k
                              AVG_X_E(U_s(IMJMK, L), U_s(IJMK, L), I),&
                              AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I), JM)
 
                         Usl_T = AVG_Z(&                    !i, j, k+1/2
                              AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I),&
                              AVG_X_E(U_s(IMJKP, L), U_s(IJKP, L), I), K)

                         Usl_B = AVG_Z(&                    !i, j, k-1/2
                              AVG_X_E(U_s(IMJKM, L), U_s(IJKM, L), I),&
                              AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I), KM)

                         IF (SHEAR)  THEN
                 
                              Vsl_E = AVG_X(&               !i+1/2, j, k
                                   AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)),&
                                   AVG_Y_N((V_s(IPJMK, L)-VSH(IPJMK)+&
                                   VSH(IJMK)+SRT*1d0/oDX_E(I)),&
                                   (V_s(IPJK, L)-VSH(IPJK)+VSH(IJK)&
                                   +SRT*1d0/oDX_E(I))), I)

                              Vsl_W = AVG_X(&               !i-1/2, j, k
                                   AVG_Y_N((V_s(IMJMK, L)-VSH(IMJMK)+&
                                   VSH(IJMK)-SRT*1d0/oDX_E(IM1(I))),&
                                   (V_s(IMJK, L)-VSH(IMJK)+VSH(IJK)&
                                   -SRT*1d0/oDX_E(IM1(I)))),&
                                   AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)), IM)

                         ELSE
                              Vsl_E = AVG_X(&               !i+1/2, j, k
                                   AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)),&
                                   AVG_Y_N(V_s(IPJMK, L), V_s(IPJK, L)), I )

                              Vsl_W = AVG_X(&               !i-1/2, j, k
                                   AVG_Y_N(V_s(IMJMK, L), V_s(IMJK, L)),&
                                   AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)), IM )
                         ENDIF

                         Vsl_T = AVG_Z(&                    !i, j, k+1/2
                              AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)),&
                              AVG_Y_N(V_s(IJMKP, L), V_s(IJKP, L)), K )
 
                         Vsl_B = AVG_Z(&                    !i, j, k-1/2
                              AVG_Y_N(V_s(IJMKM, L), V_s(IJKM, L)),&
                              AVG_Y_N(V_s(IJMK, L), V_s(IJK, L)), KM )

                         Wsl_N = AVG_Y(&                    !i, j+1/2, k
                              AVG_Z_T(W_s(IJKM, L), W_s(IJK, L)),&
                              AVG_Z_T(W_s(IJPKM, L), W_s(IJPK, L)), J )
 
                         Wsl_S = AVG_Y(&                    !i, j-1/2, k
                              AVG_Z_T(W_s(IJMKM, L), W_s(IJMK, L)),&
                              AVG_Z_T(W_s(IJKM, L), W_s(IJK, L)), JM )

                         Wsl_E = AVG_X(&                    !i+1/2, j, k
                              AVG_Z_T(W_s(IJKM, L), W_s(IJK, L)),&
                              AVG_Z_T(W_s(IPJKM, L), W_s(IPJK, L)), I)

                         Wsl_W = AVG_X(&                    !i-1/2, j, k
                              AVG_Z_T(W_s(IMJKM, L), W_s(IMJK, L)),&
                              AVG_Z_T(W_s(IJKM, L), W_s(IJK, L)), IM )
     
                         IF(CYLINDRICAL) THEN
                              Usl_C = AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I) !i, j, k
                              Wsl_C = AVG_Z_T(W_s(IJKM, L), W_s(IJK, L))    !i, j, k
                         ELSE
                              Usl_C = ZERO
                              Wsl_C = ZERO
                         ENDIF

!     Check for IS surfaces and modify solids velocity-comp accordingly
            IF(ANY_IS_DEFINED) THEN
              IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) Usl_N = AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I)
              IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) Usl_S = AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I)
              IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) Usl_T = AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I)
              IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) Usl_B = AVG_X_E(U_s(IMJK, L), U_s(IJK, L), I)
              IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) Vsl_E = AVG_Y_N(V_s(IJMK, L), V_s(IJK, L))
              IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) Vsl_W = AVG_Y_N(V_s(IJMK, L), V_s(IJK, L))
              IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) Vsl_T = AVG_Y_N(V_s(IJMK, L), V_s(IJK, L))
              IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) Vsl_B = AVG_Y_N(V_s(IJMK, L), V_s(IJK, L))
              IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) Wsl_N = AVG_Z_T(W_s(IJKM, L), W_s(IJK, L))
              IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) Wsl_S = AVG_Z_T(W_s(IJKM, L), W_s(IJK, L))
              IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) Wsl_E = AVG_Z_T(W_s(IJKM, L), W_s(IJK, L))
              IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) Wsl_W = AVG_Z_T(W_s(IJKM, L), W_s(IJK, L))
            ENDIF

!     Find components of Lth solids phase continuum strain rate
!     tensor, D_sl, at center of THETA cell-(i, j, k) 
                         D_sl(1,1) = ( U_s(IJK,L) - U_s(IMJK,L) ) * oDX(I)
                         D_sl(1,2) = HALF * ( (Usl_N - Usl_S) * oDY(J) +&
                              (Vsl_E - Vsl_W) * oDX(I) )
                         D_sl(1,3) = HALF * ( (Wsl_E - Wsl_W) * oDX(I) +&
                              (Usl_T - Usl_B) * (oX(I)*oDZ(K)) - Wsl_C * oX(I) )
                         D_sl(2,1) = D_sl(1,2)
                         D_sl(2,2) = ( V_s(IJK,L) - V_s(IJMK,L) ) * oDY(J)
                         D_sl(2,3) = HALF * ( (Vsl_T - Vsl_B) * (oX(I)*oDZ(K)) +&
                              (Wsl_N - Wsl_S) * oDY(J) )
                         D_sl(3,1) = D_sl(1,3)
                         D_sl(3,2) = D_sl(2,3)
                         D_sl(3,3) = ( W_s(IJK,L) - W_s(IJKM,L) ) * (oX(I)*oDZ(K)) +&
                               Usl_C * oX(I)
     
!     Calculate trace of the D_sl dot D_sm 
!     (normal matrix multiplication)
                         trD_s2_ip(IJK,M,L) = 0.D0 !Initialize the totalizer
                         DO 50 I1 = 1,3
                         DO 60 I2 = 1,3
                              trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,M,L)+&
                                   D_sl(I1,I2)*D_s(I1,I2)  
 60                      CONTINUE
 50                      CONTINUE
                     ELSE
                        trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,L,M)
                     ENDIF ! for L > M
                  ELSE
                       trD_s2_ip(IJK,M,M) = trD_s2(IJK,M)
                  ENDIF !for m NE L
               ENDDO
             ENDIF ! if kt_type = IA theory


!     Start definition of Relative Velocity
!     Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            USCM = AVG_X_E(U_S(IMJK,1),U_S(IJK,1),I) 
            VSCM = AVG_Y_N(V_S(IJMK,1),V_S(IJK,1)) 
            WSCM = AVG_Z_T(W_S(IJKM,1),W_S(IJK,1)) 
     
!     magnitude of gas-solids relative velocity
            VREL_array(IJK) = SQRT((UGC - USCM)**2 + &
            (VGC - VSCM)**2 + &
            (WGC - WSCM)**2)

!     Frictional-flow stress tensor
!     Calculate the second invariant of the deviator of D_s
            I2_devD_s(IJK) = ( (D_s(1,1)-D_s(2,2))**2&
            +(D_s(2,2)-D_s(3,3))**2&
            +(D_s(3,3)-D_s(1,1))**2 )/6.&
            + D_s(1,2)**2 + D_s(2,3)**2 + D_s(3,1)**2
     
            IF(SIMONIN) THEN
     
!     parameters for defining Tau_12: time-scale of the fluid turbulent motion
!     viewed by the particles (crossing trajectory effect)
               IF(SQRT(USCM**2+VSCM**2+WSCM**2) .GT. zero) THEN
                  Cos_Theta(IJK) = ((UGC-USCM)*USCM+(VGC-VSCM)*VSCM+(WGC-WSCM)*WSCM)/ &
                  (SQRT((UGC-USCM)**2+(VGC-VSCM)**2+(WGC-WSCM)**2)*  &
                  SQRT(USCM**2+VSCM**2+WSCM**2))
               ELSE
                  Cos_Theta(IJK) = ONE
               ENDIF
            ENDIF

            IF(.NOT.GRANULAR_ENERGY) THEN !algebraic granular energy equation
               IF(EP_g(IJK) .GE. EP_star_array(IJK)) THEN
    
!     Boyle-Massoudi Stress term
                  IF(V_ex .NE. ZERO) THEN
                     DEP_soDX  = ( EP_s(IJKE, M) - EP_s(IJK, M) ) * oDX_E(I)&
                     * ( ONE / ( (oDX_E(IM)/oDX_E(I)) + ONE ) ) +&
                     ( EP_s(IJK, M) - EP_s(IJKW, M) ) * oDX_E(IM)&
                     * ( ONE / ( (oDX_E(I)/oDX_E(IM)) + ONE ) )
                     DEP_soDY  = ( EP_s(IJKN, M) - EP_s(IJK, M) ) * oDY_N(J)&
                     * ( ONE / ( (oDY_N(JM)/oDY_N(J)) + ONE ) ) +&
                     ( EP_s(IJK, M) - EP_s(IJKS, M) ) * oDY_N(JM)&
                     * ( ONE / ( (oDY_N(J)/oDY_N(JM)) + ONE ) )
                     DEP_soXDZ  = (( EP_s(IJKT, M) - EP_s(IJK, M) )&
                     * oX(I)*oDZ_T(K)&
                     * ( ONE / ( (oDZ_T(KM)/oDZ_T(K)) + ONE ) ) +&
                     ( EP_s(IJK, M) - EP_s(IJKB, M) ) * oX(I)*oDZ_T(KM)&
                     * ( ONE / ( (oDZ_T(K)/oDZ_T(KM)) + ONE ) ) ) /&
                     X(I)
                     M_s(1,1) = DEP_soDX * DEP_soDX
                     M_s(1,2) = DEP_soDX * DEP_soDY
                     M_s(1,3) = DEP_soDX * DEP_soXDZ
                     M_s(2,1) = DEP_soDX * DEP_soDY
                     M_s(2,2) = DEP_soDY * DEP_soDY
                     M_s(2,3) = DEP_soDY * DEP_soXDZ
                     M_s(3,1) = DEP_soDX * DEP_soXDZ
                     M_s(3,2) = DEP_soDY * DEP_soXDZ
                     M_s(3,3) = DEP_soXDZ * DEP_soXDZ
                     trM_s(IJK)    = M_s(1,1) + M_s(2,2) + M_s(3,3)
                     trDM_s(IJK) = ZERO
                     DO 40 I1 = 1,3
                        DO 30 I2 = 1,3
                           trDM_s(IJK) = trDM_s(IJK) + D_s(I1,I2)*M_s(I1,I2)
 30                     CONTINUE
 40                  CONTINUE

                  ENDIF
               ENDIF
            ENDIF

         Endif                  ! Fluid_at
 200  Continue                  ! outer IJK loop
      
      Return
      End 

