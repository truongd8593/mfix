!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_s                                               C
!  Purpose: Calculate granular stress terms (granular viscosity        C
!     bulk viscosity, solids pressure) & granular conductivity         C
!                                                                      C
!  Author: W. Rogers                                Date: 04-mar-92    C
!  Reviewer: M. Syamlal                             Date: 16-MAR-92    C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for cylindrical geometry                     C
!  Author: M. Syamlal                               Date: 15-MAY-92    C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add volume-weighted averaging statement functions for      C
!     variable grid capability                                         C
!  Author:  W. Rogers                               Date: 21-JUL-92    C
!  Reviewer: P. Nicoletti                           Date: 11-DEC-92    C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: Add Boyle-Massoudi stress terms                            C
!  Author: M. Syamlal                               Date: 2-NOV-95     C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: MFIX 2.0 mods  (old name CALC_THETA)                       C
!  Author: M. Syamlal                               Date: 24-APR-96    C
!                                                                      C
!  Author: Sreekanth Pannala, ORNL                  Date: 10-08-05     C
!  Revision Number:X                                                   C
!  Purpose: Rewrite different modules to increase modularity           C
!                                                                      C
!  Author: QX, Iowa State                                              C
!  Revision Number:X                                                   C
!  Purpose: Variable density- changed RO_S(M) to RO_S(IJK,M)          C
!                                                                      C
!  Author: Handan Liu                                                  C
!  Revision Number:X                                                   C
!  Purpose: Added OpenMP to some routines - gt_algebraic,              C
!     init_mu_s                                                        C
!                                                                      C
!  Comments:                                                           C
!     GRANULAR_ENERGY = .FALSE.                                        C
!        EP_g < EP_star   -->    friction_schaeffer                    C
!        EP_g >= EP_star  -->    viscous (algebraic)                   C
!                                                                      C
!     GRANULAR_ENERGY = .TRUE.                                         C
!        FRICTION = .TRUE.                                             C
!           EP_s(IJK,M) > EPS_f_min  -->  friction + viscous(pde)      C
!           EP_s(IJK,M) < EP_f_min   -->  viscous (pde)                C
!        FRICTION = .FALSE.                                            C
!           EP_g < EP_star  -->  friction_schaeffer + viscous(pde)     C
!           EP_g >= EP_star -->  viscous (pde)                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
     
      SUBROUTINE CALC_MU_s(M, IER)

!-----------------------------------------------
! Modules 
!----------------------------------------------- 
      USE run
      USE vshear
      USE visc_s
      USE physprop
      USE constant
      USE fldvar
      USE compar
      USE indices
      USE geometry
      USE qmom_kinetic_equation
      Implicit NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell index
      INTEGER :: IJK      
! blend factor
      DOUBLE PRECISION :: BLEND
!----------------------------------------------- 
! Functions
!-----------------------------------------------       
      DOUBLE PRECISION, EXTERNAL :: BLEND_FUNCTION
!----------------------------------------------- 
! Include statement functions
!-----------------------------------------------      
      Include 'function.inc'
!-----------------------------------------------

! GHD Theory is called only for the mixture granular energy, i.e. for m == mmax
      IF (TRIM(KT_TYPE) == 'GHD' .AND. M /= MMAX) RETURN

      IF (SHEAR) CALL add_shear(M)

! Initialize/calculate all the quantities needed for various options
      CALL INIT_MU_S(M, IER)    
     
      IF (SHEAR) call remove_shear(M)

      IF(MU_s0 /= UNDEFINED) RETURN ! constant solids viscosity case 

! Viscous-flow stress tensor
      IF (.NOT. QMOMK) THEN
! if QMOMK then do not solve algebraic or PDE form of granular 
! temperature governing equation
         IF(.NOT.GRANULAR_ENERGY) then
            IF(SUBGRID_TYPE /= UNDEFINED_C) THEN
               IF (TRIM(SUBGRID_TYPE) .EQ. 'IGCI') THEN
                  CALL SUBGRID_STRESS_IGCI(M, IER)
               ELSEIF (TRIM(SUBGRID_TYPE) .EQ. 'MILIOLI') THEN
                  CALL SUBGRID_STRESS_MILIOLI(M, IER)
               ENDIF                    
            ELSE
              call gt_algebraic(M,IER)   ! algebraic granular energy equation
            ENDIF
         ELSE   ! granular energy transport equation
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               CALL gt_pde_ia_nonep(M,IER) ! complete polydisperse IA theory
            ELSEIF (TRIM(KT_TYPE) .EQ. 'GD_99' ) THEN
               CALL gt_pde_gd_99(M,IER) ! monodisperse GD theory
            ELSEIF(TRIM(KT_TYPE) .EQ. 'GTSH') THEN
               CALL gt_pde_gtsh(M,IER) ! GTSH theory
            ELSEIF (TRIM(KT_TYPE) == 'GHD') THEN
               CALL TRANSPORT_COEFF_GHD(M,IER) ! GHD theory for mixture temperature
            ELSE
               CALL gt_pde(M,IER) ! This is also used whith Simonin or Ahmadi models
            ENDIF
         ENDIF
      ENDIF
    
! Frictional stress tensors
! Schaeffer's frictional formulation      
      IF (SCHAEFFER .AND. CLOSE_PACKED(M)) call friction_schaeffer(M,IER)
! Princeton's frictional implementation
      IF (FRICTION .AND. CLOSE_PACKED(M)) call friction_princeton(M,IER) 
      
      IF(BLENDING_STRESS) THEN
         DO 200 IJK = ijkstart3, ijkend3
            blend =  blend_function(IJK)
            Mu_s_c(IJK,M) = Mu_s_v(IJK)
            Mu_s(IJK,M) = (1.0d0-blend)*Mu_s_p(IJK) &
                + blend*Mu_s_v(IJK) + Mu_s_f(IJK)

! Bulk viscosity in Mth solids phase   
            LAMBDA_s_c(IJK,M)= Lambda_s_v(IJK)
            LAMBDA_s(IJK,M) = (1.0d0-blend)*LAMBDA_s_p(IJK) &
                + blend*Lambda_s_v(IJK) + Lambda_s_f(IJK)

! Solids pressure in the Mth solids phase 
! Note that the plastic pressure component (represented here by P_s_p)
! is calculated in a separate routine (see calc_p_star) which is then
! directly incorporated into each of thhe solids momentum equations 
! (see source_u_s, source_v_s and source_w_s).
            P_s_c(IJK,M) = P_s_v(IJK)
            P_s(IJK,M) = (1.0d0-blend)*P_s_p(IJK) + blend*P_s_v(IJK) &
                + P_s_f(IJK) 

! Boyle-Massoudi stress coefficient
! Alpha_s is only calculated in the algebraic granular energy
! subroutine.  No other values of alpha_s are set (i.e. no
! plastic, viscous or frictional components)
            ALPHA_s(IJK,M) = (1.0d0-blend)*ALPHA_s_p(IJK) &
                + blend*ALPHA_s_v(IJK) + ALPHA_s_f(IJK)
 200     ENDDO
      ELSE   ! else branch of if(blending_stress), no blending stress then we have
         Mu_s_c(:,M) = Mu_s_v(:)
         Mu_s(:,M) = Mu_s_p(:) + Mu_s_v(:) + Mu_s_f(:)
     
         LAMBDA_s_c(:,M)= Lambda_s_v(:)
         LAMBDA_s(:,M) = LAMBDA_s_p(:) + Lambda_s_v(:) + Lambda_s_f(:)
         
         P_s_c(:,M) = P_s_v(:)
         P_s(:,M) = P_s_p(:) + P_s_v(:) + P_s_f(:)
      ENDIF  ! end if/else (blending_stress)
      
      RETURN
      END SUBROUTINE CALC_MU_s




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FRICTION_SCHAEFFER                                      C
!  Purpose: Add frictional-flow stress terms                           C
!                                                                      C
!  Author: M. Syamlal                               Date: 10-FEB-93    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine friction_schaeffer(M, IER)

!-----------------------------------------------
! Modules 
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
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! maximum value of solids viscosity in poise
      DOUBLE PRECISION :: MAX_MU_s
      PARAMETER (MAX_MU_s = 1000.D0)
!-----------------------------------------------
! Local variables 
!-----------------------------------------------
! cell index
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! sum of all solids volume fractions
      DOUBLE PRECISION :: SUM_EPS_CP
! factor in frictional-flow stress terms
      DOUBLE PRECISION :: qxP_s
!----------------------------------------------- 
! Include statement functions
!-----------------------------------------------      
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 's_pr2.inc'
!----------------------------------------------- 

      DO 200 IJK = ijkstart3, ijkend3       
         IF ( FLUID_AT(IJK) ) THEN
            
! added closed pack, this has to be consistent with the normal frictional force
! see for example source_v_s.f. Tardos Powder Tech. 92 (1997) 61-74 explains in
! his equation (3) that solids normal and shear frictional stresses have to be 
! treated consistently. --> sof May 24 2005. 
     
            IF(EP_g(IJK) .LT. EP_g_blend_end(IJK)) THEN
! part copied from source_v_s.f (sof)
               SUM_EPS_CP=0.0
               DO MM=1,MMAX
                  IF (CLOSE_PACKED(MM)) SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
               END DO
! end of part copied
     
!            P_star(IJK) = Neg_H(EP_g(IJK),EP_star_array(IJK))
     
! Frictional-flow stress tensor
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
! Schaeffer (1987)
     
               qxP_s           = SQRT( (4.D0 * Sin2_Phi) * I2_devD_s(IJK))
               MU_s_p(IJK)     = P_star(IJK) * Sin2_Phi&
                   / (qxP_s + SMALL_NUMBER) &
                   *(EP_S(IJK,M)/SUM_EPS_CP) ! added by sof for consistency
                                             ! with solids pressure treatment
               MU_s_p(IJK)     = MIN(MU_s_p(IJK), to_SI*MAX_MU_s)
               
               LAMBDA_s_p(IJK) = ZERO
     
! when solving for the granular energy equation (PDE) setting theta = 0 is done 
! in solve_granular_energy.f to avoid convergence problems. (sof)
               IF(.NOT.GRANULAR_ENERGY) THETA_m(IJK, M) = ZERO
            ENDIF
         ENDIF   ! end if (fluid_at(ijk)


 200  CONTINUE   ! outer IJK loop
      
      RETURN
      END SUBROUTINE FRICTION_SCHAEFFER



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine Gt_algebraic (M, IER)

!-----------------------------------------------
! Modules 
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
      USE trace
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables 
!----------------------------------------------- 
! cell index
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! Sum of all solids volume fractions
      DOUBLE PRECISION :: SUM_EPS_CP
! Coefficients of quadratic equation
      DOUBLE PRECISION :: aq, bq, cq
! Constant in equation for mth solids phase pressure
      DOUBLE PRECISION :: K_1m
! Constant in equation for mth solids phase bulk viscosity
      DOUBLE PRECISION :: K_2m
! Constant in equation for mth solids phase viscosity
      DOUBLE PRECISION :: K_3m
! Constant in equation for mth solids phase dissipation
      DOUBLE PRECISION :: K_4m
! Constant in Boyle-Massoudi stress term
      DOUBLE PRECISION :: K_5m
! Factor in frictional-flow stress terms
      DOUBLE PRECISION :: qxP_s
! Value of EP_s * SQRT( THETA )for Mth solids phase continuum
      DOUBLE PRECISION :: EP_sxSQRTHETA
! Value of EP_s * EP_s * THETA for Mth solids phase continuum
      DOUBLE PRECISION :: EP_s2xTHETA, temp_local
!----------------------------------------------- 
! Functions
!----------------------------------------------- 
! radial distribution function      
      DOUBLE PRECISION, EXTERNAL :: G_0
!----------------------------------------------- 
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'  
!----------------------------------------------- 

!$omp parallel do default(shared)                                    &
!$omp private(IJK, K_1m, K_2m, K_3m, K_4m, K_5m, temp_local,         &
!$omp         aq, bq, cq, EP_sxSQRTHETA, EP_s2xTHETA)  
     
       DO IJK = ijkstart3, ijkend3
   
         IF ( FLUID_AT(IJK) ) THEN

            IF(EP_g(IJK) .GE. EP_g_blend_start(IJK)) THEN

! Calculate K_1m, K_2m, K_3m, K_4m
               K_1m = 2.D0 * (ONE + C_e) * RO_S(IJK,M) * G_0(IJK, M, M)
               K_3m = HALF * D_p(IJK,M) * RO_S(IJK,M) * (&
                   ( (SQRT_PI / (3.D0*(3.D0 - C_e))) *&
                   (HALF*(3d0*C_e+ONE) + 0.4D0*(ONE + C_e)*(3.D0*C_e - ONE)*&
                   EP_s(IJK,M)*G_0(IJK, M,M)) ) + 8.D0*EP_s(IJK,M)&
                   *G_0(IJK, M,M)*(ONE + C_e)/ (5.D0*SQRT_PI) )
               K_2m = 4.D0 * D_p(IJK,M) * RO_S(IJK,M) * (ONE + C_e) *&
                   EP_s(IJK,M) * G_0(IJK, M,M) / (3.D0 * SQRT_PI) - 2.D0/3.D0 * K_3m
               K_4m = 12.D0 * (ONE - C_e*C_e) *&
                   RO_S(IJK,M) * G_0(IJK, M,M) / (D_p(IJK,M) * SQRT_PI)
               aq   = K_4m*EP_s(IJK,M)
               bq   = K_1m*EP_s(IJK,M)*trD_s_C(IJK,M)
               cq   = -(K_2m*trD_s_C(IJK,M)*trD_s_C(IJK,M)&
                   + 2.D0*K_3m*trD_s2(IJK,M))
     
! Boyle-Massoudi Stress term
               IF(V_ex .NE. ZERO) THEN
                  K_5m = 0.4 * (ONE + C_e) * G_0(IJK, M,M) * RO_S(IJK,M) *&
                  ( (V_ex * D_p(IJK,M)) / (ONE - EP_s(IJK,M) * V_ex) )**2
                  bq   = bq + EP_s(IJK,M) * K_5m * (trM_s(IJK) + 2.D0 * trDM_s(IJK))
               ELSE
                  K_5m = ZERO
               ENDIF
     
! Calculate EP_sxSQRTHETA and EP_s2xTHETA
               temp_local = bq**2 - 4.D0 * aq * cq
               EP_sxSQRTHETA = (-bq + SQRT(temp_local))&
                  / ( 2.D0 * K_4m )
               EP_s2xTHETA = EP_sxSQRTHETA * EP_sxSQRTHETA
               
               IF(EP_s(IJK,M) > SMALL_NUMBER)THEN
! Find pseudo-thermal temperature in the Mth solids phase
                  THETA_m(IJK,M) = EP_s2xTHETA/(EP_s(IJK,M)*EP_s(IJK,M))
               ELSE
                  THETA_m(IJK,M) = ZERO
               ENDIF
     
! Find pressure in the Mth solids phase
               P_s_v(IJK) = K_1m * EP_s2xTHETA
     
! bulk viscosity in Mth solids phase
               LAMBDA_s_v(IJK) = K_2m * EP_sxSQRTHETA
     
! shear viscosity in Mth solids phase
               MU_s_v(IJK) = K_3m * EP_sxSQRTHETA
     
! Boyle-Massoudi stress coefficient
               ALPHA_s(IJK, M) = -K_5m * EP_s2xTHETA
            ENDIF

         ENDIF   ! Fluid_at
      ENDDO
!$omp end parallel do

      RETURN
      END SUBROUTINE GT_ALGEBRAIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE                                                  C
!  Purpose: Calculate granular stress terms (viscosity, bulk viscosity C
!     solids pressure) & granular conductivity                         C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University      Date: 6-FEB-98     C
!                                                                      C
!  Revision Number:1                                                   C
!  Purpose: Add Simonin and Ahmadi models                              C
!  Author: Sofiane Benyahia, Fluent Inc.            Date: 02-01-05     C
!                                                                      C
!  Literature/Document References:                                     C
!     Lun, C.K.K., S.B. Savage, D.J. Jeffrey, and N. Chepurniy,        C
!        Kinetic theories for granular flow - inelastic particles in   C
!        Couette-flow and slightly inelastic particles in a general    C
!        flow field. Journal of Fluid Mechanics, 1984. 140(MAR):       C
!        p. 223-256                                                    C
!                                                                      C
!     Simonin, O., 1996. Combustion and turbulence in two-phase flows, C
!        Von Karman institute for fluid dynamics, lecture series,      C
!        1996-02                                                       C
!     Balzer, G., Simonin, O., Boelle, A., and Lavieville, J., 1996,   C
!        A unifying modelling approach for the numerical prediction    C
!        of dilute and dense gas-solid two phase flow. CFB5, 5th int.  C
!        conf. on circulating fluidized beds, Beijing, China.          C
!     Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent   C
!        flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,    C
!        No. 6, pp. 1203-1228.                                         C
!                                                                      C      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine gt_pde (M, IER)

!-----------------------------------------------
! Modules 
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
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell indices
      INTEGER :: IJK, I, J, K
! solids phase index
      INTEGER :: MM
! use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Kth, Kth_star
! defining parametrs for Simonin and Ahmadi models
      DOUBLE PRECISION :: Tau_12_st, Tau_2_c, Tau_2, Zeta_r, C_Beta
      DOUBLE PRECISION :: Sigma_c, Zeta_c, Omega_c, Zeta_c_2, C_mu, X_21, Nu_t
      DOUBLE PRECISION :: MU_2_T_Kin, Mu_2_Col, Kappa_kin, Kappa_Col
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
!
      DOUBLE PRECISION :: DGA, C_d, Re
! sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
! (Same for conductivity)

!-----------------------------------------------      
! Functions
!-----------------------------------------------     
! radial distribution function
      DOUBLE PRECISION, EXTERNAL :: G_0
! dg0/dep
      DOUBLE PRECISION, EXTERNAL :: DG_0DNU

!-----------------------------------------------
! Include statement functions
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
     
! Defining a single particle drag coefficient (similar to one defined in drag_gs)
            RE = D_p(IJK,M)*VREL_array(IJK)*ROP_G(IJK)/&
               (MU_G(IJK) + SMALL_NUMBER)
            IF(RE .LE. 1000D0)THEN
               C_d = (24.D0/(Re+SMALL_NUMBER)) * (ONE + 0.15D0 * Re**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF

! This is from Wen-Yu correlation, you can put here your own single particle drag
            DgA = 0.75D0 * C_d * VREL_array(IJK) * ROP_g(IJK) / D_p(IJK,M)
! set value for 1st iteration and 1st time step
            IF(VREL_array(IJK) == ZERO) DgA = LARGE_NUMBER 

! Define some time scales and constants related to Simonin and Ahmadi models
            IF(SIMONIN .OR. AHMADI) THEN
               C_mu = 9.0D-02
! particle relaxation time. For very dilute flows avoid singularity by
! redefining the drag as single partilce drag
               IF(Ep_s(IJK,M) > DIL_EP_S .AND. F_GS(IJK,1) > small_number) THEN
                  Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
               ELSE             !for dilute flows, drag equals single particle drag law
                  Tau_12_st = RO_S(IJK,M)/DgA
               ENDIF            !for dilute flows
! time scale of turbulent eddies
               Tau_1(ijk) = 3.d0/2.d0*C_MU*K_Turb_G(IJK)/(E_Turb_G(IJK)+small_number)
            ENDIF


! Define some time scales and constants and K_12 related to Simonin model only	  
            IF(SIMONIN) THEN
! This is Zeta_r**2 as defined by Simonin
               Zeta_r = 3.0d0 * VREL_array(IJK)**2 / (2.0d0*K_Turb_G(IJK)+small_number)
! parameters for defining Tau_12: time-scale of the fluid turbulent motion
! viewed by the particles (crossing trajectory effect)
               C_Beta = 1.8d0 - 1.35d0*Cos_Theta(IJK)**2
! Lagrangian Integral time scale: Tau_12	    
               Tau_12(ijk) = Tau_1(ijk)/sqrt(ONE+C_Beta*Zeta_r)
! Defining the inter-particle collision time
               IF(Ep_s(IJK,M) > DIL_EP_S) THEN
                  Tau_2_c = D_p(IJK,M)/(6.d0*Ep_s(IJK,M)*G_0(IJK,M,M) &
                  *DSQRT(16.d0*(Theta_m(ijk,m)+Small_number)/PI))
               ELSE             ! assign it a large number
                  Tau_2_c = LARGE_NUMBER
               ENDIF
     
               Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
! Zeta_c: const. to be used in the K_2 Diffusion coefficient.
               Zeta_c  = (ONE+ C_e)*(49.d0-33.d0*C_e)/100.d0
               Omega_c = 3.d0*(ONE+ C_e)**2 *(2.d0*C_e-ONE)/5.d0
               Zeta_c_2= 2./5.*(ONE+ C_e)*(3.d0*C_e-ONE)

! mixed time scale in the generalized Simonin theory (switch between dilute
! and kinetic theory formulation of the stresses)
               Tau_2 = ONE/(2./Tau_12_st+Sigma_c/Tau_2_c)
! The ratio of densities
               X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))
! The ratio of these two time scales.
               Nu_t =  Tau_12(ijk)/Tau_12_st
     
! Definition of an "algebraic" form of of Simonin K_12 PDE. This is obtained
! by equating the dissipation term to the exchange terms in the PDE and 
! neglecting all other terms, i.e. production, convection and diffusion.
! This works because Tau_12 is very small for heavy particles
               K_12(ijk) = Nu_t / (ONE+Nu_t*(ONE+X_21)) * &
                   (2.d+0 *K_Turb_G(IJK) + 3.d+0 *X_21*theta_m(ijk,m))
! Realizability Criteria         
               IF(K_12(ijk) > DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))) THEN
                  K_12(ijk) = DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))
               ENDIF
            ENDIF               ! for Simonin


! This is added for consistency of multi-particles kinetic theory. Solids pressure,
! viscosity and conductivity must be additive. Thus non-linear terms (eps^2) are 
! corrected so the stresses of two identical solids phases are equal to those
! of a single solids phase. sof June 15 2005.
            SUM_EpsGo = ZERO
            DO MM = 1, MMAX
               SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,M,MM)
            ENDDO 
     
! Find pressure in the Mth solids phase
            P_s_v(IJK) = ROP_s(IJK,M)*(1d0+ 4.D0 * Eta *&
                SUM_EpsGo)*Theta_m(IJK,M)
     
! implement Simonin (same as granular) and Ahmadi solids pressures
            IF(SIMONIN) THEN
               P_s_v(IJK) = P_s_v(IJK) ! no changes to solids pressure
            ELSE IF(AHMADI) THEN
               P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M) * ( (ONE + 4.0D0* &
                   SUM_EpsGo ) + HALF*(ONE - C_e*C_e) )
            ENDIF

! find bulk and shear viscosity
            Mu_s_v(IJK) = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(IJK,M)*RO_S(IJK,M))/96d0
            Mu_b_v(IJK) = (256d0*Mu_s_v(IJK)*EP_s(IJK,M)*SUM_EpsGo)&
                /(5d0*Pi)

! added Ro_g = 0 for granular flows (no gas). sof Aug-02-2005 
            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN !sof modifications (May 20 2005)
               Mu_star = Mu_s_v(IJK)
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Mu_star = ZERO
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               Mu_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)* Mu_s_v(IJK)/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M) &
                   + 2.0d0*SWITCH*DgA/RO_S(IJK,M)* Mu_s_v(IJK))
            ELSE
               Mu_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)*Mu_s_v(IJK)/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M)+ &
                   (2d0*SWITCH*F_gs(IJK,M)*Mu_s_v(IJK)/(RO_S(IJK,M)*EP_s(IJK,M))) )
            ENDIF
     
! shear viscosity in Mth solids phase  (add to frictional part)
            Mu_s_v(IJK) =&
                ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
                G_0(IJK,M,M)))*(ONE+1.6d0*Eta*SUM_EpsGo)&
                *(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
                SUM_EpsGo)+(0.6d0*Mu_b_v(IJK)*Eta))


! implement Simonin and Ahmadi solids viscosities
            IF(SIMONIN) THEN
! Defining Simonin solids turbulent Kinetic (MU_2_T_Kin) and collisional (Mu_2_Col)
! viscosities
               MU_2_T_Kin = (2.0d0/3.0d0*K_12(ijk)*Nu_t + Theta_m(IJK,M) * &
                   (ONE+ zeta_c_2*EP_s(IJK,M)*G_0(IJK,M,M)))*Tau_2
               Mu_2_Col = 8.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* (MU_2_T_Kin+ &
                   D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
               Mu_b_v(IJK) = 5.d0/3.d0*EP_s(IJK,M)*RO_S(IJK,M)*Mu_2_Col
               Mu_s_v(IJK) = EP_s(IJK,M)*RO_S(IJK,M)*(MU_2_T_Kin + Mu_2_Col)
               
            ELSE IF(AHMADI) THEN
              IF(EP_s(IJK,M) < (ONE-EP_star_array(ijk))) THEN
                 Tmp_Ahmadi_Const = &
                    ONE/(ONE+ Tau_1(ijk)/Tau_12_st * &
                    (ONE-EP_s(IJK,M)/(ONE-EP_star_array(ijk)))**3)
              ELSE
                 Tmp_Ahmadi_Const = ONE
              ENDIF
! Defining Ahmadi shear and bulk viscosities. Ahmadi coefficient 0.0853 in C_mu
! was replaced by 0.1567 to include 3/2*sqrt(3/2) because K = 3/2 Theta_m
               Mu_s_v(IJK) = Tmp_Ahmadi_Const &
                   *0.1045d0*(ONE/G_0(IJK,M,M)+3.2d0*EP_s(IJK,M)+12.1824d0*   &
                   G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M))*D_p(IJK,M)*RO_S(IJK,M)*  &
                   DSQRT(Theta_m(IJK,M))
! This is a guess of what Mu_b might be by taking 5/3 of the collisional viscosity
! contribution. In this case col. visc. is the eps^2 contribution to Mu_s_v(IJK). This
! might be changed later if communications with Ahmadi reveals a diffrent appoach
               Mu_b_v(IJK) = 5.d0/3.d0* Tmp_Ahmadi_Const                  &
                   *0.1045d0*(12.1824d0*G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M)) &
                   *D_p(IJK,M)*RO_S(IJK,M)* DSQRT(Theta_m(IJK,M))
            ENDIF               !for simonin or ahmadi viscosity
            
            
            Kth=75d0*RO_S(IJK,M)*D_p(IJK,M)*DSQRT(Pi*Theta_m(IJK,M))/&
                (48d0*Eta*(41d0-33d0*Eta))
            
            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN ! sof modifications (May 20 2005)
               Kth_star=Kth
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Kth_star = ZERO
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               Kth_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)* Kth/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M) &
                   + 1.2d0*SWITCH*DgA/RO_S(IJK,M)* Kth)
            ELSE
               Kth_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*Theta_m(IJK,M)*Kth/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M)+ &
                   (1.2d0*SWITCH*F_gs(IJK,M)*Kth/(RO_S(IJK,M)*EP_s(IJK,M))) )
            ENDIF
     
! granular conductivity in Mth solids phase
            Kth_s(IJK,M) = Kth_star/G_0(IJK,M,M)*(&
                ( ONE + (12d0/5.d0)*Eta*SUM_EpsGo )&
                * ( ONE + (12d0/5.d0)*Eta*Eta*(4d0*Eta-3d0)* SUM_EpsGo )&
                + (64d0/(25d0*Pi)) * (41d0-33d0*Eta) * (Eta*SUM_EpsGo)**2 )


! implement Simonin and Ahmadi solids conductivities
            IF(SIMONIN) THEN
! Defining Simonin's Solids Turbulent Kinetic diffusivity: Kappa
               Kappa_kin = (9.d0/10.d0*K_12(ijk)*Nu_t + 3.0D0/2.0D0 * &
                   Theta_m(IJK,M)*(ONE+ Omega_c*EP_s(IJK,M)*G_0(IJK,M,M)))/&
                   (9.d0/(5.d0*Tau_12_st) + zeta_c/Tau_2_c)
               Kappa_Col = 18.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* & 
                   (Kappa_kin+ 5.d0/9.d0*D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
               Kth_s(IJK,M) =  EP_s(IJK,M)*RO_S(IJK,M)*(Kappa_kin + Kappa_Col)
     
            ELSEIF(AHMADI) THEN
! Defining Ahmadi conductivity from his equation 42 in Cao and Ahmadi 1995 paper
! note the constant 0.0711 is now 0.1306 because K = 3/2 theta_m
               Kth_s(IJK,M) = 0.1306D0*RO_S(IJK,M)*D_p(IJK,M)*(ONE+C_e**2)* &
                   (ONE/G_0(IJK,M,M)+4.8D0*EP_s(IJK,M)+12.1184D0 &
                   *EP_s(IJK,M)*EP_s(IJK,M)*G_0(IJK,M,M) )  &
                   *DSQRT(Theta_m(IJK,M))
            ENDIF

            LAMBDA_S_V(IJK) = Eta*Mu_b_v(IJK) - (2d0*Mu_s_v(IJK))/3d0
     
! granular 'conductivity' in the Mth solids phase associated
! with gradient in volume fraction
            
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

         ENDIF   ! Fluid_at
 200  CONTINUE   ! outer IJK loop
      

      RETURN
      END SUBROUTINE GT_PDE




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_GD_99                                            C
!  Purpose: Implement kinetic theory of Garzo and Dufty (1999) for     C
!     calculation of granular stress terms and granular conductivity   C
!                                                                      C
!  Author: Janine E. Galvin                                            C
!                                                                      C
!  Literature/Document References:                                     C
!    Garzo, V., and Dufty, J., Homogeneous cooling state for a         C
!    granular mixture, Physical Review E, 1999, Vol 60 (5), 5706-      C
!    5713                                                              C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      Subroutine gt_pde_gd_99 (M, IER)

!-----------------------------------------------
! Modules 
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
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell Indices
      INTEGER :: IJK, I, J, K
! solids phase index
      INTEGER :: L 
! Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Kth_star
!
      DOUBLE PRECISION :: DGA, C_d, Re
!
      DOUBLE PRECISION :: D_PM, M_PM, NU_PM, EP_SM, RO_SM, ROP_SM
!
      DOUBLE PRECISION :: chi, dChiOdphi
!
      DOUBLE PRECISION :: c_star, zeta0_star, nu_eta_star, &
                          gamma_star, eta_k_star, eta_star, eta0, &
                          kappa0, nu_kappa_star, kappa_k_star, &
                          qmu_k_star, qmu_star, kappa_star, press_star

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
! (Same for conductivity)

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
! radial distribution function
      DOUBLE PRECISION, EXTERNAL :: G_0
! dg0/dep
      DOUBLE PRECISION, EXTERNAL :: DG_0DNU
!-----------------------------------------------
! Include statement functions
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

! local aliases
             D_PM = D_P(IJK,M)
             RO_SM = RO_S(IJK,M)
             ROP_SM = ROP_S(IJK,M)
             EP_SM = EP_S(IJK,M)
             M_PM = (PI/6.d0)*D_PM**3 * RO_SM
             NU_PM = ROP_SM/M_PM
             Chi = G_0(IJK,M,M)
             dChiOdphi = DG_0DNU(EP_SM)

! Defining a single particle drag coefficient (similar to one defined in drag_gs)
             Re = D_pM*VREL_array(IJK)*ROP_G(IJK)/&
                (MU_G(IJK) + small_number)
             IF(Re .LE. 1000D0)THEN
                C_d = (24.D0/(Re+SMALL_NUMBER)) * &
                   (ONE + 0.15D0 * Re**0.687D0)
             ELSE
                C_d = 0.44D0
             ENDIF

! This is from Wen-Yu correlation, you can put here your own single particle drag
             DgA = 0.75D0 * C_d * VREL_array(IJK) * ROP_g(IJK) / D_PM
             IF (VREL_array(IJK) == ZERO) THEN
                DgA = LARGE_NUMBER     ! for 1st iteration and 1st time step
             ENDIF
    
! Pressure/Viscosity/Bulk Viscosity
! Note: k_boltz = M_PM
!-----------------------------------
! Find pressure in the Mth solids phase
             press_star = 1.d0 + 2.d0*(1.d0+C_E)*EP_SM*Chi

! n*k_boltz = n*m = ep_s*ro_s
             P_s_v(IJK) = ROP_sM*Theta_m(IJK,M)*press_star
    
! find bulk and shear viscosity
             c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
                / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

             zeta0_star = (5.d0/12.d0)*Chi*(1.d0 - C_E*C_E) &
                * (1.d0 + (3.d0/32.d0)*c_star)

             nu_eta_star = Chi*(1.d0 - 0.25d0*(1.d0-C_E)*(1.d0-C_E)) &
                * (1.d0-(c_star/64.d0))

             gamma_star = (4.d0/5.d0)*(32.d0/PI)*EP_SM*EP_SM &
                * Chi*(1.d0+C_E) * (1.d0 - (c_star/32.d0))

             eta_k_star = (1.d0 - (2.d0/5.d0)*(1.d0+C_E)*(1.d0-3.d0*C_E) &
                * EP_SM*Chi ) / (nu_eta_star - 0.5d0*zeta0_star)

             eta_star = eta_k_star*(1.d0 + (4.d0/5.d0)*EP_SM*Chi &
                * (1.d0+C_E) ) + (3.d0/5.d0)*gamma_star

             eta0 = 5.0d0*M_PM*DSQRT(Theta_m(IJK,M)/PI) / (16.d0*D_PM*D_PM)
    
! added Ro_g = 0 for granular flows (no gas). 
             IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN 
                Mu_star = eta0
             ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                Mu_star = ZERO
             ELSEIF(EP_SM < DIL_EP_S) THEN
                Mu_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*eta0 / &
                   ( RO_S(IJK,M)*EP_SM*Chi*Theta_m(IJK,M) + &
                   2.d0*DgA*eta0/RO_S(IJK,M) )
             ELSE
                Mu_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*eta0 / &
                   ( RO_SM*EP_SM*Chi*Theta_m(IJK,M) + &
                   (2.d0*F_gs(IJK,M)*eta0/(RO_SM*EP_SM)) )
             ENDIF

! shear viscosity in Mth solids phase  (add to frictional part)
             Mu_s_v(IJK) = Mu_star * eta_star
             Mu_b_v(IJK) = Mu_star * gamma_star 

! second viscosity
             LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)

 
! Granular Conductivity/Dufour Coefficient
!-----------------------------------
             kappa0 = (15.d0/4.d0)*eta0

             nu_kappa_star = (Chi/3.d0)*(1.d0+C_E) * ( 1.d0 + &
                (33.d0/16.d0)*(1.d0-C_E) + ((19.d0-3.d0*C_E)/1024.d0)*&
                c_star)
!             nu_mu_star = nu_kappa_star

             kappa_k_star = (2.d0/3.d0)*(1.d0 +0.5d0*(1.d0+press_star)*&
                c_star + (3.d0/5.d0)*EP_SM*Chi*(1.d0+C_E)*(1.d0+C_E) * &
                (2.d0*C_E - 1.d0 + ( 0.5d0*(1.d0+C_E) - 5.d0/&
                (3*(1.d0+C_E))) * c_star ) ) / (nu_kappa_star - &
                2.d0*zeta0_star)

             kappa_star = kappa_k_star * (1.d0 + (6.d0/5.d0)*EP_SM* &
                Chi*(1.d0+C_E) ) + (256.d0/25.d0)*(EP_SM* &
                EP_SM/PI)*Chi*(1.d0+C_E)*(1.d0+(7.d0/32.d0)* &
                c_star)

             IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN ! sof modifications (May 20 2005)
                Kth_star= kappa0
             ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                Kth_star = ZERO
             ELSEIF(EP_SM < DIL_EP_S) THEN
                Kth_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*kappa0/ &
                   (RO_SM*EP_SM*Chi*Theta_m(IJK,M) + 1.2d0*DgA*kappa0/RO_SM)
             ELSE
                Kth_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*kappa0/ &
                   (RO_SM*EP_SM*Chi*Theta_m(IJK,M)+ &
                   (1.2d0*F_gs(IJK,M)*kappa0/(RO_SM*EP_SM)) )
             ENDIF

! granular conductivity in Mth solids phase
             Kth_s(IJK,M) = Kth_star * kappa_star
   
! transport coefficient of the Mth solids phase associated
! with gradient in volume fraction in heat flux
             qmu_k_star = 2.d0*( (1.d0+EP_SM*dChiOdphi)* &
                zeta0_star*kappa_k_star + ( (press_star/3.d0) + &
                (2.d0/3.d0)* EP_SM*(1.d0+C_E)*(Chi+EP_SM* dChiOdphi) )*&
                c_star - (4.d0/5.d0)*EP_SM*Chi* (1.d0+(EP_SM/2.d0)*&
                dChiOdphi)* (1.d0+C_E) * ( C_E*(1.d0-C_E)+0.25d0*&
                ((4.d0/3.d0)+C_E* (1.d0-C_E))*c_star ) ) / &
                (2.d0*nu_kappa_star-3.d0*zeta0_star)

             qmu_star = qmu_k_star*(1.d0+(6.d0/5.d0)*EP_SM*Chi*&
                (1.d0+C_E) )

             IF (EP_SM .LT. SMALL_NUMBER) THEN
                Kphi_s(IJK,M) = ZERO
             ELSE   
                Kphi_s(IJK,M) = (Theta_m(IJK,M)*Kth_star/NU_PM)*qmu_star
             ENDIF

          ENDIF   ! Fluid_at
 200  CONTINUE   ! outer IJK loop

      RETURN
      END SUBROUTINE GT_PDE_GD_99



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_GTSH                                             C
!  Purpose: Implement kinetic theory of Garzo, Tenneti, Subramaniam    C
!     Hrenya (2012) for calculation of granular stress terms and       C
!                                                                      C
!  Author: Sofiane Benyahia                                            C
!                                                                      C
!  Literature/Document References:                                     C
!     Garzo, Tenneti, Subramaniam, Hrenya (2012) J. Fluid Mech.        C
!     712, pp 129-404                                                  C
!                                                                      C
!  Comments:                                                           C
!     And also based on C.M. Hrenya hand-notes dated Sep 2013          C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      Subroutine gt_pde_gtsh (M, IER)

!-----------------------------------------------
! Modules 
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
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell Indices
      INTEGER :: IJK, I, J, K
!
      DOUBLE PRECISION :: D_PM, M_PM, NU_PM, EP_SM
!
      DOUBLE PRECISION :: chi, dChiOdphi, eta0
!
      DOUBLE PRECISION :: nu0, nuN, etaK
      DOUBLE PRECISION :: dZeta_dT, dGama_dT, NuK, Kth0, KthK
      DOUBLE PRECISION :: Rdissdphi, Kphidphi, Re_T, dGamadn, dRdphi
      DOUBLE PRECISION :: denom
      DOUBLE PRECISION :: dSdphi, R_dphi, Tau_st, dPsidn, MuK

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
! (Same for conductivity)

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
! radial distribution function
      DOUBLE PRECISION, EXTERNAL :: G_0
! dg0/dep
      DOUBLE PRECISION, EXTERNAL :: DG_0DNU
! function gamma: eq. (8.1), S_star and K_phi in GTSH theory
      DOUBLE PRECISION, EXTERNAL :: G_gtsh
      DOUBLE PRECISION, EXTERNAL :: S_star
      DOUBLE PRECISION, EXTERNAL :: K_phi
      DOUBLE PRECISION, EXTERNAL :: R_d
!-----------------------------------------------
! Include statement functions
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

! local aliases
             D_PM = D_P(IJK,M)
             EP_SM = EP_S(IJK,M)
             M_PM = (PI/6.d0)*D_PM**3 * RO_S(IJK,M)
             NU_PM = ROP_S(IJK,M)/M_PM
             Chi = G_0(IJK,M,M)
             dChiOdphi = DG_0DNU(EP_SM)

       
! note that T = (m_pm*theta_m)
! Pressure/Viscosity/Bulk Viscosity
!-----------------------------------            
! solids pressure, eq (6.14) of GTSH theory
             P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M)* &
                (one+2d0*(one+C_E)*Chi*EP_SM)

! evaluating shear viscosity, eq (7.3) of GTSH theory
! starting with nu_0 equ (7.6-7.8)
             eta0 = 0.3125d0/(dsqrt(pi)*D_PM**2)*M_pm*&
                dsqrt(theta_m(ijk,m))

             nu0 = (96.d0/5.d0)*(EP_SM/D_PM)*DSQRT(Theta_m(IJK,M)/PI)
!             nu0 = NU_PM*M_pm*theta_m(ijk,m)/eta0

             nuN = 0.25d0*nu0*Chi*(3d0-C_E)*(one+C_E) * &
                    (one+0.7375d0*A2_gtsh(ijk))
! defining kinetic part of shear viscosity nuK  equ (7.7)
             etaK = rop_s(ijk,m)*theta_m(ijk,m) / (nuN-0.5d0*( &
                EDT_s_ip(ijk,M,M)-xsi_gtsh(ijk)/theta_m(ijk,m) - &
                2d0*G_gtsh(EP_SM, Chi, IJK, M)/M_PM)) * (one -0.4d0 * &
                (one+C_E)*(one-3d0*C_E)*EP_SM*Chi)

! bulk viscosity lambda eq. (7.5)
             Mu_b_v(IJK) = 25.6d0/pi * EP_SM**2 * Chi *(one+C_E) * &
                (one - A2_gtsh(ijk)/16d0)*eta0

! Finally shear viscosity, eq (7.9) of GTSH theory
             Mu_s_v(IJK) = etaK*(one+0.8d0*EP_SM*Chi*(one+C_E)) + &
                0.6d0*Mu_b_v(IJK)

! Now let's define the true bulk viscosity as defined in MFIX 
             LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)


! Conductivity/Dufour coefficient
!-----------------------------------
! Calculate conductivity Kth_s(IJK,M), eq. 7.12 GTSH theory. 
! Start with calculating dZeta/dT and dGama/dT
! note that 1/Tau = (3d0*pi*mu_g(ijk)*D_PM/M_p)**2 defined
! under eq. 8.2 GTSH
             dZeta_dT = -0.5d0*xsi_gtsh(ijk)/(M_pm*theta_m(ijk,m))

             dGama_dT = 3d0*pi*D_PM**2*RO_g(ijk)*K_phi(EP_SM)/ &
                (2d0*M_pm*dsqrt(theta_m(ijk,m)))  
             dGama_dT = zero  ! this is giving neg. KthK for dilute 
                              ! flows, set it to zero for now.

! evaluating eq (7.16) in GTSH
             NuK = nu0*(one+C_E)/3d0*Chi*( one+2.0625d0*(one-C_E)+ &
                ((947d0-579*C_E)/256d0*A2_gtsh(ijk)) )

! evaluating eq. (7.13)
             Kth0 = 3.75d0*eta0/M_pm

! evaluating kinetic conductivity Kk eq. (7.14)
! note that 1/2m/T Psi and m dZeta_dT cancel out.
             KthK = zero
             IF(EP_SM > SMALL_NUMBER) KthK = 2d0/3d0*Kth0*nu0 / (NuK - &
                2d0*EDT_s_ip(ijk,M,M) - 2d0*theta_m(ijk,m)*dGama_dT) * &
                (one+2d0*A2_gtsh(ijk)+0.6d0*EP_SM*Chi* &
                (one+C_E)**2*(2*C_E-one+A2_gtsh(ijk)*(one+C_E)))

! the conductivity Kth from eq (7.17) in GTSH theory:
             Kth_s(IJK,M) = KthK*(one+1.2d0*EP_SM*Chi*(one+C_E)) + &
                (10.24d0/pi* EP_SM**2*Chi*(one+C_E)*(one+0.4375d0* &
                A2_gtsh(ijk))*Kth0)

! Finaly notice that conductivity K in eq (7.10) must be 
! multiplied by m because of grad(T)
             Kth_s(IJK,M) = M_pm * Kth_s(IJK,M)

! Calculate the Dufour coefficient Kphi_s(IJK,M) in equation (7.18) of
! GTSH theory.
! First, calculate terms in 2 n/m x Gama_n, dRdiss/dphi and dK_phi/dphi.
! Notice that 2 n/m Gama_n = 2 phi/m Gama_phi, so multiply the
! deriviatives of Rdiss and K_phi by phi to avoid possible division by
! phi.
             Rdissdphi = ZERO
             IF(EP_SM > SMALL_NUMBER) Rdissdphi = &
                1.5d0*dsqrt(EP_SM/2d0)+135d0/64d0*EP_SM*(dlog(EP_SM)+one) +&
                11.26d0*EP_SM*(one-10.2*EP_SM+49.71d0*EP_SM**2-87.08d0* &
                EP_SM**3) - EP_SM*dlog(epM)*(Chi+EP_SM*dChiOdphi)

             Kphidphi = EP_SM*(0.212d0*0.142d0*EP_SM**0.788d0/&
                (one-EP_SM)**4.454d0 - 4.454d0*K_phi(EP_SM)/(one-EP_SM))

             Re_T = ro_g(ijk)*d_p(ijk,m)*dsqrt(theta_m(ijk,m)) / &
                mu_g(ijk)

! The term phi x Gama_phi becomes
             dGamadn = 3d0*pi*D_pm*Mu_g(ijk)*(Rdissdphi+Re_T*Kphidphi)

! Second, calculate terms in 2 rho x Psi_n, which is same as ro_s x 
! EP_SM x Psi_n
! Take EP_SM inside the derivative dS_star/dphi to avoid singularities 
! as EP_SM -> 0

! calculating the term phi*dRd/dphi
             dRdphi = zero
             IF((EP_SM > SMALL_NUMBER) .AND. (EP_SM < 0.4d0)) THEN
                denom = one+0.681d0*EP_SM-8.48d0*EP_SM**2+8.16d0*EP_SM**3
                dRdphi = (1.5d0*dsqrt(EP_SM/2d0)+135d0/64d0*EP_SM*&
                   (dlog(EP_SM)+one)+ 17.14d0*EP_SM)/denom - EP_SM*&
                   (one+3d0*dsqrt(EP_SM/2d0) + 135d0/64d0*EP_SM*&
                   dlog(EP_SM)+17.14*EP_SM)/denom**2 * &
                   (0.681d0-16.96d0*EP_SM+24.48d0*EP_SM**2)
             ELSEIF(EP_SM > 0.4d0) THEN
                dRdphi = 10d0*(one+2d0*EP_SM)/(one-EP_SM)**4
             ENDIF

! calculating the term phi*dS_star/dphi
             dSdphi = zero
             IF(EP_SM >= 0.1d0) THEN
                R_dphi = R_d(EP_SM)
                denom = one+3.5d0*dsqrt(EP_SM)+5.9d0*EP_SM
                dSdphi = 2d0*R_dphi*dRdphi/(Chi*denom) - &
                   EP_SM*R_dphi**2 * (dChiOdphi/(Chi**2*denom) + &
                   (1.75d0/dsqrt(EP_SM)+5.9d0)/(Chi*denom**2))
             ENDIF

! defining the relaxation time Tau_st
             Tau_st = M_pm/(3d0*pi*mu_g(ijk)*D_pm)

! The term phi x Psi_n becomes
             dPsidn = dsqrt(pi)*D_pm**4*VREL_array(IJK)**2 / &
                (36d0*Tau_st**2*dsqrt(theta_m(ijk,m))) * dSdphi

! Now compute the kinetic contribution to Dufour coef. Muk eq (7.20) GSTH
             Muk = ZERO  ! This is assumed to avoid /0 for EP_SM = 0
             IF(EP_SM > SMALL_NUMBER) Muk = Kth0*Nu0*M_pm*&
                theta_m(ijk,m)/NU_PM / (NuK-1.5d0*(EDT_s_ip(ijk,M,M)-&
                xsi_gtsh(ijk)/theta_m(ijk,m))) * ( KthK/(Kth0*Nu0)*&
                (2d0/M_pm*dGamadn-ro_s(IJK,m)/(M_pm* theta_m(ijk,m))*&
                dPsidn + EDT_s_ip(ijk,M,M)*(one+EP_SM/Chi*dChiOdphi)) +&
                2d0/3d0*A2_gtsh(ijk) + 0.8d0*EP_SM*Chi* (one+C_E)*&
                (one+0.5d0*EP_SM/Chi*dChiOdphi)*(C_E*(C_E-one)+ &
                A2_gtsh(ijk)/6d0*(16d0-3d0*C_E+3d0*C_E**2)))

! Finaly compute the Dufour coefficient Mu (Kphi_s(IJK,M)) from eq (7.22) GTSH
             Kphi_s(IJK,M) = Muk*(one+1.2d0*EP_SM*Chi*(one+C_E))

          ENDIF   ! Fluid_at
 200  CONTINUE   ! outer IJK loop

      RETURN
      END SUBROUTINE GT_PDE_GTSH



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_IA_NONEP                                         C
!  Purpose: Implement kinetic theory of Iddir & Arastoopour (2005) for C
!     calculation of granular stress terms and granular conductivity   C
!                                                                      C
!  Author: Janine E. Galvin, Univeristy of Colorado                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Iddir, Y.H., PhD Modeling of the multiphase mixture of particles C
!        using the kinetic theory approach, PhD Dissertation in        C
!        Chemical Engineering, Illinois Institute of Technology, 2004  C
!     Iddir, Y.H., H. Arastoopour, and C.M. Hrenya, Analysis of binary C
!        and ternary granular mixtures behavior using the kinetic      C
!        theory approach. Powder Technology, 2005, p. 117-125.         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine gt_pde_ia_nonep (M, IER)

!-----------------------------------------------
! Modules 
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
      USE kintheory
      USE ur_facs
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------  
! cell indices
      INTEGER :: IJK, I, J, K
! solids phase index
      INTEGER :: L 
! use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Mu_s_dil, Kth_star, K_s_dil, XI_star
! variables for Iddir equipartition model
      DOUBLE PRECISION :: P_s_sum, P_s_MM, P_s_LM
      DOUBLE PRECISION :: MU_common_term, K_common_term
      DOUBLE PRECISION :: Mu_sM_sum, MU_s_MM, MU_s_LM, MU_sM_LM, MU_sL_LM
      DOUBLE PRECISION :: XI_sM_sum, XI_s_v
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R0p_lm, R1p_lm, R8p_lm, R9p_lm, Bp_lm,&
                          R5p_lm, R6p_lm, R7p_lm
      DOUBLE PRECISION :: K_s_sum, K_s_MM, K_s_LM
!      
      DOUBLE PRECISION :: Re, C_d, DgA
! Sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo
! Current value of Kth_sl_ip (i.e., without underrelaxation)
      DOUBLE PRECISION :: Kth_sL_iptmp 
!----------------------------------------------- 
! Function subroutines
!----------------------------------------------- 
! radial distribution function
      DOUBLE PRECISION, EXTERNAL :: G_0
!-----------------------------------------------
! Include statement functions
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
             M_PM = (PI/6.d0)*D_PM**3 * RO_S(IJK,M)
             NU_PM = ROP_S(IJK,M)/M_PM

             P_s_MM = NU_PM*Theta_m(IJK,M)

             MU_s_dil = (5.d0/96.d0)*D_PM* RO_S(IJK,M)*&
                 DSQRT(PI*Theta_m(IJK,M)/M_PM)

             IF(.NOT.SWITCH_IA .OR. RO_G(IJK) == ZERO) THEN 
                Mu_star = MU_s_dil ! do nothing... granular flow
             ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
                Mu_star = ZERO
             ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
                Mu_star = MU_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                    (SUM_EpsGo + 2.0d0*DgA*MU_s_dil &
                    / (RO_S(IJK,M)**2 *(Theta_m(IJK,M)/M_PM)))
             ELSE
                Mu_star = MU_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                    (SUM_EpsGo + 2.0d0*F_gs(IJK,M)*MU_s_dil &
                    / (RO_S(IJK,M)**2 *EP_s(IJK,M)*(Theta_m(IJK,M)/M_PM)))
             ENDIF
   
             MU_s_MM = (Mu_star/G_0(IJK,M,M))*&
                 (1.d0+(4.d0/5.d0)*(1.d0+C_E)*SUM_EpsGo)**2

             DO L = 1, MMAX
                D_PL = D_P(IJK,L)
                M_PL = (PI/6.d0)*D_PL**3 * RO_S(IJK,L)
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
                   R0p_lm = (1.d0/(Ap_lm**1.5 * Dp_lm**2.5))+ &
                       ((15.d0*Bp_lm*Bp_lm)/(2.d0* Ap_lm**2.5 *&
                       Dp_lm**3.5))+&
                       ((175.d0*(Bp_lm**4))/(8.d0*Ap_lm**3.5 * &
                       Dp_lm**4.5))
                   R1p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3)))+ &
                       ((9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4))+&
                       ((30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 * &
                       Dp_lm**5))
  
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
             ENDDO

! Find the term proportional to the identity matrix
! (pressure in the Mth solids phase)
             P_s_v(IJK) = P_s_sum + P_S_MM

! Find the term proportional to the gradient in velocity
! of phase M  (shear viscosity in the Mth solids phase)
             MU_s_v(IJK) = MU_sM_sum + MU_sL_ip(IJK,M,M)
             XI_s_v = XI_sM_sum + XI_sL_ip(IJK,M,M)
     
! bulk viscosity in the Mth solids phase
             LAMBDA_s_v(IJK) = -(2.d0/3.d0)*Mu_s_v(IJK) + XI_s_v


! find the granular conductivity in Mth solids phase      
!----------------------------------- 
            K_s_sum = ZERO

            K_s_dil = (75.d0/384.d0)*D_PM* RO_S(IJK,M)*&
                 DSQRT(PI*Theta_m(IJK,M)/M_PM)
               
             IF(.NOT.SWITCH_IA .OR. RO_G(IJK) == ZERO) THEN 
                Kth_star = K_s_dil ! do nothing... granular flow
             ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
                Kth_star = ZERO
        
             ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
                Kth_star = K_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                    (SUM_EpsGo+ 1.2d0*DgA*K_s_dil &
                    / (RO_S(IJK,M)**2 *(Theta_m(IJK,M)/M_PM)))
             ELSE
                Kth_star = K_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                    (SUM_EpsGo+ 1.2d0*F_gs(IJK,M)*K_s_dil &
                    / (RO_S(IJK,M)**2 *EP_s(IJK,M)*(Theta_m(IJK,M)/M_PM)))
             ENDIF

! Kth doesn't include the mass.      
             K_s_MM = (Kth_star/(M_PM*G_0(IJK,M,M)))*&  
                 (1.d0+(3.d0/5.d0)*(1.d0+C_E)*(1.d0+C_E)*SUM_EpsGo)**2

             DO L = 1, MMAX
                D_PL = D_P(IJK,L)
                M_PL = (PI/6.d0)*D_PL**3 *RO_S(IJK,L)
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
                   Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/2.d0
                   Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                       Theta_m(IJK,M) ))/(2.d0*MPSUM)
                   Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+&
                       M_PL*Theta_m(IJK,L) ))/(2.d0*MPSUM*MPSUM)

                   R0p_lm = (1.d0/(Ap_lm**1.5 * Dp_lm**2.5))+&
                       ((15.d0*Bp_lm*Bp_lm)/(2.d0* Ap_lm**2.5 * Dp_lm**3.5))+&
                       ((175.d0*(Bp_lm**4))/(8.d0*Ap_lm**3.5 * Dp_lm**4.5))

                   R1p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3)))+ &
                       ((9.d0*Bp_lm*Bp_lm)/(Ap_lm**2.5 * Dp_lm**4))+&
                       ((30.d0*Bp_lm**4)/(2.d0*Ap_lm**3.5 * Dp_lm**5))

                   R5p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**3 ) )+ &
                       ((5.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**4))+&
                       ((14.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**5))

                   R6p_lm = (1.d0/(Ap_lm**3.5 * Dp_lm**3))+ &
                       ((7.d0*Bp_lm*Bp_lm)/(Ap_lm**4.5 * Dp_lm**4))+&
                       ((126.d0*Bp_lm**4)/(5.d0*Ap_lm**5.5 * Dp_lm**5))

                   R7p_lm = (3.d0/(2.d0*Ap_lm**2.5 * Dp_lm**4))+ &
                       ((10.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**5))+&
                       ((35.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**6))

                   R8p_lm = (1.d0/(2.d0*Ap_lm**1.5 * Dp_lm**4))+ &
                       ((6.d0*Bp_lm*Bp_lm)/(Ap_lm**2.5 * Dp_lm**5))+&
                       ((25.d0*Bp_lm**4)/(Ap_lm**3.5 * Dp_lm**6))

                   R9p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**3))+ &
                       ((15.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**4))+&
                       ((70.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**5))
                        
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
                    !Kth_sL_ip(IJK,M,L) = K_common_term*NU_PM*NU_PL*(&
                   Kth_sl_iptmp = K_common_term*NU_PM*NU_PL*(&
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

                   Kth_sL_ip(IJK,M,L) = (ONE-UR_Kth_sml)*Kth_sL_ip(IJK,M,L) +&
                       UR_Kth_sml * Kth_sl_iptmp

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

! granular conductivity in Mth solids phase
             Kth_s(IJK,M) = K_s_sum

          ENDIF   ! Fluid_at
 200  CONTINUE   ! outer IJK loop

      RETURN
      END SUBROUTINE GT_PDE_IA_NONEP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FRICTION_PRINCETON                                      C
!  Purpose: Calculate frictional contributions to granular stress      C
!     terms based on model by Srivastava & Sundaresan (2003)           C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University    Date: 20-APR-98    C
!                                                                      C
!  Literature/Document References:                                     C
!     Srivastava, A., and Sundaresan, S., Analysis of a frictional-    C
!        kinetic model for gas-particle flow, Powder Technology,       C
!        2003, 129, 72-85.                                             C
!     Benyahia, S., Validation study of two continuum granular         C
!        frictional flow theories, Industrial & Engineering Chemistry  C
!        Research, 2008, 47, 8926-8932.                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      
      Subroutine Friction_princeton(M, IER)

!-----------------------------------------------
! Modules 
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
      USE trace
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!----------------------------------------------- 
! Local variables
!-----------------------------------------------
! cell index
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! used to compute frictional terms
      DOUBLE PRECISION :: Chi, Pc, Mu_zeta,Phin,PfoPc, N_Pff
      DOUBLE PRECISION :: ZETA
! sum of all solids volume fractions
      DOUBLE PRECISION :: SUM_EPS_CP
!     parameters in pressure linearization; simple averaged Dp
      DOUBLE PRECISION :: dpc_dphi, dp_avg
!-----------------------------------------------
! Function subroutines
!-----------------------------------------------
! radial distribution function      
      DOUBLE PRECISION, EXTERNAL :: G_0
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!----------------------------------------------- 

      DO 200 IJK = ijkstart3, ijkend3       
     
         IF ( FLUID_AT(IJK) ) THEN

! close_packed was added for consistency with the Schaeffer model
! I'm also extending this model in the case where more than 1 solids
! phase are used (not completed yet). sof May24 2005.

               IF (EP_g(IJK) .LT. (ONE-eps_f_min)) THEN
     
! part copied from source_v_s.f (sof)
                  SUM_EPS_CP = ZERO
                  dp_avg = ZERO
                  DO MM=1,SMAX
                  dp_avg = dp_avg + D_p(IJK,MM)
                     IF (CLOSE_PACKED(MM)) SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
                  END DO
                  dp_avg = dp_avg/DFLOAT(SMAX)
! end of part copied
     
                  IF (SAVAGE.EQ.1) THEN !form of Savage (not to be used with GHD theory)
                     Mu_zeta = &
                        ((2d0+ALPHA)/3d0)*((Mu_s_v(IJK)/(Eta*(2d0-Eta)*&
                        G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                        G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                        EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b_v(IJK)*Eta))
                     ZETA = &
                        ((48d0*Eta*(1d0-Eta)*RO_S(IJK,M)*EP_s(IJK,M)*&
                        EP_s(IJK,M)*G_0(IJK,M,M)*&
                        (Theta_m(IJK,M)**1.5d0))/&
                        (SQRT_Pi*D_p(IJK,M)*2d0*Mu_zeta))**0.5d0
                     
                  ELSEIF (SAVAGE.EQ.0) THEN !S:S form
                     ZETA = (SMALL_NUMBER +&
                        trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                        trD_s_C(IJK,M))/3.d0))**0.5d0
                     
                  ELSE          !combined form
                     IF(TRIM(KT_TYPE) == 'GHD') THEN
                       ZETA = ((Theta_m(IJK,M)/dp_avg**2) +&
                          (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                          trD_s_C(IJK,M))/3.d0)))**0.5d0
                     ELSE
                       ZETA = ((Theta_m(IJK,M)/D_p(IJK,M)**2) +&
                          (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                          trD_s_C(IJK,M))/3.d0)))**0.5d0
                     ENDIF                     
                  ENDIF

                  
                  IF ((ONE-EP_G(IJK)) .GT. ((ONE-ep_star_array(ijk))-delta)) THEN
! Linearized form of Pc; this is more stable and provides continuous function.
                     dpc_dphi = (to_SI*Fr)*((delta**5)*(2d0*(ONE-&
                         ep_star_array(IJK)-delta) - 2d0*eps_f_min)+&
                         ((ONE-ep_star_array(ijk)-delta)-eps_f_min)*&
                         (5*delta**4))/(delta**10)

                     Pc = (to_SI*Fr)*(((ONE-ep_star_array(IJK)-delta) -&
                        EPS_f_min)**N_Pc)/(delta**D_Pc)
                     Pc = Pc + dpc_dphi*((ONE-EP_G(IJK))+delta-(ONE-&
                        ep_star_array(IJK)))

                   ! Pc = 1d25*(((ONE-EP_G(IJK))- (ONE-ep_star_array(ijk)))**10d0) ! old commented Pc
                  ELSE
                     Pc = (to_SI*Fr)*(((ONE-EP_G(IJK)) - EPS_f_min)**N_Pc) / &
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
                  Lambda_s_f(IJK) = -2d0*Mu_s_f(IJK)/3d0
     
! modification of the stresses in case of more than one solids phase are used (sof)
! This is NOT done when mixture mom. eq. are solved (i.e. for GHD theory)
                  IF(TRIM(KT_TYPE) /= 'GHD') THEN
                     P_s_f(IJK) = P_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                     Mu_s_f(IJK) = Mu_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                     Lambda_s_f(IJK) = Lambda_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                  ENDIF
               
            ENDIF   ! end if ep_g < 1-eps_f_min
         ENDIF   ! Fluid_at
 200  CONTINUE   ! outer IJK loop
      

      RETURN
      END SUBROUTINE FRICTION_PRINCETON



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_IGCI                                     C
!  Purpose: Calculate solids viscosity and pressure using subgrid      C
!     model                                                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes & make consistent with variable density      C
!     feature.                                                         C
!  Author: Janine Galvin, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., Pannala, S., Benyahia, S., & Sundaresan S.,            C
!        Validation studies on filtered model equations for gas-       C
!        particle flows in risers, Industrial & Engineering Chemistry  C
!        Research, 2012, 51(4), 2094-2103                              C
!                                                                      C
!  Comments:                                                           C
!     Still needs to be reviewed for accuracy with source material     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine subgrid_stress_igci(M, IER)

!-----------------------------------------------
! Modules 
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
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!----------------------------------------------- 
! Local variables
!-----------------------------------------------  
! cell indices
      INTEGER :: IJK, I, J, K
! Igci models
      DOUBLE PRECISION :: Mu_sub,ps_sub
      DOUBLE PRECISION :: pressurefac
      DOUBLE PRECISION :: viscosityfac,ps_kinetic,ps_total, extra_visfactor
      DOUBLE PRECISION :: mu_kinetic,mu_total
! factors to correct subgrid effects from the wall
      DOUBLE PRECISION :: wfactor_Ps, wfactor_mus
! the inverse Froude number, or dimensionless Filtersize
      DOUBLE PRECISION :: Inv_Froude
! one particle terminal settling velocity
      DOUBLE PRECISION :: vt
! the filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
!----------------------------------------------- 
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!----------------------------------------------- 
  
      DO IJK = ijkstart3, ijkend3

         IF ( FLUID_AT(IJK) ) THEN
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 

! initialize
            wfactor_Ps = ONE   ! for P_S
            wfactor_mus = ONE   ! for MU_s

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
            vt = GRAVITY*D_p0(M)*D_p0(M)*(RO_S(IJK,M) - RO_g(IJK)) / &
               (18.0d0*MU_G(IJK))

! FilterSIZE calculation for each specific gridcell volume
            IF(DO_K) THEN  
               filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
            ELSE
               filtersize = filter_size_ratio * DSQRT(AXY(IJK))
            ENDIF


! dimensionless inverse of Froude number
            Inv_Froude =  filtersize * GRAVITY / vt**2     

! various factor needed:
            pressurefac = 0.48d0*(Inv_Froude**0.86)*&
               (ONE-EXP(-Inv_Froude/1.4))
            viscosityfac = 0.37d0*(Inv_Froude**1.22)
            Extra_visfactor = ONE/(0.28d0*(Inv_Froude**0.43)+ONE)

            IF (EP_s(IJK,M) .LE. 0.0131) THEN
               Ps_kinetic = -10.4d0*(EP_s(IJK,m)**2)+0.31d0*EP_s(IJK,m)
            ELSEIF (EP_s(IJK,M) .LE. 0.290) THEN
               Ps_kinetic = -0.185d0*(EP_s(IJK,m)**3)+&
                  0.066d0*(EP_s(IJK,m)**2)-0.000183d0*EP_s(IJK,m)+&
                  0.00232d0
            ELSEIF (EP_s(IJK,M) .LE. 0.595) THEN
               Ps_kinetic = -0.00978d0*EP_s(IJK,m)+0.00615d0
            ELSE
               Ps_kinetic = -6.62d0*(EP_s(IJK,m)**3)+&
                  49.5d0*(EP_s(IJK,m)**2)-50.3d0*EP_s(IJK,m)+13.8d0
            ENDIF

            Ps_sub = pressurefac*(EP_s(IJK,M)-0.59d0)*&
               (-1.69d0*EP_s(IJK,M)-4.61d0*(EP_s(IJK,M)**2)+&
               11.d0*(EP_s(IJK,M)**3))

            IF (Ps_sub .GE. ZERO) THEN
               Ps_total=Ps_kinetic+Ps_sub 
            ELSE
               Ps_total=Ps_kinetic  
            ENDIF

            IF (EP_s(IJK,M) .LE. 0.02) THEN
               Mu_kinetic = 1720.d0*(EP_s(IJK,m)**4)-&
                  215.d0*(EP_s(IJK,m)**3) + 9.81d0*(EP_s(IJK,m)**2)-&
                  0.207d0*EP_s(IJK,m)+0.00254d0
            ELSEIF (EP_s(IJK,M) .LE. 0.2) THEN
               Mu_kinetic = 2.72d0*(EP_s(IJK,m)**4)-&
                  1.55d0*(EP_s(IJK,m)**3)+0.329d0*(EP_s(IJK,m)**2)-&
                  0.0296d0*EP_s(IJK,m)+0.00136d0
            ELSEIF (EP_s(IJK,M) .LE. 0.6095) THEN
               Mu_kinetic = -0.0128d0*(EP_s(IJK,m)**3)+&
                  0.0107d0*(EP_s(IJK,m)**2)-0.0005d0*EP_s(IJK,m)+&
                  0.000335d0
            ELSE
               Mu_kinetic = 23.6d0*(EP_s(IJK,m)**2)-&
                  28.0d0*EP_s(IJK,m)+8.30d0
            ENDIF

            Mu_sub = Extra_visfactor*viscosityfac*&
               (EP_s(IJK,M)-0.59d0)*(-1.22d0*EP_s(IJK,M)-&
               0.7d0*(EP_s(IJK,M)**2)-2.d0*(EP_s(IJK,M)**3))
   
            IF (Mu_sub .GE. ZERO) THEN
               Mu_total = Mu_kinetic+Mu_sub 
            ELSE
               Mu_total = Mu_kinetic  
            ENDIF

            IF (SUBGRID_WALL) THEN
               CALL SUBGRID_STRESS_WALL(wfactor_Ps,wfactor_Mus,vt,IJK)
            ENDIF

! pressure
            P_s_v(IJK) = Ps_total * wfactor_Ps * (vt**2) * &
               RO_S(IJK,M)

! shear viscosity
            Mu_s_v(IJK) = Mu_total * wfactor_mus * (vt**3) * &
               RO_S(IJK,M)/GRAVITY
         
! set an arbitrary value in case value gets negative (this should 
! not happen unless filtersize becomes unrelastic w.r.t. gridsize)
            IF (P_s_v(IJK) .LE. SMALL_NUMBER) P_s_v(IJK) = SMALL_NUMBER
            IF (Mu_s_v(IJK) .LE. SMALL_NUMBER) Mu_s_v(IJK)= SMALL_NUMBER

! solid second viscosity, assuming the bulk viscosity is ZERO
            lambda_s_v(IJK) = (-2.0d0/3.0d0)*Mu_s_v(IJK)

! granular temperature is zeroed in all LES/Subgrid model
            THETA_m(IJK, M) = ZERO
            
         ENDIF   ! endif (fluid_at(IJK))
         
      ENDDO   ! outer IJK loop

      RETURN
      END SUBROUTINE SUBGRID_STRESS_IGCI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_MILIOLI                                  C
!  Purpose: Calculate solids viscosity and pressure using subgrid      C
!     model                                                            C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: Minor changes & make consistent with variable density      C
!     feature.                                                         C
!  Author: Janine Galvin, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Milioli, C. C., et al., Filtered two-fluid models of fluidized   C
!        gas-particle flows: new constitutive relations, AICHE J,      C
!        doi: 10.1002/aic.14130                                        C
!                                                                      C
!  Comments:                                                           C
!     Still needs to be reviewed for accuracy with source material     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      Subroutine subgrid_stress_MILIOLI(M, IER)

!-----------------------------------------------
! Modules 
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
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!----------------------------------------------- 
! Local variables
!-----------------------------------------------  
! cell indices
      INTEGER :: IJK, I, J, K
! Milioli model
      DOUBLE PRECISION :: cvisc_pot,cvisc_num,cvisc_den,Cvisc
      DOUBLE PRECISION :: cpress_pot,cpress_num,cpress_den,Cpress
! factor functions to correct subgrid effects from the wall
      DOUBLE PRECISION :: wfactor_Ps, wfactor_mus
! the inverse Froude number, or dimensionless Filtersize
      DOUBLE PRECISION :: Inv_Froude
! one particle terminal settling velocity
      DOUBLE PRECISION :: vt
! the filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
!----------------------------------------------- 
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!----------------------------------------------- 
  
      DO IJK = ijkstart3, ijkend3

         IF ( FLUID_AT(IJK) ) THEN
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 

! initialize
            wfactor_Ps = ONE   ! for P_S
            wfactor_mus = ONE   ! for MU_s

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
            vt = GRAVITY*D_p0(M)*D_p0(M)*(RO_S(IJK,M) - RO_g(IJK)) / &
               (18.0d0*MU_G(IJK))

! FilterSIZE calculation for each specific gridcell volume
            IF(DO_K) THEN  
               filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
            ELSE
               filtersize = filter_size_ratio * DSQRT(AXY(IJK))
            ENDIF

! dimensionless inverse of Froude number
            Inv_Froude =  filtersize * GRAVITY / vt**2     

! Cvisc:
            cvisc_pot = (0.59d0-(ONE-EP_g(IJK)))
            cvisc_num = (0.7d0*(ONE-EP_g(IJK))*cvisc_pot)
            cvisc_den = (0.8d0+(17.d0*cvisc_pot*cvisc_pot*cvisc_pot))
            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
                (ONE-EP_g(IJK)) .LE. 0.59) THEN
               cvisc=(cvisc_num/cvisc_den)
            ELSE
               cvisc=ZERO
            ENDIF

! aCvisc:
!            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
!                (ONE-EP_g(IJK)) .LE. 0.59) THEN
!               acvisc(IJK) = (0.7d0*(ONE-EP_g(IJK))*(0.59d0-&
!                  (ONE-EP_g(IJK))))/(0.8d0+17.d0*(0.59d0-&
!                  (ONE-EP_g(IJK)))*(0.59d0-(ONE-EP_g(IJK)))*&
!                  (0.59d0-(ONE-EP_g(IJK))))
!            ELSE
!               acvisc(IJK)=ZERO
!            ENDIF

! Cpress:
            cpress_pot = (0.59d0-(ONE-EP_g(IJK)))
            cpress_num = (0.4d0*(ONE-EP_g(IJK))*cpress_pot)
            cpress_den = (0.5d0+(13.d0*cpress_pot*cpress_pot*cpress_pot))
            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
                (ONE-EP_g(IJK)) .LE. 0.59) THEN
               cpress = (cpress_num/cpress_den)
            ELSE
               cpress = ZERO
            ENDIF

! aCpress
!            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
!                (ONE-EP_g(IJK)) .LE. 0.59) THEN
!                acpress(IJK) = (0.4d0*(ONE-EP_g(IJK))*(0.59d0-&
!                   (ONE-EP_g(IJK))))/(0.5d0+13.d0*(0.59d0-&
!                   (ONE-EP_g(IJK)))*(0.59d0-(ONE-EP_g(IJK)))*&
!                   (0.59d0-(ONE-EP_g(IJK))))
!            ELSE
!            acpress(IJK)=ZERO
!            ENDIF

            IF (SUBGRID_WALL) THEN
               CALL SUBGRID_STRESS_WALL(wfactor_Ps,wfactor_Mus,vt,IJK)
            ENDIF

! solid filtered pressure
            P_s_v(IJK) = RO_S(IJK,M) * Inv_froude**(2/7) * &
               filtersize**2 * DSQRT( I2_devD_s(IJK) )**2 * &
               cpress * wfactor_Ps    !16/7-2=2/7 in [Pa or kg/m.s2]

! solids filtered shear viscosity
            Mu_s_v(IJK) = RO_S(IJK,M) * filtersize**2 * &
               DSQRT( I2_devD_s(IJK) ) * cvisc * wfactor_mus  ! [kg/m.s]

! set an arbitrary value in case value gets negative (this should 
! not happen unless filtersize becomes unrelastic w.r.t. gridsize)
            IF (P_s_v(IJK) .LE. SMALL_NUMBER) P_s_v(IJK) = SMALL_NUMBER
            IF (Mu_s_v(IJK) .LE. SMALL_NUMBER) Mu_s_v(IJK)= SMALL_NUMBER

! solid second viscosity, assuming the bulk viscosity is ZERO
            lambda_s_v(IJK) = (-2.0d0/3.0d0)*Mu_s_v(IJK)

! granular temperature is zeroed in all LES/Subgrid model
            THETA_m(IJK, M) = ZERO
            
         ENDIF   ! endif (fluid_at(IJK))
         

      ENDDO   ! outer IJK loop

      RETURN
      END SUBROUTINE SUBGRID_STRESS_MILIOLI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_WALL                                     C
!  Purpose: Calculate subgrid corrections arising from wall to solids  C
!     viscosity and pressure.                                          C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Revision: 1                                                         C
!  Author: Janine Galvin, June 2013                                    C
!                                                                      C
!  Literature/Document References:                                     C
!     Igci, Y., and Sundaresan, S., Verification of filtered two-      C
!        fluid models for gas-particle flows in risers, AICHE J.,      C
!        2011, 57 (10), 2691-2707.                                     C
!                                                                      C
!  Comments: Currently only valid for free-slip wall but no checks     C
!     are made to ensure user has selected free-slip wall when this    C
!     option is invoked                                                C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine subgrid_stress_wall(lfactor_ps, lfactor_mus, vt, &
                 IJK)

!-----------------------------------------------
! Modules 
!-----------------------------------------------
      USE param
      USE param1
      USE constant, only : GRAVITY
      USE cutcell, only : DWALL
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! factor to correct the solids pressure 
      DOUBLE PRECISION, INTENT(OUT) :: lfactor_ps
! factor to correct the solids viscosity
      DOUBLE PRECISION, INTENT(OUT) :: lfactor_mus
! one particle terminal settling velocity
      DOUBLE PRECISION, INTENT(IN) :: vt
! current ijk index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! values are only correct for free-slip walls
      DOUBLE PRECISION, PARAMETER :: aps=9.14d0, bps=0.345d0,&
                                     amus=5.69d0, bmus=0.228d0      
!----------------------------------------------- 
! Local variables
!-----------------------------------------------  
! dimensionless distance to the wall
      DOUBLE PRECISION :: x_d
!----------------------------------------------- 

! initialize
      lfactor_ps = ONE
      lfactor_mus = ONE

! dimensionless distance to the Wall
      x_d = DWALL(IJK) * GRAVITY / vt**2  
! wall function for pressure
      lfactor_Ps = ONE / ( ONE + aps * (EXP(-bps*x_d)) )
! wall function for viscosity
      lfactor_mus = ONE / ( ONE + amus * (EXP(-bmus*x_d)) )

      RETURN
      END SUBROUTINE SUBGRID_STRESS_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine add_shear(M) 

!-----------------------------------------------
! Modules 
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell index
      INTEGER :: IJK
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3         
         IF (FLUID_AT(IJK)) THEN  
            V_s(ijk,m)=V_s(IJK,m)+VSH(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD_SHEAR



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine remove_shear(M)

!-----------------------------------------------
! Modules 
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE vshear
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! cell index
      INTEGER :: IJK
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
            V_s(IJK,m)=V_s(IJK,m)-VSH(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE REMOVE_SHEAR



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Subroutine: INIT_MU_S                                                C
!                                                                      C
! Revision:                                                            C
! Handan Liu added OpenMP, 2012-2013                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      Subroutine init_mu_s (M,IER)

!-----------------------------------------------
! Modules 
!-----------------------------------------------     
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
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------  
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION :: D_s(3,3), D_sl(3,3)
! U_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION :: U_s_N, Usl_N
! U_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION :: U_s_S, Usl_S
! U_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION :: U_s_T, Usl_T
! U_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION :: U_s_B, Usl_B
! U_s at the center of the THETA cell-(i, j, k)
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: U_s_C, Usl_C
! V_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION :: V_s_E, Vsl_E
! V_s at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION :: V_s_W, Vsl_W
! V_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION :: V_s_T, Vsl_T
! V_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION :: V_s_B, Vsl_B
! W_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION :: W_s_E, Wsl_E
! W_s at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION :: W_s_W, Wsl_W
! W_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION :: W_s_N, Wsl_N
! W_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION :: W_s_S, Wsl_S
! W_s at the center of the THETA cell-(i, j, k).
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: W_s_C, Wsl_C
! Cell center value of solids and gas velocities 
      DOUBLE PRECISION :: USCM, UGC, VSCM, VGC, WSCM, WGC,&
                          SqrtVs, SqrtVgMinusVs 
! Local DO-LOOP counters and phase index
      INTEGER :: I1, I2, MM
! d(EP_sm)/dX
      DOUBLE PRECISION :: DEP_soDX
! d(EP_sm)/dY
      DOUBLE PRECISION :: DEP_soDY
! d(EP_sm)/XdZ
      DOUBLE PRECISION :: DEP_soXDZ
! Solids volume fraction gradient tensor
      DOUBLE PRECISION :: M_s(3,3)
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,&
                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,&
                 IM, JM, KM
      INTEGER :: IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP,&
                 IJMKM, IJPKM
! solids phase index
      INTEGER :: L
! shear related reciprocal time scale
      DOUBLE PRECISION :: SRT

!-----------------------------------------------          
!     Functions
!-----------------------------------------------          
      DOUBLE PRECISION, EXTERNAL :: G_0
!----------------------------------------------- 
!     Include statement functions
!----------------------------------------------- 
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
!----------------------------------------------- 

      IF(MU_s0 == UNDEFINED) THEN ! fixes a bug noted by VTech
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
      ENDIF

      IF (SHEAR) SRT=(2d0*V_sh/XLENGTH)

!$omp  parallel do default(shared)                                             &
!$omp  private( I, J, K, IJK, IM, JM, KM, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,  &
!$omp           IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, IMJPK, IMJMK, IMJKM,IMJKP, &
!$omp           IJPKM, IJMKM, IJMKP, IPJMK, IPJKM, U_s_N, U_s_S, U_s_T, U_s_B, &
!$omp           V_s_E, V_s_W, V_s_T, V_s_B, W_s_N, W_s_S, W_s_E, W_s_W, U_s_C, &
!$omp           W_s_C, D_s, L, Usl_N, Usl_S, Usl_T, Usl_B, Vsl_E, Vsl_W, Vsl_T,&
!$omp           Vsl_B, Wsl_n, Wsl_S, Wsl_E, Wsl_W, Usl_C, Wsl_C, D_sl,         &
!$omp           UGC, VGC, WGC, USCM, VSCM, WSCM, SqrtVs, SqrtVgMinusVs,        &
!$omp           DEP_soDX, DEP_soDY, DEP_soXDZ, M_s, I1, I2)       

      DO IJK = ijkstart3, ijkend3 

         IF ( FLUID_AT(IJK) ) THEN
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
            ENDIF

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
     
! Check for IS surfaces and modify solids velocity-comp accordingly
            IF(ANY_IS_DEFINED) THEN
               IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) U_s_N = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
               IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) U_s_S = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
               IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) U_s_T = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
               IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) U_s_B = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
               IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) V_s_E = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
               IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) V_s_W = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
               IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) V_s_T = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
               IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) V_s_B = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
               IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) W_s_N = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
               IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) W_s_S = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
               IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) W_s_E = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
               IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) W_s_W = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
            ENDIF

! Find components of Mth solids phase continuum strain rate
! tensor, D_s, at center of THETA cell-(i, j, k)
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

            IF(CUT_CELL_AT(IJK))  CALL CG_CALC_VEL_S_GRAD(IJK,M,D_s, IER)

! Calculate the trace of D_s
            trD_s_C(IJK,M) = D_s(1,1) + D_s(2,2) + D_s(3,3)

! Calculate trace of the square of D_s
            trD_s2(IJK,M) = 0.D0 !Initialize the totalizer
            DO I1 = 1,3
               DO I2 = 1,3
                  trD_s2(IJK,M) = trD_s2(IJK,M) + D_s(I1,I2)*D_s(I1,I2)
               ENDDO
            ENDDO
 
! use this fact to prevent underflow during theta calculation
            IF (trD_s2(IJK,M) == zero)trD_s_C(IJK,M) = zero 

! The trace of D_sm dot D_sl is required in the implementation of
! Iddir's (2004) kinetic theory 
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               DO L = 1,MMAX
                  IF (L .NE. M) THEN
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

! Check for IS surfaces and modify solids velocity-comp accordingly
                        IF(ANY_IS_DEFINED) THEN
                           IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) &
                              Usl_N = AVG_X_E(U_s(IMJK,L),U_s(IJK,L), I)
                           IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) &
                              Usl_S = AVG_X_E(U_s(IMJK,L),U_s(IJK,L), I)
                           IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) &
                              Usl_T = AVG_X_E(U_s(IMJK,L),U_s(IJK,L), I)
                           IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) &
                              Usl_B = AVG_X_E(U_s(IMJK,L),U_s(IJK,L), I)
                           IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) &
                              Vsl_E = AVG_Y_N(V_s(IJMK,L),V_s(IJK,L))
                           IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) &
                              Vsl_W = AVG_Y_N(V_s(IJMK,L),V_s(IJK,L))
                           IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) &
                              Vsl_T = AVG_Y_N(V_s(IJMK,L),V_s(IJK,L))
                           IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) &
                              Vsl_B = AVG_Y_N(V_s(IJMK,L),V_s(IJK,L))
                           IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) &
                              Wsl_N = AVG_Z_T(W_s(IJKM,L),W_s(IJK,L))
                           IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) &
                              Wsl_S = AVG_Z_T(W_s(IJKM,L),W_s(IJK,L))
                           IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) &
                              Wsl_E = AVG_Z_T(W_s(IJKM,L),W_s(IJK,L))
                           IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) &
                              Wsl_W = AVG_Z_T(W_s(IJKM,L),W_s(IJK,L))
                        ENDIF

! Find components of Lth solids phase continuum strain rate
! tensor, D_sl, at center of THETA cell-(i, j, k) 
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
     
! Calculate trace of the D_sl dot D_sm 
! (normal matrix multiplication)
                        trD_s2_ip(IJK,M,L) = 0.D0 ! initialize 
                        DO I1 = 1,3
                           DO I2 = 1,3
                              trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,M,L)+&
                                 D_sl(I1,I2)*D_s(I1,I2)  
                           ENDDO
                        ENDDO

                     ELSE   ! elseif (m<=m)
                        trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,L,M)
                     ENDIF ! end if/else (l>m)
                  ELSE   ! elseif (L=M)   
                     trD_s2_ip(IJK,M,M) = trD_s2(IJK,M)
                  ENDIF   ! end if/else (m.NE.l)
               ENDDO   ! end do (l=1,mmax)
            ENDIF   ! endif (kt_type = IA theory)


! Start definition of Relative Velocity
! Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            USCM = AVG_X_E(U_S(IMJK,1),U_S(IJK,1),I) 
            VSCM = AVG_Y_N(V_S(IJMK,1),V_S(IJK,1)) 
            WSCM = AVG_Z_T(W_S(IJKM,1),W_S(IJK,1)) 
     
! magnitude of gas-solids relative velocity
            VREL_array(IJK) = SQRT((UGC - USCM)**2 + &
                (VGC - VSCM)**2 + (WGC - WSCM)**2)

! Frictional-flow stress tensor
! Calculate the second invariant of the deviator of D_s
            I2_devD_s(IJK) = ( (D_s(1,1)-D_s(2,2))**2&
                +(D_s(2,2)-D_s(3,3))**2&
                +(D_s(3,3)-D_s(1,1))**2 )/6.&
                + D_s(1,2)**2 + D_s(2,3)**2 + D_s(3,1)**2


            IF(SIMONIN) THEN
! parameters for defining Tau_12: time-scale of the fluid turbulent motion
! viewed by the particles (crossing trajectory effect)
               SqrtVs = SQRT(USCM**2+VSCM**2+WSCM**2)
               SqrtVgMinusVs = SQRT((UGC-USCM)**2+(VGC-VSCM)**2+(WGC-WSCM)**2)
               IF(SqrtVs > Small_Number .AND. SqrtVgMinusVs > &
                  Small_Number .AND. EP_S(IJK,1) > ZERO_EP_S) THEN
                  Cos_Theta(IJK) = ((UGC-USCM)*USCM+(VGC-VSCM)*VSCM+&
                     (WGC-WSCM)*WSCM)/ (SqrtVgMinusVs * SqrtVs)
               ELSE
                  Cos_Theta(IJK) = ZERO ! no solids -> tau_12 = tau_1
               ENDIF
            ENDIF

            IF(.NOT.GRANULAR_ENERGY) THEN 
! algebraic granular energy equation                    
               IF(EP_g(IJK) .GE. EP_star_array(IJK)) THEN
! Boyle-Massoudi Stress term
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
                      DO I1 = 1,3
                         DO I2 = 1,3
                            trDM_s(IJK) = trDM_s(IJK) + D_s(I1,I2)*M_s(I1,I2)
                         ENDDO
                      ENDDO 
                  ENDIF   ! end if(v_ex .ne. zero)
               ENDIF   ! end if (ep_g >=ep_star_array)
            ENDIF   ! end if (.not.granular_energy)


         ENDIF   ! end if (fluid_at)
      ENDDO   ! end outer IJK loop
!$omp end parallel do
      
      RETURN
      END SUBROUTINE INIT_MU_S 

