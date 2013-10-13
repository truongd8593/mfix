!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_CONSTProp                                           C
!  Purpose: This module sets all the constant physical properties      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-MAY-97  C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_CONSTPROP 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE visc_s
      USE visc_g
      USE energy
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE run
      USE funits 
      USE drag
      USE compar 
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices      
      INTEGER :: IJK, M, N, I, J
      DOUBLE PRECISION old_value, DP_TMP(MMAX)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

! Initialize transport coefficients to zero everywhere
      MU_gt = ZERO
      LAMBDA_GT = ZERO
      MU_s = ZERO
      LAMBDA_s_c = ZERO
      LAMBDA_s = ZERO
      K_g = ZERO
      K_s = ZERO
      DIF_g = ZERO
      DIF_S = ZERO
      F_GS = ZERO
      F_SS = ZERO

! Set the flag for recalculating gas viscosity.
      RECALC_VISC_G = (ENERGY_EQ .OR. L_SCALE0/=ZERO .OR. K_EPSILON)

! Set default value for virtual mass coefficient
      Cv = HALF

! Variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) 
! kinetic theory
      IF (TRIM(KT_TYPE) == 'IA_NONEP') THEN
         MU_sM_ip = ZERO
         MU_sL_ip = ZERO
         XI_sM_ip = ZERO
         XI_sL_ip = ZERO
         Fnu_s_ip = ZERO
         FT_sM_ip = ZERO
         FT_sL_ip = ZERO
         Kth_sL_ip = ZERO
         Knu_sM_ip = ZERO
         Knu_sL_ip = ZERO
         Kvel_s_ip = ZERO
         ED_ss_ip = ZERO
         EDvel_sL_ip = ZERO
      ENDIF
      IF (TRIM(KT_TYPE) == 'IA_NONEP' .OR. TRIM(KT_TYPE) == 'GD_99') THEN
         EDT_s_ip = ZERO
         EDvel_sM_ip = ZERO
      ENDIF
     
! Set specified constant physical properties values
      DO IJK = ijkstart3, ijkend3

! All wall cells: FLAG >= 100
         IF (WALL_AT(IJK)) THEN
            RO_G(IJK) = ZERO 
            MU_G(IJK) = ZERO 
            K_G(IJK) = ZERO 
            C_PG(IJK) = ZERO 
            MW_MIX_G(IJK) = ZERO 
         ELSE
! Fluid and inflow/outlfow cells: FLAG < 100
            IF (RO_G0 /= UNDEFINED) RO_G(IJK) = RO_G0 
            IF (C_PG0 /= UNDEFINED) C_PG(IJK) = C_PG0 
            IF (MW_AVG /= UNDEFINED) MW_MIX_G(IJK) = MW_AVG 
! Strictly fluid cells: FLAG = 1
            IF(FLUID_AT(IJK)) THEN
               IF (MU_G0 /= UNDEFINED) THEN
                  MU_G(IJK) = MU_G0
                  MU_GT(IJK) = MU_G0
                  LAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
               ENDIF
               IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0 
               IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = DIF_G0 
            ENDIF
         ENDIF 

      ENDDO 


      DO M = 1, MMAX 
         DO IJK = ijkstart3, ijkend3
! All wall cells: FLAG >= 100
            IF (WALL_AT(IJK)) THEN 
               P_S(IJK,M) = ZERO 
               MU_S(IJK,M) = ZERO 
               LAMBDA_S(IJK,M) = ZERO 
               ALPHA_S(IJK,M) = ZERO 
               K_S(IJK,M) = ZERO 
               C_PS(IJK,M) = ZERO 
               D_p(IJK,M) = ZERO
               RO_S(IJK,M) = ZERO
            ELSE
! Fluid and inflow/outlfow cells: FLAG < 100
               IF (RO_S0(M) /= UNDEFINED) RO_S(IJK,M) = RO_S0(M)
               IF (C_PS0 /= UNDEFINED) C_PS(IJK,M) = C_PS0 
               IF (D_P0(M) /= UNDEFINED) D_P(IJK,M) = D_P0(M)
! Strictly fluid cells: FLAG = 1
               IF (FLUID_AT(IJK)) THEN
                  IF (MU_S0 /= UNDEFINED) THEN 
                     P_S(IJK,M) = ZERO 
                     MU_S(IJK,M) = MU_S0 
                     LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M) 
                     ALPHA_S(IJK,M) = ZERO 
                  ENDIF 
                  IF (K_S0 /= UNDEFINED) K_S(IJK,M) = K_S0 
                  IF (DIF_S0 /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = DIF_S0
               ENDIF
            ENDIF 

! set ep_star_array to user input ep_star in all cells.
            EP_star_array(ijk) = ep_star
            IF(EP_S_MAX(M) == UNDEFINED) EP_S_MAX(M) = ONE-EP_STAR
! this probably should not be used anymore            
            EP_S_CP = 1.D0 - EP_STAR

! initializing Sreekanth blending stress parameters (sof)
! changed blend_start to 0.99*ep_star from 0.97*ep_star [ceaf 2006-03-17]
! changed blend_end to 1.01*ep_star from 1.03*ep_star [ceaf 2006-03-17]
! added option for sigmoid function [sp 2006-10-24]
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

         ENDDO   ! end loop over ijk
      ENDDO   ! end loop over MMAX


      IF (RO_G0 == ZERO .AND. MMAX > 0) THEN
         IF(allocated(F_GS)) F_GS = ZERO
      ENDIF 
     

! Initializing the indexing system. This doesn't need to be done if no
! correlation is used to compute ep_star.
      IF(YU_STANDISH .OR. FEDORS_LANDEL .AND. .NOT. CALL_DQMOM) THEN

! refer to Syam's dissertation
         IF (SMAX == 2) THEN
            ep_s_max_ratio(1,2) = ep_s_max(1)/                         &
               (ep_s_max(1)+(1.-ep_s_max(1))*ep_s_max(2)) 
         ENDIF

         DO I = 1, MMAX
! Why is there a bailout condition for undefined D_P0?
            IF(D_P0(I) == UNDEFINED) RETURN
            DP_TMP(I) = D_P0(I)
            M_MAX(I) = I
         ENDDO

! Rearrange the indices from coarsest particles to finest to be 
! used in CALC_ep_star. Done here because it may need to be done
! for auto_restart
         DO I = 1, MMAX
            DO J = I, MMAX                  
               IF(DP_TMP(I) < DP_TMP(J)) THEN
                  old_value = DP_TMP(I)
                  DP_TMP(I) = DP_TMP(J)
                  DP_TMP(J) = old_value
               ENDIF                  
            ENDDO
         ENDDO
         DO I = 1, MMAX
            DO J = 1, MMAX
               IF(DP_TMP(I) == D_P0(J) .AND. D_P0(I) .NE. D_P0(J)) THEN
                  M_MAX(I) = J 
               ENDIF
            ENDDO
         ENDDO
      ENDIF    ! if Yu-standish or Fedors-Landel and .not. call_dqmom

      RETURN  
      END SUBROUTINE SET_CONSTPROP
