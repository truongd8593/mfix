!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_constprop                                           C
!  Purpose: This routine sets various constant physical properties     C
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
      USE param1, only: zero, half, one, undefined
      USE fldvar
      USE visc_s
      USE visc_g
      USE energy
      USE geometry
      USE indices
      USE physprop
      USE constant, only: ep_s_max_ratio, d_p_ratio, ep_s_max, m_max
      use constant, only: ep_star, l_scale0
      USE run
      USE drag, only: f_gs, f_ss
      USE compar 
      use kintheory
      use mms, only: use_mms

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices      
      INTEGER :: IJK, M, N, I, J
      DOUBLE PRECISION :: old_value, DP_TMP(SMAX)
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

! Variables specific to various kinetic theory models
      IF (KT_TYPE_ENUM == IA_2005) THEN
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
      IF (KT_TYPE_ENUM == IA_2005 .OR. KT_TYPE_ENUM == GD_1999 .OR.  &
          KT_TYPE_ENUM == GTSH_2012) THEN
         EDT_s_ip = ZERO
         EDvel_sM_ip = ZERO
      ENDIF
      IF(KT_TYPE_ENUM == GTSH_2012) THEN
         A2_gtsh = ZERO
         xsi_gtsh = zero
      ENDIF
     
! Set specified constant physical properties values
      DO IJK = ijkstart3, ijkend3

! All wall cells: FLAG >= 100
         IF (WALL_AT(IJK) .AND. .NOT.USE_MMS) THEN
            RO_G(IJK) = ZERO 
            MU_G(IJK) = ZERO 
            K_G(IJK) = ZERO 
            C_PG(IJK) = ZERO 
            MW_MIX_G(IJK) = ZERO 
         ELSE
! Fluid and inflow/outflow cells: FLAG < 100
            IF (RO_G0 /= UNDEFINED) RO_G(IJK) = RO_G0 
            IF (C_PG0 /= UNDEFINED) C_PG(IJK) = C_PG0 
            IF (MW_AVG /= UNDEFINED) MW_MIX_G(IJK) = MW_AVG 

! Strictly fluid cells: FLAG = 1
            IF(FLUID_AT(IJK) .OR. USE_MMS) THEN
               IF (MU_G0 /= UNDEFINED) THEN
                  MU_G(IJK) = MU_G0
                  MU_GT(IJK) = MU_G0
                  LAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
               ENDIF
               IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0 
               IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = DIF_G0
            ELSE
! ONLY inflow/outflow cells: FLAG .NE. 1 and FLAG < 100
! initialize transport coefficients to zero in inflow/outflow cells
               IF (MU_G0 /= UNDEFINED) THEN
                  MU_G(IJK) = ZERO
                  MU_GT(IJK) = ZERO
                  LAMBDA_GT(IJK) = ZERO
               ENDIF
               IF (K_G0 /= UNDEFINED) K_G(IJK) = ZERO
               IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = ZERO
            ENDIF
         ENDIF 

      ENDDO 


      DO M = 1, MMAX 
         DO IJK = ijkstart3, ijkend3
! All wall cells: FLAG >= 100
            IF (WALL_AT(IJK) .AND. .NOT.USE_MMS) THEN 
               P_S(IJK,M) = ZERO 
               MU_S(IJK,M) = ZERO 
               LAMBDA_S(IJK,M) = ZERO 
               ALPHA_S(IJK,M) = ZERO 
               K_S(IJK,M) = ZERO 
               C_PS(IJK,M) = ZERO 
               D_p(IJK,M) = ZERO
               RO_S(IJK,M) = ZERO
            ELSE
! Fluid and inflow/outflow cells: FLAG < 100
               IF (RO_S0(M) /= UNDEFINED) RO_S(IJK,M) = RO_S0(M)
               IF (C_PS0(M) /= UNDEFINED) C_PS(IJK,M) = C_PS0(M)
               IF (D_P0(M) /= UNDEFINED) D_P(IJK,M) = D_P0(M)

! Strictly fluid cells: FLAG = 1
               IF(FLUID_AT(IJK) .OR. USE_MMS) THEN
                  IF (MU_S0 /= UNDEFINED) THEN 
                     P_S(IJK,M) = ZERO 
                     MU_S(IJK,M) = MU_S0 
                     LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M) 
                     ALPHA_S(IJK,M) = ZERO 
                  ENDIF 
                  IF (K_S0(M) /= UNDEFINED) K_S(IJK,M) = K_S0(M)
                  IF (DIF_S0 /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = DIF_S0
               ELSE
! ONLY inflow/outflow cells: FLAG .NE. 1 and FLAG < 100
! initialize transport coefficients to zero in inflow/outflow cells
                  IF (MU_S0 /= UNDEFINED) THEN 
                     P_S(IJK,M) = ZERO 
                     MU_S(IJK,M) = ZERO
                     LAMBDA_S(IJK,M) = ZERO
                     ALPHA_S(IJK,M) = ZERO
                  ENDIF
                  IF (K_S0(M) /= UNDEFINED) K_S(IJK,M) = ZERO
                  IF (DIF_S0 /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = ZERO 
               ENDIF
            ENDIF 


! set ep_star_array to user input ep_star in all cells.
            EP_star_array(ijk) = ep_star
! initializing blending stress parameters 
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
     

! Initializing parameters needed if a correlation is used to compute
! ep_star: initializing the indexing system.
      IF(YU_STANDISH .OR. FEDORS_LANDEL) THEN
         DO M = 1, SMAX
            IF(EP_S_MAX(M) == UNDEFINED) EP_S_MAX(M) = ONE-EP_STAR
         ENDDO 

         IF (.NOT.CALL_DQMOM) THEN

! refer to Syam's dissertation
            IF (SMAX == 2) THEN
               ep_s_max_ratio(1,2) = ep_s_max(1)/ &
                  (ep_s_max(1)+(1.-ep_s_max(1))*ep_s_max(2)) 
            ENDIF

! initialize local variables            
            DO I = 1, SMAX
               DP_TMP(I) = D_P0(I)
               M_MAX(I) = I
            ENDDO
            
! Rearrange the indices from coarsest particles to finest to be 
! used in CALC_ep_star. Done here because it may need to be done
! for auto_restart
            DO I = 1, SMAX
               DO J = I, SMAX                  
                  IF(DP_TMP(I) < DP_TMP(J)) THEN
                     old_value = DP_TMP(I)
                     DP_TMP(I) = DP_TMP(J)
                     DP_TMP(J) = old_value
                  ENDIF                  
               ENDDO
            ENDDO

            DO I = 1, SMAX
               DO J = 1, SMAX
                  IF(DP_TMP(I) == D_P0(J) .AND. D_P0(I) .NE. D_P0(J)) THEN
                     M_MAX(I) = J 
                  ENDIF
               ENDDO
            ENDDO
         ENDIF    ! if .not. call_dqmom
      ELSE   ! if .not. Yu-standish or Fedors-Landel
         EP_S_MAX(:) = ZERO
         EP_S_MAX_RATIO(:,:) = ZERO
         D_P_RATIO(:,:) = ZERO
         M_MAX(:) = ZERO
      ENDIF

      RETURN  
      END SUBROUTINE SET_CONSTPROP
