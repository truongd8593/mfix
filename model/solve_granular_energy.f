!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_GRANULAR_ENERGY(IER)                             C
!  Purpose: Solve granular energy equations                            C
!                                                                      C
!                                                                      C
!  Author: K. Agrawal                                 Date:            C
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
      SUBROUTINE SOLVE_GRANULAR_ENERGY(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run
      USE physprop
      USE visc_s
      USE geometry
      USE fldvar
      USE constant
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs 
      USE pgcor
      USE pscor
      USE leqsol 
      USE bc
      USE energy
      USE rxns
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3
      USE compar      
      USE mflux     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      phase index 
      INTEGER          M, L 
! 
! 
!                       Cp * Flux 
      DOUBLE PRECISION CpxFlux_E(DIMENSION_3), CpxFlux_N(DIMENSION_3), CpxFlux_T(DIMENSION_3) 
! 
! 
      DOUBLE PRECISION apo, sourcelhs, sourcerhs 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!
!                      particle mass and diameter
      DOUBLE PRECISION M_PM,D_PM 
!
!                      total solids volume fraction
      DOUBLE PRECISION TOT_EPS(DIMENSION_3)      
!      
!                      net production of all solids phase
      DOUBLE PRECISION TOT_SUM_RS(DIMENSION_3)
!
!                      old value of total solids number density
      DOUBLE PRECISION TOT_NO(DIMENSION_3)
!
!                      small value of theta_m, 1 cm2/s2 = 1e-4 m2/s2
      DOUBLE PRECISION smallTheta
!
!                      temporary variable for volume x interphase transfer 
!                      coefficient between solids phases
      DOUBLE PRECISION VxTC_ss(DIMENSION_3,DIMENSION_LM)
!-----------------------------------------------
      INCLUDE 'radtn1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'radtn2.inc'
      
      call lock_ambm
      call lock_tmp_array

      smallTheta = (to_SI)**4 * ZERO_EP_S

      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 

      IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN

          DO M = 1, MMAX 
               DO IJK = ijkstart3, ijkend3

                    D_PM = D_P(IJK,M)
                    M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)

! In Iddir & Arastoopour (2005) the granular temperature includes 
! mass of the particle in the definition.
!
                    IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                       CpxFlux_E(IJK) = (1.5D0/M_PM) * Flux_sE(IJK,M)
                       CpxFlux_N(IJK) = (1.5D0/M_PM) * Flux_sN(IJK,M)
                       CpxFlux_T(IJK) = (1.5D0/M_PM) * Flux_sT(IJK,M)
                    ELSE ! in case added mass is used.
                       CpxFlux_E(IJK) = (1.5D0/M_PM) * Flux_sSE(IJK)
                       CpxFlux_N(IJK) = (1.5D0/M_PM) * Flux_sSN(IJK)
                       CpxFlux_T(IJK) = (1.5D0/M_PM) * Flux_sST(IJK)
                    ENDIF

                    IF (FLUID_AT(IJK)) THEN 
                         CALL SOURCE_IA_NONEP_GRANULAR_ENERGY (SOURCELHS, SOURCERHS, IJK, M, IER) 
                         APO = (1.5D0/M_PM)*ROP_SO(IJK,M)*VOL(IJK)*ODT 
                         S_P(IJK) = APO + SOURCELHS + (1.5d0/M_PM)*ZMAX(SUM_R_S(IJK,M)) * VOL(IJK) 
                         S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
                              (1.50d0/M_PM)*THETA_M(IJK,M)*ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
                         EPS(IJK) = EP_S(IJK,M)
                    ELSE 
                         EPS(IJK) = ZERO 
                         S_P(IJK) = ZERO 
                         S_C(IJK) = ZERO 
                    ENDIF 
               END DO 
 
               CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8), U_S(1,M), &
                    V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M, IER)
           
               CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), BC_HW_THETA_M(1,M), &
                    BC_C_THETA_M(1,M), M, A_M, B_M, IER) 

! override bc settings if Johnson-Jackson bcs are specified
               CALL BC_THETA (M, A_M, B_M, IER)
 
               CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M, IER)
          ENDDO ! for M = 1, mmax

! use partial elimination on collisional dissipation term: SUM(Nip)
!     SUM( ED_s_ip* (Theta_p-Theta_i))
          IF (MMAX > 1) THEN 
              CALL CALC_VTC_SS (VXTC_SS, IER)   
              CALL PARTIAL_ELIM_IA (THETA_M, VXTC_SS, A_M, B_M, IER)
          ENDIF

! Adjusting the values of theta_m to zero when Ep_g < EP_star (Shaeffer, 1987)
! This is done here instead of calc_mu_s.f to avoid convergence problems. (sof)
          IF (SCHAEFFER) THEN
               DO M = 1, MMAX
                    DO IJK = ijkstart3, ijkend3
                         IF (FLUID_AT(IJK) .AND. EP_g(IJK) .LT. EP_g_blend_start(ijk)) THEN 

                              D_PM = D_P(IJK,M)
                              M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
                              A_M(IJK,1,M) = ZERO 
                              A_M(IJK,-1,M) = ZERO 
                              A_M(IJK,2,M) = ZERO 
                              A_M(IJK,-2,M) = ZERO 
                              A_M(IJK,3,M) = ZERO 
                              A_M(IJK,-3,M) = ZERO 
                              A_M(IJK,0,M) = -ONE
! In Iddir & Arastoopour (2005) the granular temperature includes
! mass of the particle in the definition.  Systems with small
! particles can give rise to smaller temperatures than the standard
! zero_ep_s		  
                              B_M(IJK,M) = -smallTheta*M_PM
                         ENDIF
                    ENDDO
               ENDDO ! for M
          ENDIF
! End of Shaeffer adjustments, sof.

          DO M = 1, MMAX 

               CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, &
                    NUM_RESID(RESID_TH,M), DEN_RESID(RESID_TH,M), RESID(RESID_TH,M),&
                    MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO, IER) 

               CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8), IER) 

!              call check_ab_m(a_m, b_m, m, .true., ier)
!              write(*,*)      resid(resid_th, m), max_resid(resid_th, m),&
!                              I_OF(ijk_resid(resid_th, m)), J_OF(ijk_resid(resid_th, m))
!              call write_ab_m(a_m, b_m, ijkmax2, m, ier)

!              call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!                   &    ier)

               CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8), LEQI, &
                    LEQM, IER) 

               CALL SOLVE_LIN_EQ ('Theta_m', 7, THETA_M(1,M), A_M, B_M, M, LEQI, LEQM, &
                    LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER) 
!              call out_array(Theta_m(1,m), 'Theta_m')

! Remove very small negative values of theta caused by leq solvers
               CALL ADJUST_THETA (M, IER) 
               IF (IER /= 0) RETURN                    !large negative granular temp -> divergence 

          ENDDO
 

      ELSEIF (TRIM(KT_TYPE) .EQ. 'GHD') THEN

          TOT_NO(:) = ZERO
          TOT_SUM_RS(:) = ZERO
          TOT_EPS(:) = ZERO
          CALL CALC_NFLUX (IER)  

          DO IJK = ijkstart3, ijkend3

              DO L = 1,SMAX
                  M_PM = (PI/6.d0)*(D_P(IJK,L)**3)*RO_S(L)
                  TOT_SUM_RS(IJK) = TOT_SUM_RS(IJK) + SUM_R_S(IJK,L)/M_PM
                  TOT_NO(IJK) = TOT_NO(IJK) + ROP_SO(IJK,L)/M_PM
                  IF (FLUID_AT(IJK)) TOT_EPS(IJK) = TOT_EPS(IJK) + EP_S(IJK,L)
              ENDDO
              M = MMAX 

! total number density is used for GHD theory
              CpxFlux_E(IJK) = 1.5D0 * Flux_nE(IJK)
              CpxFlux_N(IJK) = 1.5D0 * Flux_nN(IJK)
              CpxFlux_T(IJK) = 1.5D0 * Flux_nT(IJK)

              IF (FLUID_AT(IJK)) THEN
                  CALL SOURCE_GHD_GRANULAR_ENERGY (SOURCELHS, SOURCERHS, IJK, IER) 
                  APO = 1.5D0 * TOT_NO(IJK)*VOL(IJK)*ODT 
                  S_P(IJK) = APO + SOURCELHS + 1.5d0 *ZMAX(TOT_SUM_RS(IJK)) * VOL(IJK) 
                  S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
                              1.50d0 * THETA_M(IJK,M)*ZMAX(-TOT_SUM_RS(IJK)) * VOL(IJK)
                  EPS(IJK) = TOT_EPS(IJK)
              ELSE 
                  EPS(IJK) = ZERO 
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
              ENDIF 
          ENDDO 
 
          CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8), U_S(1,M), &
                    V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M, IER)
           
          CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), BC_HW_THETA_M(1,M), &
                    BC_C_THETA_M(1,M), M, A_M, B_M, IER) 

! override bc settings if Johnson-Jackson bcs are specified
          CALL BC_THETA (M, A_M, B_M, IER)
 
          CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M, IER)

          CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, &
                    NUM_RESID(RESID_TH,M), DEN_RESID(RESID_TH,M), RESID(RESID_TH,M),&
                    MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO, IER) 

          CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8), IER)

!              call check_ab_m(a_m, b_m, m, .true., ier)
!              write(*,*)      resid(resid_th, m), max_resid(resid_th, m),&
!                              I_OF(ijk_resid(resid_th, m)), J_OF(ijk_resid(resid_th, m))
!              call write_ab_m(a_m, b_m, ijkmax2, m, ier)

!              call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!                   &    ier)

          CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8), LEQI, &
                    LEQM, IER) 

          CALL SOLVE_LIN_EQ ('Theta_m', 7, THETA_M(1,M), A_M, B_M, M, LEQI, LEQM, &
                    LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER) 
!              call out_array(Theta_m(1,m), 'Theta_m')

! Remove very small negative values of theta caused by leq solvers
          CALL ADJUST_THETA (M, IER) 
          IF (IER /= 0) RETURN                    !large negative granular temp -> divergence

      ELSE     ! default KT in MFIX or GD_99 theory
        DO M = 1, MMAX 


         DO IJK = ijkstart3, ijkend3

            IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
               CpxFlux_E(IJK) = 1.5D0 * Flux_sE(IJK,M)
               CpxFlux_N(IJK) = 1.5D0 * Flux_sN(IJK,M)
               CpxFlux_T(IJK) = 1.5D0 * Flux_sT(IJK,M)
            ELSE
               CpxFlux_E(IJK) = 1.5D0 * Flux_sSE(IJK)
               CpxFlux_N(IJK) = 1.5D0 * Flux_sSN(IJK)
               CpxFlux_T(IJK) = 1.5D0 * Flux_sST(IJK)
            ENDIF
 
            IF (FLUID_AT(IJK)) THEN 

                 IF (TRIM(KT_TYPE) .EQ. 'GD_99') THEN
                      CALL SOURCE_GD_99_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER)
                 ELSE
                      CALL SOURCE_GRANULAR_ENERGY (SOURCELHS, SOURCERHS, IJK, M, IER) 
                 ENDIF

                 APO = 1.5D0*ROP_SO(IJK,M)*VOL(IJK)*ODT 
                 S_P(IJK) = APO + SOURCELHS + 1.5D0 * ZMAX(SUM_R_S(IJK,M)) * VOL(IJK) 
                 S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
                      1.5D0 * THETA_M(IJK,M)*ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
                 EPS(IJK) = EP_S(IJK,M) 

            ELSE 

               EPS(IJK) = ZERO 
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 

            ENDIF 
         END DO 

         CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8), U_S(1,M), &
            V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M, IER) 

         CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), BC_HW_THETA_M(1,M), &
            BC_C_THETA_M(1,M), M, A_M, B_M, IER) 

         CALL BC_THETA (M, A_M, B_M, IER)     !override bc settings if 
                                              !Johnson-Jackson bcs are specified

         CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M, IER)


! Adjusting the values of theta_m to zero when Ep_g < EP_star (Shaeffer, 1987)
! This is done here instead of calc_mu_s.f to avoid convergence problems. (sof)

         IF (SCHAEFFER) THEN
           DO IJK = ijkstart3, ijkend3

              IF (FLUID_AT(IJK) .AND. EP_g(IJK) .LT. EP_g_blend_start(ijk)) THEN 

                 A_M(IJK,1,M) = ZERO 
                 A_M(IJK,-1,M) = ZERO 
                 A_M(IJK,2,M) = ZERO 
                 A_M(IJK,-2,M) = ZERO 
                 A_M(IJK,3,M) = ZERO 
                 A_M(IJK,-3,M) = ZERO 
                 A_M(IJK,0,M) = -ONE   
                 B_M(IJK,M) = -smallTheta
              ENDIF
           END DO
         ENDIF
! End of Shaeffer adjustments, sof.

         CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, NUM_RESID(RESID_TH,M), &
            DEN_RESID(RESID_TH,M), RESID(RESID_TH,M), &
            MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO, IER) 

         CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8), IER) 

!        call check_ab_m(a_m, b_m, m, .true., ier)
!          write(*,*)
!     &      resid(resid_th, m), max_resid(resid_th, m),
!     &      ijk_resid(resid_th, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)

!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)

         CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8), LEQI, &
            LEQM, IER) 

         CALL SOLVE_LIN_EQ ('Theta_m', 7, THETA_M(1,M), A_M, B_M, M, LEQI, LEQM, &
                             LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER) 
!          call out_array(Theta_m(1,m), 'Theta_m')

! Remove very small negative values of theta caused by leq solvers
         CALL ADJUST_THETA (M, IER) 
         IF (IER /= 0) RETURN                    !large negative granular temp -> divergence 
        END DO 
      ENDIF     ! for kt_type
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_GRANULAR_ENERGY

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
