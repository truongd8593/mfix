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
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, TxCp => Array4
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
      INTEGER          m 
! 
!                      Septadiagonal matrix A_m 
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source term on LHS.  Must be positive. 
!      DOUBLE PRECISION S_p(DIMENSION_3) 
! 
!                      Source term on RHS 
!      DOUBLE PRECISION S_C(DIMENSION_3) 
! 
!                      Solids volume fraction 
!      DOUBLE PRECISION EPs(DIMENSION_3) 
! 
!                      ROP * Cp 
!      DOUBLE PRECISION ROPxCp(DIMENSION_3) 
! 
! 
      DOUBLE PRECISION apo, sourcelhs, sourcerhs 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
      INCLUDE 'radtn1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'radtn2.inc'
      
      call lock_ambm
      call lock_tmp_array

!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      DO M = 1, MMAX 
!
!
         DO IJK = ijkstart3, ijkend3
!
            IF (FLUID_AT(IJK)) THEN 
!
               CALL SOURCE_GRANULAR_ENERGY (SOURCELHS, SOURCERHS, IJK, M, IER) 
               TXCP(IJK) = 1.5D0*THETA_M(IJK,M) 
               APO = 1.5D0*ROP_SO(IJK,M)*VOL(IJK)*ODT 
               S_P(IJK) = APO + SOURCELHS + ZMAX(SUM_R_S(IJK,M)) * VOL(IJK) 
               S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
	                  THETA_M(IJK,M)*ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
               EPS(IJK) = EP_S(IJK,M) 
!
            ELSE 
!
               EPS(IJK) = ZERO 
               TXCP(IJK) = ZERO 
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 
!
            ENDIF 
         END DO 
         CALL CONV_DIF_PHI (TXCP, KTH_S(1,M), DISCRETIZE(8), U_S(1,M), &
            V_S(1,M), W_S(1,M), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M, IER) 
!
         CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), BC_HW_THETA_M(1,M), &
            BC_C_THETA_M(1,M), M, A_M, B_M, IER) 
!
         CALL BC_THETA (M, A_M, B_M, IER)        !override bc settings if 
!                                        Johnson-Jackson bcs are specified
!
         CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M, IER) 
!
!
! Adjusting the values of theta_m to zero when Ep_g < EP_star (Shaeffer, 1987)
! This is done here instead of calc_mu_s.f to avoid convergence problems. (sof)
!
         IF (SCHAEFFER) THEN
           DO IJK = ijkstart3, ijkend3
!
              IF (FLUID_AT(IJK) .AND. EP_g(IJK) .LT. EP_star) THEN 
!

                 A_M(IJK,1,M) = ZERO 
                 A_M(IJK,-1,M) = ZERO 
                 A_M(IJK,2,M) = ZERO 
                 A_M(IJK,-2,M) = ZERO 
                 A_M(IJK,3,M) = ZERO 
                 A_M(IJK,-3,M) = ZERO 
                 A_M(IJK,0,M) = -ONE 		  
                 B_M(IJK,M) = ZERO
	      ENDIF
	   END DO
	 ENDIF	 
! End of Shaeffer adjustments, sof.
!
         CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, RESID(RESID_TH,M), &
            MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO, IER) 
!
         CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8), IER) 
!
!        call check_ab_m(a_m, b_m, m, .true., ier)
!          write(*,*)
!     &      resid(resid_th, m), max_resid(resid_th, m),
!     &      ijk_resid(resid_th, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('Theta_m', THETA_M(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(8), LEQ_TOL(8),IER) 
!          call out_array(Theta_m(1,m), 'Theta_m')
!
!
!        Remove very small negative values of theta caused by leq solvers
         CALL ADJUST_THETA (M, IER) 
         IF (IER /= 0) RETURN                    !large negative granular temp -> divergence 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_GRANULAR_ENERGY 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
