!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
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
      SUBROUTINE SOLVE_VEL_STAR(IER) 
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
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs 
      USE pgcor
      USE pscor
      USE leqsol
      Use ambm 
      Use tmp_array1, VxF_gs => Arraym1
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
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
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
!                      Volume x average at momentum cell centers 
!      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M) 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
!
      call lock_ambm
      call lock_tmp_array1


      IF (MMAX == 0) CALL ZERO_ARRAY (VXF_GS(1,1), IJKMAX2, IER) 
!
!  2.1 Calculate U_m_star and residuals
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      CALL CONV_DIF_U_G (A_M, B_M, IER) 
      CALL CONV_DIF_U_S (A_M, B_M, IER) 
!
      CALL SOURCE_U_G (A_M, B_M, IER) 
      CALL SOURCE_U_S (A_M, B_M, IER) 
!
      IF (MMAX > 0) CALL VF_GS_X (F_GS, VXF_GS, IER) 
!
      CALL CALC_D_E (A_M, VXF_GS, D_E, IER) 
      IF (MMAX > 0) CALL CALC_E_E (A_M, MCP, E_E, IER) 
!
      IF (MMAX > 0) CALL PARTIAL_ELIM_U (U_G, U_S, VXF_GS, A_M, B_M, IER) 
!
      CALL ADJUST_A_U_G (A_M, B_M, IER) 
      CALL ADJUST_A_U_S (A_M, B_M, IER) 
!
      IF (MOMENTUM_X_EQ(0)) THEN 
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_U,0), &
            MAX_RESID(RESID_U,0), IJK_RESID(RESID_U,0), IER) 
         CALL UNDER_RELAX_U (U_G, A_M, B_M, 0, UR_FAC(3), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_u, 0), max_resid(resid_u, 0),
!     &    ijk_resid(resid_u, 0)
!
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
            CALL CALC_RESID_U (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_U,M), MAX_RESID(RESID_U,M), IJK_RESID(RESID_U,M), IER) 
            CALL UNDER_RELAX_U (U_S(1,M), A_M, B_M, M, UR_FAC(3), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_u, m), max_resid(resid_u, m),
!     &      ijk_resid(resid_u, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
      END DO 
      IF (MOMENTUM_X_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_U,0), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('U_g', U_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(3), LEQ_TOL(3),IER) 
!        call out_array(u_g, 'u_g')
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
            CALL ADJUST_LEQ (RESID(RESID_U,M), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('U_s', U_S(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(3), LEQ_TOL(3),IER) 
!          call out_array(u_s(1,m), 'u_s')
         ENDIF 
      END DO 
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      CALL CONV_DIF_V_G (A_M, B_M, IER) 
      CALL CONV_DIF_V_S (A_M, B_M, IER) 
!
      CALL SOURCE_V_G (A_M, B_M, IER) 
      CALL SOURCE_V_S (A_M, B_M, IER) 
!
      IF (MMAX > 0) CALL VF_GS_Y (F_GS, VXF_GS, IER) 
!
      CALL CALC_D_N (A_M, VXF_GS, D_N, IER) 
      IF (MMAX > 0) CALL CALC_E_N (A_M, MCP, E_N, IER) 
!
      IF (MMAX > 0) CALL PARTIAL_ELIM_V (V_G, V_S, VXF_GS, A_M, B_M, IER) 
!
      CALL ADJUST_A_V_G (A_M, B_M, IER) 
      CALL ADJUST_A_V_S (A_M, B_M, IER) 
!
      IF (MOMENTUM_Y_EQ(0)) THEN 
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_V,0), &
            MAX_RESID(RESID_V,0), IJK_RESID(RESID_V,0), IER) 
         CALL UNDER_RELAX_V (V_G, A_M, B_M, 0, UR_FAC(4), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_v, 0), max_resid(resid_v, 0),
!     &    ijk_resid(resid_v, 0)
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Y_EQ(M)) THEN 
            CALL CALC_RESID_V (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_V,M), MAX_RESID(RESID_V,M), IJK_RESID(RESID_V,M), IER) 
            CALL UNDER_RELAX_V (V_S(1,M), A_M, B_M, M, UR_FAC(4), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_v, m), max_resid(resid_v, m),
!     &      ijk_resid(resid_v, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
      END DO 
      IF (MOMENTUM_Y_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_V,0), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('V_g', V_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(4), LEQ_TOL(4),IER) 
!        call out_array(v_g, 'v_g')
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Y_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
            CALL ADJUST_LEQ (RESID(RESID_V,M), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('V_s', V_S(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(4), LEQ_TOL(4),IER) 
!          call out_array(v_s(1,m), 'v_s')
         ENDIF 
      END DO 
      IF (NO_K)THEN
        call unlock_ambm
        call unlock_tmp_array1
        RETURN  
      ENDIF
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      CALL CONV_DIF_W_G (A_M, B_M, IER) 
      CALL CONV_DIF_W_S (A_M, B_M, IER) 
!
      CALL SOURCE_W_G (A_M, B_M, IER) 
      CALL SOURCE_W_S (A_M, B_M, IER) 
!
      IF (MMAX > 0) CALL VF_GS_Z (F_GS, VXF_GS, IER) 
!
      CALL CALC_D_T (A_M, VXF_GS, D_T, IER) 
      IF (MMAX > 0) CALL CALC_E_T (A_M, MCP, E_T, IER) 
!
      IF (MMAX > 0) CALL PARTIAL_ELIM_W (W_G, W_S, VXF_GS, A_M, B_M, IER) 
!
      CALL ADJUST_A_W_G (A_M, B_M, IER) 
      CALL ADJUST_A_W_S (A_M, B_M, IER) 
!
      IF (MOMENTUM_Z_EQ(0)) THEN 
         CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_W,0), &
            MAX_RESID(RESID_W,0), IJK_RESID(RESID_W,0), IER) 
         CALL UNDER_RELAX_W (W_G, A_M, B_M, 0, UR_FAC(5), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &      resid(resid_w, 0), max_resid(resid_w, 0),
!     &      ijk_resid(resid_w, 0)
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Z_EQ(M)) THEN 
            CALL CALC_RESID_W (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_W,M), MAX_RESID(RESID_W,M), IJK_RESID(RESID_W,M), IER) 
            CALL UNDER_RELAX_W (W_S(1,M), A_M, B_M, M, UR_FAC(5), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_w, m), max_resid(resid_w, m),
!     &      ijk_resid(resid_w, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
      END DO 
      IF (MOMENTUM_Z_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_W,0), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('W_g', W_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(5), LEQ_TOL(5),IER) 
!        call out_array(w_g, 'w_g')
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Z_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
            CALL ADJUST_LEQ (RESID(RESID_W,M), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('W_s', W_S(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(5), LEQ_TOL(5),IER) 
!          call out_array(w_s(1,m), 'w_s')
         ENDIF 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array1

      RETURN  
      END SUBROUTINE SOLVE_VEL_STAR 
