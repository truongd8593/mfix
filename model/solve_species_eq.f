!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_SPECIES_EQ(IER)                                  C
!  Purpose: Solve species mass balance equations                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 11-FEB-98  C
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
      SUBROUTINE SOLVE_SPECIES_EQ(IER) 
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
      USE bc
      USE energy
      USE rxns
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
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
!                      species index 
      INTEGER          n 
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
!                      Volume x average gama at cell centers 
!      DOUBLE PRECISION VxGama(DIMENSION_3, DIMENSION_M) 
! 
! 
      DOUBLE PRECISION apo 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      
      call lock_ambm
      call lock_tmp_array

!//SP
      write(*,*) 'entered solve_species',myPE
 
!
!     Fluid phase species mass balance equations
!
      IF (SPECIES_EQ(0)) THEN 
         DO N = 1, NMAX(0) 
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!// 350 1206 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3	    
            DO IJK = ijkstart3, ijkend3 
!
               IF (FLUID_AT(IJK)) THEN 
                  APO = ROP_G(IJK)*VOL(IJK)*ODT 
                  S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK))+ROX_GC(IJK,N))*VOL(IJK) 
                  S_C(IJK) = APO*X_GO(IJK,N) + X_G(IJK,N)*ZMAX((-SUM_R_G(IJK)))&
                     *VOL(IJK) + R_GP(IJK,N)*VOL(IJK) 
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 
            END DO 
            CALL CONV_DIF_PHI (X_G(1,N), DIF_G(1,N), DISCRETIZE(7), U_G, V_G, &
               W_G, ROP_G, 0, A_M, B_M, IER) 
!
!//SP
	    write(*,*) 'aft CONV_DIF_PHI in solve_species', myPE
!
            CALL BC_PHI (BC_X_G(1,N), BC_XW_G(1,N), BC_HW_X_G(1,N), BC_C_X_G(1,&
               N), 0, A_M, B_M, IER) 
!
!//SP
	    write(*,*) 'aft BC_PHI in solve_species', myPE
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, X_G(1,N), 0, A_M, B_M, IER) 
!
!//SP
	    write(*,*) 'aft source_phi in solve_species', myPE
!//SP
!           call mfix_exit(myPE)   !//SP

!
            CALL CALC_RESID_S (X_G(1,N), A_M, B_M, 0, RESID(RESID_X+(N-1),0), &
               MAX_RESID(RESID_X+(N-1),0), IJK_RESID(RESID_X+(N-1),0), &
	       ZERO_X_GS, IER) 
!//SP
	    write(*,*) 'aft CALC_RESID_S in solve_species', myPE
!
            CALL UNDER_RELAX_S (X_G(1,N), A_M, B_M, 0, UR_FAC(7), IER) 
!//SP
            write(*,*) 'aft UNDER_RELAX_S in solve_species', myPE
!
!          call check_ab_m(a_m, b_m, 0, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!          write(*,*)
!     &    resid(resid_x+(N-1), 0), max_resid(resid_x+(N-1), 0),
!     &    ijk_resid(resid_x+(N-1), 0)
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
            CALL ADJUST_LEQ (RESID(RESID_X+(N-1),0), LEQ_IT(7), LEQ_METHOD(7), &
               LEQI, LEQM, IER) 
!//SP
            write(*,*) 'aft ADJUST_LEQ in solve_species', myPE
!
            CALL SOLVE_LIN_EQ ('X_g', X_G(1,N), A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 
!//SP
            write(*,*) 'aft SOLVE_LIN_EQ in solve_species', myPE
!	    call mfix_exit
!          call out_array(X_g(1, N), 'X_g')
            CALL BOUND_X (X_G(1,N), IJKMAX2, IER) 
!//SP
            write(*,*) 'aft BOUND_X in solve_species', myPE
!
         END DO 
      ENDIF 
!
!
!     Granular phase species balance equations
!
      DO M = 1, MMAX 
!
         IF (SPECIES_EQ(M)) THEN 
!
!
            DO N = 1, NMAX(M) 
!
               CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
!
!// 350 1025 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
               DO IJK = ijkstart3, ijkend3 
!
                  IF (FLUID_AT(IJK)) THEN 
!
                     APO = ROP_S(IJK,M)*VOL(IJK)*ODT 
                     S_P(IJK) = APO + (ZMAX(SUM_R_S(IJK,M))+ROX_SC(IJK,M,N))*&
                        VOL(IJK) 
                     S_C(IJK) = APO*X_SO(IJK,M,N) + X_S(IJK,M,N)*ZMAX((-SUM_R_S&
                        (IJK,M)))*VOL(IJK) + R_SP(IJK,M,N)*VOL(IJK) 
!
                     EPS(IJK) = EP_S(IJK,M) 
!
                  ELSE 
!
                     S_P(IJK) = ZERO 
                     S_C(IJK) = ZERO 
                     EPS(IJK) = ZERO 
!
                  ENDIF 
               END DO 
               CALL CONV_DIF_PHI (X_S(1,M,N), DIF_S(1,M,N), DISCRETIZE(7), U_S(&
                  1,M), V_S(1,M), W_S(1,M), ROP_S(1,M), M, A_M, B_M, IER) 
!
               CALL BC_PHI (BC_X_S(1,M,N), BC_XW_S(1,M,N), BC_HW_X_S(1,M,N), &
                  BC_C_X_S(1,M,N), M, A_M, B_M, IER) 
!
!
               CALL SOURCE_PHI (S_P, S_C, EPS, X_S(1,M,N), M, A_M, B_M, IER) 
!
               CALL CALC_RESID_S (X_S(1,M,N), A_M, B_M, M, RESID(RESID_X+(N-1),&
                  M), MAX_RESID(RESID_X+(N-1),M), IJK_RESID(RESID_X+(N-1),M), &
                  ZERO_X_GS, IER) 
!
!
               CALL UNDER_RELAX_S (X_S(1,M,N), A_M, B_M, M, UR_FAC(7), IER) 
!           call check_ab_m(a_m, b_m, m, .false., ier)
!           write(*,*)
!     &      resid(resid_x+(N-1), m), max_resid(resid_x+(N-1), m),
!     &      ijk_resid(resid_x+(N-1), m)
!           call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!           call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &     ier)
!
               CALL ADJUST_LEQ (RESID(RESID_X+(N-1),M), LEQ_IT(7), LEQ_METHOD(7&
                  ), LEQI, LEQM, IER) 
!
               CALL SOLVE_LIN_EQ ('X_s', X_S(1,M,N), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 
!            call out_array(X_s(1,m,n), 'X_s')
               CALL BOUND_X (X_S(1,M,N), IJKMAX2, IER) 
!
            END DO 
         ENDIF 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_SPECIES_EQ 
