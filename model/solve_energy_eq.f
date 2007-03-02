!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_ENERGY_EQ(IER)                                   C
!  Purpose: Solve energy equations                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To eliminate kinetic solids calculations when doing DES    C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
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
      SUBROUTINE SOLVE_ENERGY_EQ(IER) 
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
      Use tmp_array, S_p => ARRAY1, S_C => ARRAY2, EPs => ARRAY3
      Use tmp_array1, VxGama => ARRAYm1
      USE compar   
      USE discretelement 
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
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      phase index 
      INTEGER          m 
      INTEGER          TEMP_MMAX
! 
!                       Cp * Flux 
      DOUBLE PRECISION CpxFlux_E(DIMENSION_3), CpxFlux_N(DIMENSION_3), CpxFlux_T(DIMENSION_3) 
! 
! 
      DOUBLE PRECISION apo 
! 
!                      Indices 
      INTEGER          IJK, I, J, K 
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
      call lock_tmp_array1

      TEMP_MMAX = MMAX
      IF(DISCRETE_ELEMENT) THEN
         MMAX = 0   ! Only the gas calculations are needed
      END IF            
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 

      DO IJK = IJKSTART3, IJKEND3
!
         I = I_OF(IJK)
	 J = J_OF(IJK)
	 K = K_OF(IJK) 

         IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
	   CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gE(IJK)
	 ENDIF
	 
         IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
	   CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gN(IJK)
	 ENDIF
	 
         IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
	   CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gT(IJK)
	 ENDIF

         IF (FLUID_AT(IJK)) THEN 
            APO = ROP_GO(IJK)*C_PG(IJK)*VOL(IJK)*ODT 
            S_P(IJK) = APO + S_RPG(IJK)*VOL(IJK) 
            S_C(IJK)=APO*T_GO(IJK)-HOR_G(IJK)*VOL(IJK)+S_RCG(IJK)*VOL(IJK) 
         ELSE 
!
            S_P(IJK) = ZERO 
            S_C(IJK) = ZERO 
!
         ENDIF 
      END DO 
      CALL CONV_DIF_PHI (T_g, K_G, DISCRETIZE(6), U_G, V_G, W_G, CpxFlux_E, CpxFlux_N, CpxFlux_T, 0, A_M&
         , B_M, IER) 
!
      CALL BC_PHI (T_g, BC_T_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, 0, A_M, B_M, IER) 
!
      CALL SOURCE_PHI (S_P, S_C, EP_G, T_G, 0, A_M, B_M, IER) 
!
      DO M = 1, MMAX 
!
         DO IJK = IJKSTART3, IJKEND3
!
           I = I_OF(IJK)
	   J = J_OF(IJK)
	   K = K_OF(IJK) 

           IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
             CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sE(IJK,M)
	   ENDIF
	 
           IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
	     CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sN(IJK,M)
	   ENDIF
	 
           IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
             CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sT(IJK,M)
	   ENDIF
	   
            IF (FLUID_AT(IJK)) THEN 
               APO = ROP_SO(IJK,M)*C_PS(IJK,M)*VOL(IJK)*ODT 
               S_P(IJK) = APO + S_RPS(IJK,M)*VOL(IJK) 
               S_C(IJK) = APO*T_SO(IJK,M) - HOR_S(IJK,M)*VOL(IJK) + S_RCS(IJK,M&
                  )*VOL(IJK) 
!
               VXGAMA(IJK,M) = GAMA_GS(IJK,M)*VOL(IJK) 
               EPS(IJK) = EP_S(IJK,M) 
!
            ELSE 
!
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 
               VXGAMA(IJK,M) = ZERO 
               EPS(IJK) = ZERO 
!
            ENDIF 
         END DO 
         CALL CONV_DIF_PHI (T_s(1,M), K_S(1,M), DISCRETIZE(6), U_S(1,M), V_S(1,&
            M), W_S(1,M), CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M, IER) 
!
         CALL BC_PHI (T_s(1,M), BC_T_S(1,M), BC_TW_S(1,M), BC_HW_T_S(1,M), BC_C_T_S(1,M)&
            , M, A_M, B_M, IER) 
!
         CALL SOURCE_PHI (S_P, S_C, EPS, T_S(1,M), M, A_M, B_M, IER) 
      END DO 
      IF (MMAX > 0) CALL PARTIAL_ELIM_S (T_G, T_S, VXGAMA, A_M, B_M, IER) 
!
      CALL CALC_RESID_S (T_G, A_M, B_M, 0, RESID(RESID_T,0), MAX_RESID(RESID_T,&
         0), IJK_RESID(RESID_T,0), ZERO, IER) 
!
      CALL UNDER_RELAX_S (T_G, A_M, B_M, 0, UR_FAC(6), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_t, 0), max_resid(resid_t, 0),
!     &    ijk_resid(resid_t, 0)
!
!
      DO M = 1, MMAX 
!
         CALL CALC_RESID_S (T_S(1,M), A_M, B_M, M, RESID(RESID_T,M), MAX_RESID(&
            RESID_T,M), IJK_RESID(RESID_T,M), ZERO, IER) 
!
         CALL UNDER_RELAX_S (T_S(1,M), A_M, B_M, M, UR_FAC(6), IER) 
      END DO 
      CALL ADJUST_LEQ(RESID(RESID_T,0),LEQ_IT(6),LEQ_METHOD(6),LEQI,LEQM,IER) 
!         call test_lin_eq(a_m(1, -3, 0), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)
!
      CALL SOLVE_LIN_EQ ('T_g', 6, T_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)  

!       bound the temperature
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.WALL_AT(IJK))&
              T_g(IJK) = MIN(TMAX, MAX(TMIN, T_g(IJK)))
	 ENDDO

!        call out_array(T_g, 'T_g')
!
      DO M = 1, MMAX 
!
         CALL ADJUST_LEQ (RESID(RESID_T,M), LEQ_IT(6), LEQ_METHOD(6), LEQI, &
            LEQM, IER) 
!         call test_lin_eq(a_m(1, -3, M), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)
!
         CALL SOLVE_LIN_EQ ('T_s', 6, T_S(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER) 

!       bound the temperature
        DO IJK = IJKSTART3, IJKEND3
          IF(.NOT.WALL_AT(IJK))&
            T_s(IJK, M) = MIN(TMAX, MAX(TMIN, T_s(IJK, M))) 
        ENDDO

      END DO 
      
      call unlock_ambm
      call unlock_tmp_array
      call unlock_tmp_array1
      
      MMAX = TEMP_MMAX
      
      RETURN  
      END SUBROUTINE SOLVE_ENERGY_EQ 


!// Comments on the modifications for DMP version implementation      
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
