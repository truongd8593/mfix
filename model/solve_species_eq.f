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
      USE matrix 
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
      USE compar       
      USE mpi_utility      
      USE sendrecv      
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
      INTEGER          ln 
! 
! 
      DOUBLE PRECISION apo
      DOUBLE PRECISION errorpercent(0:MMAX)
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM 
 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
      
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
      DOUBLE PRECISION , EXTERNAL :: Check_conservation 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
     
     
      call lock_ambm
      call lock_tmp_array
 
!
!     Fluid phase species mass balance equations
!
      IF (SPECIES_EQ(0)) THEN 
         DO LN = 1, NMAX(0) 
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!$omp    parallel do private(IJK, APO)
            DO IJK = ijkstart3, ijkend3 
!
               IF (FLUID_AT(IJK)) THEN
	          if(strong_conservation(7)) then 
                     S_P(IJK) = (ROX_GC(IJK,LN))*VOL(IJK) 
		     S_C(IJK) = ROP_GO(IJK)*VOL(IJK)*ODT*X_GO(IJK,LN) + R_GP(IJK,LN)*VOL(IJK) 
		   else 
                     APO = ROP_G(IJK)*VOL(IJK)*ODT 
                     S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK))+ROX_GC(IJK,LN))*VOL(IJK) 
                     S_C(IJK) = APO*X_GO(IJK,LN) + X_G(IJK,LN)*ZMAX((-SUM_R_G(IJK)))&
                      *VOL(IJK) + R_GP(IJK,LN)*VOL(IJK)
		   endif
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 
            END DO
	    
	     
	    if(strong_conservation(7)) then
              CALL CONV_DIF_PHI_Strong (X_G(1,LN), DIF_G(1,LN), DISCRETIZE(7), U_G, V_G, &
               W_G, ROP_G, 0, A_M, B_M, IER)
	    else 
              CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN), DISCRETIZE(7), U_G, V_G, &
               W_G, ROP_G, 0, A_M, B_M, IER)
	    endif 

!
            CALL BC_PHI (BC_X_G(1,LN), BC_XW_G(1,LN), BC_HW_X_G(1,LN), BC_C_X_G(1,&
               LN), 0, A_M, B_M, IER) 
!
	    if(strong_conservation(7)) then
              CALL SOURCE_PHI_Strong (S_P, S_C, EP_G, X_G(1,LN), 0, A_M, B_M, IER)
	    else 
              CALL SOURCE_PHI (S_P, S_C, EP_G, X_G(1,LN), 0, A_M, B_M, IER)
	    endif 
!
            CALL CALC_RESID_S (X_G(1,LN), A_M, B_M, 0, RESID(RESID_X+(LN-1),0), &
               MAX_RESID(RESID_X+(LN-1),0), IJK_RESID(RESID_X+(LN-1),0), &
	       ZERO_X_GS, IER) 

            CALL UNDER_RELAX_S (X_G(1,LN), A_M, B_M, 0, UR_FAC(7), IER) 

            CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),0), LEQ_IT(7), LEQ_METHOD(7), &
                LEQI, LEQM, IER)
!
            CALL SOLVE_LIN_EQ ('X_g', X_G(1,LN), A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 

            CALL BOUND_X (X_G(1,LN), IJKMAX2, IER) 

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
            DO LN = 1, NMAX(M) 
!
               CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
!
!$omp    parallel do private(IJK, APO)
               DO IJK = ijkstart3, ijkend3 
!
                  IF (FLUID_AT(IJK)) THEN 
!
	            if(strong_conservation(7)) then 
                      S_P(IJK) = (ROX_SC(IJK,M,LN))*VOL(IJK) 
		      S_C(IJK) = ROP_SO(IJK,M)*VOL(IJK)*ODT*X_SO(IJK,M,LN) + R_SP(IJK,M,LN)*VOL(IJK)
		    else 
                      APO = ROP_S(IJK,M)*VOL(IJK)*ODT 
                      S_P(IJK) = APO + (ZMAX(SUM_R_S(IJK,M))+ROX_SC(IJK,M,LN))*&
                        VOL(IJK) 
                      S_C(IJK) = APO*X_SO(IJK,M,LN) + X_S(IJK,M,LN)*ZMAX((-SUM_R_S&
                        (IJK,M)))*VOL(IJK) + R_SP(IJK,M,LN)*VOL(IJK)
		    endif 
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
	       if(strong_conservation(7)) then
                 CALL CONV_DIF_PHI_Strong (X_S(1,M,LN), DIF_S(1,M,LN), DISCRETIZE(7), U_S(&
                  1,M), V_S(1,M), W_S(1,M), ROP_S(1,M), M, A_M, B_M, IER)
	       else 
                 CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), DISCRETIZE(7), U_S(&
                  1,M), V_S(1,M), W_S(1,M), ROP_S(1,M), M, A_M, B_M, IER)
	       endif 
!
               CALL BC_PHI (BC_X_S(1,M,LN), BC_XW_S(1,M,LN), BC_HW_X_S(1,M,LN), &
                  BC_C_X_S(1,M,LN), M, A_M, B_M, IER) 
!
!
	       if(strong_conservation(7)) then
                 CALL SOURCE_PHI_Strong (S_P, S_C, EPS, X_S(1,M,LN), M, A_M, B_M, IER)
	       else 
                 CALL SOURCE_PHI (S_P, S_C, EPS, X_S(1,M,LN), M, A_M, B_M, IER)
	       endif 
!
               CALL CALC_RESID_S (X_S(1,M,LN), A_M, B_M, M, RESID(RESID_X+(LN-1),&
                  M), MAX_RESID(RESID_X+(LN-1),M), IJK_RESID(RESID_X+(LN-1),M), &
                  ZERO_X_GS, IER) 
!
!
               CALL UNDER_RELAX_S (X_S(1,M,LN), A_M, B_M, M, UR_FAC(7), IER) 
!           call check_ab_m(a_m, b_m, m, .false., ier)
!           write(*,*)
!     &      resid(resid_x+(LN-1), m), max_resid(resid_x+(LN-1), m),
!     &      ijk_resid(resid_x+(LN-1), m)
!           call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!           call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &     ier)
!
               CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),M), LEQ_IT(7), LEQ_METHOD(7&
                  ), LEQI, LEQM, IER) 
!
               CALL SOLVE_LIN_EQ ('X_s', X_S(1,M,LN), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 
!            call out_array(X_s(1,m,LN), 'X_s')
!               CALL BOUND_X (X_S(1,M,LN), IJKMAX2, IER) 
!
            END DO 
         ENDIF 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_SPECIES_EQ 
      

