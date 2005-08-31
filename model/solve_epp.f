!TO DO
! check the formulation based on MCp
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_Epp(NORMs, RESs, IER)                            C
!  Purpose: Solve solids volume fraction correction equation           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-SEP-96  C
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
      SUBROUTINE SOLVE_EPP(NORMS, RESS, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE geometry
      USE pscor
      USE residual
      USE leqsol 
      USE physprop
      Use ambm
      Use tmp_array1, B_mMAX => ARRAYm1
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!Parameter to make tolerance for residual scaled with max value compatible with
!residual scaled with first iteration residual.  Increase it, to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1 !5.0D2  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Normalization factor for solids volume fraction residual
      DOUBLE PRECISION NORMs, NORMSloc
!
!                      solids volume fraction residual
      DOUBLE PRECISION RESs
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
! 

      INTEGER          M
!
!                      linear equation solver method and iterations
      INTEGER          LEQM, LEQI
!
!-----------------------------------------------
!
      call lock_ambm
      call lock_tmp_array1

!
!  Form the sparse matrix equation
!
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!//
      CALL ZERO_ARRAY (EPP, IER) 
!
      CALL CONV_SOURCE_EPP (A_M, B_M, B_mmax, IER) 
!
!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &  ier)
!
!  Find average residual, maximum residual and location
!
      NORMSloc = NORMS
      IF(NORMS == ZERO) THEN
        CALL CALC_RESID_PP (B_MMAX, ONE, RESID(RESID_P,1), MAX_RESID(RESID_P,1), &
         IJK_RESID(RESID_P,1), IER) 
	 NORMSloc = RESID(RESID_P,1)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMSloc, RESID(RESID_P,1), MAX_RESID(RESID_P,1), &
         IJK_RESID(RESID_P,1), IER) 
      RESS = RESID(RESID_P,1) 
!      write(*,*)
!     &    resid(resid_p, 1), max_resid(resid_p, 1),
!     &    ijk_resid(resid_p, 1)
!
!  Solve EP_s_prime equation
!
!
      CALL ADJUST_LEQ(RESID(RESID_P,1),LEQ_IT(1),LEQ_METHOD(1),LEQI,LEQM,IER) 
!
      CALL SOLVE_LIN_EQ ('EPp', EPP, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER) 
!      call out_array(EPp, 'EPp')
!
      call unlock_tmp_array1
      call unlock_ambm

      RETURN  
      END SUBROUTINE SOLVE_EPP 
