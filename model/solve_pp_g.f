!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_Pp_g(NORMg, RESg, IER)
!  Purpose: Solve fluid pressure correction equation                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
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
      SUBROUTINE SOLVE_PP_G(NORMG, RESG, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE geometry
      USE pgcor
      USE residual
      USE leqsol 
      USE run
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
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1   !5.0D2  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Normalization factor for gas pressure residual
      DOUBLE PRECISION NORMg, NORMGloc
!
!                      gas pressure residual
      DOUBLE PRECISION RESg
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

      INTEGER          M 
!
!                      linear equation solver method and iterations
      INTEGER          LEQM, LEQI

!-----------------------------------------------
!
      call lock_ambm
      call lock_tmp_array1

!
!     If gas momentum equations in x and y directions are not solved return
!
!
      CALL ZERO_ARRAY (PP_G, IER) 
      IF (.NOT.(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0)) .AND. RO_G0 .NE. UNDEFINED) THEN
        call unlock_ambm
        call unlock_tmp_array1
        RETURN  
      ENDIF
!
!  Form the sparse matrix equation
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      CALL CONV_PP_G (A_M, B_M, IER) 
!
      CALL SOURCE_PP_G (A_M, B_M, B_MMAX, IER) 
!
!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &  ier)
!
!  Find average residual, maximum residual and location
!
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
        CALL CALC_RESID_PP (B_MMAX, ONE, RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0), IER)
         NORMGloc = RESID(RESID_P,0)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMGloc, RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0), IER) 
      RESG = RESID(RESID_P,0) 
!      write(*,*)
!     &    resid(resid_p, 0), max_resid(resid_p, 0),
!     &    ijk_resid(resid_p, 0)
!
!  Solve P_g_prime equation
!
!
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)
!      CALL ADJUST_LEQ(RESID(RESID_P,0),LEQ_IT(1),LEQ_METHOD(1),LEQI,LEQM,IER) 
!
      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER) 
!      call out_array(Pp_g, 'Pp_g')
!
      call unlock_tmp_array1
      call unlock_ambm

      RETURN  
      END SUBROUTINE SOLVE_PP_G 

!// Comments on the modifications for DMP version implementation      
!//    Initialize PP_G with call zero_array
