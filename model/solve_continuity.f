!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_CONTINUITY(IER)                                  C
!  Purpose: Solve for volume fractions                                 C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2-JUL-96   C
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
      SUBROUTINE SOLVE_CONTINUITY(IER) 
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
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE residual
      USE cont
      USE leqsol
      Use ambm 
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
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
! 
!-----------------------------------------------
!

      call lock_ambm

      IF (DO_CONT(0)) THEN 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!
         CALL CONV_ROP_G (A_M, B_M, IER) 
!
         CALL SOURCE_ROP_G (A_M, B_M, IER) 
!
         CALL CALC_RESID_C (ROP_G, A_M, B_M, 0, RESID(RESID_RO,0), MAX_RESID(&
            RESID_RO,0), IJK_RESID(RESID_RO,0), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .true., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_ro, 0), max_resid(resid_ro, 0),
!     &    ijk_resid(resid_ro, 0)
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_RO,0), LEQ_IT(2), LEQ_METHOD(2), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('ROP_g', ROP_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(2), LEQ_TOL(2),IER) 
         CALL ADJUST_ROP (ROP_G, IER) 
!        call out_array(ROP_g, 'rop_g')
      ENDIF 
!
!
      DO M = 1, MMAX 
         IF (DO_CONT(M)) THEN 
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
            CALL CONV_ROP_S (A_M, B_M, M, IER) 
            CALL SOURCE_ROP_S (A_M, B_M, M, IER) 
            CALL CALC_RESID_C (ROP_S(1,M), A_M, B_M, M, RESID(RESID_RO,M), &
               MAX_RESID(RESID_RO,M), IJK_RESID(RESID_RO,M), IER) 
!
!          call check_ab_m(a_m, b_m, m, .true., ier)
!          write(*,*)
!     &    resid(resid_ro, m), max_resid(resid_ro, m),
!     &    ijk_resid(resid_ro, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1,
!     &      DO_K, ier)
!
            CALL ADJUST_LEQ (RESID(RESID_RO,M), LEQ_IT(2), LEQ_METHOD(2), LEQI&
               , LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('ROP_s', ROP_S(1,M), A_M, B_M, M, LEQI, LEQM,&
	                     LEQ_SWEEP(2), LEQ_TOL(2),IER) 
            CALL ADJUST_ROP (ROP_S(1,M), IER) 
!          call out_array(rop_s(1,m), 'rop_s')
         ENDIF 
      END DO 
      
      call unlock_ambm

      
      RETURN  
      END SUBROUTINE SOLVE_CONTINUITY 
