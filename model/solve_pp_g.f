!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_Pp_g
!  Purpose: Solve fluid pressure correction equation                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_PP_G(NORMG, RESG, IER) 

!-----------------------------------------------
! Modules
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
! Local parameters
!-----------------------------------------------
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.  
! Increase it to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1   !5.0D2  
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Normalization factor for gas pressure correction residual.
! At start of the iterate loop normg will either be 1 (i.e. not
! normalized) or a user defined value given by norm_g.  If norm_g
! was set to zero then the normalization is based on dominate 
! term in the equation
      DOUBLE PRECISION, INTENT(IN) :: NORMg
! gas pressure correction residual
      DOUBLE PRECISION, INTENT(OUT) :: RESg
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M
! Normalization factor for gas pressure correction residual
      DOUBLE PRECISION :: NORMGloc
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! temporary use of global arrays:
! arraym1 (locally b_mmax) 
! vector B_M based on dominate term in correction equation
!      DOUBLE PRECISION :: B_MMAX(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m, vector B_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      call lock_ambm
      call lock_tmp_array1

! initializing 
      CALL ZERO_ARRAY (PP_G, IER) 
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      ENDDO 

! If gas momentum equations in x and y directions are not solved return
      IF (.NOT.(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0)) .AND.&
          RO_G0 .NE. UNDEFINED) THEN
        call unlock_ambm
        call unlock_tmp_array1
        RETURN  
      ENDIF

! Forming the sparse matrix equation.      
      CALL CONV_PP_G (A_M, B_M, IER) 
      CALL SOURCE_PP_G (A_M, B_M, B_MMAX, IER) 

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)


! Find average residual, maximum residual and location
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P,0), &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0), IER)
         NORMGloc = RESID(RESID_P,0)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMGloc, NUM_RESID(RESID_P,0),  &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0), IER) 
      RESG = RESID(RESID_P,0) 
!      write(*,*) resid(resid_p, 0), max_resid(resid_p, 0), &
!         ijk_resid(resid_p, 0)


! Solve P_g_prime equation
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)
!      CALL ADJUST_LEQ(RESID(RESID_P,0),LEQ_IT(1),LEQ_METHOD(1),LEQI,LEQM,IER) 

!     call check_symmetry(A_m, 0, IER)
!     call test_lin_eq(A_M, LEQ_IT(1),LEQ_METHOD(1), LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1),0,IER)
      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_M, B_M, 0, LEQI, LEQM, &
                         LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER) 

!      call out_array(Pp_g, 'Pp_g')

      call unlock_tmp_array1
      call unlock_ambm

      RETURN  
      END SUBROUTINE SOLVE_PP_G 

