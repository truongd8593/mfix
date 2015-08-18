!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_SPECIES_EQ                                        C
!  Purpose: Solve species mass balance equations in matrix equation    C
!     form Ax=b. The center coefficient (ap) and source vector (b)     C
!     are negative.  The off-diagonal coefficients are positive.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 11-FEB-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_SPECIES_EQ(IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE ChiScheme
      USE bc
      USE compar
      USE drag
      USE energy
      USE fldvar
      USE geometry
      USE indices
      USE leqsol
      USE matrix
      USE mflux
      USE mpi_utility
      USE output
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE residual
      USE run
      USE rxns
      USE sendrecv
      USE toleranc
      USE ur_facs
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
      use functions
      use ps
      use toleranc

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M
! species index
      INTEGER :: LN
! previous time step term
      DOUBLE PRECISION :: apo
! Indices
      INTEGER :: IJK
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI
! tmp array to pass to set_chi
      DOUBLE PRECISION :: X_s_temp(DIMENSION_3, DIMENSION_N_s)

! Arrays for storing errors:
! 130 - Gas phase species equation diverged
! 131 - Solids phase species equation diverged
! 13x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! temporary use of global arrays:
! array1 (locally s_p)
! source lhs: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
!      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c)
! source rhs vector: constant part becomes part of b_m vector
!      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
! alias for solids volume fraction
!      DOUBLE PRECISION :: eps(DIMENSION_3)
! array4 (locally vxgama)
!      DOUBLE PRECISION :: vxgama(DIMENSION_3)
! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: Check_conservation
!-----------------------------------------------

      call lock_ambm       ! locks arrys a_m and b_m
      call lock_tmp_array  ! locks array1,array2,array3,array4
                           ! (locally s_p, s_c, eps, vxgama)

! Initialize error flag.
      Err_l = 0

! Fluid phase species mass balance equations
! ---------------------------------------------------------------->>>
      IF (SPECIES_EQ(0)) THEN
         IF(chi_scheme) call set_chi(DISCRETIZE(7), X_g, NMAX(0), &
                                     U_g, V_g, W_g)

! looping over species
         DO LN = 1, NMAX(0)
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)
!!$omp    parallel do private(IJK, APO)
            DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN
! calculate the source terms to be used in the a matrix and b vector
                   APO = ROP_GO(IJK)*VOL(IJK)*ODT
                   S_P(IJK) = APO+&
                      (ZMAX(SUM_R_G(IJK))+ROX_GC(IJK,LN))*VOL(IJK)
                   S_C(IJK) = APO*X_GO(IJK,LN) + &
                      X_G(IJK,LN)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) +&
                      R_GP(IJK,LN)*VOL(IJK)
               ELSE
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
               ENDIF
            ENDDO

! calculate the convection-diffusion terms
            IF(.NOT.ADDED_MASS) THEN
               CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN), &
                  DISCRETIZE(7), U_G, V_G, W_G, &
                  Flux_gE, Flux_gN, Flux_gT, 0, A_M, B_M)
            ELSE
               CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN),&
                  DISCRETIZE(7), U_G, V_G, W_G, &
                  Flux_gSE, Flux_gSN, Flux_gST, 0, A_M, B_M)
            ENDIF

! calculate standard bc
            CALL BC_PHI (X_G(1,LN), BC_X_G(1,LN), BC_XW_G(1,LN), &
               BC_HW_X_G(1,LN), BC_C_X_G(1,LN), 0, A_M, B_M)

! set the source terms in a and b matrix form
            CALL SOURCE_PHI (S_P, S_C, EP_G, X_G(1,LN), 0, A_M, B_M)

! Add point sources.
            IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (X_G(1,LN), &
               PS_X_G(:,LN), PS_MASSFLOW_G, 0, A_M, B_M)

            CALL CALC_RESID_S (X_G(1,LN), A_M, B_M, 0, &
               NUM_RESID(RESID_X+(LN-1),0), &
               DEN_RESID(RESID_X+(LN-1),0), RESID(RESID_X+(LN-1),0), &
               MAX_RESID(RESID_X+(LN-1),0), IJK_RESID(RESID_X+(LN-1),0), &
               ZERO_X_GS)

            CALL UNDER_RELAX_S (X_G(1,LN), A_M, B_M, 0, UR_FAC(7))

            CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),0), LEQ_IT(7), &
               LEQ_METHOD(7), LEQI, LEQM)

! solve the phi equation
            CALL SOLVE_LIN_EQ ('X_g', 7, X_G(1,LN), A_M, B_M, 0, &
               LEQI, LEQM, LEQ_SWEEP(7), LEQ_TOL(7), LEQ_PC(7), IER)

! Check for linear solver divergence.
               IF(ier == -2) Err_l(myPE) = 130

            CALL BOUND_X (X_G(1,LN), IJKMAX2)

         ENDDO    ! end do loop (ln = 1, nmax(0)
         IF(chi_scheme) call unset_chi()
      ENDIF
! end fluid phase species equations
! ----------------------------------------------------------------<<<

! Solids phase species balance equations
! ---------------------------------------------------------------->>>
      DO M = 1, SMAX
         IF (SPECIES_EQ(M)) THEN
            IF(chi_scheme) THEN
               DO LN = 1, NMAX(M)
                  DO IJK = ijkstart3, ijkend3
                     X_S_temp(IJK, LN) = X_S(IJK,M,LN)
                  ENDDO
                ENDDO
              call set_chi(DISCRETIZE(7), X_S_temp, NMAX(M), &
                           U_S(1,M), V_S(1,M), W_S(1,M))
            ENDIF ! for chi_scheme

            DO LN = 1, NMAX(M)
               CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

!!$omp    parallel do private(IJK, APO)
               DO IJK = ijkstart3, ijkend3
                  IF (FLUID_AT(IJK)) THEN
                    APO = ROP_SO(IJK,M)*VOL(IJK)*ODT
                    S_P(IJK) = APO + &
                       (ZMAX(SUM_R_S(IJK,M))+ROX_SC(IJK,M,LN))*VOL(IJK)
                    S_C(IJK) = APO*X_SO(IJK,M,LN) + &
                       X_S(IJK,M,LN)*ZMAX((-SUM_R_S(IJK,M)))*VOL(IJK) + &
                       R_SP(IJK,M,LN)*VOL(IJK)
                    EPS(IJK) = EP_S(IJK,M)
                  ELSE
                     S_P(IJK) = ZERO
                     S_C(IJK) = ZERO
                     EPS(IJK) = ZERO
                  ENDIF
               ENDDO

               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), &
                    DISCRETIZE(7), U_S(1,M), V_S(1,M), W_S(1,M), &
                    Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M)
               ELSE
                  CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), &
                    DISCRETIZE(7), U_S(1,M), V_S(1,M), W_S(1,M), &
                    Flux_sSE, Flux_sSN, Flux_sST, M, A_M, B_M)
               ENDIF

               CALL BC_PHI (X_S(1,M,LN), BC_X_S(1,M,LN), &
                  BC_XW_S(1,M,LN), BC_HW_X_S(1,M,LN), &
                  BC_C_X_S(1,M,LN), M, A_M, B_M)

               CALL SOURCE_PHI (S_P, S_C, EPS, X_S(1,M,LN), M, A_M, B_M)

! Add point sources.
               IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (X_S(1,M,LN), &
                  PS_X_S(:,M,LN), PS_MASSFLOW_S(:,M), M, A_M, B_M)

               CALL CALC_RESID_S (X_S(1,M,LN), A_M, B_M, M, &
                  NUM_RESID(RESID_X+(LN-1),M), &
                  DEN_RESID(RESID_X+(LN-1),M), RESID(RESID_X+(LN-1),&
                  M), MAX_RESID(RESID_X+(LN-1),M), IJK_RESID(RESID_X+(LN-1),M), &
                  ZERO_X_GS)

               CALL UNDER_RELAX_S (X_S(1,M,LN), A_M, B_M, M, UR_FAC(7))

!               call check_ab_m(a_m, b_m, m, .false., ier)
!               write(*,*) resid(resid_x+(LN-1), m), &
!                  max_resid(resid_x+(LN-1), m), &
!                  ijk_resid(resid_x+(LN-1), m)
!               call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!               call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1,-3,M),&
!                  1, DO_K, ier)

               CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),M), LEQ_IT(7), &
                  LEQ_METHOD(7), LEQI, LEQM)

               CALL SOLVE_LIN_EQ ('X_s', 7, X_S(1,M,LN), A_M, B_M, M,&
                  LEQI, LEQM, LEQ_SWEEP(7), LEQ_TOL(7), LEQ_PC(7), IER)

! Check for linear solver divergence.
               IF(ier == -2) Err_l(myPE) = 131

               CALL BOUND_X (X_S(1,M,LN), IJKMAX2)
!               call out_array(X_s(1,m,LN), 'X_s')

            END DO

            if(chi_scheme) call unset_chi()
         ENDIF ! check for any species in phase m
      END DO ! for m = 1, mmax
! end solids phases species equations
! ----------------------------------------------------------------<<<

      call unlock_ambm
      call unlock_tmp_array


! If the linear solver diverged, species mass fractions may take on
! unphysical values. To prevent them from propogating through the domain
! or causing failure in other routines, force an exit from iterate and
! reduce the time step.
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)

      RETURN
      END SUBROUTINE SOLVE_SPECIES_EQ
