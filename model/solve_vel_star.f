!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E, D_N and D_T claculations       C
!           and account for the averaged Solid-Solid drag              C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEb-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To flag the solids calculations: Continuum or DES          C
!  And to change the gas arrays incorporating the drag, when doing DES C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Author: Rahul Garg                                 Date: 01-Aug-07  C
!  Purpose: For Discrete element simulation, no need to call           C
!   DES_CALC_D_E to caculate pressure correction coeff's. CALC_D_E     C
!   already contains this spl case. Same for North and Top cases.      C

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
      USE ghdtheory
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs 
      USE pgcor
      USE pscor
      USE leqsol
      Use ambm
      Use tmp_array1,  VxF_gs => Arraym1
      Use tmp_array,  VxF_ss => ArrayLM     !S. Dartevelle, LANL, MARCH 2004
      USE compar
      USE discretelement
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
      INTEGER          m, UVEL, VVEL, WVEL 
      INTEGER          IJK
! 
!                      temporary velocity arrays 
      DOUBLE PRECISION U_gtmp(DIMENSION_3),  V_gtmp(DIMENSION_3), W_gtmp(DIMENSION_3)
      DOUBLE PRECISION U_stmp(DIMENSION_3, DIMENSION_M),  V_stmp(DIMENSION_3, DIMENSION_M), W_stmp(DIMENSION_3, DIMENSION_M)
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

      call lock_ambm
      call lock_tmp_array1
      call lock_tmp_array2
      
!     Store the velocities so that the order of solving the momentum equations does not matter
      DO IJK = ijkstart3, ijkend3
        U_gtmp(IJK) = U_g(IJK)
        V_gtmp(IJK) = V_g(IJK)
        W_gtmp(IJK) = W_g(IJK)
      ENDDO
      DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
          U_stmp(IJK, M) = U_s(IJK, M)
          V_stmp(IJK, M) = V_s(IJK, M)
          W_stmp(IJK, M) = W_s(IJK, M)
        ENDDO
      ENDDO
      
      IF(DISCRETE_ELEMENT) THEN
         UVEL = 0
         VVEL = 0
         WVEL = 0
      END IF

      IF (MMAX == 0) CALL ZERO_ARRAY (VXF_GS(1,1), IER)
      IF (MMAX == 1) CALL ZERO_ARRAY (VXF_SS(1,1), IER)

! Calculate U_m_star and residuals
! ----------------------------------------
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      
      CALL CONV_DIF_U_G (A_M, B_M, IER) 
      
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CONV_DIF_U_S (A_M, B_M, IER) 
      END IF

      CALL SOURCE_U_G (A_M, B_M, IER) 
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL SOURCE_U_S (A_M, B_M, IER) 
      END IF

      IF (MMAX > 0) CALL VF_GS_X (VXF_GS, IER)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL VF_SS_X (VXF_SS, IER)   !S. Dartevelle, LANL, Feb.2004
      END IF

      CALL CALC_D_E (A_M, VXF_GS, VXF_SS, D_E, IER)  !S. Dartevelle, LANL, Feb.2004
  
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL CALC_E_E (A_M, MCP, E_E, IER) 
         IF (MMAX > 0) CALL PARTIAL_ELIM_U (U_G, U_S, VXF_GS, A_M, B_M, IER) 
      END IF

      CALL ADJUST_A_U_G (A_M, B_M, IER) 
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL ADJUST_A_U_S (A_M, B_M, IER)
      END IF 
      
      
      IF(DES_CONTINUUM_COUPLED) THEN
         UVEL = 1
         CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
         UVEL = 0
      END IF

      IF (MOMENTUM_X_EQ(0)) THEN 
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, NUM_RESID(RESID_U,0), &
            DEN_RESID(RESID_U,0), RESID(RESID_U,0), &
            MAX_RESID(RESID_U,0), IJK_RESID(RESID_U,0), IER) 
         CALL UNDER_RELAX_U (U_G, A_M, B_M, 0, UR_FAC(3), IER) 

!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_u, 0), max_resid(resid_u, 0),
!     &    ijk_resid(resid_u, 0)
!
      ENDIF 
      

      DO M = 1, MMAX   
       IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
         IF (MOMENTUM_X_EQ(M)) THEN 
            CALL CALC_RESID_U (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M,  &
               NUM_RESID(RESID_U,M), DEN_RESID(RESID_U,M), RESID&
               (RESID_U,M), MAX_RESID(RESID_U,M), IJK_RESID(RESID_U,M), IER) 
            CALL UNDER_RELAX_U (U_S(1,M), A_M, B_M, M, UR_FAC(3), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_u, m), max_resid(resid_u, m),
!     &      ijk_resid(resid_u, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
       ENDIF ! for GHD Theory
      END DO 

      
      IF (MOMENTUM_X_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)

         CALL ADJUST_LEQ (RESID(RESID_U,0), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
            LEQM, IER) 

         CALL SOLVE_LIN_EQ ('U_g', 3, U_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(3), LEQ_TOL(3),  LEQ_PC(3), IER) 
!        call out_array(u_g, 'u_g')
      ENDIF 

      DO M = 1, MMAX 
       IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
         IF (MOMENTUM_X_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
            CALL ADJUST_LEQ (RESID(RESID_U,M), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
               LEQM, IER) 

            CALL SOLVE_LIN_EQ ('U_s', 3, U_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
               LEQ_SWEEP(3), LEQ_TOL(3),  LEQ_PC(3), IER) 
!          call out_array(u_s(1,m), 'u_s')
         ENDIF 
       ENDIF ! for GHD Theory
      ENDDO 
      
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 


! Calculate V_m_star and residuals
! ----------------------------------------
      CALL CONV_DIF_V_G (A_M, B_M, IER) 
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CONV_DIF_V_S (A_M, B_M, IER) 
      END IF

      CALL SOURCE_V_G (A_M, B_M, IER) 
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL SOURCE_V_S (A_M, B_M, IER)
      END IF 

      IF (MMAX > 0) CALL VF_GS_Y (VXF_GS, IER)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL VF_SS_Y (VXF_SS, IER)    !S. Dartevelle, LANL, Feb.2004
      END IF

      CALL CALC_D_N (A_M, VXF_GS, VXF_SS, D_N, IER) !S. Dartevelle, LANL, Feb.2004

      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL CALC_E_N (A_M, MCP, E_N, IER) 
         IF (MMAX > 0) CALL PARTIAL_ELIM_V (V_G, V_S, VXF_GS, A_M, B_M, IER) 
      END IF

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      CALL ADJUST_A_V_G (A_M, B_M, IER)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL ADJUST_A_V_S (A_M, B_M, IER) 
      END IF

      IF(DES_CONTINUUM_COUPLED) THEN
         VVEL = 1     
         CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
         VVEL = 0
      END IF
      
      IF (MOMENTUM_Y_EQ(0)) THEN 
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, NUM_RESID(RESID_V,0), &
            DEN_RESID(RESID_V,0), RESID(RESID_V,0), &
            MAX_RESID(RESID_V,0), IJK_RESID(RESID_V,0), IER) 
         CALL UNDER_RELAX_V (V_G, A_M, B_M, 0, UR_FAC(4), IER) 

!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_v, 0), max_resid(resid_v, 0),
!     &    ijk_resid(resid_v, 0)
      ENDIF 

      DO M = 1, MMAX 
       IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
         IF (MOMENTUM_Y_EQ(M)) THEN 
            CALL CALC_RESID_V (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, &
                NUM_RESID(RESID_V,M), DEN_RESID(RESID_V,M), RESID&
               (RESID_V,M), MAX_RESID(RESID_V,M), IJK_RESID(RESID_V,M), IER) 
            CALL UNDER_RELAX_V (V_S(1,M), A_M, B_M, M, UR_FAC(4), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_v, m), max_resid(resid_v, m),
!     &      ijk_resid(resid_v, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
       ENDIF ! for GHD Theory
      ENDDO 

      IF (MOMENTUM_Y_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)

         CALL ADJUST_LEQ (RESID(RESID_V,0), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
            LEQM, IER) 

         CALL SOLVE_LIN_EQ ('V_g', 4, V_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(4), LEQ_TOL(4),  LEQ_PC(4), IER) 
!        call out_array(v_g, 'v_g')
      ENDIF 

      DO M = 1, MMAX 
       IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
         IF (MOMENTUM_Y_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)

            CALL ADJUST_LEQ (RESID(RESID_V,M), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
               LEQM, IER) 

            CALL SOLVE_LIN_EQ ('V_s', 4, V_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
               LEQ_SWEEP(4), LEQ_TOL(4), LEQ_PC(4), IER) 
!          call out_array(v_s(1,m), 'v_s')
         ENDIF 
       ENDIF ! for GHD Theory
      END DO 



! Calculate W_m_star and residuals
! ----------------------------------------
      IF (DO_K)THEN

        DO M = 0, MMAX 
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
        END DO 
        CALL CONV_DIF_W_G (A_M, B_M, IER) 

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
        IF(.NOT.DISCRETE_ELEMENT) THEN
          CALL CONV_DIF_W_S (A_M, B_M, IER) 
        END IF

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
         
        CALL SOURCE_W_G (A_M, B_M, IER) 

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
    
        IF(.NOT.DISCRETE_ELEMENT) THEN    
          CALL SOURCE_W_S (A_M, B_M, IER) 
        END IF

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
! call mfix_exit(myPE)

        IF (MMAX > 0) CALL VF_GS_Z (VXF_GS, IER)
        IF(.NOT.DISCRETE_ELEMENT) THEN
          IF (MMAX > 0) CALL VF_SS_Z (VXF_SS, IER)   !S. Dartevelle, LANL, Feb.2004
        END IF
        
        CALL CALC_D_T (A_M, VXF_GS, VXF_SS, D_T, IER) !S. Dartevelle, LANL, Feb.2004

        IF(.NOT.DISCRETE_ELEMENT) THEN
          IF (MMAX > 0) CALL CALC_E_T (A_M, MCP, E_T, IER) 
          IF (MMAX > 0) CALL PARTIAL_ELIM_W (W_G, W_S, VXF_GS, A_M, B_M, IER)
        END IF 

        CALL ADJUST_A_W_G (A_M, B_M, IER)
        IF(.NOT.DISCRETE_ELEMENT) THEN 
          CALL ADJUST_A_W_S (A_M, B_M, IER)
        END IF 

        IF(DIMN.EQ.3) THEN
          IF(DES_CONTINUUM_COUPLED) THEN
            WVEL = 1
            CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
            WVEL = 0
          END IF
        END IF
                                              
        IF (MOMENTUM_Z_EQ(0)) THEN 
          CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, NUM_RESID(RESID_W,0), &
            DEN_RESID(RESID_W,0), RESID(RESID_W,0), &
            MAX_RESID(RESID_W,0), IJK_RESID(RESID_W,0), IER) 
          CALL UNDER_RELAX_W (W_G, A_M, B_M, 0, UR_FAC(5), IER) 

!        call check_ab_m(a_m, b_m, 0, .false., ier)
!     &      resid(resid_w, 0), max_resid(resid_w, 0),
!     &      ijk_resid(resid_w, 0)
        ENDIF 

        DO M = 1, MMAX 
         IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
          IF (MOMENTUM_Z_EQ(M)) THEN 
            CALL CALC_RESID_W (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, &
                NUM_RESID(RESID_W,M), DEN_RESID(RESID_W,M), RESID&
               (RESID_W,M), MAX_RESID(RESID_W,M), IJK_RESID(RESID_W,M), IER) 
            CALL UNDER_RELAX_W (W_S(1,M), A_M, B_M, M, UR_FAC(5), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_w, m), max_resid(resid_w, m),
!     &      ijk_resid(resid_w, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
          ENDIF 
         ENDIF ! for GHD Theory
        ENDDO 

        IF (MOMENTUM_Z_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)

          CALL ADJUST_LEQ (RESID(RESID_W,0), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
            LEQM, IER) 

          CALL SOLVE_LIN_EQ ('W_g', 5, W_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER) 
!        call out_array(w_g, 'w_g')
        ENDIF 

        DO M = 1, MMAX 
         IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
          IF (MOMENTUM_Z_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)

            CALL ADJUST_LEQ (RESID(RESID_W,M), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
               LEQM, IER) 

            CALL SOLVE_LIN_EQ ('W_s', 5, W_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
               LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER) 
!          call out_array(w_s(1,m), 'w_s')
          ENDIF 
         ENDIF ! for GHD Theory
        ENDDO 
      ENDIF
      
!     Now update all velocity components
      DO IJK = ijkstart3, ijkend3
        U_g(IJK) = U_gtmp(IJK)
        V_g(IJK) = V_gtmp(IJK)
        W_g(IJK) = W_gtmp(IJK)
      ENDDO
      DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
        U_s(IJK, M) = U_stmp(IJK, M)
        V_s(IJK, M) = V_stmp(IJK, M)
        W_s(IJK, M) = W_stmp(IJK, M)
        ENDDO
      ENDDO

! modification for GHD theory to compute species velocity: Ui = Joi/(mi ni) + U.
      IF(TRIM(KT_TYPE) == 'GHD') THEN
	DO M = 1, SMAX
          CALL GHDMassFlux(M,IER) ! to compute solid species mass flux
	  CALL UpdateSpeciesVelocities(M,IER) ! located at end of ghdMassFlux.f file
        ENDDO 
      ENDIF
! end of modification for GHD theory

      call unlock_ambm
      call unlock_tmp_array1
      call unlock_tmp_array2

      RETURN  
      END SUBROUTINE SOLVE_VEL_STAR 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
