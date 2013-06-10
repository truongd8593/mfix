!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_ENERGY_EQ                                         C
!  Purpose: Solve energy equations                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Eliminate energy calculations when doing DEM               C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_ENERGY_EQ(IER) 

!-----------------------------------------------
! Modules
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
      USE des_thermo
      USE mflux     
      USE mpi_utility      
      USE sendrecv 
      USE ps

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
      INTEGER :: TMP_SMAX
!  Cp * Flux 
      DOUBLE PRECISION :: CpxFlux_E(DIMENSION_3), &
                          CpxFlux_N(DIMENSION_3), &
                          CpxFlux_T(DIMENSION_3) 
! previous time step term  
      DOUBLE PRECISION :: apo 
! Indices 
      INTEGER :: IJK, I, J, K 
! linear equation solver method and iterations 
      INTEGER :: LEQM, LEQI 

! Local/Global error flags.
      LOGICAL :: lDiverged, gDiverged

! temporary use of global arrays:
! arraym1 (locally vxgama) 
! the volume x average gas-solids heat transfer at cell centers
!      DOUBLE PRECISION :: VXGAMA(DIMENSION_3, DIMENSION_M)
! array1 (locally s_p) 
! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive      
!      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c)
! source vector: constant part becomes part of b_m vector
!      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
!      DOUBLE PRECISION :: EPS(DIMENSION_3)      

! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'radtn1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'radtn2.inc'
!-----------------------------------------------

! Initialize error flag.
      lDiverged = .FALSE.

      call lock_ambm         ! locks arrys a_m and b_m
      call lock_tmp_array    ! locks arraym1 (locally vxgama) 
      call lock_tmp_array1   ! locks array1, array2, array3 
                             ! (locally s_p, s_c, eps) 

      TMP_SMAX = SMAX
      IF(DISCRETE_ELEMENT) THEN
         TMP_SMAX = 0   ! Only the gas calculations are needed
      ENDIF            

! initializing      
      DO M = 0, TMP_SMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      ENDDO 

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK) 

         IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gE(IJK)
            ELSE
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gSE(IJK)
            ENDIF
         ENDIF
         
         IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gN(IJK)
            ELSE
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gSN(IJK)
            ENDIF
         ENDIF
         
         IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gT(IJK)
            ELSE
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gST(IJK)
            ENDIF
         ENDIF

         IF (FLUID_AT(IJK)) THEN 
            APO = ROP_GO(IJK)*C_PG(IJK)*VOL(IJK)*ODT 
            S_P(IJK) = APO + S_RPG(IJK)*VOL(IJK) 
            S_C(IJK) = APO*T_GO(IJK)-HOR_G(IJK)*VOL(IJK)+S_RCG(IJK)*VOL(IJK)
         ELSE 
            S_P(IJK) = ZERO 
            S_C(IJK) = ZERO 
         ENDIF 
      ENDDO 

! Account for heat transfer between the discrete particles and the gas phase.
      IF( DES_CONV_EQ ) CALL DES_Hgm(S_C, S_P)

! calculate the convection-diffusion terms
      CALL CONV_DIF_PHI (T_g, K_G, DISCRETIZE(6), U_G, V_G, W_G, &
         CpxFlux_E, CpxFlux_N, CpxFlux_T, 0, A_M, B_M, IER) 

! calculate standard bc
      CALL BC_PHI (T_g, BC_T_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, 0, A_M, B_M, IER) 

! set the source terms in a and b matrix equation form
      CALL SOURCE_PHI (S_P, S_C, EP_G, T_G, 0, A_M, B_M, IER) 

! add point sources
      IF(POINT_SOURCE) CALL POINT_SOURCE_ENERGY_EQ (T_g, DIM_N_g, &
         Thigh_g, Tlow_g, Tcom_g, Ahigh_g, Alow_g, MW_g, BC_T_g,  &
         BC_X_g, BC_MASSFLOW_g, 0, A_M, B_M, IER)


      DO M = 1, TMP_SMAX 
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK) 

            IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sE(IJK,M)
               ELSE   ! M=M_AM is the only phase for which virtual mass is added
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sSE(IJK)
               ENDIF
            ENDIF
         
            IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sN(IJK,M)
               ELSE
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sSN(IJK)
               ENDIF
            ENDIF

            IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sT(IJK,M)
               ELSE
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sST(IJK)
               ENDIF
            ENDIF

            IF (FLUID_AT(IJK)) THEN 
               APO = ROP_SO(IJK,M)*C_PS(IJK,M)*VOL(IJK)*ODT 
               S_P(IJK) = APO + S_RPS(IJK,M)*VOL(IJK) 
               S_C(IJK) = APO*T_SO(IJK,M) - HOR_S(IJK,M)*VOL(IJK) + &
                  S_RCS(IJK,M)*VOL(IJK) 
               VXGAMA(IJK,M) = GAMA_GS(IJK,M)*VOL(IJK) 
               EPS(IJK) = EP_S(IJK,M) 
            ELSE 
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 
               VXGAMA(IJK,M) = ZERO 
               EPS(IJK) = ZERO 
            ENDIF 
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
         CALL CONV_DIF_PHI (T_s(1,M), K_S(1,M), DISCRETIZE(6), &
            U_S(1,M), V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, &
            CpxFlux_T, M, A_M, B_M, IER) 

! calculate standard bc
         CALL BC_PHI (T_s(1,M), BC_T_S(1,M), BC_TW_S(1,M), &
            BC_HW_T_S(1,M), BC_C_T_S(1,M), M, A_M, B_M, IER) 

! set the source terms in a and b matrix equation form
         CALL SOURCE_PHI (S_P, S_C, EPS, T_S(1,M), M, A_M, B_M, IER) 

! Add point sources.
         IF(POINT_SOURCE) CALL POINT_SOURCE_ENERGY_EQ (T_s(:,M),       &
            DIM_N_s, Thigh_s(M,:), Tlow_s(M,:), Tcom_s(M,:),           &
            Ahigh_s(:,M,:), Alow_s(:,M,:), MW_s, BC_T_s(:,M),          &
            BC_X_s(:,M,:), BC_MASSFLOW_s(:,M), M, A_M, B_M, IER)

      ENDDO   ! end do (m=1,tmp_smax)

! use partial elimination on interphase heat transfer term
      IF (TMP_SMAX > 0) CALL PARTIAL_ELIM_S (T_G, T_S, VXGAMA, A_M, B_M, IER) 

      CALL CALC_RESID_S (T_G, A_M, B_M, 0, NUM_RESID(RESID_T,0),& 
         DEN_RESID(RESID_T,0), RESID(RESID_T,0), MAX_RESID(RESID_T,&
         0), IJK_RESID(RESID_T,0), ZERO, IER) 

      CALL UNDER_RELAX_S (T_G, A_M, B_M, 0, UR_FAC(6), IER) 

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      write(*,*) &
!         resid(resid_t, 0), max_resid(resid_t, 0), &
!         ijk_resid(resid_t, 0)


      DO M = 1, TMP_SMAX 
         CALL CALC_RESID_S (T_S(1,M), A_M, B_M, M, NUM_RESID(RESID_T,M), &
            DEN_RESID(RESID_T,M), RESID(RESID_T,M), MAX_RESID(&
            RESID_T,M), IJK_RESID(RESID_T,M), ZERO, IER) 

         CALL UNDER_RELAX_S (T_S(1,M), A_M, B_M, M, UR_FAC(6), IER) 
      ENDDO 

! set/adjust linear equation solver method and iterations 
      CALL ADJUST_LEQ(RESID(RESID_T,0), LEQ_IT(6), LEQ_METHOD(6), &
         LEQI, LEQM, IER) 
!      call test_lin_eq(a_m(1, -3, 0), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)

      CALL SOLVE_LIN_EQ ('T_g', 6, T_G, A_M, B_M, 0, LEQI, LEQM, &
         LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)  
! Check for linear solver divergence.
      IF(ier == -2) lDiverged = .TRUE.

! bound temperature in any fluid or flow boundary cells
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.WALL_AT(IJK))&
            T_g(IJK) = MIN(TMAX, MAX(TMIN, T_g(IJK)))
      ENDDO

!      call out_array(T_g, 'T_g')

      DO M = 1, TMP_SMAX 
         CALL ADJUST_LEQ (RESID(RESID_T,M), LEQ_IT(6), LEQ_METHOD(6), &
            LEQI, LEQM, IER) 
!         call test_lin_eq(a_m(1, -3, M), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)
         CALL SOLVE_LIN_EQ ('T_s', 6, T_S(1,M), A_M, B_M, M, LEQI, &
            LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER) 

! Check for linear solver divergence.
         IF(ier == -2) lDiverged = .TRUE.

! bound temperature in any fluid or flow boundary cells
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.WALL_AT(IJK))&
               T_s(IJK, M) = MIN(TMAX, MAX(TMIN, T_s(IJK, M))) 
         ENDDO
      ENDDO   ! end do (m=1, tmp_smax)
      
      call unlock_ambm
      call unlock_tmp_array
      call unlock_tmp_array1

! If the linear solver diverged, temperatures may take on unphysical
! values. To prevent them from propogating through the domain or
! causing failure in other routines, force an exit from iterate and
! reduce the time step.

      CALL GLOBAL_ALL_OR(lDiverged, gDiverged)
      if(gDiverged) IER = -100

      
      RETURN  
      END SUBROUTINE SOLVE_ENERGY_EQ 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_ENERGY_EQ                                  C

!  Purpose: Adds point sources to the energy equations.                C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_ENERGY_EQ(T_x, lDIM_N, Thigh, Tlow, &
         Tcom, Ahigh, Alow, lMW, BC_T, BC_X, BC_FLOW, M, A_M, B_M, IER) 

      use compar    
      use constant
      use geometry
      use indices
      use physprop
      use ps
      use run

! To be removed upon complete integration of point source routines.
      use bc
      use usr

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------

      INTEGER, intent(in) :: lDIM_N

      DOUBLE PRECISION, intent(in) :: Thigh(lDIM_N)
      DOUBLE PRECISION, intent(in) :: Tlow( lDIM_N)
      DOUBLE PRECISION, intent(in) :: Tcom(lDIM_N)
      DOUBLE PRECISION, intent(in) :: Ahigh(7, lDIM_N)
      DOUBLE PRECISION, intent(in) :: Alow( 7, lDIM_N)
      DOUBLE PRECISION, intent(in) :: lMW(lDIM_N)

! Vector b_m 
      DOUBLE PRECISION, INTENT(IN) :: T_x(DIMENSION_3) 


! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: BC_T(DIMENSION_BC)
      DOUBLE PRECISION, INTENT(IN) :: BC_X(DIMENSION_BC, lDIM_N)

! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: BC_FLOW(DIMENSION_BC) 

      INTEGER, intent(in) :: M

! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 

! Vector b_m 
      DOUBLE PRECISION, intent(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M) 

! Error index 
      INTEGER, intent(INOUT) :: IER 

!----------------------------------------------- 
! Local Variables
!----------------------------------------------- 

! Indices 
      INTEGER :: IJK, I, J, K
      INTEGER :: BCV, N

! terms of bm expression
      DOUBLE PRECISION pSource, lMass

      DOUBLE PRECISION intCp_Tijk, intCp_Tref


!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

! External function. Integrates the temperature-dependent specific 
! heat from zero to T.
      DOUBLE PRECISION, EXTERNAL :: calc_ICpoR


      BC_LP: do BCV = 50, DIMENSION_BC
         if(POINT_SOURCES(BCV) == 0) cycle BC_LP

         do k = BC_K_B(BCV), BC_K_T(BCV)
         do j = BC_J_S(BCV), BC_J_N(BCV)
         do i = BC_I_W(BCV), BC_I_E(BCV)
            ijk = funijk(i,j,k)
            if(fluid_at(ijk)) then


               if(A_M(IJK,0,M) == -ONE .AND. B_M(IJK,M) == -T_x(IJK)) then
                  B_M(IJK,M) = -BC_T(BCV)
               else

                  lMass = BC_FLOW(BCV) * (VOL(IJK)/PS_VOLUME(BCV))

                  pSource = 0.0d0
                  do N=1,NMAX(M)

! Part of the thermal source: 0K --> T_BC
                     intCp_Tref = calc_ICpOR(BC_T(BCV), Thigh(N),  &
                        Tlow(N), Tcom(N), Ahigh(1,N), Alow(1,N)) * &
                        (GAS_CONST_cal / lMW(N))

! Part of the thermal source: 0K --> T_x(IJK)
                     intCp_Tijk = calc_ICpOR(T_x(IJK), Thigh(N),    &
                        Tlow(N), Tcom(N), Ahigh(1,N), Alow(1,N)) *  &
                        (GAS_CONST_cal / lMW(N))

! Total Thermal Source: T_BC --> T_x
                     pSource = pSource + lMass * BC_X(BCV,N) * &
                        (intCp_Tijk - intCp_Tref)

                  enddo

                  if(UNITS == 'SI') pSource = pSource * 4.183925d3

               endif


            endif

         enddo
         enddo
         enddo

      enddo BC_LP


      RETURN
      END SUBROUTINE POINT_SOURCE_ENERGY_EQ
