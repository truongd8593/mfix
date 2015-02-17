!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_V_s                                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative;                                              C
!           See source_v_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_DIF_V_S(A_M, B_M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
! maximum number of computational cells, number of solids phases
      USE param, only: dimension_3, dimension_m

! kinetic theories
      USE run, only: kt_type_enum
      USE run, only: ghd_2007
! run time flag for deferred correction
      USE run, only: def_cor
! run time flag to solve y momentum equation
      USE run, only: momentum_y_eq
! discretization scheme for indicated equation
      USE run, only: discretize

! number of solids phases
      USE physprop, only: mmax

! solids phase viscosity
      USE visc_s, only: mu_s

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Solids phase index
      INTEGER :: M
!-----------------------------------------------

      DO M = 1, MMAX
        IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
           (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

          IF  (MOMENTUM_Y_EQ(M)) THEN

! IF DEFERRED CORRECTION IS USED TO SOLVE V_S
             IF (DEF_COR) THEN
                CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER)
                IF (DISCRETIZE(4) > 1)CALL STORE_A_V_SDC (A_M(1,-3,M), M, B_M)
             ELSE

! NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR V_S
                IF (DISCRETIZE(4) == 0) THEN         ! 0 & 1 => FOUP
                   CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER)
                ELSE
                   CALL STORE_A_V_S1 (A_M(1,-3,M), M)
                ENDIF
             ENDIF

             CALL DIF_V_IS (MU_S(1,M), A_M, B_M, M, IER)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_V_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_V_s0                                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. FOUP                                         C
!           See source_v_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_V_S0(A_V_S, M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc
      USE fldvar
      USE output
      USE compar
      USE mflux
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_V_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_s(DIMENSION_3, -3:3, M:M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                 IJKNE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER :: IJMK, JM, IJKS
      INTEGER :: IJKM, KM, IJKB, IJKBN, IJPKM
! Face mass flux
      DOUBLE PRECISION :: Flux
! Diffusion parameter
      DOUBLE PRECISION :: D_f
! for cartesian grid:
      DOUBLE PRECISION :: AW,HW,VELW
!-----------------------------------------------

! Calculate convection-diffusion fluxes through each of the faces


!!$omp      parallel do                                                  &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,      &
!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, D_f,   &
!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,       &
!!$omp&             IJMK, JM, IJKS,     &
!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_N(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKN
            ELSE
               IJKC = IJK
            ENDIF
            JP = JP1(J)
            IJKE = EAST_OF(IJK)
            IJKNE = EAST_OF(IJKN)

! East face (i+1/2, j+1/2, k)
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IJK) * Flux_sE(IJK,M) + &
                       Theta_V_ne(IJK) * Flux_sE(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Ve_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                     ONEoDX_E_V(IJK)*AYZ_V(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                     ODX_E(I)*AYZ_V(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_V_S(IJK,E,M) = D_F
               A_V_S(IPJK,W,M) = D_F + Flux
            ELSE
               A_V_S(IJK,E,M) = D_F - Flux
               A_V_S(IPJK,W,M) = D_F
            ENDIF

! North face (i, j+1, k)
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJK) * Flux_sN(IJK,M) + &
                       Theta_Vn(IJK) * Flux_sN(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       alpha_Vn_c(IJK) ,AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKN,M)*ONEoDY_N_V(IJK)*AXZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJPK,M))
               D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_V_S(IJK,N,M) = D_F
               A_V_S(IJPK,S,M) = D_F + Flux
            ELSE
               A_V_S(IJK,N,M) = D_F - Flux
               A_V_S(IJPK,S,M) = D_F
            ENDIF

! Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJK) * Flux_sT(IJK,M) + &
                          Theta_V_st(IJK) * Flux_sT(IJPK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Vt_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                        OX(I)*ONEoDZ_T_V(IJK)*AXY_V(IJK)
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                        OX(I)*ODZ_T(K)*AXY_V(IJK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_V_S(IJK,T,M) = D_F
                  A_V_S(IJKP,B,M) = D_F + Flux
               ELSE
                  A_V_S(IJK,T,M) = D_F - Flux
                  A_V_S(IJKP,B,M) = D_F
               ENDIF
            ENDIF

! West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWN = NORTH_OF(IJKW)
               IMJPK = JP_OF(IMJK)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_se(IMJK) * Flux_sE(IMJK,M) +&
                          Theta_V_ne(IMJK) * Flux_sE(IMJPK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Ve_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)*&
                        ONEoDX_E_V(IMJK)*AYZ_V(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
                  D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)*&
                        ODX_E(IM)*AYZ_V(IMJK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_V_S(IJK,W,M) = D_F + Flux
               ELSE
                  A_V_S(IJK,W,M) = D_F
               ENDIF
            ENDIF

! South face (i, j, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               JM = JM1(J)
               IJKS = SOUTH_OF(IJK)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Vn_bar(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_Vn(IJMK) * Flux_sN(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          alpha_Vn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJKC,M)*ONEoDY_N_V(IJMK)*AXZ_V(IJMK)
               ELSE
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
                  D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_V_S(IJK,S,M) = D_F + Flux
               ELSE
                  A_V_S(IJK,S,M) = D_F
               ENDIF
            ENDIF

! Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  KM = KM1(K)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBN = NORTH_OF(IJKB)
                  IJPKM = JP_OF(IJKM)
                  IF(CUT_V_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_V_nt(IJKM) * Flux_sT(IJKM,M) + &
                          Theta_V_st(IJKM) * Flux_sT(IJPKM,M))
                     CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                             ALPHA_Vt_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                           OX(I)*ONEoDZ_T_V(IJKM)*AXY_V(IJKM)
                  ELSE
                     Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
                     D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                           OX(I)*ODZ_T(KM)*AXY_V(IJKM)
                  ENDIF
                  IF (Flux >= ZERO) THEN
                     A_V_S(IJK,B,M) = D_F + Flux
                  ELSE
                     A_V_S(IJK,B,M) = D_F
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_n)
      ENDDO   ! end do ijk


      RETURN
      END SUBROUTINE STORE_A_V_S0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_V_sdc                                           C
!  Purpose: To use deferred correction method to solve the v-momentum  C
!           equation. This method combines first order upwind and a    C
!           user specified higher order method                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date: 8-APR-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_V_SDC(A_V_S, M, B_M)
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE cutcell
      USE discretization, ONLY: fpfoi_of
      USE fldvar
      USE fun_avg
      USE function3
      USE functions
      USE geometry
      USE indices
      USE matrix
      USE mflux
      USE output
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE sendrecv
      USE sendrecv3
      USE toleranc
      USE visc_s
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_V_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_s(DIMENSION_3, -3:3, M:M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                 IJKNE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER :: IJMK, JM, IJKS
      INTEGER :: IJKM, KM, IJKB, IJKBN, IJPKM
      INTEGER :: IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER :: IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! indicator for shear
      INTEGER :: incr
! Deferred correction contribution from high order method
      DOUBLE PRECISION :: MOM_HO
! low order approximation
      DOUBLE PRECISION :: MOM_LO
! convection factor at the face
      DOUBLE PRECISION :: Flux
! deferred correction contributions from each face
      DOUBLE PRECISION :: EAST_DC
      DOUBLE PRECISION :: WEST_DC
      DOUBLE PRECISION :: NORTH_DC
      DOUBLE PRECISION :: SOUTH_DC
      DOUBLE PRECISION :: TOP_DC
      DOUBLE PRECISION :: BOTTOM_DC

! for cartersian grid:
        DOUBLE PRECISION :: AW,HW,VELW

! temporary use of global arrays:
! array1 (locally u)
! the x directional velocity
!      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)
! the y directional velocity
!      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww)
! the z directional velocity
!      DOUBLE PRECISION :: WW(DIMENSION_3)

!-----------------------------------------------

      call lock_tmp4_array
      call lock_tmp_array   ! locks array1, array2, array3 (locally u, v, ww)
      call lock_xsi_array


! Calculate convection factors
! ---------------------------------------------------------------->>>
! Send recv the third ghost layer
      IF (FPFOI) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = V_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKN = NORTH_OF(IJK)

! East face (i+1/2, j+1/2, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_V_se(IJK) * U_S(IJK,M) +&
                      Theta_V_ne(IJK) * U_S(IJPK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                    ALPHA_Ve_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J)
         ENDIF

! North face (i, j+1, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_Vn_bar(IJK) * V_S(IJK,M) + &
                      Theta_Vn(IJK) * V_S(IJPK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                    alpha_Vn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M))
         ENDIF

! Top face (i, j+1/2, k+1/2)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_V_nt(IJK) * W_S(IJK,M) + &
                          Theta_V_st(IJK) * W_S(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Vt_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE
            IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J)
         ENDIF
      ENDDO   ! end do ijk

! shear
      incr=2

      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, &
                     XSI_T,incr)


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>

!!$omp      parallel do         &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,      &
!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,  &
!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,       &
!!$omp&             IJMK, JM, IJKS,     &
!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM, &
!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC )
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_N(IJK)) THEN
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKN
            ELSE
               IJKC = IJK
            ENDIF
            JP = JP1(J)
            IJKE = EAST_OF(IJK)
            IJKNE = EAST_OF(IJKN)

! Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))


! DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j+1/2, k)
            IF(U(IJK) >= ZERO)THEN
               MOM_LO = V_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IPJK,M), V_S(IJK,M), &
                                   V_S(IMJK,M), V_S(IM_OF(IMJK),M))
            ELSE
                MOM_LO = V_S(IPJK,M)
                IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IPJK,M), &
                                    V_S(IP_OF(IPJK),M), TMP4(IPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IJK)*V_S(IPJK,M)+ &
                                      (1.0-XSI_E(IJK))*V_S(IJK,M)

            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IJK) * Flux_sE(IJK,M) +&
                       Theta_V_ne(IJK) * Flux_sE(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Ve_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
            ENDIF
            EAST_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
            IF(V(IJK) >= ZERO)THEN
               MOM_LO = V_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJPK,M), V_S(IJK,M), &
                                   V_S(IJMK,M), V_S(JM_OF(IJMK),M))
            ELSE
               MOM_LO = V_S(IJPK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJPK,M), &
                                   V_S(JP_OF(IJPK),M), TMP4(JPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJK)*V_S(IJPK,M)+ &
                                      (1.0-XSI_N(IJK))*V_S(IJK,M)

            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJK) * Flux_sN(IJK,M) + &
                       Theta_Vn(IJK) * Flux_sN(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       alpha_Vn_c(IJK) ,AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJPK,M))
            ENDIF
            NORTH_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = V_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJKP,M), V_S(IJK,M), &
                                      V_S(IJKM,M), V_S(KM_OF(IJKM),M))
               ELSE
                  MOM_LO = V_S(IJKP,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJKP,M), &
                                      V_S(KP_OF(IJKP),M), TMP4(KPPP4))
               ENDIF

               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJK)*V_S(IJKP,M)+ &
                                         (1.0-XSI_T(IJK))*V_S(IJK,M)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJK) * Flux_sT(IJK,M) + &
                          Theta_V_st(IJK) * Flux_sT(IJPK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Vt_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
               ENDIF
               TOP_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
               TOP_DC = ZERO
            ENDIF   ! end if do_k


! DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK)
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            IJKWN = NORTH_OF(IJKW)
            IMJPK = JP_OF(IMJK)
            IF(U(IMJK) >= ZERO)THEN
              MOM_LO = V_S(IMJK,M)
              IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IMJK,M), &
                                  V_S(IM_OF(IMJK),M), TMP4(IMMM4))
            ELSE
               MOM_LO = V_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IMJK,M), V_S(IJK,M), &
                                   V_S(IPJK,M), V_S(IP_OF(IPJK),M))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IMJK)*V_S(IJK,M)+ &
                                      (1.0-XSI_E(IMJK))*V_S(IMJK,M)

            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IMJK) * Flux_sE(IMJK,M) +&
                       Theta_V_ne(IMJK) * Flux_sE(IMJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Ve_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
            ENDIF
            WEST_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
            IJMK = JM_OF(IJK)
            JM = JM1(J)
            IJKS = SOUTH_OF(IJK)
            IF(V(IJMK) >= ZERO)THEN
               MOM_LO = V_S(IJMK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJMK,M), &
                                   V_S(JM_OF(IJMK),M), TMP4(JMMM4))
            ELSE
               MOM_LO = V_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJMK,M), V_S(IJK,M), &
                                   V_S(IJPK,M), V_S(JP_OF(IJPK),M))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJMK)*V_S(IJK,M)+ &
                                      (1.0-XSI_N(IJMK))*V_S(IJMK,M)

            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJMK) * Flux_sN(IJMK,M) + &
                       Theta_Vn(IJMK) * Flux_sN(IJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       alpha_Vn_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
            ENDIF
            SOUTH_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               IJKBN = NORTH_OF(IJKB)
               IJPKM = JP_OF(IJKM)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = V_S(IJKM,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJKM,M), &
                                      V_S(KM_OF(IJKM),M), TMP4(KMMM4))
               ELSE
                  MOM_LO = V_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(V_S(IJKM,M), V_S(IJK,M), &
                                      V_S(IJKP,M), V_S(KP_OF(IJKP),M))
               ENDIF

               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJKM)*V_S(IJK,M)+ &
                                         (1.0-XSI_T(IJKM))*V_S(IJKM,M)

               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJKM) * Flux_sT(IJKM,M) + &
                          Theta_V_st(IJKM) * Flux_sT(IJPKM,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Vt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
               ELSE
                  Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
               ENDIF
               BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
               BOTTOM_DC = ZERO
            ENDIF   ! end if do_k

! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC+&
                         BOTTOM_DC-TOP_DC

         ENDIF   ! end if flow_at_n
      ENDDO   ! end do ijk

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_V_SDC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_V_s1                                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. Higher order                                 C
!  See source_v_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_V_S1(A_V_S, M)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc
      USE fldvar
      USE output
      USE vshear
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE mflux
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_V_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_V_s(DIMENSION_3, -3:3, M:M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                 IJKNE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER :: IJMK, JM, IJKS
      INTEGER :: IJKM, KM, IJKB, IJKBN, IJPKM
! indicator for shear
      INTEGER incr
! Face mass flux
      DOUBLE PRECISION :: Flux
! Diffusion parameter
      DOUBLE PRECISION :: D_f
! for cartesian grid:
      DOUBLE PRECISION :: AW,HW,VELW

! temporary use of global arrays:
! array1 (locally u)
! the x directional velocity
!      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)
! the y directional velocity
!      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww)
! the z directional velocity
!      DOUBLE PRECISION :: WW(DIMENSION_3)

!-----------------------------------------------

      call lock_tmp_array
      call lock_xsi_array

! Calculate convection factors
! ---------------------------------------------------------------->>>
!!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKN = NORTH_OF(IJK)

! East face (i+1/2, j+1/2, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_V_se(IJK) * U_S(IJK,M) +&
                      Theta_V_ne(IJK) * U_S(IJPK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                    ALPHA_Ve_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J)
         ENDIF

! North face (i, j+1, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_Vn_bar(IJK) * V_S(IJK,M) + &
                      Theta_Vn(IJK) * V_S(IJPK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                    alpha_Vn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M))
         ENDIF

! Top face (i, j+1/2, k+1/2)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_V_nt(IJK) * W_S(IJK,M) + &
                          Theta_V_st(IJK) * W_S(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Vt_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE
            IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J)
         ENDIF
      ENDDO   ! end do ijk


! shear indicator
      incr=2

      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, &
                     XSI_T,incr)

! shear: update to true velocity
      IF (SHEAR) THEN
!!$omp      parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               V(IJK)=V(IJK)+VSH(IJK)
            ENDIF
         ENDDO
      ENDIF


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>

!!$omp      parallel do         &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,      &
!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,  &
!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,       &
!!$omp&             IJMK, JM, IJKS,     &
!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_N(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKN
            ELSE
               IJKC = IJK
            ENDIF
            JP = JP1(J)
            IJKE = EAST_OF(IJK)
            IJKNE = EAST_OF(IJKN)


! East face (i+1/2, j+1/2, k)
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IJK) * Flux_sE(IJK,M) +&
                       Theta_V_ne(IJK) * Flux_sE(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       ALPHA_Ve_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                     ONEoDX_E_V(IJK)*AYZ_V(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
               ODX_E(I)*AYZ_V(IJK)
            ENDIF
            A_V_S(IJK,E,M) = D_F - XSI_E(IJK)*Flux
            A_V_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*Flux


! North face (i, j+1, k)
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJK) * Flux_sN(IJK,M) + &
                       Theta_Vn(IJK) * Flux_sN(IJPK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                       alpha_Vn_c(IJK) ,AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKN,M)*ONEoDY_N_V(IJK)*AXZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJPK,M))
               D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK)
            ENDIF
            A_V_S(IJK,N,M) = D_F - XSI_N(IJK)*Flux
            A_V_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*Flux


! Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJK) * Flux_sT(IJK,M) + &
                          Theta_V_st(IJK) * Flux_sT(IJPK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Vt_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                        OX(I)*ONEoDZ_T_V(IJK)*AXY_V(IJK)
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                        OX(I)*ODZ_T(K)*AXY_V(IJK)
               ENDIF
               A_V_S(IJK,T,M) = D_F - XSI_T(IJK)*Flux
               A_V_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*Flux
            ENDIF


! West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWN = NORTH_OF(IJKW)
               IMJPK = JP_OF(IMJK)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_se(IMJK) * Flux_sE(IMJK,M) +&
                          Theta_V_ne(IMJK) * Flux_sE(IMJPK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          ALPHA_Ve_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)*&
                        ONEoDX_E_V(IMJK)*AYZ_V(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
                  D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)*&
                        ODX_E(IM)*AYZ_V(IMJK)
               ENDIF
               A_V_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF


! South face (i, j, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               JM = JM1(J)
               IJKS = SOUTH_OF(IJK)
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Vn_bar(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_Vn(IJMK) * Flux_sN(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                          alpha_Vn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJKC,M)*ONEoDY_N_V(IJMK)*AXZ_V(IJMK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
                  D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK)
               ENDIF
               A_V_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF


! Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  KM = KM1(K)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBN = NORTH_OF(IJKB)
                  IJPKM = JP_OF(IJKM)
                  IF(CUT_V_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_V_nt(IJKM) * Flux_sT(IJKM,M) + &
                             Theta_V_st(IJKM) * Flux_sT(IJPKM,M))
                     CALL GET_INTERPOLATION_TERMS_S(IJK,M,'V_MOMENTUM',&
                             ALPHA_Vt_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                           OX(I)*ONEoDZ_T_V(IJKM)*AXY_V(IJKM)
                  ELSE
                     Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
                     D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                           OX(I)*ODZ_T(KM)*AXY_V(IJKM)
                  ENDIF
                  A_V_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*Flux
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if flow_at_n
      ENDDO   ! end do ijk

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_V_S1
