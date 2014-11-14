!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_U_s                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative;                                              C
!           See source_u_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_DIF_U_S(A_M, B_M, IER)

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
! run time flag to solve x momentum equation
      USE run, only: momentum_x_eq
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

          IF (MOMENTUM_X_EQ(M)) THEN

! IF DEFERRED CORRECTION IS USED TO SOLVE U_S
             IF(DEF_COR)THEN
                CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER)
                IF (DISCRETIZE(3) > 1) CALL STORE_A_U_SDC (A_M(1,-3,M), M, B_M, IER)
             ELSE

! NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR U_S
                IF (DISCRETIZE(3) == 0) THEN         ! 0 & 1 => FOUP
                   CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER)
                ELSE
                   CALL STORE_A_U_S1 (A_M(1,-3,M), M, IER)
                ENDIF
             ENDIF

            CALL DIF_U_IS (MU_S(1,M), A_M, B_M, M, IER)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_U_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_s0                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. FOUP                                         C
!           See source_u_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_U_S0(A_U_S, M, IER)

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
! Septadiagonal matrix A_U_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_s(DIMENSION_3, -3:3, M:M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                 IJKNE, IJKP, IJKT, IJKTE
      INTEGER :: IMJK, IM, IJKW
      INTEGER :: IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER :: IJKM, KM, IPJKM, IJKB, IJKBE
! Face mass flux
      DOUBLE PRECISION :: Flux
! Diffusion parameter
      DOUBLE PRECISION :: D_f
! for cartesian grid:
      DOUBLE PRECISION :: AW,HW,VELW
!-----------------------------------------------

! Calculate convection-diffusion fluxes through each of the faces

!!$omp      parallel do         &
!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  D_f,      &
!!$omp&                    IMJK, IM, IJKW,      &
!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,        &
!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKE = EAST_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKE
            ELSE
               IJKC = IJK
            ENDIF
            IP = IP1(I)
            IJKN = NORTH_OF(IJK)
            IJKNE = EAST_OF(IJKN)

! East face (i+1, j, k)
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IJK) * Flux_sE(IJK,M) + &
                       Theta_Ue(IJK) * Flux_sE(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       alpha_Ue_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKE,M)*ONEoDX_E_U(IJK)*AYZ_U(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
               D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_U_S(IJK,E,M) = D_F
               A_U_S(IPJK,W,M) = D_F + Flux
            ELSE
               A_U_S(IJK,E,M) = D_F - Flux
               A_U_S(IPJK,W,M) = D_F
            ENDIF

! North face (i+1/2, j+1/2, k)
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJK) * Flux_sN(IJK,M) + &
                       Theta_U_ne(IJK) * Flux_sN(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Un_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)*&
                     ONEoDY_N_U(IJK)*AXZ_U(IJK)
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)*&
                     ODY_N(J)*AXZ_U(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_U_S(IJK,N,M) = D_F
               A_U_S(IJPK,S,M) = D_F + Flux
            ELSE
               A_U_S(IJK,N,M) = D_F - Flux
               A_U_S(IJPK,S,M) = D_F
            ENDIF

! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTE = EAST_OF(IJKT)

               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJK) * Flux_sT(IJK,M) + &
                          Theta_U_te(IJK) * Flux_sT(IPJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Ut_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)*&
                        OX_E(I)*ONEoDZ_T_U(IJK)*AXY_U(IJK)
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)*&
                        OX_E(I)*ODZ_T(K)*AXY_U(IJK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_U_S(IJK,T,M) = D_F
                  A_U_S(IJKP,B,M) = D_F + Flux
               ELSE
                  A_U_S(IJK,T,M) = D_F - Flux
                  A_U_S(IJKP,B,M) = D_F
               ENDIF
            ENDIF

! West face (i, j, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)

               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Ue_bar(IMJK) * Flux_sE(IMJK,M) + &
                          Theta_Ue(IMJK) * Flux_sE(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          alpha_Ue_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJKC,M)*ONEoDX_E_U(IMJK)*AYZ_U(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
                  D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_U_S(IJK,W,M) = D_F + Flux
               ELSE
                  A_U_S(IJK,W,M) = D_F
               ENDIF
            ENDIF

! South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               JM = JM1(J)
               IPJMK = IP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)

               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_nw(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_U_ne(IJMK) * Flux_sN(IPJMK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Un_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*&
                        ONEoDY_N_U(IJMK)*AXZ_U(IJMK)
               ELSE
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
                  D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*&
                        ODY_N(JM)*AXZ_U(IJMK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_U_S(IJK,S,M) = D_F + Flux
               ELSE
                  A_U_S(IJK,S,M) = D_F
               ENDIF
            ENDIF

! Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  KM = KM1(K)
                  IPJKM = IP_OF(IJKM)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBE = EAST_OF(IJKB)
                  IF(CUT_U_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_U_tw(IJKM) * Flux_sT(IJKM,M) + &
                             Theta_U_te(IJKM) * Flux_sT(IPJKM,M))
                     CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                             ALPHA_Ut_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*&
                           OX_E(I)*ONEoDZ_T_U(IJKM)*AXY_U(IJKM)
                  ELSE
                     Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
                     D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*&
                           OX_E(I)*ODZ_T(KM)*AXY_U(IJKM)
                  ENDIF
                  IF (Flux >= ZERO) THEN
                     A_U_S(IJK,B,M) = D_F + Flux
                  ELSE
                     A_U_S(IJK,B,M) = D_F
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_e)
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE STORE_A_U_S0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_sdc                                           C
!  Purpose: To use deferred correction method to solve the u-momentum  C
!           equation. This method combines first order upwind and a    C
!           user specified higher order method                         C
!                                                                      C
!  Author: C. GUENTHER                                Date: 8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_U_SDC(A_U_S, M, B_M, IER)

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
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE sendrecv
      USE sendrecv3
      USE mflux
      USE cutcell
      USE fun_avg
      USE function3
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_U_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_s(DIMENSION_3, -3:3, M:M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                 IJKNE, IJKP, IJKT, IJKTE
      INTEGER :: IMJK, IM, IJKW
      INTEGER :: IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER :: IJKM, KM, IPJKM, IJKB, IJKBE
      INTEGER :: IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER :: IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: D_f
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
! External functions
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
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
            TMP4(IJK4) = U_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         IP = IP1(I)
         IPJK = IP_OF(IJK)
         IJKE = EAST_OF(IJK)

! East face (i+1, j, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_Ue_bar(IJK) * U_S(IJK,M) + &
                      Theta_Ue(IJK) * U_S(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                    alpha_Ue_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP)
         ENDIF

! North face (i+1/2, j+1/2, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_U_nw(IJK) * V_S(IJK,M) + &
                      Theta_U_ne(IJK) * V_S(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                    ALPHA_Un_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I)
         ENDIF

! Top face (i+1/2, j, k+1/2)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_S(IJK,M) + &
                          Theta_U_te(IJK) * W_S(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Ut_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE
            IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I)
         ENDIF
      ENDDO   ! end do ijk


! shear
      incr=1

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N,&
                     XSI_T,incr)


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>

!!$omp      parallel do         &
!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  D_f,      &
!!$omp&                    IMJK, IM, IJKW,      &
!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,        &
!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_E(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKE
            ELSE
               IJKC = IJK
            ENDIF
            IP = IP1(I)
            IJKN = NORTH_OF(IJK)
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


! DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
            IF(U(IJK) >= ZERO)THEN
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IPJK,M), U_S(IJK,M), &
                                   U_S(IMJK,M), U_S(IM_OF(IMJK),M))
            ELSE
               MOM_LO = U_S(IPJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IPJK,M), &
                                   U_S(IP_OF(IPJK),M), TMP4(IPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IJK)*U_S(IPJK,M)+ &
                                       (1.0-XSI_E(IJK))*U_S(IJK,M)

            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IJK) * Flux_sE(IJK,M) + &
                       Theta_Ue(IJK) * Flux_sE(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       alpha_Ue_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
                Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
            ENDIF
            EAST_DC = Flux *(MOM_LO - MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
            IF(V(IJK) >= ZERO)THEN
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJPK,M), U_S(IJK,M), &
                                   U_S(IJMK,M), U_S(JM_OF(IJMK),M))
            ELSE
               MOM_LO = U_S(IJPK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJPK,M), &
                                   U_S(JP_OF(IJPK),M), TMP4(JPPP4))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJK)*U_S(IJPK,M)+ &
                                      (1.0-XSI_N(IJK))*U_S(IJK,M)

            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJK) * Flux_sN(IJK,M) + &
                       Theta_U_ne(IJK) * Flux_sN(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Un_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
            ENDIF
            NORTH_DC = Flux * (MOM_LO - MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTE = EAST_OF(IJKT)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = U_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJKP,M), U_S(IJK,M), &
                                      U_S(IJKM,M), U_S(KM_OF(IJKM),M))
               ELSE
                  MOM_LO = U_S(IJKP,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKP,M), &
                                      U_S(KP_OF(IJKP),M), TMP4(KPPP4))
               ENDIF
               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJK)*U_S(IJKP,M)+ &
                                         (1.0-XSI_T(IJK))*U_S(IJK,M)

               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJK) * Flux_sT(IJK,M) + &
                          Theta_U_te(IJK) * Flux_sT(IPJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Ut_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
               ENDIF
               TOP_DC = Flux *(MOM_LO - MOM_HO)
            ELSE
               TOP_DC = ZERO
            ENDIF   ! end if do_k


! DEFERRED CORRECTION CONTRIBUTION AT THE West face (i, j, k)
            IMJK = IM_OF(IJK)
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            IF(U(IMJK) >= ZERO)THEN
               MOM_LO = U_S(IMJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IMJK,M), &
                                   U_S(IM_OF(IMJK),M), TMP4(IMMM4))
            ELSE
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IMJK,M), U_S(IJK,M), &
                          U_S(IPJK,M), U_S(IP_OF(IPJK),M))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_E(IMJK)*U_S(IJK,M)+ &
                                      (1.0-XSI_E(IMJK))*U_S(IMJK,M)

            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IMJK) * Flux_sE(IMJK,M) + &
                       Theta_Ue(IMJK) * Flux_sE(IJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       alpha_Ue_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
            ENDIF
            WEST_DC = Flux * (MOM_LO - MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK)
            JM = JM1(J)
            IPJMK = IP_OF(IJMK)
            IJKS = SOUTH_OF(IJK)
            IJKSE = EAST_OF(IJKS)
            IF(V(IJMK) >= ZERO)THEN
               MOM_LO = U_S(IJMK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJMK,M), &
                                   U_S(JM_OF(IJMK),M), TMP4(JMMM4))
            ELSE
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJMK,M), U_S(IJK,M), &
                                   U_S(IJPK,M), U_S(JP_OF(IJPK),M))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJMK)*U_S(IJK,M)+ &
                                      (1.0-XSI_N(IJMK))*U_S(IJMK,M)

            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJMK) * Flux_sN(IJMK,M) + &
                       Theta_U_ne(IJMK) * Flux_sN(IPJMK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Un_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
            ENDIF
            SOUTH_DC = Flux * (MOM_LO - MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               KM = KM1(K)
               IPJKM = IP_OF(IJKM)
               IJKB = BOTTOM_OF(IJK)
               IJKBE = EAST_OF(IJKB)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = U_S(IJKM,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKM,M), &
                                      U_S(KM_OF(IJKM),M), TMP4(KMMM4))
               ELSE
                  MOM_LO = U_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJKM,M), U_S(IJK,M), &
                                      U_S(IJKP,M), U_S(KP_OF(IJKP),M))
               ENDIF
               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJKM)*U_S(IJK,M)+ &
                                         (1.0-XSI_T(IJKM))*U_S(IJKM,M)

               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJKM) * Flux_sT(IJKM,M) + &
                          Theta_U_te(IJKM) * Flux_sT(IPJKM,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Ut_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
               ELSE
                  Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
               ENDIF
               BOTTOM_DC = Flux * (MOM_LO - MOM_HO)
            ELSE
               BOTTOM_DC = ZERO
            ENDIF   ! end if do_k

! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC+&
                         BOTTOM_DC-TOP_DC

         ENDIF   ! end if flow_at_e
      ENDDO   ! end do ijk

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_U_SDC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_s1                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. Higher order                                 C
!  See source_u_s                                                      C
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

      SUBROUTINE STORE_A_U_S1(A_U_S, M, IER)

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
! Septadiagonal matrix A_U_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_s(DIMENSION_3, -3:3, M:M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                 IJKNE, IJKP, IJKT, IJKTE
      INTEGER :: IMJK, IM, IJKW
      INTEGER :: IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER :: IJKM, KM, IPJKM, IJKB, IJKBE
! indicator for shear
      INTEGER :: incr
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
!!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         IP = IP1(I)
         IPJK = IP_OF(IJK)
         IJKE = EAST_OF(IJK)

! East face (i+1, j, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_Ue_bar(IJK) * U_S(IJK,M) + &
                      Theta_Ue(IJK) * U_S(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                    alpha_Ue_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP)
         ENDIF

! North face (i+1/2, j+1/2, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_U_nw(IJK) * V_S(IJK,M) + &
                      Theta_U_ne(IJK) * V_S(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                    ALPHA_Un_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I)
         ENDIF

! Top face (i+1/2, j, k+1/2)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_S(IJK,M) + &
                          Theta_U_te(IJK) * W_S(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Ut_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE
            IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I)
         ENDIF
      ENDDO   ! end do ijk


! shear indicator
      incr=1

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N, &
                     XSI_T,incr)

! shear: update to true velocity
      IF (SHEAR) THEN
!!$omp  parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               V(IJK)=V(IJK)+VSHE(IJK)
            ENDIF
         ENDDO
      ENDIF


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>
!!$omp      parallel do         &
!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,     &
!!$omp&                    IMJK, IM, IJKW,      &
!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,        &
!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_E(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKE = EAST_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKE
            ELSE
               IJKC = IJK
            ENDIF
            IP = IP1(I)
            IJKN = NORTH_OF(IJK)
            IJKNE = EAST_OF(IJKN)


! East face (i+1, j, k)
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IJK) * Flux_sE(IJK,M) + &
                       Theta_Ue(IJK) * Flux_sE(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       alpha_Ue_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKE,M)*ONEoDX_E_U(IJK)*AYZ_U(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
               D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK)
            ENDIF
            A_U_S(IJK,E,M) = D_F - XSI_E(IJK) * Flux
            A_U_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK)) * Flux


! North face (i+1/2, j+1/2, k)
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJK) * Flux_sN(IJK,M) + &
                       Theta_U_ne(IJK) * Flux_sN(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                       ALPHA_Un_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)*&
                     ONEoDY_N_U(IJK)*AXZ_U(IJK)
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)*&
                     ODY_N(J)*AXZ_U(IJK)
            ENDIF
            A_U_S(IJK,N,M) = D_F - XSI_N(IJK) * Flux
            A_U_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK)) * Flux


! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTE = EAST_OF(IJKT)
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJK) * Flux_sT(IJK,M) + &
                          Theta_U_te(IJK) * Flux_sT(IPJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Ut_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)*&
                        OX_E(I)*ONEoDZ_T_U(IJK)*AXY_U(IJK)
               ELSE
                  Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)*&
                        OX_E(I)*ODZ_T(K)*AXY_U(IJK)
               ENDIF
               A_U_S(IJK,T,M) = D_F - XSI_T(IJK) * Flux
               A_U_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK)) * Flux
            ENDIF


! West face (i, j, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Ue_bar(IMJK) * Flux_sE(IMJK,M) + &
                          Theta_Ue(IMJK) * Flux_sE(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          alpha_Ue_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJKC,M)*ONEoDX_E_U(IMJK)*AYZ_U(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
                  D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK)
               ENDIF
               A_U_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK)) * Flux
            ENDIF


! South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               JM = JM1(J)
               IPJMK = IP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_nw(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_U_ne(IJMK) * Flux_sN(IPJMK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                          ALPHA_Un_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*&
                        ONEoDY_N_U(IJMK)*AXZ_U(IJMK)
               ELSE
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
                  D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*&
                        ODY_N(JM)*AXZ_U(IJMK)
               ENDIF
               A_U_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK)) * Flux
            ENDIF


! Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  KM = KM1(K)
                  IPJKM = IP_OF(IJKM)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBE = EAST_OF(IJKB)
                  IF(CUT_U_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_U_tw(IJKM) * Flux_sT(IJKM,M) + &
                             Theta_U_te(IJKM) * Flux_sT(IPJKM,M))
                     CALL GET_INTERPOLATION_TERMS_S(IJK,M,'U_MOMENTUM',&
                             ALPHA_Ut_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*&
                           OX_E(I)*ONEoDZ_T_U(IJKM)*AXY_U(IJKM)
                  ELSE
                     Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
                     D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                   AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*&
                           OX_E(I)*ODZ_T(KM)*AXY_U(IJKM)
                  ENDIF
                  A_U_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM)) * Flux
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if flow_at_e
      ENDDO   ! end do ijk

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_U_S1
