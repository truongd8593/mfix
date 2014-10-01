!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_W_s                                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative;                                              C
!           See source_w_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^CC

      SUBROUTINE CONV_DIF_W_S(A_M, B_M, IER)

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
! run time flag to solve z momentum equation
      USE run, only: momentum_z_eq
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
!------------------------------------------------


      DO M = 1, MMAX
        IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
           (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

          IF (MOMENTUM_Z_EQ(M)) THEN

! IF DEFERRED CORRECTION IS TO BE USED TO SOLVE W_S
             IF (DEF_COR) THEN
                CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER)
                IF (DISCRETIZE(5) > 1)CALL STORE_A_W_SDC (A_M(1,-3,M), M, B_M, IER)
             ELSE

! NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR W_S
                IF (DISCRETIZE(5) == 0) THEN         ! 0 & 1 => FOUP
                   CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER)
                ELSE
                   CALL STORE_A_W_S1 (A_M(1,-3,M), M, IER)
                ENDIF
             ENDIF

            CALL DIF_W_IS (MU_S(1,M), A_M, B_M, M, IER)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_W_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_s0                                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. FOUP                                         C
!           See source_w_s                                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-Jun-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_W_S0(A_W_S, M, IER)

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
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_W_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_s(DIMENSION_3, -3:3, M:M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                 IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER :: IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER :: IJKM, KM, IJKB
! Face mass flux
      DOUBLE PRECISION :: Flux
! Diffusion parameter
      DOUBLE PRECISION :: D_f
! for cartesian grid:
      DOUBLE PRECISION :: AW,HW,VELW
!-----------------------------------------------
! Include statment functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

! Calculate convection-diffusion fluxes through each of the faces

!!$omp      parallel do         &
!!$omp&     private(IJK,  I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, &
!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, D_f,        &
!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,       &
!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,       &
!!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_T(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKT
            ELSE
               IJKC = IJK
            ENDIF
            KP = KP1(K)
            IJKE = EAST_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKTE = EAST_OF(IJKT)

! East face (i+1/2, j, k+1/2)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_sE(IJK,M) + &
                       Theta_W_te(IJK) * Flux_sE(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                     ONEoDX_E_W(IJK)*AYZ_W(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                     ODX_E(I)*AYZ_W(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_W_S(IJK,E,M) = D_F
               A_W_S(IPJK,W,M) = D_F + Flux
            ELSE
               A_W_S(IJK,E,M) = D_F - Flux
               A_W_S(IPJK,W,M) = D_F
            ENDIF

! North face (i, j+1/2, k+1/2)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_sN(IJK,M) + &
                       Theta_W_tn(IJK) * Flux_sN(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)*&
                     ONEoDY_N_W(IJK)*AXZ_W(IJK)
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)*&
                     ODY_N(J)*AXZ_W(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_W_S(IJK,N,M) = D_F
               A_W_S(IJPK,S,M) = D_F + Flux
            ELSE
               A_W_S(IJK,N,M) = D_F - Flux
               A_W_S(IJPK,S,M) = D_F
            ENDIF

! Top face (i, j, k+1)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_sT(IJK,M) + &
                       Theta_Wt(IJK) * Flux_sT(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKT,M)*ONEoDZ_T_W(IJK)*AXY_W(IJK)
            ELSE
               Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
               D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK)
            ENDIF
            IF (Flux >= ZERO) THEN
               A_W_S(IJK,T,M) = D_F
               A_W_S(IJKP,B,M) = D_F + Flux
            ELSE
               A_W_S(IJK,T,M) = D_F - Flux
               A_W_S(IJKP,B,M) = D_F
            ENDIF

! West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IMJKP = KP_OF(IMJK)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_be(IMJK) * Flux_sE(IMJK,M) + &
                          Theta_W_te(IMJK) * Flux_sE(IMJKP,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          ALPHA_We_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)*&
                        ONEoDX_E_W(IMJK)*AYZ_W(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
                  D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)*&
                        ODX_E(IM)*AYZ_W(IMJK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_W_S(IJK,W,M) = D_F + Flux
               ELSE
                  A_W_S(IJK,W,M) = D_F
               ENDIF
            ENDIF

! South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               JM = JM1(J)
               IJMKP = KP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_bn(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_W_tn(IJMK) * Flux_sN(IJMKP,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          ALPHA_Wn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                        ONEoDY_N_W(IJMK)*AXZ_W(IJMK)
               ELSE
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
                  D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                        ODY_N(JM)*AXZ_W(IJMK)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_W_S(IJK,S,M) = D_F + Flux
               ELSE
                  A_W_S(IJK,S,M) = D_F
               ENDIF
            ENDIF

! Bottom face (i, j, k)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Wt_bar(IJKM) * Flux_sT(IJKM,M) + &
                          Theta_Wt(IJKM) * Flux_sT(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          alpha_Wt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJK,M)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM)
               ELSE
                  Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
                  D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM)
               ENDIF
               IF (Flux >= ZERO) THEN
                  A_W_S(IJK,B,M) = D_F + Flux
               ELSE
                  A_W_S(IJK,B,M) = D_F
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_t)
      ENDDO   ! end do ijk


      RETURN
      END SUBROUTINE STORE_A_W_S0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_sdc                                           C
!  Purpose: To use deferred correction method to solve the w-momentum  C
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

      SUBROUTINE STORE_A_W_SDC(A_W_S, M, B_M, IER)

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
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE sendrecv
      USE sendrecv3
      USE mflux
      USE cutcell
      USE fun_avg
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_W_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_s(DIMENSION_3, -3:3, M:M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                 IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER :: IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER :: IJKM, KM, IJKB
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
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'function3.inc'
!-----------------------------------------------

      call lock_tmp4_array
      call lock_tmp_array   ! locks array1, array2, array3 (locally u, v, ww)
      call lock_xsi_array

! Calculate convection factors
! ---------------------------------------------------------------->>>
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = W_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF


!!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         IJKP = KP_OF(IJK)

! East face (i+1/2, j, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_W_be(IJK) * U_S(IJK,M) + &
                      Theta_W_te(IJK) * U_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    ALPHA_We_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K)
         ENDIF

! North face (i, j+1/2, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_W_bn(IJK) * V_S(IJK,M) + &
                      Theta_W_tn(IJK) * V_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    ALPHA_Wn_c(IJK) ,AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K)
         ENDIF

! Top face (i, j, k+1)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            WW(IJK) = (Theta_Wt_bar(IJK) * W_S(IJK,M) + &
                       Theta_Wt(IJK) * W_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    alpha_Wt_c(IJK),AW,HW,VELW)
            WW(IJK) = WW(IJK) * AW
         ELSE
            WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M))
         ENDIF
      ENDDO

! shear
      incr=0

      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N,&
                     XSI_T,incr)


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>

!!!$omp      parallel do        &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f, &
!!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,      &
!!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,      &
!!!$omp&             IJKM, KM, IJKB, &
!!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_T(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKT
            ELSE
               IJKC = IJK
            ENDIF
            KP = KP1(K)
            IJKE = EAST_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKTE = EAST_OF(IJKT)

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


! DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k+1/2)
            IF(U(IJK) >= ZERO)THEN
               MOM_LO = W_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_S(IPJK,M), W_S(IJK,M), &
                                   W_S(IMJK,M), W_S(IM_OF(IMJK),M))
            ELSE
               MOM_LO = W_S(IPJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IPJK,M), &
                                   W_S(IP_OF(IPJK),M), TMP4(IPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IJK)*W_S(IPJK,M)+ &
                                      (1.0-XSI_E(IJK))*W_S(IJK,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_sE(IJK,M) + &
                       Theta_W_te(IJK) * Flux_sE(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
            ENDIF
            EAST_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
            IF(V(IJK) >= ZERO)THEN
               MOM_LO = W_S(IJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJPK,M), W_S(IJK,M), &
                                  W_S(IJMK,M), W_S(JM_OF(IJMK),M))
            ELSE
               MOM_LO = W_S(IJPK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJPK,M), &
                                  W_S(JP_OF(IJPK),M), TMP4(JPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJK)*W_S(IJPK,M)+ &
                                      (1.0-XSI_N(IJK))*W_S(IJK,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_sN(IJK,M) + &
                       Theta_W_tn(IJK) * Flux_sN(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
            ENDIF
            NORTH_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
            IF(WW(IJK) >= ZERO)THEN
               MOM_LO = W_S(IJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJKP,M), W_S(IJK,M), &
                                  W_S(IJKM,M), W_S(KM_OF(IJKM),M))
            ELSE
               MOM_LO = W_S(IJKP,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJKP,M), &
                                  W_S(KP_OF(IJKP),M), TMP4(KPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_T(IJK)*W_S(IJKP,M)+ &
                                      (1.0-XSI_T(IJK))*W_S(IJK,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_sT(IJK,M) + &
                       Theta_Wt(IJK) * Flux_sT(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
            ENDIF
            TOP_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK)
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            IJKWT = TOP_OF(IJKW)
            IMJKP = KP_OF(IMJK)
            IF(U(IMJK) >= ZERO)THEN
               MOM_LO = W_S(IMJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IMJK,M), &
                                  W_S(IM_OF(IMJK),M), TMP4(IMMM4))
            ELSE
               MOM_LO = W_S(IJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IMJK,M), W_S(IJK,M), &
                                  W_S(IPJK,M), W_S(IP_OF(IPJK),M))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IMJK)*W_S(IJK,M)+ &
                                      (1.0-XSI_E(IMJK))*W_S(IMJK,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IMJK) * Flux_sE(IMJK,M) + &
                       Theta_W_te(IMJK) * Flux_sE(IMJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_We_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
            ENDIF
            WEST_DC = Flux*(MOM_LO-MOM_HO)


! CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK)
            JM = JM1(J)
            IJMKP = KP_OF(IJMK)
            IJKS = SOUTH_OF(IJK)
            IJKST = TOP_OF(IJKS)
            IF(V(IJMK) >= ZERO)THEN
               MOM_LO = W_S(IJMK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJMK,M), &
                                  W_S(JM_OF(IJMK),M), TMP4(JMMM4))
            ELSE
               MOM_LO = W_S(IJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJMK,M), W_S(IJK,M), &
                                  W_S(IJPK,M), W_S(JP_OF(IJPK),M))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJMK)*W_S(IJK,M)+ &
                                      (1.0-XSI_N(IJMK))*W_S(IJMK,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJMK) * Flux_sN(IJMK,M) + &
                       Theta_W_tn(IJMK) * Flux_sN(IJMKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_Wn_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
            ENDIF
            SOUTH_DC = Flux*(MOM_LO-MOM_HO)


! DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
            IJKM = KM_OF(IJK)
            KM = KM1(K)
            IJKB = BOTTOM_OF(IJK)
            IF(WW(IJK) >= ZERO)THEN
               MOM_LO = W_S(IJKM,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJKM,M), &
                                  W_S(KM_OF(IJKM),M), TMP4(KMMM4))
            ELSE
               MOM_LO = W_S(IJK,M)
               IF(FPFOI) MOM_HO = FPFOI_OF(W_S(IJKM,M), W_S(IJK,M), &
                                  W_S(IJKP,M), W_S(KP_OF(IJKP),M))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_T(IJKM)*W_S(IJK,M)+ &
                                      (1.0-XSI_T(IJKM))*W_S(IJKM,M)

            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJKM) * Flux_sT(IJKM,M) + &
                       Theta_Wt(IJKM) * Flux_sT(IJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       alpha_Wt_c(IJKM),AW,HW,VELW)
               Flux = Flux * AW
            ELSE
               Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
            ENDIF
            BOTTOM_DC = Flux*(MOM_LO-MOM_HO)


! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC+&
                         BOTTOM_DC-TOP_DC

         ENDIF   ! end if flow_at_t
      ENDDO   ! end do ijk

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_SDC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_s1                                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!           The off-diagonal coefficients calculated here must be      C
!           positive. The center coefficient and the source vector     C
!           are negative. Higher order                                 C
!  See source_w_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE STORE_A_W_S1(A_W_S, M, IER)

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
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE mflux
      USE cutcell
      USE fun_avg
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_W_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_s(DIMENSION_3, -3:3, M:M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                 IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER :: IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER :: IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER :: IJKM, KM, IJKB
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
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------'

      call lock_tmp_array
      call lock_xsi_array

! Calculate convection factors
! ---------------------------------------------------------------->>>
!!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         IJKP = KP_OF(IJK)

! East face (i+1/2, j, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_W_be(IJK) * U_S(IJK,M) + &
                      Theta_W_te(IJK) * U_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    ALPHA_We_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE
            U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K)
         ENDIF

! North face (i, j+1/2, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_W_bn(IJK) * V_S(IJK,M) + &
                      Theta_W_tn(IJK) * V_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    ALPHA_Wn_c(IJK) ,AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE
            V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K)
         ENDIF

! Top face (i, j, k+1)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            WW(IJK) = (Theta_Wt_bar(IJK) * W_S(IJK,M) + &
                       Theta_Wt(IJK) * W_S(IJKP,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                    alpha_Wt_c(IJK),AW,HW,VELW)
            WW(IJK) = WW(IJK) * AW
         ELSE
            WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M))
         ENDIF
      ENDDO   ! end do ijk

! shear indicator
      incr=0

      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N,&
                     XSI_T,incr)


! shear: update to true velocity
      IF (SHEAR) THEN
!!!$omp      parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               V(IJK)=V(IJK)+VSH(IJK)
            ENDIF
         ENDDO
      ENDIF


! Calculate convection-diffusion fluxes through each of the faces
! ---------------------------------------------------------------->>>
!!$omp      parallel do         &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,      &
!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,  &
!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,       &
!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,       &
!!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3
         IF (FLOW_AT_T(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IF (WALL_AT(IJK)) THEN
               IJKC = IJKT
            ELSE
               IJKC = IJK
            ENDIF
            KP = KP1(K)
            IJKE = EAST_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKTE = EAST_OF(IJKT)


! East face (i+1/2, j, k+1/2)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_sE(IJK,M) + &
                       Theta_W_te(IJK) * Flux_sE(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                     ONEoDX_E_W(IJK)*AYZ_W(IJK)
            ELSE
               Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                             AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                     ODX_E(I)*AYZ_W(IJK)
            ENDIF
            A_W_S(IJK,E,M) = D_F - XSI_E(IJK)*Flux
            A_W_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*Flux


! North face (i, j+1/2, k+1/2)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_sN(IJK,M) + &
                       Theta_W_tn(IJK) * Flux_sN(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)*&
                     ONEoDY_N_W(IJK)*AXZ_W(IJK)
            ELSE
               Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                             AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)*&
                     ODY_N(J)*AXZ_W(IJK)
            ENDIF
            A_W_S(IJK,N,M) = D_F - XSI_N(IJK)*Flux
            A_W_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*Flux


! Top face (i, j, k+1)
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_sT(IJK,M) + &
                       Theta_Wt(IJK) * Flux_sT(IJKP,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                       alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_S(IJKT,M)*ONEoDZ_T_W(IJK)*AXY_W(IJK)
            ELSE
               Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
               D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK)
            ENDIF
            A_W_S(IJK,T,M) = D_F - XSI_T(IJK)*Flux
            A_W_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*Flux


! West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IMJKP = KP_OF(IMJK)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_be(IMJK) * Flux_sE(IMJK,M) + &
                          Theta_W_te(IMJK) * Flux_sE(IMJKP,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          ALPHA_We_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)*&
                        ONEoDX_E_W(IMJK)*AYZ_W(IMJK)
               ELSE
                  Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
                  D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                                AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)*&
                        ODX_E(IM)*AYZ_W(IMJK)
               ENDIF
               A_W_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF


! South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               JM = JM1(J)
               IJMKP = KP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_bn(IJMK) * Flux_sN(IJMK,M) + &
                          Theta_W_tn(IJMK) * Flux_sN(IJMKP,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          ALPHA_Wn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                        ONEoDY_N_W(IJMK)*AXZ_W(IJMK)
               ELSE
                  Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
                  D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                                AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                        ODY_N(JM)*AXZ_W(IJMK)
               ENDIF
               A_W_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF


! Bottom face (i, j, k)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Wt_bar(IJKM) * Flux_sT(IJKM,M) + &
                          Theta_Wt(IJKM) * Flux_sT(IJK,M))
                  CALL GET_INTERPOLATION_TERMS_S(IJK,M,'W_MOMENTUM',&
                          alpha_Wt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_S(IJK,M)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM)
               ELSE
                  Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
                  D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM)
               ENDIF
               A_W_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*Flux
            ENDIF

         ENDIF   ! end if flow_at_t
      ENDDO   ! end do ijk


      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_S1
