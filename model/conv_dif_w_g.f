!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_W_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
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
      SUBROUTINE CONV_DIF_W_G(A_M, B_M, IER)
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
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE compar

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!
      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       IF DEFERRED CORRECTION IS USED TO SOLVE W_G
      IF (DEF_COR) THEN
        CALL STORE_A_W_G0 (A_M(1,-3,0), IER)
        IF (DISCRETIZE(5) > 1)CALL STORE_A_W_GDC (A_M(1,-3,0), B_M(1,0), IER)
      ELSE
!
        IF (DISCRETIZE(5) == 0) THEN               ! 0 & 1 => FOUP
          CALL STORE_A_W_G0 (A_M(1,-3,0), IER)
        ELSE
          CALL STORE_A_W_G1 (A_M(1,-3,0), IER)
        ENDIF
      ENDIF
!
      CALL DIF_W_IS (MU_GT, A_M, B_M, 0, IER)
!

      RETURN
      END SUBROUTINE CONV_DIF_W_G
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_g0(A_W_g, IER)                               C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
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
      SUBROUTINE STORE_A_W_G0(A_W_G, IER)
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
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc
      USE physprop
      USE fldvar
      USE output
      USE compar
      USE mflux
      USE fun_avg
      USE functions
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                       IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER          IJKM, KM, IJKB
!
!                      Solids phase
      INTEGER          M
!
!                      Face mass flux
      DOUBLE PRECISION Flux
!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_W_g
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: AW,HW,VELW
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0

!!!$omp      parallel do                                               &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK, D_f,    &
!!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!!!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_T(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK )
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
!
!           East face (i+1/2, j, k+1/2)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_gE(IJK) + Theta_W_te(IJK) * Flux_gE(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
                  ),MU_GT(IJKTE),I),K)*ONEoDX_E_W(IJK)*AYZ_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
                  ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN
               A_W_G(IJK,E) = D_F
               A_W_G(IPJK,W) = D_F + Flux
            ELSE
               A_W_G(IJK,E) = D_F - Flux
               A_W_G(IPJK,W) = D_F
            ENDIF
!
!           North face (i, j+1/2, k+1/2)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_gN(IJK) + Theta_W_tn(IJK) * Flux_gN(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
                  ),MU_GT(IJKTN),J),K)*ONEoDY_N_W(IJK)*AXZ_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
                  ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN
               A_W_G(IJK,N) = D_F
               A_W_G(IJPK,S) = D_F + Flux
            ELSE
               A_W_G(IJK,N) = D_F - Flux
               A_W_G(IJPK,S) = D_F
            ENDIF
!
!           Top face (i, j, k+1)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_gT(IJK) + Theta_Wt(IJK) * Flux_gT(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_GT(IJKT)*ONEoDZ_T_W(IJK)*AXY_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
               D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN
               A_W_G(IJK,T) = D_F
               A_W_G(IJKP,B) = D_F + Flux
            ELSE
               A_W_G(IJK,T) = D_F - Flux
               A_W_G(IJKP,B) = D_F
            ENDIF
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IMJKP = KP_OF(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_be(IMJK) * Flux_gE(IMJK) + Theta_W_te(IMJK) * Flux_gE(IMJKP))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                      IJKWT),MU_GT(IJKT),IM),K)*ONEoDX_E_W(IMJK)*AYZ_W(IMJK)
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
                   D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                   IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_W_G(IJK,W) = D_F + Flux
               ELSE
                  A_W_G(IJK,W) = D_F
               ENDIF
            ENDIF
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               JM = JM1(J)
               IJMKP = KP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_bn(IJMK) * Flux_gN(IJMK) + Theta_W_tn(IJMK) * Flux_gN(IJMKP))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKST),MU_GT(IJKT),JM),K)*ONEoDY_N_W(IJMK)*AXZ_W(IJMK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
                  D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_W_G(IJK,S) = D_F + Flux
               ELSE
                  A_W_G(IJK,S) = D_F
               ENDIF
            ENDIF
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Wt_bar(IJKM) * Flux_gT(IJKM) + Theta_Wt(IJKM) * Flux_gT(IJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_GT(IJK)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
                  D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_W_G(IJK,B) = D_F + Flux
               ELSE
                  A_W_G(IJK,B) = D_F
               ENDIF
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE STORE_A_W_G0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_gdc(A_W_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE W-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date:8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
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
      SUBROUTINE STORE_A_W_GDC(A_W_G, B_M, IER)
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
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc
      USE physprop
      USE fldvar
      USE output
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE sendrecv
      USE sendrecv3
      USE mflux
      USE fun_avg
      USE function3
      USE functions
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                       IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER          IJKM, KM, IJKB
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
!
! loezos
        INTEGER incr
! loezos

!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_W_g
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!       FACE VELOCITY
        DOUBLE PRECISION V_F
!
!       DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
        DOUBLE PRECISION MOM_HO
!
!       LOW ORDER APPROXIMATION
        DOUBLE PRECISION MOM_LO
!
!       CONVECTION FACTOR AT THE FACE
        DOUBLE PRECISION Flux
!
!       DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
        DOUBLE PRECISION        EAST_DC
        DOUBLE PRECISION        WEST_DC
        DOUBLE PRECISION        NORTH_DC
        DOUBLE PRECISION        SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: AW,HW,VELW
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!
!-----------------------------------------------
!
!---------------------------------------------------------------
!       EXTERNAL FUNCTIONS
!---------------------------------------------------------------
        DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      call lock_tmp4_array
      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = W_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF


!!!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         IJKP = KP_OF(IJK)
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           East face (i+1/2, j, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_W_be(IJK) * U_G(IJK) + Theta_W_te(IJK) * U_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE   ! Original terms
            U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K)
         ENDIF
!
!
!           North face (i, j+1/2, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_W_bn(IJK) * V_G(IJK) + Theta_W_tn(IJK) * V_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE   ! Original terms
            V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K)
         ENDIF
!
!
!           Top face (i, j, k+1)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            WW(IJK) = (Theta_Wt_bar(IJK) * W_G(IJK) + Theta_Wt(IJK) * W_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJK),AW,HW,VELW)
            WW(IJK) = WW(IJK) * AW
         ELSE   ! Original terms
            WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP))
         ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO

! loezos
        incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
         XSI_T,incr)
!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!!!$omp      parallel do                                               &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!!!$omp&             IJKM, KM, IJKB, &
!!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
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
            IJKP = KP_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKTE = EAST_OF(IJKT)
!
!           Third Ghost layer information
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
!
!

!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k+1/2)
!
            IF(U(IJK) >= ZERO)THEN
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IPJK), W_G(IJK), &
                            W_G(IMJK), W_G(IM_OF(IMJK)))
            ELSE
              MOM_LO = W_G(IPJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IPJK), &
                            W_G(IP_OF(IPJK)), TMP4(IPPP4))
            ENDIF
            IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_E(IJK)*W_G(IPJK)+ &
                               (1.0-XSI_E(IJK))*W_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_gE(IJK) + Theta_W_te(IJK) * Flux_gE(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
!
            IF(V(IJK) >= ZERO)THEN
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJPK), W_G(IJK), &
                            W_G(IJMK), W_G(JM_OF(IJMK)))
            ELSE
              MOM_LO = W_G(IJPK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJPK), &
                            W_G(JP_OF(IJPK)), TMP4(JPPP4))
            ENDIF
            IF (.NOT. FPFOI ) &
                       MOM_HO = XSI_N(IJK)*W_G(IJPK)+ &
                                (1.0-XSI_N(IJK))*W_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_gN(IJK) + Theta_W_tn(IJK) * Flux_gN(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
!
            IF(WW(IJK) >= ZERO)THEN
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJKP), W_G(IJK), &
                            W_G(IJKM), W_G(KM_OF(IJKM)))
            ELSE
              MOM_LO = W_G(IJKP)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKP), &
                            W_G(KP_OF(IJKP)), TMP4(KPPP4))
            ENDIF
            IF (.NOT. FPFOI ) &
                       MOM_HO = XSI_T(IJK)*W_G(IJKP)+ &
                                (1.0-XSI_T(IJK))*W_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_gT(IJK) + Theta_Wt(IJK) * Flux_gT(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            TOP_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k+1/2)
!
            IMJK = IM_OF(IJK)
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            IJKWT = TOP_OF(IJKW)
            IMJKP = KP_OF(IMJK)
            IF(U(IMJK) >= ZERO)THEN
              MOM_LO = W_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IMJK), &
                            W_G(IM_OF(IMJK)), TMP4(IMMM4))
            ELSE
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IMJK), W_G(IJK), &
                            W_G(IPJK), W_G(IP_OF(IPJK)))
            ENDIF
            IF (.NOT. FPFOI ) &
                       MOM_HO = XSI_E(IMJK)*W_G(IJK)+ &
                                (1.0-XSI_E(IMJK))*W_G(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IMJK) * Flux_gE(IMJK) + Theta_W_te(IMJK) * Flux_gE(IMJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k+1/2)
!
            IJMK = JM_OF(IJK)
            JM = JM1(J)
            IJMKP = KP_OF(IJMK)
            IJKS = SOUTH_OF(IJK)
            IJKST = TOP_OF(IJKS)
            IF(V(IJMK) >= ZERO)THEN
              MOM_LO = W_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJMK), &
                            W_G(JM_OF(IJMK)), TMP4(JMMM4))
            ELSE
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJMK), W_G(IJK), &
                            W_G(IJPK), W_G(JP_OF(IJPK)))
            ENDIF
            IF (.NOT. FPFOI ) &
                        MOM_HO = XSI_N(IJMK)*W_G(IJK)+ &
                                 (1.0-XSI_N(IJMK))*W_G(IJMK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJMK) * Flux_gN(IJMK) + Theta_W_tn(IJMK) * Flux_gN(IJMKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
!
            IJKM = KM_OF(IJK)
            KM = KM1(K)
            IJKB = BOTTOM_OF(IJK)
            IF(WW(IJK) >= ZERO)THEN
              MOM_LO = W_G(IJKM)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKM), &
                            W_G(KM_OF(IJKM)), TMP4(KMMM4))
            ELSE
              MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJKM), W_G(IJK), &
                            W_G(IJKP), W_G(KP_OF(IJKP)))
            ENDIF
            IF (.NOT. FPFOI ) &
                       MOM_HO = XSI_T(IJKM)*W_G(IJK)+ &
                                (1.0-XSI_T(IJKM))*W_G(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJKM) * Flux_gT(IJKM) + Theta_Wt(IJKM) * Flux_gT(IJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJKM),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
!
!               CONTRIBUTION DUE TO DEFERRED CORRECTION
!
            B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
                                +BOTTOM_DC-TOP_DC
!
         ENDIF
      END DO

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_GDC


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_g1(A_W_g, IER)                               C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. Higher order C
!  discretization.                                                     C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:19-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
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
      SUBROUTINE STORE_A_W_G1(A_W_G, IER)
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
      USE parallel
      USE matrix
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc
      USE physprop
      USE fldvar
      USE output
      USE vshear
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE mflux
      USE fun_avg
      USE functions
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,&
                       IJKTE, IJKP, IJKT, IJKTN, IJK
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST
      INTEGER          IJKM, KM, IJKB
!
! loezos
        INTEGER incr
! loezos

!                      Face mass flux
      DOUBLE PRECISION Flux

!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_W_g
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3)
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                       XSI_t(DIMENSION_3)
!      DOUBLE PRECISION U(DIMENSION_3),&
!                       V(DIMENSION_3), WW(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: AW,HW,VELW
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!

!!!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         IJKP = KP_OF(IJK)
!
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           East face (i+1/2, j, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_W_be(IJK) * U_G(IJK) + Theta_W_te(IJK) * U_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE   ! Original terms
            U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K)
         ENDIF
!
!
!           North face (i, j+1/2, k+1/2)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_W_bn(IJK) * V_G(IJK) + Theta_W_tn(IJK) * V_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE   ! Original terms
            V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K)
         ENDIF
!
!
!           Top face (i, j, k+1)
         IF(CUT_W_TREATMENT_AT(IJK)) THEN
            WW(IJK) = (Theta_Wt_bar(IJK) * W_G(IJK) + Theta_Wt(IJK) * W_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJK),AW,HW,VELW)
            WW(IJK) = WW(IJK) * AW
         ELSE   ! Original terms
            WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP))
         ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO

! loezos
        incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
               XSI_T,incr)

! loezos
! update to true velocity
      IF (SHEAR) THEN
!!!$omp parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
           V(IJK)=V(IJK)+VSH(IJK)
          END IF
        END DO
      END IF
! loezos

!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!!!!$omp      parallel do                                               &
!!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!!!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!!!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!!!!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3
!
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
!
!           East face (i+1/2, j, k+1/2)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_be(IJK) * Flux_gE(IJK) + Theta_W_te(IJK) * Flux_gE(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
                  ),MU_GT(IJKTE),I),K)*ONEoDX_E_W(IJK)*AYZ_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
                  ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_W_G(IJK,E) = D_F - XSI_E(IJK)*Flux
!
            A_W_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*Flux
!
!           North face (i, j+1/2, k+1/2)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_W_bn(IJK) * Flux_gN(IJK) + Theta_W_tn(IJK) * Flux_gN(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
                  ),MU_GT(IJKTN),J),K)*ONEoDY_N_W(IJK)*AXZ_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
                  ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_W_G(IJK,N) = D_F - XSI_N(IJK)*Flux
!
            A_W_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*Flux
!
!           Top face (i, j, k+1)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_W_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Wt_bar(IJK) * Flux_gT(IJK) + Theta_Wt(IJK) * Flux_gT(IJKP))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_GT(IJKT)*ONEoDZ_T_W(IJK)*AXY_W(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
               D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            A_W_G(IJK,T) = D_F - XSI_T(IJK)*Flux
            A_W_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*Flux
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IMJKP = KP_OF(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_be(IMJK) * Flux_gE(IMJK) + Theta_W_te(IMJK) * Flux_gE(IMJKP))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_We_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                      IJKWT),MU_GT(IJKT),IM),K)*ONEoDX_E_W(IMJK)*AYZ_W(IMJK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
                  D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                     IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               A_W_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               JM = JM1(J)
               IJMKP = KP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_W_bn(IJMK) * Flux_gN(IJMK) + Theta_W_tn(IJMK) * Flux_gN(IJMKP))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',ALPHA_Wn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKST),MU_GT(IJKT),JM),K)*ONEoDY_N_W(IJMK)*AXZ_W(IJMK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
                  D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               A_W_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_W_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Wt_bar(IJKM) * Flux_gT(IJKM) + Theta_Wt(IJKM) * Flux_gT(IJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',alpha_Wt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_GT(IJK)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
                  D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               A_W_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*Flux
            ENDIF
         ENDIF
      END DO

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_G1

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

