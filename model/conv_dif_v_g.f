!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_v_g                                                      C
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
      SUBROUTINE CONV_DIF_V_G(A_M, B_M, IER)
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
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN
!       IF DEFERRED CORRECTION IS TO BE USED TO SOLVE V_G
!
      IF (DEF_COR) THEN
        CALL STORE_A_V_G0 (A_M(1,-3,0), IER)
        IF (DISCRETIZE(4) > 1)CALL STORE_A_V_GDC (A_M(1,-3,0), B_M(1,0), IER)
      ELSE
!
        IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP
          CALL STORE_A_V_G0 (A_M(1,-3,0), IER)
        ELSE
          CALL STORE_A_V_G1 (A_M(1,-3,0))
        ENDIF
      ENDIF

      CALL DIF_V_IS (MU_GT, A_M, B_M, 0, IER)

      RETURN
      END SUBROUTINE CONV_DIF_V_G
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_g0(A_V_g, IER)                               C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-JUN-96   C
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
      SUBROUTINE STORE_A_V_G0(A_V_G, IER)
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
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                       IJKNE, IJKP, IJKT, IJKTN, IJK
!
!                      Indices
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM
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
!                      Septadiagonal matrix A_V_g
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3)
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

!!!$omp      parallel do                                            &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,  &
!!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, D_f, &
!!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                  &
!!!$omp&             IJMK, JM, IJKS,                                &
!!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN
!
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
!
!
!           East face (i+1/2, j+1/2, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IJK) * Flux_gE(IJK) +Theta_V_ne(IJK) * Flux_gE(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
                  ),MU_GT(IJKNE),I),J)*ONEoDX_E_V(IJK)*AYZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
                  ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN
               A_V_G(IJK,E) = D_F
               A_V_G(IPJK,W) = D_F + Flux
            ELSE
               A_V_G(IJK,E) = D_F - Flux
               A_V_G(IPJK,W) = D_F
            ENDIF
!
!           North face (i, j+1, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJK) * Flux_gN(IJK) + Theta_Vn(IJK) * Flux_gN(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJK) ,AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_GT(IJKN)*ONEoDY_N_V(IJK)*AXZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
               D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN
               A_V_G(IJK,N) = D_F
               A_V_G(IJPK,S) = D_F + Flux
            ELSE
               A_V_G(IJK,N) = D_F - Flux
               A_V_G(IJPK,S) = D_F
            ENDIF
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJK) * Flux_gT(IJK) + Theta_V_st(IJK) * Flux_gT(IJPK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKN),MU_GT(IJKTN),K),J)*OX(I)*ONEoDZ_T_V(IJK)*AXY_V(IJK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                    IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_V_G(IJK,T) = D_F
                  A_V_G(IJKP,B) = D_F + Flux
               ELSE
                  A_V_G(IJK,T) = D_F - Flux
                  A_V_G(IJKP,B) = D_F
               ENDIF
            ENDIF
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWN = NORTH_OF(IJKW)
               IMJPK = JP_OF(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_se(IMJK) * Flux_gE(IMJK) +Theta_V_ne(IMJK) * Flux_gE(IMJPK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                     IJKWN),MU_GT(IJKN),IM),J)*ONEoDX_E_V(IMJK)*AYZ_V(IMJK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
                  D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                     IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_V_G(IJK,W) = D_F + Flux
               ELSE
                  A_V_G(IJK,W) = D_F
               ENDIF
            ENDIF
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               JM = JM1(J)
               IJKS = SOUTH_OF(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Vn_bar(IJMK) * Flux_gN(IJMK) + Theta_Vn(IJMK) * Flux_gN(IJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_GT(IJKC)*ONEoDY_N_V(IJMK)*AXZ_V(IJMK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
                  D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN
                  A_V_G(IJK,S) = D_F + Flux
               ELSE
                  A_V_G(IJK,S) = D_F
               ENDIF
            ENDIF
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  KM = KM1(K)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBN = NORTH_OF(IJKB)
                  IJPKM = JP_OF(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_V_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_V_nt(IJKM) * Flux_gT(IJKM) + Theta_V_st(IJKM) * Flux_gT(IJPKM))
                     CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ONEoDZ_T_V(IJKM)*AXY_V(IJKM)
                  ELSE   ! Original terms
                     Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
                     D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                        MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                        )
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF (Flux >= ZERO) THEN
                     A_V_G(IJK,B) = D_F + Flux
                  ELSE
                     A_V_G(IJK,B) = D_F
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE STORE_A_V_G0

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_gdc(A_V_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE V-MOMENTUM  C
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
      SUBROUTINE STORE_A_V_GDC(A_V_G, B_M, IER)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
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
      USE visc_g
      USE vshear
      USE xsi
      USE xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3

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
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                       IJKNE, IJKP, IJKT, IJKTN, IJK
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
!
! loezos
        INTEGER incr
! loezos

!                      Septadiagonal matrix A_V_g
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
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
            TMP4(IJK4) = V_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!!!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKN = NORTH_OF(IJK)
!
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           East face (i+1/2, j+1/2, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_V_se(IJK) * U_g(IJK) +Theta_V_ne(IJK) * U_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ve_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
         ELSE   ! Original terms
            U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J)
         ENDIF
!
!
!           North face (i, j+1, k)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_Vn_bar(IJK) * V_g(IJK) + Theta_Vn(IJK) * V_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Vn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
         ELSE   ! Original terms
            V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK))
         ENDIF
!
!
!           Top face (i, j+1/2, k+1/2)
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_V_nt(IJK) * W_g(IJK) + Theta_V_st(IJK) * W_g(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Vt_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE   ! Original terms
            IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J)
         ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO

! loezos
        incr=2
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr)
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!!!$omp      parallel do                                             &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!!!$omp&             IJMK, JM, IJKS,                                 &
!!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM , &
!!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
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
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j+1/2, k)
!
                IF(U(IJK) >= ZERO)THEN
                    MOM_LO = V_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IPJK), V_G(IJK), &
                            V_G(IMJK), V_G(IM_OF(IMJK)))
                ELSE
                    MOM_LO = V_G(IPJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IPJK), &
                            V_G(IP_OF(IPJK)), TMP4(IPPP4))
                ENDIF
                IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_E(IJK)*V_G(IPJK)+(1.0-XSI_E(IJK))*V_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_V_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_V_se(IJK) * Flux_gE(IJK) +Theta_V_ne(IJK) * Flux_gE(IJPK))
                   CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IJK),AW,HW,VELW)
                   Flux = Flux * AW
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
!
                IF(V(IJK) >= ZERO)THEN
                    MOM_LO = V_G(IJK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJPK), V_G(IJK), &
                            V_G(IJMK), V_G(JM_OF(IJMK)))
                ELSE
                    MOM_LO = V_G(IJPK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJPK), &
                            V_G(JP_OF(IJPK)), TMP4(JPPP4))
                ENDIF
                IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_N(IJK)*V_G(IJPK)+(1.0-XSI_N(IJK))*V_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_V_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_Vn_bar(IJK) * Flux_gN(IJK) + Theta_Vn(IJK) * Flux_gN(IJPK))
                   CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJK) ,AW,HW,VELW)
                   Flux = Flux * AW
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
!
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IF(WW(IJK) >= ZERO)THEN
                    MOM_LO = V_G(IJK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJKP), V_G(IJK), &
                            V_G(IJKM), V_G(KM_OF(IJKM)))
                ELSE
                    MOM_LO = V_G(IJKP)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKP), &
                            V_G(KP_OF(IJKP)), TMP4(KPPP4))
                ENDIF
                IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_T(IJK)*V_G(IJKP)+(1.0-XSI_T(IJK))*V_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_V_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_V_nt(IJK) * Flux_gT(IJK) + Theta_V_st(IJK) * Flux_gT(IJPK))
                   CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJK),AW,HW,VELW)
                   Flux = Flux * AW
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                TOP_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
                TOP_DC = ZERO

            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j+1/2, k)
!
            IMJK = IM_OF(IJK)
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            IJKWN = NORTH_OF(IJKW)
            IMJPK = JP_OF(IMJK)
            IF(U(IMJK) >= ZERO)THEN
              MOM_LO = V_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IMJK), &
                            V_G(IM_OF(IMJK)), TMP4(IMMM4))
            ELSE
              MOM_LO = V_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IMJK), V_G(IJK), &
                            V_G(IPJK), V_G(IP_OF(IPJK)))
            ENDIF
            IF (.NOT. FPFOI ) &
              MOM_HO = XSI_E(IMJK)*V_G(IJK)+(1.0-XSI_E(IMJK))*V_G(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IMJK) * Flux_gE(IMJK) +Theta_V_ne(IMJK) * Flux_gE(IMJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
              Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
!
            IJMK = JM_OF(IJK)
            JM = JM1(J)
            IJKS = SOUTH_OF(IJK)
            IF(V(IJMK) >= ZERO)THEN
              MOM_LO = V_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJMK), &
                            V_G(JM_OF(IJMK)), TMP4(JMMM4))
            ELSE
              MOM_LO = V_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJMK), V_G(IJK), &
                            V_G(IJPK), V_G(JP_OF(IJPK)))
            ENDIF
            IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_N(IJMK)*V_G(IJK)+(1.0-XSI_N(IJMK))*V_G(IJMK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJMK) * Flux_gN(IJMK) + Theta_Vn(IJMK) * Flux_gN(IJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j+1/2, k-1/2)
!
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               IJKBN = NORTH_OF(IJKB)
               IJPKM = JP_OF(IJKM)
               IF(WW(IJK) >= ZERO)THEN
                 MOM_LO = V_G(IJKM)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKM), &
                            V_G(KM_OF(IJKM)), TMP4(KMMM4))
               ELSE
                 MOM_LO = V_G(IJK)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJKM), V_G(IJK), &
                            V_G(IJKP), V_G(KP_OF(IJKP)))
               ENDIF
               IF (.NOT. FPFOI ) &
                      MOM_HO = XSI_T(IJKM)*V_G(IJK)+(1.0-XSI_T(IJKM))*V_G(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJKM) * Flux_gT(IJKM) + Theta_V_st(IJKM) * Flux_gT(IJPKM))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
               BOTTOM_DC = ZERO
            ENDIF
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
      END SUBROUTINE STORE_A_V_GDC


!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_g1(A_V_g, IER)                               C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:20-MAR-97   C
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
      SUBROUTINE STORE_A_V_G1(A_V_G)
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
      USE vshear
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
!                      Indices
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,&
                       IJKNE, IJKP, IJKT, IJKTN, IJK
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM
!
! start loezos
      INTEGER incr
! end loezos

!                      Face mass flux
      DOUBLE PRECISION Flux

!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_V_g
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3)
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
!!!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKN = NORTH_OF(IJK)
!
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         IF(CUT_V_TREATMENT_AT(IJK)) THEN
!           East face (i+1/2, j+1/2, k)
            U(IJK) = (Theta_V_se(IJK) * U_g(IJK) +Theta_V_ne(IJK) * U_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ve_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW
!           North face (i, j+1, k)
            V(IJK) = (Theta_Vn_bar(IJK) * V_g(IJK) + Theta_Vn(IJK) * V_g(IJPK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Vn_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               WW(IJK) = (Theta_V_nt(IJK) * W_g(IJK) + Theta_V_st(IJK) * W_g(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Vt_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF
         ELSE   ! Original terms
            U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J)
            V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK))
            IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J)
         ENDIF

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO

! loezos
        incr=2
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr)

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
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!!!$omp      parallel do                                             &
!!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!!!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!!!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!!!$omp&             IJMK, JM, IJKS,                                 &
!!!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN
!
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
!
!           East face (i+1/2, j+1/2, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_V_se(IJK) * Flux_gE(IJK) +Theta_V_ne(IJK) * Flux_gE(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IJK),AW,HW,VELW)
               Flux = Flux * AW
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
                  ),MU_GT(IJKNE),I),J)*ONEoDX_E_V(IJK)*AYZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
                  ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_V_G(IJK,E) = D_F - XSI_E(IJK)*Flux
!
            A_V_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*Flux
!
!
!           North face (i, j+1, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_V_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Vn_bar(IJK) * Flux_gN(IJK) + Theta_Vn(IJK) * Flux_gN(IJPK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJK) ,AW,HW,VELW)
               Flux = Flux * AW
               D_F = MU_GT(IJKN)*ONEoDY_N_V(IJK)*AXZ_V(IJK)
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
               D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK)
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            A_V_G(IJK,N) = D_F - XSI_N(IJK)*Flux
!
            A_V_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*Flux
!
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_nt(IJK) * Flux_gT(IJK) + Theta_V_st(IJK) * Flux_gT(IJPK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKN),MU_GT(IJKTN),K),J)*OX(I)*ONEoDZ_T_V(IJK)*AXY_V(IJK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_V_G(IJK,T) = D_F - XSI_T(IJK)*Flux
!
               A_V_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*Flux
            ENDIF
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IMJK)) THEN
               IM = IM1(I)
               IJKW = WEST_OF(IJK)
               IJKWN = NORTH_OF(IJKW)
               IMJPK = JP_OF(IMJK)
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_V_se(IMJK) * Flux_gE(IMJK) +Theta_V_ne(IMJK) * Flux_gE(IMJPK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Ve_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                     IJKWN),MU_GT(IJKN),IM),J)*ONEoDX_E_V(IMJK)*AYZ_V(IMJK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
                  D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                     IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_V_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK)
            IF (.NOT.FLOW_AT_N(IJMK)) THEN
               JM = JM1(J)
               IJKS = SOUTH_OF(IJK)
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_V_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Vn_bar(IJMK) * Flux_gN(IJMK) + Theta_Vn(IJMK) * Flux_gN(IJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',alpha_Vn_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW
                  D_F = MU_GT(IJKC)*ONEoDY_N_V(IJMK)*AXZ_V(IJMK)
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
                  D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK)
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_V_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF (.NOT.FLOW_AT_N(IJKM)) THEN
                  KM = KM1(K)
                  IJKB = BOTTOM_OF(IJK)
                  IJKBN = NORTH_OF(IJKB)
                  IJPKM = JP_OF(IJKM)
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_V_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_V_nt(IJKM) * Flux_gT(IJKM) + Theta_V_st(IJKM) * Flux_gT(IJPKM))
                     CALL GET_INTERPOLATION_TERMS_G(IJK,'V_MOMENTUM',ALPHA_Vt_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW
                     D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ONEoDZ_T_V(IJKM)*AXY_V(IJKM)
                  ELSE   ! Original terms
                     Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
                     D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                        MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                        )
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
                  A_V_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*Flux
               ENDIF
            ENDIF
         ENDIF
      END DO



      call unlock_tmp_array
      call unlock_xsi_array


      RETURN
      END SUBROUTINE STORE_A_V_G1

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
