!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_U_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_g                                                      C
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
      SUBROUTINE CONV_DIF_U_G(A_M, B_M, IER) 
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
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
 
!-----------------------------------------------
!
!
!
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN

            
!	IF DEFERRED CORRECTION IS USED TO SOLVE U_G
!
      IF(DEF_COR)THEN
	CALL STORE_A_U_G0(A_M(1,-3,0), IER)
        IF (DISCRETIZE(3) > 1) CALL STORE_A_U_GDC (A_M(1,-3,0),B_M(1,0), IER)
      ELSE
!
!	NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR U_G   
!
        IF (DISCRETIZE(3) == 0) THEN               ! 0 & 1 => FOUP 
          CALL STORE_A_U_G0 (A_M(1,-3,0), IER) 
        ELSE 
          CALL STORE_A_U_G1 (A_M(1,-3,0), IER) 
        ENDIF 
!
      ENDIF
            
      CALL DIF_U_IS (MU_GT, A_M, B_M, 0, IER) 
            
      RETURN  
      END SUBROUTINE CONV_DIF_U_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_g0(A_U_g, IER)                               C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
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
      SUBROUTINE STORE_A_U_G0(A_U_G, IER) 
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
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                      IJKNE, IJKP, IJKT, IJKTE
 
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE
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
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: AW,HW,VELW
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

            

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0 
      
!!!$omp      parallel do                                                  &
!!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,     &
!!!$omp&                    IJKNE, IJKP, IJKT, IJKTE, D_f,   &
!!!$omp&                    IMJK, IM, IJKW,                               &
!!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                 &
!!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)                 
      DO IJK = ijkstart3, ijkend3
!
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
!
!
!           East face (i+1, j, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IJK) * Flux_gE(IJK) + Theta_Ue(IJK) * Flux_gE(IPJK)) 
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IJK),AW,HW,VELW)
               Flux = Flux * AW 
               D_F = MU_GT(IJKE)*ONEoDX_E_U(IJK)*AYZ_U(IJK) 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
               D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN 
               A_U_G(IJK,E) = D_F 
               A_U_G(IPJK,W) = D_F + Flux
            ELSE 
               A_U_G(IJK,E) = D_F - Flux
               A_U_G(IPJK,W) = D_F 
            ENDIF 
!
!
!           North face (i+1/2, j+1/2, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJK) * Flux_gN(IJK) +Theta_U_ne(IJK) * Flux_gN(IPJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJK),AW,HW,VELW)
               Flux = Flux * AW 
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
                  ),MU_GT(IJKNE),J),I)*ONEoDY_N_U(IJK)*AXZ_U(IJK) 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
                 ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (Flux >= ZERO) THEN 
               A_U_G(IJK,N) = D_F 
               A_U_G(IJPK,S) = D_F + Flux
            ELSE 
               A_U_G(IJK,N) = D_F - Flux
               A_U_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJK) * Flux_gT(IJK) + Theta_U_te(IJK) * Flux_gT(IPJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ONEoDZ_T_U(IJK)*AXY_U(IJK) 
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,T) = D_F 
                  A_U_G(IJKP,B) = D_F + Flux
               ELSE 
                  A_U_G(IJK,T) = D_F - Flux
                  A_U_G(IJKP,B) = D_F 
               ENDIF 
            ENDIF 
!
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Ue_bar(IMJK) * Flux_gE(IMJK) + Theta_Ue(IMJK) * Flux_gE(IJK)) 
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = MU_GT(IJKC)*ONEoDX_E_U(IMJK)*AYZ_U(IMJK)  
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
                  D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,W) = D_F + Flux
               ELSE 
                  A_U_G(IJK,W) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_nw(IJMK) * Flux_gN(IJMK) +Theta_U_ne(IJMK) * Flux_gN(IPJMK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKSE),MU_GT(IJKE),JM),I)*ONEoDY_N_U(IJMK)*AXZ_U(IJMK) 
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
                  D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,S) = D_F + Flux
               ELSE 
                  A_U_G(IJK,S) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_U_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_U_tw(IJKM) * Flux_gT(IJKM) + Theta_U_te(IJKM) * Flux_gT(IPJKM))
                     CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW 
                     D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                        MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ONEoDZ_T_U(IJKM)*AXY_U(&
                        IJKM) 
                  ELSE   ! Original terms
                     Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
                     D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                       MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                       IJKM) 
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF (Flux >= ZERO) THEN 
                     A_U_G(IJK,B) = D_F + Flux
                  ELSE 
                     A_U_G(IJK,B) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
    
      RETURN  
      END SUBROUTINE STORE_A_U_G0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_GDC(A_U_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                Date: 8-APR-99   C
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
      SUBROUTINE STORE_A_U_GDC(A_U_G, B_M, IER) 
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
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar   
      USE sendrecv
      USE sendrecv3
      USE mflux 
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
!                      Error index
      INTEGER          IER
!
!                      Indices
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                       IJKNE, IJKP, IJKT, IJKTE
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
!
! loezos
	INTEGER incr
! loezos

!                      Diffusion parameter
      DOUBLE PRECISION D_f

!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
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
!
!-----------------------------------------------
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'function3.inc'


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
            TMP4(IJK4) = U_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!!!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
       DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           East face (i+1, j, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_Ue_bar(IJK) * U_g(IJK) + Theta_Ue(IJK) * U_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW 
         ELSE   ! Original terms
            U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
         ENDIF
!
!
!           North face (i+1/2, j+1/2, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_U_nw(IJK) * V_g(IJK) + Theta_U_ne(IJK) * V_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW 
         ELSE   ! Original terms
            V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
         ENDIF
!
!
!           Top face (i+1/2, j, k+1/2)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_g(IJK) + Theta_U_te(IJK) * W_g(IPJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW 
            ENDIF
         ELSE   ! Original terms
            IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
         ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!!!$omp      parallel do                                                 &
!!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!!!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!!!$omp&                    IMJK, IM, IJKW,                              &
!!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!!!$omp&              MOM_HO, MOM_LO, EAST_DC,WEST_DC,NORTH_DC,&
!!!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
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
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
!
		IF(U(IJK) >= ZERO)THEN
		    MOM_LO = U_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IPJK), U_G(IJK), & 
                            U_G(IMJK), U_G(IM_OF(IMJK)))
		ELSE
		    MOM_LO = U_G(IPJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IPJK), & 
                            U_G(IP_OF(IPJK)), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_E(IJK)*U_G(IPJK)+ &
                               (1.0-XSI_E(IJK))*U_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_U_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_Ue_bar(IJK) * Flux_gE(IJK) + Theta_Ue(IJK) * Flux_gE(IPJK)) 
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW 
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
		EAST_DC = Flux *(MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
!            
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = U_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJPK), U_G(IJK), & 
                            U_G(IJMK), U_G(JM_OF(IJMK)))
		ELSE
		    MOM_LO = U_G(IJPK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJPK), & 
                            U_G(JP_OF(IJPK)), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_N(IJK)*U_G(IJPK)+ &
                                 (1.0-XSI_N(IJK))*U_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_U_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_U_nw(IJK) * Flux_gN(IJK) +Theta_U_ne(IJK) * Flux_gN(IPJK))
                   CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJK),AW,HW,VELW)
                   Flux = Flux * AW 
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
		NORTH_DC = Flux *(MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
!
            IF (DO_K) THEN 
                IJKP = KP_OF(IJK) 
                IJKT = TOP_OF(IJK) 
                IJKTE = EAST_OF(IJKT) 
		IF(WW(IJK) >= ZERO)THEN
		   MOM_LO = U_G(IJK)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJKP), U_G(IJK), & 
                            U_G(IJKM), U_G(KM_OF(IJKM)))
		ELSE
		   MOM_LO = U_G(IJKP)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJKP), & 
                            U_G(KP_OF(IJKP)), TMP4(KPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_T(IJK)*U_G(IJKP)+ &
                               (1.0-XSI_T(IJK))*U_G(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                IF(CUT_U_TREATMENT_AT(IJK)) THEN
                   Flux = (Theta_U_tw(IJK) * Flux_gT(IJK) + Theta_U_te(IJK) * Flux_gT(IPJK))
                   CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJK),AW,HW,VELW)
                   Flux = Flux * AW 
                ELSE   ! Original terms
                   Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
                ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
		TOP_DC = Flux *(MOM_LO - MOM_HO)
		
	    ELSE
		TOP_DC = ZERO
		
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i, j, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK)
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = U_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IMJK), & 
                            U_G(IM_OF(IMJK)), TMP4(IMMM4))
 	    ELSE
	      MOM_LO = U_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IMJK), U_G(IJK), & 
                            U_G(IPJK), U_G(IP_OF(IPJK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_E(IMJK)*U_G(IJK)+ &
                               (1.0-XSI_E(IMJK))*U_G(IMJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IMJK) * Flux_gE(IMJK) + Theta_Ue(IMJK) * Flux_gE(IJK)) 
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IMJK),AW,HW,VELW)
               Flux = Flux * AW 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
	    WEST_DC = Flux * (MOM_LO - MOM_HO)

!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IPJMK = IP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKSE = EAST_OF(IJKS) 
 	    IF(V(IJMK) >= ZERO)THEN
	      MOM_LO = U_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJMK), & 
                            U_G(JM_OF(IJMK)), TMP4(JMMM4))
	    ELSE
	      MOM_LO = U_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJMK), U_G(IJK), & 
                            U_G(IJPK), U_G(JP_OF(IJPK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_N(IJMK)*U_G(IJK)+ &
                              (1.0-XSI_N(IJMK))*U_G(IJMK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJMK) * Flux_gN(IJMK) +Theta_U_ne(IJMK) * Flux_gN(IPJMK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJMK),AW,HW,VELW)
               Flux = Flux * AW 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
	    SOUTH_DC = Flux *(MOM_LO - MOM_HO)
	    
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i+1/2, j, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IPJKM = IP_OF(IJKM) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB)
	       IF(WW(IJK) >= ZERO)THEN
		 MOM_LO = U_G(IJKM)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJKM), & 
                            U_G(KM_OF(IJKM)), TMP4(KMMM4))
	       ELSE
		 MOM_LO = U_G(IJK)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJKM), U_G(IJK), & 
                            U_G(IJKP), U_G(KP_OF(IJKP)))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJKM)*U_G(IJK)+ &
                               (1.0-XSI_T(IJKM))*U_G(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJKM) * Flux_gT(IJKM) + Theta_U_te(IJKM) * Flux_gT(IPJKM))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJKM),AW,HW,VELW)
                  Flux = Flux * AW 
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
	       BOTTOM_DC = Flux * (MOM_LO - MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
	      
            ENDIF
!
!	    CONTRIBUTION DUE TO DEFERRED CORRECTION
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
      END SUBROUTINE STORE_A_U_GDC 


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_g1(A_U_g, IER)                               C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
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
      SUBROUTINE STORE_A_U_G1(A_U_G, IER) 
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
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      USE mflux 
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
!                      Error index
      INTEGER          IER
!
!                      Indices
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                       IJKNE, IJKP, IJKT, IJKTE
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE

! loezos                     SHEAR VELOCITY
      INTEGER incr    
!loezos

!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Face mass flux
      DOUBLE PRECISION Flux
!
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                      XSI_t(DIMENSION_3)
!      DOUBLE PRECISION U(DIMENSION_3),&
!                      V(DIMENSION_3), WW(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: AW,HW,VELW
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_tmp_array
      call lock_xsi_array


!
!  Calculate convection factors
!
!!!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
	 J=J_OF(IJK)
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!           East face (i+1, j, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            U(IJK) = (Theta_Ue_bar(IJK) * U_g(IJK) + Theta_Ue(IJK) * U_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IJK),AW,HW,VELW)
            U(IJK) = U(IJK) * AW 
         ELSE   ! Original terms
            U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
         ENDIF
!
!
!           North face (i+1/2, j+1/2, k)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            V(IJK) = (Theta_U_nw(IJK) * V_g(IJK) + Theta_U_ne(IJK) * V_g(IPJK))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJK),AW,HW,VELW)
            V(IJK) = V(IJK) * AW 
         ELSE   ! Original terms
            V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
         ENDIF
!
!
!           Top face (i+1/2, j, k+1/2)
         IF(CUT_U_TREATMENT_AT(IJK)) THEN
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_g(IJK) + Theta_U_te(IJK) * W_g(IPJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJK),AW,HW,VELW)
               WW(IJK) = WW(IJK) * AW 
            ENDIF
         ELSE   ! Original terms
            IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
         ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 

! loezos      
!update V to true velocity
      IF (SHEAR) THEN
!!!$omp parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSHE(IJK)	
          END IF
        END DO
      END IF
! loezos


!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!!!$omp      parallel do                                                 &
!!!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!!!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!!!$omp&                    IMJK, IM, IJKW,                              &
!!!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!!!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3
!
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
!
!           East face (i+1, j, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_Ue_bar(IJK) * Flux_gE(IJK) + Theta_Ue(IJK) * Flux_gE(IPJK)) 
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IJK),AW,HW,VELW)
               Flux = Flux * AW 
               D_F = MU_GT(IJKE)*ONEoDX_E_U(IJK)*AYZ_U(IJK) 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
               D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_U_G(IJK,E) = D_F - XSI_E(IJK) * Flux
!
            A_U_G(IPJK,W) = D_F + (ONE - XSI_E(IJK)) * Flux
!
!
!           North face (i+1/2, j+1/2, k)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_U_TREATMENT_AT(IJK)) THEN
               Flux = (Theta_U_nw(IJK) * Flux_gN(IJK) +Theta_U_ne(IJK) * Flux_gN(IPJK))
               CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJK),AW,HW,VELW)
               Flux = Flux * AW 
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
                   ),MU_GT(IJKNE),J),I)*ONEoDY_N_U(IJK)*AXZ_U(IJK) 
            ELSE   ! Original terms
               Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
                   ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_U_G(IJK,N) = D_F - XSI_N(IJK) * Flux 
!
            A_U_G(IJPK,S) = D_F + (ONE - XSI_N(IJK)) * Flux 
!
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_tw(IJK) * Flux_gT(IJK) + Theta_U_te(IJK) * Flux_gT(IPJK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ONEoDZ_T_U(IJK)*AXY_U(IJK) 
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                   IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_U_G(IJK,T) = D_F - XSI_T(IJK) * Flux
!
               A_U_G(IJKP,B) = D_F + (ONE - XSI_T(IJK)) * Flux 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_Ue_bar(IMJK) * Flux_gE(IMJK) + Theta_Ue(IMJK) * Flux_gE(IJK)) 
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',alpha_Ue_c(IMJK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = MU_GT(IJKC)*ONEoDX_E_U(IMJK)*AYZ_U(IMJK)  
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
                  D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_U_G(IJK,W) = D_F + (ONE - XSI_E(IMJK)) * Flux
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_U_TREATMENT_AT(IJK)) THEN
                  Flux = (Theta_U_nw(IJMK) * Flux_gN(IJMK) +Theta_U_ne(IJMK) * Flux_gN(IPJMK))
                  CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Un_c(IJMK),AW,HW,VELW)
                  Flux = Flux * AW 
                  D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                        IJKSE),MU_GT(IJKE),JM),I)*ONEoDY_N_U(IJMK)*AXZ_U(IJMK) 
               ELSE   ! Original terms
                  Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
                  D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                     IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_U_G(IJK,S) = D_F + (ONE - XSI_N(IJMK)) * Flux
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_U_TREATMENT_AT(IJK)) THEN
                     Flux = (Theta_U_tw(IJKM) * Flux_gT(IJKM) + Theta_U_te(IJKM) * Flux_gT(IPJKM))
                     CALL GET_INTERPOLATION_TERMS_G(IJK,'U_MOMENTUM',ALPHA_Ut_c(IJKM),AW,HW,VELW)
                     Flux = Flux * AW 
                     D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                       MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ONEoDZ_T_U(IJKM)*AXY_U(&
                       IJKM) 
                  ELSE   ! Original terms
                     Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
                     D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                        MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                        IJKM) 
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
                  A_U_G(IJK,B) = D_F + (ONE - XSI_T(IJKM)) * Flux
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 

      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_G1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
