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
!	IF DEFERRED CORRECTION IS USED TO SOLVE W_G
      IF (DEF_COR) THEN
	CALL STORE_A_W_G0 (A_M(1,-3,0), IER) 
	CALL STORE_A_W_GDC (A_M(1,-3,0), B_M(1,0), IER)
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
!                      Face velocity 
      DOUBLE PRECISION V_f 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_g 
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3) 
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
!$omp      parallel do                                               &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK, V_f, D_f,    &
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!$omp&             IJKM, KM, IJKB)
      DO IJK = 1, IJKMAX2 
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
            V_F = AVG_Z(U_G(IJK),U_G(IJKP),K) 
            D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
               ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_G(IJK,E) = D_F 
               A_W_G(IPJK,W) = D_F + AVG_Z(ROP_G(IJKC),ROP_G(IJKT),K)*V_F*AYZ_W&
                  (IJK) 
            ELSE 
               A_W_G(IJK,E) = D_F - AVG_Z(ROP_G(IJKE),ROP_G(IJKTE),K)*V_F*AYZ_W&
                  (IJK) 
               A_W_G(IPJK,W) = D_F 
            ENDIF 
!
!           North face (i, j+1/2, k+1/2)
            V_F = AVG_Z(V_G(IJK),V_G(IJKP),K) 
            D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
               ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_G(IJK,N) = D_F 
               A_W_G(IJPK,S) = D_F + AVG_Z(ROP_G(IJKC),ROP_G(IJKT),K)*V_F*AXZ_W&
                  (IJK) 
            ELSE 
               A_W_G(IJK,N) = D_F - AVG_Z(ROP_G(IJKN),ROP_G(IJKTN),K)*V_F*AXZ_W&
                  (IJK) 
               A_W_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1)
            V_F = AVG_Z_T(W_G(IJK),W_G(IJKP)) 
            D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_G(IJK,T) = D_F 
               A_W_G(IJKP,B) = D_F + AVG_Z(ROP_G(IJKC),ROP_G(IJKT),K)*V_F*AXY_W&
                  (IJK) 
            ELSE 
               A_W_G(IJK,T) = D_F - AVG_Z(ROP_G(IJKT),ROP_G(TOP_OF(IJKT)),KP)*&
                  V_F*AXY_W(IJK) 
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
               V_F = AVG_Z(U_G(IMJK),U_G(IMJKP),K) 
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_W_G(IJK,W) = D_F + AVG_Z(ROP_G(IJKW),ROP_G(IJKWT),K)*V_F*&
                     AYZ_W(IMJK) 
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
               V_F = AVG_Z(V_G(IJMK),V_G(IJMKP),K) 
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_W_G(IJK,S) = D_F + AVG_Z(ROP_G(IJKS),ROP_G(IJKST),K)*V_F*&
                     AXZ_W(IJMK) 
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
               V_F = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
               D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               IF (V_F >= ZERO) THEN 
                  A_W_G(IJK,B) = D_F + AVG_Z(ROP_G(IJKB),ROP_G(IJKC),KM)*V_F*&
                     AXY_W(IJKM) 
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
      Use xsi_array
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
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
!
!                      Septadiagonal matrix A_W_g 
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!	FACE VELOCITY
	DOUBLE PRECISION V_F
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION CONV_FAC
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
! 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      
      call lock_tmp_array
      call lock_xsi_array
!
!  Calculate convection factors
!
!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = 1, IJKMAX2 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP)) 
      END DO 
      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N, XSI_T) 
!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                               &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!$omp&             IJKM, KM, IJKB, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = 1, IJKMAX2 
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
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k+1/2)
!            
	    IF(U(IJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)*U(IJK)*AYZ_W(IJK) 
	      MOM_LO = W_G(IJK)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJKE),ROP_G(IJKTE),K)*U(IJK)*AYZ_W(IJK)
	      MOM_LO = W_G(IPJK)
	    ENDIF
	    MOM_HO = XSI_E(IJK)*W_G(IPJK)+(1.0-XSI_E(IJK))*W_G(IJK)
	    EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
!            
	    IF(V(IJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJKC),ROP_G(IJKT),K)*V(IJK)*AXZ_W(IJK) 
	      MOM_LO = W_G(IJK)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJKN),ROP_G(IJKTN),K)*V(IJK)*AXZ_W(IJK)
	      MOM_LO = W_G(IJPK)
	    ENDIF
	    MOM_HO = XSI_N(IJK)*W_G(IJPK)+(1.0-XSI_N(IJK))*W_G(IJK)
	    NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
!
	    IF(WW(IJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJKC),ROP_G(IJKT),K)*WW(IJK)*AXY_W(IJK) 
	      MOM_LO = W_G(IJK)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJKT),ROP_G(TOP_OF(IJKT)),KP)&
	                *WW(IJK)*AXY_W(IJK)
	      MOM_LO = W_G(IJKP)
	    ENDIF
	    MOM_HO = XSI_T(IJK)*W_G(IJKP)+(1.0-XSI_T(IJK))*W_G(IJK)
	    TOP_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k+1/2)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
            IJKWT = TOP_OF(IJKW) 
            IMJKP = KP_OF(IMJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJKW),ROP_G(IJKWT),K)*U(IMJK)*AYZ_W(IMJK) 
	      MOM_LO = W_G(IMJK)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)*U(IMJK)*AYZ_W(IMJK)
	      MOM_LO = W_G(IJK)
	    ENDIF
	    MOM_HO = XSI_E(IMJK)*W_G(IJK)+(1.0-XSI_E(IMJK))*W_G(IMJK)
	    WEST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k+1/2)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJMKP = KP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKST = TOP_OF(IJKS) 
            IF(V(IJMK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJKS),ROP_G(IJKST),K)*V(IJMK)*AXZ_W(IJMK) 
	      MOM_LO = W_G(IJMK)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)*V(IJMK)*AXZ_W(IJMK)
	      MOM_LO = W_G(IJK)
	    ENDIF
	    MOM_HO = XSI_N(IJMK)*W_G(IJK)+(1.0-XSI_N(IJMK))*W_G(IJMK)
	    SOUTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
!
            IJKM = KM_OF(IJK) 
            KM = KM1(K) 
            IJKB = BOTTOM_OF(IJK) 
	    IF(WW(IJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_G(IJKB),ROP_G(IJKC),KM)&
	                 *WW(IJKM)*AXY_W(IJKM) 
	      MOM_LO = W_G(IJKM)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)*WW(IJKM)*AXY_W(IJKM)
	      MOM_LO = W_G(IJK)
	    ENDIF
	    MOM_HO = XSI_T(IJKM)*W_G(IJK)+(1.0-XSI_T(IJKM))*W_G(IJKM)
	    BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
            B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 
      
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
      Use xsi_array
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
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
! 
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
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      
      call lock_tmp_array
      call lock_xsi_array
!
!  Calculate convection factors
!
!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = 1, IJKMAX2 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP)) 
      END DO 
      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N, XSI_T) 
!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!!$omp      parallel do                                               &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!!$omp&             IJKM, KM, IJKB)
      DO IJK = 1, IJKMAX2 
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
            D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
               ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK) 
!
            A_W_G(IJK,E) = D_F - XSI_E(IJK)*AVG_Z(ROP_G(IJKE),ROP_G(IJKTE),K)*U&
               (IJK)*AYZ_W(IJK) 
!
            A_W_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*AVG_Z(ROP_G(IJKC),ROP_G(&
               IJKT),K)*U(IJK)*AYZ_W(IJK) 
!
!           North face (i, j+1/2, k+1/2)
            D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
               ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK) 
!
            A_W_G(IJK,N) = D_F - XSI_N(IJK)*AVG_Z(ROP_G(IJKN),ROP_G(IJKTN),K)*V&
               (IJK)*AXZ_W(IJK) 
!
            A_W_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*AVG_Z(ROP_G(IJKC),ROP_G(&
               IJKT),K)*V(IJK)*AXZ_W(IJK) 
!
!           Top face (i, j, k+1)
            D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            A_W_G(IJK,T) = D_F - XSI_T(IJK)*AVG_Z(ROP_G(IJKT),ROP_G(TOP_OF(IJKT&
               )),KP)*WW(IJK)*AXY_W(IJK) 
            A_W_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*AVG_Z(ROP_G(IJKC),ROP_G(&
               IJKT),K)*WW(IJK)*AXY_W(IJK) 
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               A_W_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*AVG_Z(ROP_G(IJKW),ROP_G&
                  (IJKWT),K)*U(IMJK)*AYZ_W(IMJK) 
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               A_W_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*AVG_Z(ROP_G(IJKS),ROP_G&
                  (IJKST),K)*V(IJMK)*AXZ_W(IJMK) 
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
               D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               A_W_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*AVG_Z(ROP_G(IJKB),ROP_G&
                  (IJKC),KM)*WW(IJKM)*AXY_W(IJKM) 
            ENDIF 
         ENDIF 
      END DO 
      
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_G1 
