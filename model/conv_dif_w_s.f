!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_W_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_W_s                                                      C
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
      SUBROUTINE CONV_DIF_W_S(A_M, B_M, IER) 
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
      USE physprop
      USE visc_s
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
!                      Solids phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!
      DO M = 1, MMAX 
        IF (MOMENTUM_Z_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS TO BE USED TO SOLVE W_S
!
          IF (DEF_COR) THEN
	    CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER) 
	    CALL STORE_A_W_SDC (A_M(1,-3,M), M, B_M, IER)
          ELSE
!	   NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR W_S
!
            IF (DISCRETIZE(5) == 0) THEN         ! 0 & 1 => FOUP 
              CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER) 
            ELSE 
              CALL STORE_A_W_S1 (A_M(1,-3,M), M, IER) 
            ENDIF 
          ENDIF
!
          CALL DIF_W_IS (MU_S(1,M), A_M, B_M, M, IER) 
!
        ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CONV_DIF_W_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_s0(A_W_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. FOUP         C
!  See source_w_s                                                      C
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
      SUBROUTINE STORE_A_W_S0(A_W_S, M, IER) 
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
      USE physprop
      USE visc_s
      USE toleranc 
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
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, DIMENSION_M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private(IJK,  I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK, V_f, D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
!$omp&             IJKM, KM, IJKB)
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
            V_F = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K) 
            D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKT,M),MU_S(IJKTE,M),I),K)*ODX_E(I)*AYZ_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_S(IJK,E,M) = D_F 
               A_W_S(IPJK,W,M) = D_F + AVG_Z(ROP_S(IJKC,M),ROP_S(IJKT,M),K)*V_F&
                  *AYZ_W(IJK) 
            ELSE 
               A_W_S(IJK,E,M) = D_F - AVG_Z(ROP_S(IJKE,M),ROP_S(IJKTE,M),K)*V_F&
                  *AYZ_W(IJK) 
               A_W_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i, j+1/2, k+1/2)
            V_F = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K) 
            D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKT,M),MU_S(IJKTN,M),J),K)*ODY_N(J)*AXZ_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_S(IJK,N,M) = D_F 
               A_W_S(IJPK,S,M) = D_F + AVG_Z(ROP_S(IJKC,M),ROP_S(IJKT,M),K)*V_F&
                  *AXZ_W(IJK) 
            ELSE 
               A_W_S(IJK,N,M) = D_F - AVG_Z(ROP_S(IJKN,M),ROP_S(IJKTN,M),K)*V_F&
                  *AXZ_W(IJK) 
               A_W_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1)
            V_F = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M)) 
            D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            IF (V_F >= ZERO) THEN 
               A_W_S(IJK,T,M) = D_F 
               A_W_S(IJKP,B,M) = D_F + AVG_Z(ROP_S(IJKC,M),ROP_S(IJKT,M),K)*V_F&
                  *AXY_W(IJK) 
            ELSE 
               A_W_S(IJK,T,M) = D_F - AVG_Z(ROP_S(IJKT,M),ROP_S(TOP_OF(IJKT),M)&
                  ,KP)*V_F*AXY_W(IJK) 
               A_W_S(IJKP,B,M) = D_F 
            ENDIF 
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
               V_F = AVG_Z(U_S(IMJK,M),U_S(IMJKP,M),K) 
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWT,M),MU_S(IJKT,M),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_W_S(IJK,W,M) = D_F + AVG_Z(ROP_S(IJKW,M),ROP_S(IJKWT,M),K)*&
                     V_F*AYZ_W(IMJK) 
               ELSE 
                  A_W_S(IJK,W,M) = D_F 
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
               V_F = AVG_Z(V_S(IJMK,M),V_S(IJMKP,M),K) 
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKST,M),MU_S(IJKT,M),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_W_S(IJK,S,M) = D_F + AVG_Z(ROP_S(IJKS,M),ROP_S(IJKST,M),K)*&
                     V_F*AXZ_W(IJMK) 
               ELSE 
                  A_W_S(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
               V_F = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 
               D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               IF (V_F >= ZERO) THEN 
                  A_W_S(IJK,B,M) = D_F + AVG_Z(ROP_S(IJKB,M),ROP_S(IJKC,M),KM)*&
                     V_F*AXY_W(IJKM) 
               ELSE 
                  A_W_S(IJK,B,M) = D_F 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE STORE_A_W_S0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_sdc(A_W_s, M, B_M, IER)                      C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM  C
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
      SUBROUTINE STORE_A_W_SDC(A_W_S, M, B_M, IER) 
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
      USE physprop
      USE visc_s
      USE toleranc 
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
!                      Solids phase 
      INTEGER          M 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
! 
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
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

      call lock_tmp_array
      call lock_xsi_array
!
!  Calculate convection factors
!
!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = 1, IJKMAX2 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M)) 
      END DO 
      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
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
	            CONV_FAC = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)&
		               *U(IJK)*AYZ_W(IJK) 
		    MOM_LO = W_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Z(ROP_S(IJKE,M),ROP_S(IJKTE,M),K)&
		               *U(IJK)*AYZ_W(IJK)
		    MOM_LO = W_S(IPJK,M)
		ENDIF
		MOM_HO = XSI_E(IJK)*W_S(IPJK,M)+(1.0-XSI_E(IJK))*W_S(IJK,M)
		EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
!            
		IF(V(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Z(ROP_S(IJKC,M),ROP_S(IJKT,M),K)&
		               *V(IJK)*AXZ_W(IJK) 
		    MOM_LO = W_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Z(ROP_S(IJKN,M),ROP_S(IJKTN,M),K)&
		               *V(IJK)*AXZ_W(IJK)
		    MOM_LO = W_S(IJPK,M)
		ENDIF
		MOM_HO = XSI_N(IJK)*W_S(IJPK,M)+(1.0-XSI_N(IJK))*W_S(IJK,M)
		NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
!            
		IF(WW(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Z(ROP_S(IJKC,M),ROP_S(IJKT,M),K)&
		               *WW(IJK)*AXY_W(IJK) 
		    MOM_LO = W_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Z(ROP_S(IJKT,M),ROP_S(TOP_OF(IJKT),M),KP)&
		              *WW(IJK)*AXY_W(IJK)
		    MOM_LO = W_S(IJKP,M)
		ENDIF
		MOM_HO = XSI_T(IJK)*W_S(IJKP,M)+(1.0-XSI_T(IJK))*W_S(IJK,M)
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
	      CONV_FAC = AVG_Z(ROP_S(IJKW,M),ROP_S(IJKWT,M),K)&
	                 *U(IMJK)*AYZ_W(IMJK) 
	      MOM_LO = W_S(IMJK,M)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)*U(IMJK)*AYZ_W(IMJK)
	      MOM_LO = W_S(IJK,M)
	    ENDIF
	    MOM_HO = XSI_E(IMJK)*W_S(IJK,M)+(1.0-XSI_E(IMJK))*W_S(IMJK,M)
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
	      CONV_FAC = AVG_Z(ROP_S(IJKS,M),ROP_S(IJKST,M),K)&
	                 *V(IJMK)*AXZ_W(IJMK) 
	      MOM_LO = W_S(IJMK,M)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)*V(IJMK)*AXZ_W(IJMK)
	      MOM_LO = W_S(IJK,M)
	    ENDIF
	    MOM_HO = XSI_N(IJMK)*W_S(IJK,M)+(1.0-XSI_N(IJMK))*W_S(IJMK,M)
	    SOUTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
!
            IJKM = KM_OF(IJK) 
            KM = KM1(K) 
            IJKB = BOTTOM_OF(IJK) 
            IF(WW(IJK) >= ZERO)THEN
	      CONV_FAC = AVG_Z(ROP_S(IJKB,M),ROP_S(IJKC,M),KM)&
	                 *WW(IJKM)*AXY_W(IJKM) 
	      MOM_LO = W_S(IJKM,M)
	    ELSE
	      CONV_FAC = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)&
	                 *WW(IJKM)*AXY_W(IJKM)
	      MOM_LO = W_S(IJK,M)
	    ENDIF
	    MOM_HO = XSI_T(IJKM)*W_S(IJK,M)+(1.0-XSI_T(IJKM))*W_S(IJKM,M)
	    BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 

      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_SDC 

 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_s1(A_W_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. Higher order C
!  discretization.                                                     C
!  See source_w_s                                                      C
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
      SUBROUTINE STORE_A_W_S1(A_W_S, M, IER) 
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
      USE physprop
      USE visc_s
      USE toleranc 
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
!                      Solids phase 
      INTEGER          M 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, DIMENSION_M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

      call lock_tmp_array
      call lock_xsi_array
!
!  Calculate convection factors
!
!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = 1, IJKMAX2 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M)) 
      END DO 
      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
!$omp&             IJKM, KM, IJKB)
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
            D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKT,M),MU_S(IJKTE,M),I),K)*ODX_E(I)*AYZ_W(IJK) 
!
            A_W_S(IJK,E,M) = D_F - XSI_E(IJK)*AVG_Z(ROP_S(IJKE,M),ROP_S(IJKTE,M&
               ),K)*U(IJK)*AYZ_W(IJK) 
!
            A_W_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*AVG_Z(ROP_S(IJKC,M),&
               ROP_S(IJKT,M),K)*U(IJK)*AYZ_W(IJK) 
!
!           North face (i, j+1/2, k+1/2)
            D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKT,M),MU_S(IJKTN,M),J),K)*ODY_N(J)*AXZ_W(IJK) 
!
            A_W_S(IJK,N,M) = D_F - XSI_N(IJK)*AVG_Z(ROP_S(IJKN,M),ROP_S(IJKTN,M&
               ),K)*V(IJK)*AXZ_W(IJK) 
!
            A_W_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*AVG_Z(ROP_S(IJKC,M),&
               ROP_S(IJKT,M),K)*V(IJK)*AXZ_W(IJK) 
!
!           Top face (i, j, k+1)
            D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK) 
!
            A_W_S(IJK,T,M) = D_F - XSI_T(IJK)*AVG_Z(ROP_S(IJKT,M),ROP_S(TOP_OF(&
               IJKT),M),KP)*WW(IJK)*AXY_W(IJK) 
!
            A_W_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*AVG_Z(ROP_S(IJKC,M),&
               ROP_S(IJKT,M),K)*WW(IJK)*AXY_W(IJK) 
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWT,M),MU_S(IJKT,M),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               A_W_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*AVG_Z(ROP_S(IJKW,M),&
                  ROP_S(IJKWT,M),K)*U(IMJK)*AYZ_W(IMJK) 
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKST,M),MU_S(IJKT,M),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               A_W_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*AVG_Z(ROP_S(IJKS,M),&
                  ROP_S(IJKST,M),K)*V(IJMK)*AXZ_W(IJMK) 
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
               D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               A_W_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*AVG_Z(ROP_S(IJKB,M),&
                  ROP_S(IJKC,M),KM)*WW(IJKM)*AXY_W(IJKM) 
            ENDIF 
         ENDIF 
      END DO 

      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_S1 
