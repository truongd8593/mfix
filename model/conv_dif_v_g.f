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
!	IF DEFERRED CORRECTION IS TO BE USED TO SOLVE V_G
!
      IF (DEF_COR) THEN
	CALL STORE_A_V_G0 (A_M(1,-3,0), IER) 
        IF (DISCRETIZE(4) > 1)CALL STORE_A_V_GDC (A_M(1,-3,0), B_M(1,0), IER) 
      ELSE  
!
        IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP 
          CALL STORE_A_V_G0 (A_M(1,-3,0), IER) 
        ELSE 
          CALL STORE_A_V_G1 (A_M(1,-3,0), IER) 
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
!                      Face velocity 
      DOUBLE PRECISION V_f 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_g 
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3) 
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
      
!$omp      parallel do                                            &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,  &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, V_f, D_f, &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                  &
!$omp&             IJMK, JM, IJKS,                                &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
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
            V_F = AVG_Y(U_G(IJK),U_G(IJPK),J) 
            D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
               ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK) 
            IF (V_F >= ZERO) THEN 
               A_V_G(IJK,E) = D_F 
               A_V_G(IPJK,W) = D_F + AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)*V_F*AYZ_V&
                  (IJK) 
            ELSE 
               A_V_G(IJK,E) = D_F - AVG_Y(ROP_G(IJKE),ROP_G(IJKNE),J)*V_F*AYZ_V&
                  (IJK) 
               A_V_G(IPJK,W) = D_F 
            ENDIF 
!
!           North face (i, j+1, k)
            V_F = AVG_Y_N(V_G(IJK),V_G(IJPK)) 
            D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK) 
            IF (V_F >= ZERO) THEN 
               A_V_G(IJK,N) = D_F 
               A_V_G(IJPK,S) = D_F + AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)*V_F*AXZ_V&
                  (IJK) 
            ELSE 
               A_V_G(IJK,N) = D_F - AVG_Y(ROP_G(IJKN),ROP_G(NORTH_OF(IJKN)),JP)&
                  *V_F*AXZ_V(IJK) 
               A_V_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
               V_F = AVG_Y(W_G(IJK),W_G(IJPK),J) 
               D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_V_G(IJK,T) = D_F 
                  A_V_G(IJKP,B) = D_F + AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)*V_F*&
                     AXY_V(IJK) 
               ELSE 
                  A_V_G(IJK,T) = D_F - AVG_Y(ROP_G(IJKT),ROP_G(IJKTN),J)*V_F*&
                     AXY_V(IJK) 
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
               V_F = AVG_Y(U_G(IMJK),U_G(IMJPK),J) 
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_V_G(IJK,W) = D_F + AVG_Y(ROP_G(IJKW),ROP_G(IJKWN),J)*V_F*&
                     AYZ_V(IMJK) 
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
               V_F = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
               D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_V_G(IJK,S) = D_F + AVG_Y(ROP_G(IJKC),ROP_G(IJKS),JM)*V_F*&
                     AXZ_V(IJMK) 
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
                  V_F = AVG_Y(W_G(IJKM),W_G(IJPKM),J) 
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                     ) 
                  IF (V_F >= ZERO) THEN 
                     A_V_G(IJK,B) = D_F + AVG_Y(ROP_G(IJKB),ROP_G(IJKBN),J)*V_F&
                        *AXY_V(IJKM) 
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
      USE vshear
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar    
      USE sendrecv
      USE sendrecv3

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

!                      Diffusion parameter 
      DOUBLE PRECISION D_f
! 
!                      Septadiagonal matrix A_V_g 
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3)
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
            TMP4(IJK4) = V_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J) 
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

!$omp      parallel do                                             &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!$omp&             IJMK, JM, IJKS,                                 &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM , &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
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
	            CONV_FAC = AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)&
		              *U(IJK)*AYZ_V(IJK) 
		    MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IPJK), V_G(IJK), & 
                            V_G(IMJK), V_G(IM_OF(IMJK)))
                     ELSE
                     ENDIF
		ELSE
		    CONV_FAC = AVG_Y(ROP_G(IJKE),ROP_G(IJKNE),J)&
		              *U(IJK)*AYZ_V(IJK)
		    MOM_LO = V_G(IPJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IPJK), & 
                            V_G(IP_OF(IPJK)), TMP4(IPPP4))
                     ELSE
                     ENDIF
		ENDIF
                     IF (.NOT. FPFOI ) THEN
		      MOM_HO = XSI_E(IJK)*V_G(IPJK)+(1.0-XSI_E(IJK))*V_G(IJK)
                     ELSE
                     ENDIF
		EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
!
		IF(V(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)&
		               *V(IJK)*AXZ_V(IJK) 
		    MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJPK), V_G(IJK), & 
                            V_G(IJMK), V_G(JM_OF(IJMK)))
                     ELSE
                     ENDIF
		ELSE
		    CONV_FAC = AVG_Y(ROP_G(IJKN),ROP_G(NORTH_OF(IJKN)),JP)&
                               *V(IJK)*AXZ_V(IJK) 
		    MOM_LO = V_G(IJPK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJPK), & 
                            V_G(JP_OF(IJPK)), TMP4(JPPP4))
                     ELSE
                     ENDIF
		ENDIF
                     IF (.NOT. FPFOI ) THEN
		      MOM_HO = XSI_N(IJK)*V_G(IJPK)+(1.0-XSI_N(IJK))*V_G(IJK)
                     ELSE
                     ENDIF
		NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)&
		              *WW(IJK)*AXY_V(IJK) 
		    MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJKP), V_G(IJK), & 
                            V_G(IJKM), V_G(KM_OF(IJKM)))
                     ELSE
                     ENDIF
		ELSE
		    CONV_FAC = AVG_Y(ROP_G(IJKT),ROP_G(IJKTN),J)&
		               *WW(IJK)*AXY_V(IJK)
		    MOM_LO = V_G(IJKP)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKP), & 
                            V_G(KP_OF(IJKP)), TMP4(KPPP4))
                     ELSE
                     ENDIF
		ENDIF
                     IF (.NOT. FPFOI ) THEN
		      MOM_HO = XSI_T(IJK)*V_G(IJKP)+(1.0-XSI_T(IJK))*V_G(IJK)
                     ELSE
                     ENDIF
		TOP_DC = CONV_FAC*(MOM_LO-MOM_HO)
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
	      CONV_FAC = AVG_Y(ROP_G(IJKW),ROP_G(IJKWN),J)*U(IMJK)*AYZ_V(IMJK) 
	      MOM_LO = V_G(IMJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IMJK), & 
                            V_G(IM_OF(IMJK)), TMP4(IMMM4))
                     ELSE
                     ENDIF	 
	    ELSE
	      CONV_FAC = AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)*U(IMJK)*AYZ_V(IMJK)
	      MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IMJK), V_G(IJK), & 
                            V_G(IPJK), V_G(IP_OF(IPJK)))
                     ELSE
                     ENDIF	
	    ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_E(IMJK)*V_G(IJK)+(1.0-XSI_E(IMJK))*V_G(IMJK)
                     ELSE
                     ENDIF
	    WEST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJKS = SOUTH_OF(IJK) 
	    IF(V(IJMK) >= ZERO)THEN
	      CONV_FAC = AVG_Y(ROP_G(IJKS),ROP_G(IJKC),JM)*V(IJMK)*AXZ_U(IJMK) 
	      MOM_LO = V_G(IJMK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJMK), & 
                            V_G(JM_OF(IJMK)), TMP4(JMMM4))
                     ELSE
                     ENDIF
	    ELSE
	      CONV_FAC = AVG_Y(ROP_G(IJKC),ROP_G(IJKN),J)*V(IJMK)*AXZ_U(IJMK)
	      MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJMK), V_G(IJK), & 
                            V_G(IJPK), V_G(JP_OF(IJPK)))
                     ELSE
                     ENDIF
	    ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_N(IJMK)*V_G(IJK)+(1.0-XSI_N(IJMK))*V_G(IJMK)
                     ELSE
                     ENDIF
	    SOUTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
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
	         CONV_FAC = AVG_Y(ROP_G(IJKB),ROP_G(IJKBN),J)&
		           *WW(IJKM)*AXY_V(IJKM) 
		 MOM_LO = V_G(IJKM)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKM), & 
                            V_G(KM_OF(IJKM)), TMP4(KMMM4))
                     ELSE
                     ENDIF
	       ELSE
		 CONV_FAC=AVG_Y(ROP_G(IJK),ROP_G(IJKN),J)*WW(IJKM)*AXY_V(IJKM)
		 MOM_LO = V_G(IJK)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(V_G(IJKM), V_G(IJK), & 
                            V_G(IJKP), V_G(KP_OF(IJKP)))
                     ELSE
                     ENDIF
	       ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_T(IJKM)*V_G(IJK)+(1.0-XSI_T(IJKM))*V_G(IJKM)
                     ELSE
                     ENDIF
	       BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
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
      SUBROUTINE STORE_A_V_G1(A_V_G, IER) 
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
      USE vshear
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
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
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J) 
      END DO 

! loezos
	incr=2		
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 

! loezos    
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)  
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
!$omp      parallel do                                             &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!$omp&             IJMK, JM, IJKS,                                 &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
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
            D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
               ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK) 
!
            A_V_G(IJK,E) = D_F - XSI_E(IJK)*AVG_Y(ROP_G(IJKE),ROP_G(IJKNE),J)*U&
               (IJK)*AYZ_V(IJK) 
!
            A_V_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*AVG_Y(ROP_G(IJKC),ROP_G(&
               IJKN),J)*U(IJK)*AYZ_V(IJK) 
!
!
!           North face (i, j+1, k)
            D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK) 
            A_V_G(IJK,N) = D_F - XSI_N(IJK)*AVG_Y(ROP_G(IJKN),ROP_G(NORTH_OF(&
               IJKN)),JP)*V(IJK)*AXZ_V(IJK) 
!
            A_V_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*AVG_Y(ROP_G(IJKC),ROP_G(&
               IJKN),J)*V(IJK)*AXZ_V(IJK) 
!
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
               D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
!
               A_V_G(IJK,T) = D_F - XSI_T(IJK)*AVG_Y(ROP_G(IJKT),ROP_G(IJKTN),J&
                  )*WW(IJK)*AXY_V(IJK) 
!
               A_V_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*AVG_Y(ROP_G(IJKC),ROP_G&
                  (IJKN),J)*WW(IJK)*AXY_V(IJK) 
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
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
!
               A_V_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*AVG_Y(ROP_G(IJKW),ROP_G&
                  (IJKWN),J)*U(IMJK)*AYZ_V(IMJK) 
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
               D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK) 
!
               A_V_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*AVG_Y(ROP_G(IJKS),ROP_G&
                  (IJKC),JM)*V(IJMK)*AXZ_V(IJMK) 
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
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                     ) 
!
                  A_V_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*AVG_Y(ROP_G(IJKB),&
                     ROP_G(IJKBN),J)*WW(IJKM)*AXY_V(IJKM) 
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
