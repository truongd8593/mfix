!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
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
      SUBROUTINE CONV_DIF_V_S(A_M, B_M, IER) 
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
      USE compar    !//AIKEPARDBG
      USE sendrecv  !// 400      
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
  
 
!-----------------------------------------------
!
!



      DO M = 1, MMAX 
        IF  (MOMENTUM_Y_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED TO SOLVE V_S
           IF (DEF_COR) THEN
	     CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER)
	     CALL STORE_A_V_SDC (A_M(1,-3,M), M, B_M, IER)
           ELSE    
!
             IF (DISCRETIZE(4) == 0) THEN         ! 0 & 1 => FOUP 
               CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER) 
             ELSE 
               CALL STORE_A_V_S1 (A_M(1,-3,M), M, IER) 
             ENDIF 
           ENDIF
!
!// 400 1229 COMM A_M and B_M      
      CALL SEND_RECV(A_M, 2)
      CALL SEND_RECV(B_M, 2)

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef dif_v_is in conv_dif_v_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

           CALL DIF_V_IS (MU_S(1,M), A_M, B_M, M, IER) 

!// 400 1229 COMM A_M and B_M      
      CALL SEND_RECV(A_M, 2)
      CALL SEND_RECV(B_M, 2)

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft dif_u_is in conv_dif_v_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
	   
!
        ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CONV_DIF_V_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_s0(A_V_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_v_s                                                      C
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
      SUBROUTINE STORE_A_V_S0(A_V_S, M, IER) 
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
      USE compar    !//d
      USE sendrecv  !// 400      
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
!                      Solids phase 
      INTEGER          M 
! 
!                      Face velocity 
      DOUBLE PRECISION V_f 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, DIMENSION_M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1225 Communicate boundaries
      call send_recv(U_S,2)
      call send_recv(V_S,2)
      call send_recv(W_S,2)
      call send_recv(MU_S,2)
      call send_recv(AYZ_V,2)
      call send_recv(AXZ_V,2)
      call send_recv(AXY_V,2)      
      call send_recv(ROP_S,2)
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3 

!$omp      parallel do                                                  &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, V_f, D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM ) 
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
!// 360 1229 Check if current i,j,k resides on this PE	    
            IF(.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE
	    	    
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
            V_F = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J) 
            D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKN,M),MU_S(IJKNE,M),I),J)*ODX_E(I)*AYZ_V(IJK) 
            IF (V_F >= ZERO) THEN 
               A_V_S(IJK,E,M) = D_F 
               A_V_S(IPJK,W,M) = D_F + AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)*V_F&
                  *AYZ_V(IJK) 
            ELSE 
               A_V_S(IJK,E,M) = D_F - AVG_Y(ROP_S(IJKE,M),ROP_S(IJKNE,M),J)*V_F&
                  *AYZ_V(IJK) 
               A_V_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i, j+1, k)
            V_F = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M)) 
            D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK) 
            IF (V_F >= ZERO) THEN 
               A_V_S(IJK,N,M) = D_F 
               A_V_S(IJPK,S,M) = D_F + AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)*V_F&
                  *AXZ_V(IJK) 
            ELSE 
               A_V_S(IJK,N,M) = D_F - AVG_Y(ROP_S(IJKN,M),ROP_S(NORTH_OF(IJKN),&
                  M),JP)*V_F*AXZ_V(IJK) 
               A_V_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
               V_F = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J) 
               D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKN,M),MU_S(IJKTN,M),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_V_S(IJK,T,M) = D_F 
                  A_V_S(IJKP,B,M) = D_F + AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)*&
                     V_F*AXY_V(IJK) 
               ELSE 
                  A_V_S(IJK,T,M) = D_F - AVG_Y(ROP_S(IJKT,M),ROP_S(IJKTN,M),J)*&
                     V_F*AXY_V(IJK) 
                  A_V_S(IJKP,B,M) = D_F 
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
               V_F = AVG_Y(U_S(IMJK,M),U_S(IMJPK,M),J) 
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWN,M),MU_S(IJKN,M),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_V_S(IJK,W,M) = D_F + AVG_Y(ROP_S(IJKW,M),ROP_S(IJKWN,M),J)*&
                     V_F*AYZ_V(IMJK) 
               ELSE 
                  A_V_S(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
               V_F = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
               D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_V_S(IJK,S,M) = D_F + AVG_Y(ROP_S(IJKS,M),ROP_S(IJKC,M),JM)*&
                     V_F*AXZ_V(IJMK) 
               ELSE 
                  A_V_S(IJK,S,M) = D_F 
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
                  V_F = AVG_Y(W_S(IJKM,M),W_S(IJPKM,M),J) 
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(&
                     IJKM) 
                  IF (V_F >= ZERO) THEN 
                     A_V_S(IJK,B,M) = D_F + AVG_Y(ROP_S(IJKB,M),ROP_S(IJKBN,M),&
                        J)*V_F*AXY_V(IJKM) 
                  ELSE 
                     A_V_S(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      
!//? Check if need to COMM A_V_S ??       

      RETURN  
      END SUBROUTINE STORE_A_V_S0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_sdc(A_V_s, M, B_M, IER)                      C
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
      SUBROUTINE STORE_A_V_SDC(A_V_S, M, B_M, IER) 
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
      USE compar    !//d
      USE sendrecv  !// 400       
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
!                      Solids phase 
      INTEGER          M 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3)
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

! loezos
	INTEGER  incr
! loezos

      call lock_tmp_array
      call lock_xsi_array

!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1225 Communicate boundaries
      call send_recv(U_S,2)
      call send_recv(V_S,2)
      call send_recv(W_S,2)
      call send_recv(MU_S,2)
      call send_recv(AYZ_V,2)
      call send_recv(AXZ_V,2)
      call send_recv(AXY_V,2)      
      call send_recv(ROP_S,2)
      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)      
!
!  Calculate convection factors
!
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J) 
      END DO 

! loezos
	incr=2		
! loezos


      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,&
			incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC )
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
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j+1/2, k)
!            
		IF(U(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)&
		               *U(IJK)*AYZ_V(IJK) 
		    MOM_LO = V_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Y(ROP_S(IJKE,M),ROP_S(IJKNE,M),J)&
		              *U(IJK)*AYZ_V(IJK)
		    MOM_LO = V_S(IPJK,M)
		ENDIF
		MOM_HO = XSI_E(IJK)*V_S(IPJK,M)+(1.0-XSI_E(IJK))*V_S(IJK,M)
		EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
!
		IF(V(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)&
		               *V(IJK)*AXZ_V(IJK) 
		    MOM_LO = V_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Y(ROP_S(IJKN,M),ROP_S(NORTH_OF(IJKN),M),JP)&
                               *V(IJK)*AXZ_V(IJK) 
		    MOM_LO = V_S(IJPK,M)
		ENDIF
		MOM_HO = XSI_N(IJK)*V_S(IJPK,M)+(1.0-XSI_N(IJK))*V_S(IJK,M)
		NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)&
		               *WW(IJK)*AXY_V(IJK) 
		    MOM_LO = V_S(IJK,M)
		ELSE
		    CONV_FAC = AVG_Y(ROP_S(IJKT,M),ROP_S(IJKTN,M),J)&
		              *WW(IJK)*AXY_V(IJK)
		    MOM_LO = V_S(IJKP,M)
		ENDIF
		MOM_HO = XSI_T(IJK)*V_S(IJKP,M)+(1.0-XSI_T(IJK))*V_S(IJK,M)
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
	      CONV_FAC = AVG_Y(ROP_S(IJKW,M),ROP_S(IJKWN,M),J)&
	                 *U(IMJK)*AYZ_V(IMJK) 
	      MOM_LO = V_S(IMJK,M)
	    ELSE
	      CONV_FAC = AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)&
		           *U(IMJK)*AYZ_V(IMJK)
	      MOM_LO = V_S(IJK,M)
	    ENDIF
	    MOM_HO = XSI_E(IMJK)*V_S(IJK,M)+(1.0-XSI_E(IMJK))*V_S(IMJK,M)
	    WEST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJKS = SOUTH_OF(IJK) 
	    IF(V(IJMK) >= ZERO)THEN
	       CONV_FAC = AVG_Y(ROP_S(IJKS,M),ROP_S(IJKC,M),JM)&
	                 *V(IJMK)*AXZ_U(IJMK) 
	       MOM_LO = V_S(IJMK,M)
	    ELSE
	       CONV_FAC = AVG_Y(ROP_S(IJKC,M),ROP_S(IJKN,M),J)&
	                  *V(IJMK)*AXZ_U(IJMK)
	       MOM_LO = V_S(IJK,M)
	    ENDIF
	    MOM_HO = XSI_N(IJMK)*V_S(IJK,M)+(1.0-XSI_N(IJMK))*V_S(IJMK,M)
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
	         CONV_FAC = AVG_Y(ROP_S(IJKB,M),ROP_S(IJKBN,M),J)&
		            *WW(IJKM)*AXY_V(IJKM) 
		 MOM_LO = V_S(IJKM,M)
	       ELSE
		 CONV_FAC = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)&
		            *WW(IJKM)*AXY_V(IJKM)
		 MOM_LO = V_S(IJK,M)
	       ENDIF
	       MOM_HO = XSI_T(IJKM)*V_S(IJK,M)+(1.0-XSI_T(IJKM))*V_S(IJKM,M)
	       BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!	    CONTRIBUTION DUE TO DEFERRED CORRECTION
!
	    B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 

      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_V_SDC 

 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_s1(A_V_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_v_s                                                      C
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
      SUBROUTINE STORE_A_V_S1(A_V_S, M, IER) 
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
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar   !//d
      USE sendrecv  !// 400       
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
!                      Solids phase 
      INTEGER          M 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, DIMENSION_M) 
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
! start loezos
      INTEGER incr   
! end loezos

      call lock_tmp_array
      call lock_xsi_array

!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1229 Communicate boundaries
      call send_recv(U_S,2)
      call send_recv(V_S,2)
      call send_recv(W_S,2)
      call send_recv(MU_S,2)
      call send_recv(AYZ_V,2)
      call send_recv(AXZ_V,2)
      call send_recv(AXY_V,2)      
      call send_recv(ROP_S,2)
      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)
!
!  Calculate convection factors
!
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
	
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J) 
      END DO 

! loezos
	incr=2		
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,incr) 


! loezos      
! update to true velocity
      IF (SHEAR) THEN
!$omp      parallel do private(IJK)
        DO IJK = 1, IJKMAX2
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
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
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
            D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKN,M),MU_S(IJKNE,M),I),J)*ODX_E(I)*AYZ_V(IJK) 
!
            A_V_S(IJK,E,M) = D_F - XSI_E(IJK)*AVG_Y(ROP_S(IJKE,M),ROP_S(IJKNE,M&
               ),J)*U(IJK)*AYZ_V(IJK) 
!
            A_V_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*AVG_Y(ROP_S(IJKC,M),&
               ROP_S(IJKN,M),J)*U(IJK)*AYZ_V(IJK) 
!
!
!           North face (i, j+1, k)
            D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK) 
            A_V_S(IJK,N,M) = D_F - XSI_N(IJK)*AVG_Y(ROP_S(IJKN,M),ROP_S(&
               NORTH_OF(IJKN),M),JP)*V(IJK)*AXZ_V(IJK) 
!
            A_V_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*AVG_Y(ROP_S(IJKC,M),&
               ROP_S(IJKN,M),J)*V(IJK)*AXZ_V(IJK) 
!
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
               D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKN,M),MU_S(IJKTN,M),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
!
               A_V_S(IJK,T,M) = D_F - XSI_T(IJK)*AVG_Y(ROP_S(IJKT,M),ROP_S(&
                  IJKTN,M),J)*WW(IJK)*AXY_V(IJK) 
!
               A_V_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*AVG_Y(ROP_S(IJKC,M),&
                  ROP_S(IJKN,M),J)*WW(IJK)*AXY_V(IJK) 
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
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWN,M),MU_S(IJKN,M),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
!
               A_V_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*AVG_Y(ROP_S(IJKW,M),&
                  ROP_S(IJKWN,M),J)*U(IMJK)*AYZ_V(IMJK) 
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
               D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK) 
!
               A_V_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*AVG_Y(ROP_S(IJKS,M),&
                  ROP_S(IJKC,M),JM)*V(IJMK)*AXZ_V(IJMK) 
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
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(&
                     IJKM) 
!
                  A_V_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*AVG_Y(ROP_S(IJKB,M&
                     ),ROP_S(IJKBN,M),J)*WW(IJKM)*AXY_V(IJKM) 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 


      call unlock_tmp_array
      call unlock_xsi_array

!//? check if COMM of A_V_S is necessary?
      
      RETURN  
      END SUBROUTINE STORE_A_V_S1 
