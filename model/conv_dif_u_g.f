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
!!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): entered conv_dif_u_g')") myPE  !//AIKEPARDBG
!    write(*,"('(PE ',I2,'): DEF_COR= ',L4)") myPE,DEF_COR
!    write(*,"('(PE ',I2,'): DISCRETIZE(3)= ',I4)") myPE,DISCRETIZE(3)
!    call mfix_exit(myPE)     !//AIKEPARDBG

      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN

            
!	IF DEFERRED CORRECTION IS USED TO SOLVE U_G
!
      IF(DEF_COR)THEN
	CALL STORE_A_U_G0(A_M(1,-3,0), IER)
        CALL STORE_A_U_GDC (A_M(1,-3,0),B_M(1,0), IER)
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
      
!// 400 1213 COMM A_M and B_M      
      CALL SEND_RECV(A_M, 2)
      CALL SEND_RECV(B_M, 2)

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): bef dif_u_is in conv_dif_u_g')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
      CALL DIF_U_IS (MU_GT, A_M, B_M, 0, IER) 

!//? Check if the following COMM is extra or not?
!// 400 1213 COMM A_M and B_M      
      CALL SEND_RECV(A_M, 2)
      CALL SEND_RECV(B_M, 2)

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft dif_u_is in conv_dif_u_g')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
            
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
!                      Face velocity
      DOUBLE PRECISION V_f
!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1225 Communicate boundaries
      call send_recv(U_G,2)
      call send_recv(V_G,2)
      call send_recv(W_G,2)
      call send_recv(MU_GT,2)
      call send_recv(AYZ_U,2)
      call send_recv(AXZ_U,2)
      call send_recv(AXY_U,2)      
      call send_recv(ROP_G,2)            

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0 
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
      
!$omp      parallel do                                                  &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,     &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE, V_f, D_f,   &
!$omp&                    IMJK, IM, IJKW,                               &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                 &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)                 
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
!// 360 1225 Check if current i,j,k resides on this PE	    
            IF(.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	  
	      
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
            V_F = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
            D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
            IF (V_F >= ZERO) THEN 
               A_U_G(IJK,E) = D_F 
               A_U_G(IPJK,W) = D_F + AVG_X(ROP_G(IJKC),ROP_G(IJKE),I)*V_F*AYZ_U&
                  (IJK) 
            ELSE 
               A_U_G(IJK,E) = D_F - AVG_X(ROP_G(IJKE),ROP_G(EAST_OF(IJKE)),IP)*&
                  V_F*AYZ_U(IJK) 
               A_U_G(IPJK,W) = D_F 
            ENDIF 
!
!
!           North face (i+1/2, j+1/2, k)
            V_F = AVG_X(V_G(IJK),V_G(IPJK),I) 
            D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
               ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
            IF (V_F >= ZERO) THEN 
               A_U_G(IJK,N) = D_F 
               A_U_G(IJPK,S) = D_F + AVG_X(ROP_G(IJKC),ROP_G(IJKE),I)*V_F*AXZ_U&
                  (IJK) 
            ELSE 
               A_U_G(IJK,N) = D_F - AVG_X(ROP_G(IJKN),ROP_G(IJKNE),I)*V_F*AXZ_U&
                  (IJK) 
               A_U_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
               V_F = AVG_X(W_G(IJK),W_G(IPJK),I) 
               D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_U_G(IJK,T) = D_F 
                  A_U_G(IJKP,B) = D_F + AVG_X(ROP_G(IJKC),ROP_G(IJKE),I)*V_F*&
                     AXY_U(IJK) 
               ELSE 
                  A_U_G(IJK,T) = D_F - AVG_X(ROP_G(IJKT),ROP_G(IJKTE),I)*V_F*&
                     AXY_U(IJK) 
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
               V_F = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
               D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_U_G(IJK,W) = D_F + AVG_X(ROP_G(IJKW),ROP_G(IJKC),IM)*V_F*&
                     AYZ_U(IMJK) 
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
               V_F = AVG_X(V_G(IJMK),V_G(IPJMK),I) 
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_U_G(IJK,S) = D_F + AVG_X(ROP_G(IJKS),ROP_G(IJKSE),I)*V_F*&
                     AXZ_U(IJMK) 
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
                  V_F = AVG_X(W_G(IJKM),W_G(IPJKM),I) 
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
                  IF (V_F >= ZERO) THEN 
                     A_U_G(IJK,B) = D_F + AVG_X(ROP_G(IJKB),ROP_G(IJKBE),I)*V_F&
                        *AXY_U(IJKM) 
                  ELSE 
                     A_U_G(IJK,B) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
!//? Check if need to COMM A_U_G ??      
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
      USE compar   !//d
      USE sendrecv     !// 400
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
!
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
!
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

! loezos
	INTEGER incr
! loezos

      call lock_tmp_array
      call lock_xsi_array


!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1225 Communicate boundaries
      call send_recv(U_G,2)
      call send_recv(V_G,2)
      call send_recv(W_G,2)
      call send_recv(MU_GT,2)
      call send_recv(AYZ_U,2)
      call send_recv(AXZ_U,2)
      call send_recv(AXY_U,2)      
      call send_recv(ROP_G,2)
      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)      
!
!  Calculate convection factors
!
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
       DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp      parallel do                                                 &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!$omp&                    IMJK, IM, IJKW,                              &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
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
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
!
		IF(U(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_X(ROP_G(IJK),ROP_G(IJKE),I)*U(IJK)*AYZ_U(IJK) 
		    MOM_LO = U_G(IJK)
		ELSE
		    CONV_FAC = AVG_X(ROP_G(IJKE),ROP_G(EAST_OF(IJKE)),IP)&
                               *U(IJK)*AYZ_U(IJK)
		    MOM_LO = U_G(IPJK)
		ENDIF
		MOM_HO = XSI_E(IJK)*U_G(IPJK)+(1.0-XSI_E(IJK))*U_G(IJK)
		EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
!            
		IF(V(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_X(ROP_G(IJKC),ROP_G(IJKE),I)*V(IJK)*AXZ_U(IJK) 
		    MOM_LO = U_G(IJK)
		ELSE
		    CONV_FAC = AVG_X(ROP_G(IJKN),ROP_G(IJKNE),I)*V(IJK)*AXZ_U(IJK)
		    MOM_LO = U_G(IJPK)
		ENDIF
		MOM_HO = XSI_N(IJK)*U_G(IJPK)+(1.0-XSI_N(IJK))*U_G(IJK)
		NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
!
            IF (DO_K) THEN 
                IJKP = KP_OF(IJK) 
                IJKT = TOP_OF(IJK) 
                IJKTE = EAST_OF(IJKT) 
		IF(WW(IJK) >= ZERO)THEN
	           CONV_FAC = AVG_X(ROP_G(IJKC),ROP_G(IJKE),I)*WW(IJK)*AXY_U(IJK) 
		   MOM_LO = U_G(IJK)
		ELSE
		   CONV_FAC = AVG_X(ROP_G(IJKT),ROP_G(IJKTE),I)*WW(IJK)*AXY_U(IJK)
		   MOM_LO = U_G(IJKP)
		ENDIF
		MOM_HO = XSI_T(IJK)*U_G(IJKP)+(1.0-XSI_T(IJK))*U_G(IJK)
		TOP_DC = CONV_FAC*(MOM_LO-MOM_HO)
		
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
	      CONV_FAC = AVG_X(ROP_G(IJKW),ROP_G(IJK),IM)*U(IMJK)*AYZ_U(IMJK) 
	      MOM_LO = U_G(IMJK)
	      
	    ELSE
	      CONV_FAC = AVG_X(ROP_G(IJK),ROP_G(IJKE),I)*U(IMJK)*AYZ_U(IMJK)
	      MOM_LO = U_G(IJK)
	      
	    ENDIF
	    MOM_HO = XSI_E(IMJK)*U_G(IJK)+(1.0-XSI_E(IMJK))*U_G(IMJK)
	    WEST_DC = CONV_FAC*(MOM_LO-MOM_HO)

!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IPJMK = IP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKSE = EAST_OF(IJKS) 
 	    IF(V(IJMK) >= ZERO)THEN
	      CONV_FAC = AVG_X(ROP_G(IJKS),ROP_G(IJKSE),I)*V(IJMK)*AXZ_U(IJMK) 
	      MOM_LO = U_G(IJMK)
	    ELSE
	      CONV_FAC = AVG_X(ROP_G(IJK),ROP_G(IJKE),I)*V(IJMK)*AXZ_U(IJMK)
	      MOM_LO = U_G(IJK)
	    ENDIF
	    MOM_HO = XSI_N(IJMK)*U_G(IJK)+(1.0-XSI_N(IJMK))*U_G(IJMK)
	    SOUTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
	    
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
	         CONV_FAC = AVG_X(ROP_G(IJKB),ROP_G(IJKBE),I)*WW(IJK)*AXY_U(IJKM) 
		 MOM_LO = U_G(IJKM)
	       ELSE
		 CONV_FAC = AVG_X(ROP_G(IJK),ROP_G(IJKE),I)*WW(IJK)*AXY_U(IJKM)
		 MOM_LO = U_G(IJK)
	       ENDIF
	       MOM_HO = XSI_T(IJKM)*U_G(IJK)+(1.0-XSI_T(IJKM))*U_G(IJKM)
	       BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
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

!//? check if COMM of B_M is necessary?      

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
!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                      XSI_t(DIMENSION_3)
!      DOUBLE PRECISION U(DIMENSION_3),&
!                      V(DIMENSION_3), WW(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

! loezos                     SHEAR VELOCITY
      INTEGER incr    
!loezos

      call lock_tmp_array
      call lock_xsi_array

!//? Check if all these COMMs are necessary, added here as fool-proof approach
!// 400 1225 Communicate boundaries
      call send_recv(U_G,2)
      call send_recv(V_G,2)
      call send_recv(W_G,2)
      call send_recv(MU_GT,2)
      call send_recv(AYZ_U,2)
      call send_recv(AXZ_U,2)
      call send_recv(AXY_U,2)      
      call send_recv(ROP_G,2)
      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)      

!
!  Calculate convection factors
!
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    


!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
	 J=J_OF(IJK)
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 


!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 

! loezos      
!update V to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)
!//SP
	 DO IJK = IJKSTART3, IJKEND3
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
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp      parallel do                                                 &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!$omp&                    IMJK, IM, IJKW,                              &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
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
            D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
!
            A_U_G(IJK,E) = D_F - XSI_E(IJK)*AVG_X(ROP_G(IJKE),ROP_G(EAST_OF(&
               IJKE)),IP)*U(IJK)*AYZ_U(IJK) 
!
            A_U_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*AVG_X(ROP_G(IJKC),ROP_G(&
               IJKE),I)*U(IJK)*AYZ_U(IJK) 
!
!
!           North face (i+1/2, j+1/2, k)
            D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
               ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
!
            A_U_G(IJK,N) = D_F - XSI_N(IJK)*AVG_X(ROP_G(IJKN),ROP_G(IJKNE),I)*V&
               (IJK)*AXZ_U(IJK) 
!
            A_U_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*AVG_X(ROP_G(IJKC),ROP_G(&
               IJKE),I)*V(IJK)*AXZ_U(IJK) 
!
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!
               D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
!
               A_U_G(IJK,T) = D_F - XSI_T(IJK)*AVG_X(ROP_G(IJKT),ROP_G(IJKTE),I&
                  )*WW(IJK)*AXY_U(IJK) 
!
               A_U_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*AVG_X(ROP_G(IJKC),ROP_G&
                  (IJKE),I)*WW(IJK)*AXY_U(IJK) 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
!
               A_U_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*AVG_X(ROP_G(IJKW),ROP_G&
                  (IJKC),IM)*U(IMJK)*AYZ_U(IMJK) 
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
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
!
               A_U_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*AVG_X(ROP_G(IJKS),ROP_G&
                  (IJKSE),I)*V(IJMK)*AXZ_U(IJMK) 
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
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
!
                  A_U_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*AVG_X(ROP_G(IJKB),&
                     ROP_G(IJKBE),I)*WW(IJKM)*AXY_U(IJKM) 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 

!//? check if COMM of A_U_G is necessary?
      
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_G1 
