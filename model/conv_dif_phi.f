!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!  CONV_DIF_Phi(Phi, Dif, Disc, Uf, Vf, Wf, ROPf, M, A_m, B_m, IER)    C
!  Purpose: Determine convection diffusion terms for a sclar phi       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C

!  The diffusion at the flow boundaries is prevented by setting the 
!  diffusion coefficients at boundary cells to zero and then using a 
!  harmonic average to calculate the boundary diffusivity.  The value
!  diffusivities at the boundaries are checked in check_data_30.  Ensure
!  that harmonic avergaing is used in this routine. 
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
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
      SUBROUTINE CONV_DIF_PHI(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER) 
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
      USE run 
      USE geometry
      USE compar
      USE sendrecv
      Use xsi_array
      USE mpi_utility
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Macroscopic density
      DOUBLE PRECISION ROPf(DIMENSION_3)
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER

!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	IF(DEF_COR)THEN
	  CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER)
	  CALL CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER)
	ELSE
!
!	NO DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	  IF (DISC == 0) THEN                        
            CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER)
	  ELSE
            CALL CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER) 
          ENDIF
	ENDIF 
	
        CALL DIF_PHI_IS (DIF, A_M, B_M, M, IER)

        RETURN  
      END SUBROUTINE CONV_DIF_PHI 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!   CONV_DIF_Phi0(Phi, Dif, Disc, Uf, Vf, Wf, ROPf, M, A_m, B_m, IER)
!  Purpose: Determine convection diffusion terms for Phi balance       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
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
      SUBROUTINE CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER) 
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
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Macroscopic density
      DOUBLE PRECISION ROPf(DIMENSION_3)
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB
!
!                      Face velocity
      DOUBLE PRECISION V_f
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!!$omp      parallel do                                              &
!!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!!$omp&             IMJK, IM, IJKW,                                  &
!!$omp&             IJMK, JM, IJKS,                                  &
!!$omp&             IJKM, KM,  IJKB)                     
      DO IJK = ijkstart3, ijkend3
!
!// Determine if IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)

!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j, k)
            V_F = UF(IJK) 
            D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
            IF (V_F >= ZERO) THEN 
               A_M(IJK,E,M) = D_F 
               A_M(IPJK,W,M) = D_F + ROPF(IJK)*V_F*AYZ(IJK) 
            ELSE 
               A_M(IJK,E,M) = D_F - ROPF(IJKE)*V_F*AYZ(IJK) 
               A_M(IPJK,W,M) = D_F 
            ENDIF 
!
!
!           North face (i, j+1/2, k)
            V_F = VF(IJK) 
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
            IF (V_F >= ZERO) THEN 
               A_M(IJK,N,M) = D_F 
               A_M(IJPK,S,M) = D_F + ROPF(IJK)*V_F*AXZ(IJK) 
            ELSE 
               A_M(IJK,N,M) = D_F - ROPF(IJKN)*V_F*AXZ(IJK) 
               A_M(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               V_F = WF(IJK) 
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,T,M) = D_F 
                  A_M(IJKP,B,M) = D_F + ROPF(IJK)*V_F*AXY(IJK) 
               ELSE 
                  A_M(IJK,T,M) = D_F - ROPF(IJKT)*V_F*AXY(IJK) 
                  A_M(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               V_F = UF(IMJK) 
               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,W,M) = D_F + ROPF(IJKW)*V_F*AYZ(IMJK) 
               ELSE 
                  A_M(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
               V_F = VF(IJMK) 
               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,S,M) = D_F + ROPF(IJKS)*V_F*AXZ(IJMK) 
               ELSE 
                  A_M(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  V_F = WF(IJKM) 
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX_E(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
                  IF (V_F >= ZERO) THEN 
                     A_M(IJK,B,M) = D_F + ROPF(IJKB)*V_F*AXY(IJKM) 
                  ELSE 
                     A_M(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF
      END DO 
!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!    CONV_DIF_Phi_DC(Phi, Dif, Disc, Uf, Vf, Wf, ROPf, M, A_m, B_m, IER)
!  Purpose: TO USE DEFERRED CORRECTION IN SOLVING THE SCALAR TRANSPORT C
!  EQN. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER SPECIFIED   C
!  HIGH ORDER METHOD TO SOLVE FOR THE SCALAR PHI.
!  See source_Phi                                                  C
!                                                                      C
!  Author: C. GUENTHER                                Date: 1-ARP-99   C
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
      SUBROUTINE CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER) 
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
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      Use xsi_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Macroscopic density
      DOUBLE PRECISION ROPf(DIMENSION_3)
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IJKW
      INTEGER          IJMK, IJKS
      INTEGER          IJKM, IJKB

! loezos
      INTEGER  incr
! loezos

!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!	FACE VELOCITY
	DOUBLE PRECISION V_F
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION PHI_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION PHI_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION CONV_FAC
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION  EAST_DC
	DOUBLE PRECISION  WEST_DC
	DOUBLE PRECISION  NORTH_DC
	DOUBLE PRECISION  SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
!
!
!
!
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_xsi_array
!
!  Calculate convection factors
!
!

! loezos
	incr=0		
! loezos

      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!!$omp      parallel do                                              &
!!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!!$omp&             IMJK, IJKW,                                  &
!!$omp&             IJMK, IJKS,                                  &
!!$omp&             IJKM, IJKB, PHI_HO, PHI_LO, CONV_FAC,       &
!!$omp&             EAST_DC, WEST_DC, NORTH_DC, SOUTH_DC, TOP_DC, BOTTOM_DC)                     
!
      DO IJK = ijkstart3, ijkend3
!
!//  Determine whether IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k)
!
		V_F = UF(IJK)
		IF(V_F >= ZERO)THEN
			CONV_FAC = ROPF(IJK)*V_F*AYZ(IJK)
			PHI_LO = PHI(IJK)
		ELSE
			CONV_FAC = ROPF(IJKE)*V_F*AYZ(IJK)
			PHI_LO = PHI(IPJK)
		ENDIF
		PHI_HO = XSI_E(IJK)*PHI(IPJK)+(1.0-XSI_E(IJK))*PHI(IJK)
		EAST_DC = CONV_FAC*(PHI_LO-PHI_HO)
!
!
!           DEFERRED CORECTION COTRIBUTION AT THE North face (i, j+1/2, k)
!
		V_F = VF(IJK)
		IF(V_F >= ZERO)THEN
			CONV_FAC = ROPF(IJK)*V_F*AXZ(IJK)
			PHI_LO = PHI(IJK)
		ELSE
			CONV_FAC = ROPF(IJKN)*V_F*AXZ(IJK)
			PHI_LO = PHI(IJPK)
		ENDIF
		PHI_HO = XSI_N(IJK)*PHI(IJPK)+(1.0-XSI_N(IJK))*PHI(IJK)
		NORTH_DC = CONV_FAC*(PHI_LO-PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1/2)
!
            	IF (DO_K) THEN
                  IJKP = KP_OF(IJK) 
                  IJKT = TOP_OF(IJK)
	          V_F = WF(IJK)
		  IF(V_F >= ZERO)THEN
		  	CONV_FAC = ROPF(IJK)*V_F*AXY(IJK)
                  	PHI_LO = PHI(IJK) 
	          ELSE
		  	CONV_FAC = ROPF(IJKT)*V_F*AXY(IJK)
		  	PHI_LO = PHI(IJKP)
	          ENDIF
                  PHI_HO = XSI_T(IJK)*PHI(IJKP)+(1.0-XSI_T(IJK))*PHI(IJK)
		  TOP_DC = CONV_FAC*(PHI_LO-PHI_HO)
	    	ELSE
		  TOP_DC = ZERO
	    
            	ENDIF
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k)
!
	    	IMJK = IM_OF(IJK)
	    	IJKW = WEST_OF(IJK)
	    	V_F = UF(IMJK)
	    	IF(V_F >= ZERO)THEN
			CONV_FAC = ROPF(IJKW)*V_F*AYZ(IMJK)
			PHI_LO = PHI(IMJK)
		ELSE
			CONV_FAC = ROPF(IJK)*V_F*AYZ(IMJK)
			PHI_LO = PHI(IJK)
		ENDIF
		PHI_HO = XSI_E(IMJK)*PHI(IJK)+(ONE-XSI_E(IMJK))*PHI(IMJK)
		WEST_DC = CONV_FAC*(PHI_LO-PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k)
!
            	IJMK = JM_OF(IJK) 
            	IJKS = SOUTH_OF(IJK)
		V_F = VF(IJMK)
		IF(V_F >= ZERO)THEN
			CONV_FAC = ROPF(IJKS)*V_F*AXZ(IJMK)
			PHI_LO = PHI(IJMK)
		ELSE
			CONV_FAC = ROPF(IJK)*V_F*AXZ(IJMK)
			PHI_LO = PHI(IJK)
		ENDIF
            	PHI_HO = XSI_N(IJMK)*PHI(IJK)+(ONE-XSI_N(IJMK))*PHI(IJMK)
		SOUTH_DC = CONV_FAC*(PHI_LO-PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k-1/2)
            	IF (DO_K) THEN 
                 IJKM = KM_OF(IJK) 
                 IJKB = BOTTOM_OF(IJK)
		 V_F = WF(IJKM)
		 IF(V_F >= ZERO)THEN
			CONV_FAC = ROPF(IJKB)*V_F*AXY(IJKM)
                 	PHI_LO = PHI(IJKM)
		 ELSE
			CONV_FAC = ROPF(IJK)*V_F*AXY(IJKM)
                 	PHI_LO = PHI(IJK)
		 ENDIF
                 PHI_HO = XSI_T(IJKM)*PHI(IJK)+(1.0-XSI_T(IJKM))*PHI(IJKM)
		 BOTTOM_DC = CONV_FAC*(PHI_LO-PHI_HO)
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
      call unlock_xsi_array
!
!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI_DC 
!
!

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!    CONV_DIF_Phi1(Phi, Dif, Disc, Uf, Vf, Wf, ROPf, M, A_m, B_m, IER)
!  Purpose: Determine convection diffusion terms for gas energy eq Phi C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_Phi                                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
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
      SUBROUTINE CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,ROPF,M,A_M,B_M,IER) 
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
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE vshear
      Use xsi_array

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Macroscopic density
      DOUBLE PRECISION ROPf(DIMENSION_3)
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB
! start loezos
      INTEGER          I1, J1
      INTEGER incr
! end loezos

!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                       XSI_t(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      call lock_xsi_array
!
!  Calculate convection factors
!
!

! loezos
	incr=0		
! loezos	
 
      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T,incr) 

! loezos
!update V to true velocity      

      IF (SHEAR) THEN
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   VF(IJK)=VF(IJK)+VSH(IJK)	
          END IF
        END DO
      END IF
! loezos

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!!$omp      parallel do                                               &
!!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,          &
!!$omp&             IJKP, IJKT,    D_f,                        &
!!$omp&             IMJK, IM, IJKW,                                   &
!!$omp&             IJMK, JM,  IJKS,                                  &
!!$omp&             IJKM, KM,  IJKB )                      
!
!
!
      DO IJK = ijkstart3, ijkend3
!
!//  Determine whether IJK falls within 1 ghost layer........
       I = I_OF(IJK) 
       J = J_OF(IJK) 
       K = K_OF(IJK) 
!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j, k)
            D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
!
            A_M(IJK,E,M) = D_F - XSI_E(IJK)*ROPF(IJKE)*UF(IJK)*AYZ(IJK) 
!
            A_M(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*ROPF(IJK)*UF(IJK)*AYZ(IJK) 
!
!
!           North face (i, j+1/2, k)
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!
            A_M(IJK,N,M) = D_F - XSI_N(IJK)*ROPF(IJKN)*VF(IJK)*AXZ(IJK) 
!
            A_M(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*ROPF(IJK)*VF(IJK)*AXZ(IJK) 
!
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
!
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
!
               A_M(IJK,T,M) = D_F - XSI_T(IJK)*ROPF(IJKT)*WF(IJK)*AXY(IJK) 
!
               A_M(IJKP,B,M)=D_F+(ONE-XSI_T(IJK))*ROPF(IJK)*WF(IJK)*AXY(IJK) 
            ENDIF 
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
!
               A_M(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*ROPF(IJKW)*UF(IMJK)*AYZ&
                  (IMJK) 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
!
               A_M(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*ROPF(IJKS)*VF(IJMK)*AXZ&
                  (IJMK) 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
!
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX_E(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
!
                  A_M(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*ROPF(IJKB)*WF(IJKM)*&
                     AXY(IJKM) 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 

! loezos 
       IF (SHEAR) THEN

	 DO IJK = ijkstart3, ijkend3
          IF (FLUID_AT(IJK)) THEN  	 
	   VF(IJK)=VF(IJK)-VSH(IJK)	
	  END IF
         END DO 	
        END IF
! loezos      
      call unlock_xsi_array

!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI1 
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DIF_phi_IS(Dif, A_m, B_m, M, IER)                      C
!  Purpose: Remove diffusive fluxes across internal surfaces.          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
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
      SUBROUTINE DIF_PHI_IS(DIF, A_M, B_M, M, IER) 
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
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE output
      USE is
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Internal surface
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                       IJKE, IJKN, IJKT, IPJK, IJPK, IJKP
!
!                      Solids phase
      INTEGER          M
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
! Make user defined internal surfaces non-conducting
!
      DO L = 1, DIMENSION_IS 
         IF (IS_DEFINED(L)) THEN 
            I1 = IS_I_W(L) 
            I2 = IS_I_E(L) 
            J1 = IS_J_S(L) 
            J2 = IS_J_N(L) 
            K1 = IS_K_B(L) 
            K2 = IS_K_T(L) 

!// Limit I1, I2 and all to local processor first ghost layer
	    IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)


            DO K = K1, K2 
               DO J = J1, J2 
                  DO I = I1, I2 
                     IJK = FUNIJK(I,J,K) 
!
                     SELECT CASE (IS_PLANE(L))  
                     CASE ('E')  
                        IJKE = EAST_OF(IJK) 
                        IPJK = IP_OF(IJK) 
!
                        D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
!
                        A_M(IJK,E,M) = A_M(IJK,E,M) - D_F 
                        A_M(IPJK,W,M) = A_M(IPJK,W,M) - D_F 
!
                     CASE ('N')  
                        IJKN = NORTH_OF(IJK) 
                        IJPK = JP_OF(IJK) 
!
                        D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!
                        A_M(IJK,N,M) = A_M(IJK,N,M) - D_F 
                        A_M(IJPK,S,M) = A_M(IJPK,S,M) - D_F 
!
                     CASE ('T')  
                        IF (DO_K) THEN 
                           IJKT = TOP_OF(IJK) 
                           IJKP = KP_OF(IJK) 
!
                           D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*&
                              AXY(IJK) 
!
                           A_M(IJK,T,M) = A_M(IJK,T,M) - D_F 
                           A_M(IJKP,B,M) = A_M(IJKP,B,M) - D_F 
!
                        ENDIF 
                     CASE DEFAULT 
!
                     END SELECT 
                  END DO 
               END DO 
            END DO 
         ENDIF 
      END DO 
!
      RETURN  
      END SUBROUTINE DIF_PHI_IS 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 300 Limit I1, I2 and all to local processor first ghost layer
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
