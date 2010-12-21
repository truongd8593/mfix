!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!  CONV_DIF_Phi(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E, Flux_N, Flux_T, M, A_m, B_m, IER)    C
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
      SUBROUTINE CONV_DIF_PHI(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
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
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
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
	  CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	  if (DISC > 1) CALL CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	ELSE
!
!	NO DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	  IF (DISC == 0) THEN                        
            CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	  ELSE
            CALL CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
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
!   CONV_DIF_Phi0(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
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
      SUBROUTINE CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
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
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
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
!$omp      parallel do                                              &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!$omp&             IMJK, IM, IJKW,                                  &
!$omp&             IJMK, JM, IJKS,                                  &
!$omp&             IJKM, KM,  IJKB)                     
      DO IJK = ijkstart3, ijkend3
!
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IPJK)) THEN
                  D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*DY(J)*DZ(K)
               ENDIF
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

            IF (V_F >= ZERO) THEN 
               A_M(IJK,E,M) = D_F 
               A_M(IPJK,W,M) = D_F + FLUX_E(IJK) 
            ELSE 
               A_M(IJK,E,M) = D_F - FLUX_E(IJK) 
               A_M(IPJK,W,M) = D_F 
            ENDIF 
!
!
!           North face (i, j+1/2, k)
            V_F = VF(IJK) 
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJPK)) THEN
                  D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*DX(I)*DZ(K)
               ENDIF
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF (V_F >= ZERO) THEN 
               A_M(IJK,N,M) = D_F 
               A_M(IJPK,S,M) = D_F + FLUX_N(IJK) 
            ELSE 
               A_M(IJK,N,M) = D_F - FLUX_N(IJK) 
               A_M(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               V_F = WF(IJK) 
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJKP)) THEN
                     D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*DX(I)*DY(J)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,T,M) = D_F 
                  A_M(IJKP,B,M) = D_F + FLUX_T(IJK) 
               ELSE 
                  A_M(IJK,T,M) = D_F - FLUX_T(IJK) 
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IMJK)) THEN
                     D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*DY(J)*DZ(K)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,W,M) = D_F + FLUX_E(IMJK) 
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJMK)) THEN
                     D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*DX(I)*DZ(K)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,S,M) = D_F + FLUX_N(IJMK) 
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
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_TREATMENT_AT(IJK)) THEN
                     IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJKM)) THEN
                        D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*DX(I)*DY(J)
                     ENDIF
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF (V_F >= ZERO) THEN 
                     A_M(IJK,B,M) = D_F + FLUX_T(IJKM) 
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
!    CONV_DIF_Phi_DC(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
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
      SUBROUTINE CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
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
      USE sendrecv3
      USE indices
      Use xsi_array
      Use tmp_array
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
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
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
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4

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
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION  EAST_DC
	DOUBLE PRECISION  WEST_DC
	DOUBLE PRECISION  NORTH_DC
	DOUBLE PRECISION  SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
!
!
!---------------------------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'function3.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_xsi_array
      call lock_tmp4_array
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
            TMP4(IJK4) = PHI(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

! loezos
	incr=0		
! loezos

      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                              &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!$omp&             IMJK, IJKW,                                  &
!$omp&             IJMK, IJKS,                                  &
!$omp&             IJKM, IJKB, PHI_HO, PHI_LO,        &
!$omp&             EAST_DC, WEST_DC, NORTH_DC, SOUTH_DC, TOP_DC, BOTTOM_DC)                     
!
      DO IJK = ijkstart3, ijkend3
!
! Determine whether IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
!
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
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
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k)
!
		V_F = UF(IJK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IPJK), PHI(IJK), & 
                            PHI(IMJK), PHI(IM_OF(IMJK)))
		ELSE
		   PHI_LO = PHI(IPJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IPJK), & 
                            PHI(IP_OF(IPJK)), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
                      PHI_HO = XSI_E(IJK)*PHI(IPJK)+(1.0-XSI_E(IJK))*PHI(IJK)
		EAST_DC = FLUX_E(IJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k)
!
		V_F = VF(IJK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJPK), PHI(IJK), & 
                            PHI(IJMK), PHI(JM_OF(IJMK)))
		ELSE
		   PHI_LO = PHI(IJPK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJPK), & 
                            PHI(JP_OF(IJPK)), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		     PHI_HO = XSI_N(IJK)*PHI(IJPK)+(1.0-XSI_N(IJK))*PHI(IJK)
		NORTH_DC = FLUX_N(IJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1/2)
!
              IF (DO_K) THEN
                IJKP = KP_OF(IJK) 
                IJKT = TOP_OF(IJK)
	        V_F = WF(IJK)
		IF(V_F >= ZERO)THEN
                   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJKP),  PHI(IJK), &
                            PHI(IJKM), PHI(KM_OF(IJKM)))
	        ELSE
		   PHI_LO = PHI(IJKP)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKP),  &
                            PHI(KP_OF(IJKP)), TMP4(KPPP4))
	        ENDIF
                IF (.NOT. FPFOI ) &
                     PHI_HO = XSI_T(IJK)*PHI(IJKP)+(1.0-XSI_T(IJK))*PHI(IJK)
                TOP_DC = FLUX_T(IJK)*(PHI_LO - PHI_HO)
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
		   PHI_LO = PHI(IMJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IMJK), &
                            PHI(IM_OF(IMJK)), TMP4(IMMM4))
		ELSE
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IMJK), PHI(IJK), &
                            PHI(IPJK), PHI(IP_OF(IPJK)))
                ENDIF
                IF (.NOT. FPFOI ) &
		      PHI_HO = XSI_E(IMJK)*PHI(IJK)+(ONE-XSI_E(IMJK))*PHI(IMJK)
		WEST_DC = FLUX_E(IMJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k)
!
            	IJMK = JM_OF(IJK) 
            	IJKS = SOUTH_OF(IJK)
		V_F = VF(IJMK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJMK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJMK), & 
                            PHI(JM_OF(IJMK)), TMP4(JMMM4))
		ELSE
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJMK), PHI(IJK), & 
                            PHI(IJPK), PHI(JP_OF(IJPK)))
                ENDIF
                IF (.NOT. FPFOI ) &
            	      PHI_HO = XSI_N(IJMK)*PHI(IJK)+(ONE-XSI_N(IJMK))*PHI(IJMK)
		SOUTH_DC = FLUX_N(IJMK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k-1/2)
              IF (DO_K) THEN 
                 IJKM = KM_OF(IJK) 
                 IJKB = BOTTOM_OF(IJK)
		 V_F = WF(IJKM)
		 IF(V_F >= ZERO)THEN
                   PHI_LO = PHI(IJKM)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKM), &
                            PHI(KM_OF(IJKM)), TMP4(KMMM4))
		 ELSE
                   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJKM), PHI(IJK), &
                            PHI(IJKP), PHI(KP_OF(IJKP)))
                 ENDIF
                 IF (.NOT. FPFOI ) &
                      PHI_HO = XSI_T(IJKM)*PHI(IJK)+(1.0-XSI_T(IJKM))*PHI(IJKM)
		 BOTTOM_DC = FLUX_T(IJKM)*(PHI_LO - PHI_HO)
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
      call unlock_tmp4_array
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
!    CONV_DIF_Phi1(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
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
      SUBROUTINE CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
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
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
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
!$omp      parallel do                                               &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,          &
!$omp&             IJKP, IJKT,    D_f,                        &
!$omp&             IMJK, IM, IJKW,                                   &
!$omp&             IJMK, JM,  IJKS,                                  &
!$omp&             IJKM, KM,  IJKB )                      
!
!
!
      DO IJK = ijkstart3, ijkend3
!
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IPJK)) THEN
                  D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*DY(J)*DZ(K)
               ENDIF
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_M(IJK,E,M) = D_F - XSI_E(IJK)*FLUX_E(IJK) 
!
            A_M(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*FLUX_E(IJK) 
!
!
!           North face (i, j+1/2, k)
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJPK)) THEN
                  D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*DX(I)*DZ(K)
               ENDIF
            ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            A_M(IJK,N,M) = D_F - XSI_N(IJK)*FLUX_N(IJK) 
!
            A_M(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*FLUX_N(IJK) 
!
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
!
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJKP)) THEN
                     D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*DX(I)*DY(J)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
               A_M(IJK,T,M) = D_F - XSI_T(IJK)*FLUX_T(IJK) 
!
               A_M(IJKP,B,M)=D_F+(ONE-XSI_T(IJK))*FLUX_T(IJK) 
            ENDIF 
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IMJK)) THEN
                     D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*DY(J)*DZ(K)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
               A_M(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*FLUX_E(IMJK) 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJMK)) THEN
                     D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*DX(I)*DZ(K)
                  ENDIF
               ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
               A_M(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*FLUX_N(IJMK) 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
!
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(CUT_TREATMENT_AT(IJK)) THEN
                     IF(CUT_CELL_AT(IJK).AND.BLOCKED_CELL_AT(IJKM)) THEN
                        D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*DX(I)*DY(J)
                     ENDIF
                  ENDIF
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
                  A_M(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*FLUX_T(IJKM) 
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

!       Limit I1, I2 and all to local processor first ghost layer

	    IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)

            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)

            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)

            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)

            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)

            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

!     End of limiting to the first ghost cells of the processor....
            DO K = K1, K2 
               DO J = J1, J2 
                  DO I = I1, I2 
                     IJK = FUNIJK(I,J,K) 
!
                     SELECT CASE (TRIM(IS_PLANE(L)))  
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
