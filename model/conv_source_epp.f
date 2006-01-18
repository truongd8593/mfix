!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_SOURCE_EPp(A_m, B_m, B_mmax, IER)                 C
!  Purpose: Determine convection terms for solids volume fraction      C
!           correction equation.  Master routine                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-MAR-97   C
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
      SUBROUTINE CONV_SOURCE_EPP(A_M, B_M, B_mmax, IER) 
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
      USE fldvar
      USE run
      USE geometry
      USE compar
      USE sendrecv
      Use xsi_array

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
!                      phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
  
! 
!                      maximum term in b_m expression
      DOUBLE PRECISION B_mmax(DIMENSION_3, 0:DIMENSION_M) 
  
!
      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_SOURCE_EPP0 (A_M, B_M, B_MMAX, IER) 
      ELSE 
         CALL CONV_SOURCE_EPP1 (A_M, B_M, B_MMAX, IER) 
      ENDIF 


      RETURN  
      END SUBROUTINE CONV_SOURCE_EPP 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_SOURCE_EPp0(A_m, B_m, B_MMAX, IER)                C
!  Purpose: Determine convection terms for solids volume fraction      C
!           correction equation.  First order upwinding.               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-SEP-96  C
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
      SUBROUTINE CONV_SOURCE_EPP0(A_M, B_M, B_MMAX, IER) 
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
      USE fldvar
      USE run
      USE parallel 
      USE matrix 
      USE constant
      USE physprop
      USE rxns
      USE geometry
      USE indices
      USE pgcor
      USE pscor
      USE compar
      USE sendrecv

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
! 
!                      maximum term in b_m expression
      DOUBLE PRECISION B_mmax(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      phase index 
      INTEGER          M 
! 
!                      Indices 
      INTEGER          I, J, K, IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS 
      INTEGER          IJKB, IJKT 
! 
!                      dPodEP_s(EP_s(IJK, M)) 
      DOUBLE PRECISION K_P 
! 
!                      Mass source 
      DOUBLE PRECISION Src 
! 
!                      error message 
      CHARACTER*80     LINE(1) 
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
! 
!                      terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb, bmr
!
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'
!
!
!     Calculate convection-diffusion fluxes through each of the faces
!
      M = MCP 
!
!     CHEM & ISAT begin (nan xie)
!
! Set the source terms zero
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_S_temp = SUM_R_S
         SUM_R_S = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
!$omp    parallel do                                                   &
!$omp&   private(I, J, K, IJK, IPJK, IJPK, IJKP,                &
!$omp&           IMJK, IJMK, IJKM,                               &
!$omp&           IJKE, IJKW, IJKN, IJKS, IJKT, IJKB,                   &
!$omp&           K_P, SRC, bma, bme, bmw, bmn, bms, bmt, bmb, bmr )
      DO IJK = ijkstart3, ijkend3
!// Determine whehter IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
            IJKE = EAST_OF(IJK) 
            IJKW = WEST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKS = SOUTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IJKB = BOTTOM_OF(IJK) 
!
            A_M(IJK,0,0) = ZERO 
            B_M(IJK,0) = ZERO 
            K_P = K_CP(IJK) 
!
!
!         East face (i+1/2, j, k)
            IF (U_S(IJK,M) < ZERO) THEN 
               A_M(IJK,E,0) = (ROP_S(IJKE,M)*E_E(IJK)*K_CP(IJKE)-RO_S(M)*U_S(&
                  IJK,M))*AYZ(IJK) 
!
               A_M(IJK,0,0)=A_M(IJK,0,0)+ROP_S(IJKE,M)*E_E(IJK)*K_P*AYZ(IJK) 
!
               bme = (-ROP_S(IJKE,M)*U_S(IJK,M))*AYZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) +  bme
!
            ELSE 
               A_M(IJK,E,0) = (ROP_S(IJK,M)*E_E(IJK)*K_CP(IJKE))*AYZ(IJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_E(IJK)*K_P+RO_S(M)&
                  *U_S(IJK,M))*AYZ(IJK) 
!
               bme = (-ROP_S(IJK,M)*U_S(IJK,M))*AYZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) +  bme
!
            ENDIF 
!
!
!         West face (i-1/2, j, k)
            IF (U_S(IMJK,M) > ZERO) THEN 
               A_M(IJK,W,0) = (ROP_S(IJKW,M)*E_E(IMJK)*K_CP(IJKW)+RO_S(M)*U_S(&
                  IMJK,M))*AYZ(IMJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + ROP_S(IJKW,M)*E_E(IMJK)*K_P*AYZ(&
                  IMJK) 
!
               bmw = (ROP_S(IJKW,M)*U_S(IMJK,M))*AYZ(IMJK)
               B_M(IJK,0) = B_M(IJK,0) + bmw  
!
            ELSE 
               A_M(IJK,W,0) = (ROP_S(IJK,M)*E_E(IMJK)*K_CP(IJKW))*AYZ(IMJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_E(IMJK)*K_P-RO_S(M&
                  )*U_S(IMJK,M))*AYZ(IMJK) 
!
               bmw = (ROP_S(IJK,M)*U_S(IMJK,M))*AYZ(IMJK)
               B_M(IJK,0) = B_M(IJK,0) + bmw 
!
            ENDIF 
!
!
!         North face (i, j+1/2, k)
            IF (V_S(IJK,M) < ZERO) THEN 
               A_M(IJK,N,0) = (ROP_S(IJKN,M)*E_N(IJK)*K_CP(IJKN)-RO_S(M)*V_S(&
                  IJK,M))*AXZ(IJK) 
!
               A_M(IJK,0,0)=A_M(IJK,0,0)+ROP_S(IJKN,M)*E_N(IJK)*K_P*AXZ(IJK) 
!
               bmn = (-ROP_S(IJKN,M)*V_S(IJK,M))*AXZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmn  
!
            ELSE 
               A_M(IJK,N,0) = (ROP_S(IJK,M)*E_N(IJK)*K_CP(IJKN))*AXZ(IJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_N(IJK)*K_P+RO_S(M)&
                  *V_S(IJK,M))*AXZ(IJK) 
!
               bmn = (-ROP_S(IJK,M)*V_S(IJK,M))*AXZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmn  
!
            ENDIF 
!
!
!         South face (i, j-1/2, k)
            IF (V_S(IJMK,M) > ZERO) THEN 
               A_M(IJK,S,0) = (ROP_S(IJKS,M)*E_N(IJMK)*K_CP(IJKS)+RO_S(M)*V_S(&
                  IJMK,M))*AXZ(IJMK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + ROP_S(IJKS,M)*E_N(IJMK)*K_P*AXZ(&
                  IJMK) 
!
               bms = (ROP_S(IJKS,M)*V_S(IJMK,M))*AXZ(IJMK)
               B_M(IJK,0) = B_M(IJK,0) + bms 
!
            ELSE 
               A_M(IJK,S,0) = (ROP_S(IJK,M)*E_N(IJMK)*K_CP(IJKS))*AXZ(IJMK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_N(IJMK)*K_P-RO_S(M&
                  )*V_S(IJMK,M))*AXZ(IJMK) 
!
               bms = (ROP_S(IJK,M)*V_S(IJMK,M))*AXZ(IJMK)
               B_M(IJK,0) = B_M(IJK,0) + bms 
!
            ENDIF 
!
            IF (DO_K) THEN 
!
!           Top face (i, j, k+1/2)
               IF (W_S(IJK,M) < ZERO) THEN 
                  A_M(IJK,T,0) = (ROP_S(IJKT,M)*E_T(IJK)*K_CP(IJKT)-RO_S(M)*W_S&
                     (IJK,M))*AXY(IJK) 
!
                  A_M(IJK,0,0) = A_M(IJK,0,0) + ROP_S(IJKT,M)*E_T(IJK)*K_P*AXY(&
                     IJK) 
!
                  bmt = (-ROP_S(IJKT,M)*W_S(IJK,M))*AXY(IJK)
                  B_M(IJK,0)=B_M(IJK,0) + bmt 
!
               ELSE 
                  A_M(IJK,T,0) = (ROP_S(IJK,M)*E_T(IJK)*K_CP(IJKT))*AXY(IJK) 
!
                  A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_T(IJK)*K_P+RO_S&
                     (M)*W_S(IJK,M))*AXY(IJK) 
!
                  bmt = (-ROP_S(IJK,M)*W_S(IJK,M))*AXY(IJK)
                  B_M(IJK,0) = B_M(IJK,0) + bmt 
!
               ENDIF 
!
!
!           Bottom face (i, j, k-1/2)
               IF (W_S(IJKM,M) > ZERO) THEN 
                  A_M(IJK,B,0) = (ROP_S(IJKB,M)*E_T(IJKM)*K_CP(IJKB)+RO_S(M)*&
                     W_S(IJKM,M))*AXY(IJKM) 
!
                  A_M(IJK,0,0) = A_M(IJK,0,0) + ROP_S(IJKB,M)*E_T(IJKM)*K_P*AXY&
                     (IJKM) 
!
                  bmb = (ROP_S(IJKB,M)*W_S(IJKM,M))*AXY(IJKM) 
                  B_M(IJK,0) = B_M(IJK,0) + bmb
!
               ELSE 
                  A_M(IJK,B,0) = (ROP_S(IJK,M)*E_T(IJKM)*K_CP(IJKB))*AXY(IJKM) 
!
                  A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*E_T(IJKM)*K_P-&
                     RO_S(M)*W_S(IJKM,M))*AXY(IJKM) 
!
                  bmb = (ROP_S(IJK,M)*W_S(IJKM,M))*AXY(IJKM)
                  B_M(IJK,0)=B_M(IJK,0) + bmb
!
               ENDIF 
!
            ELSE
	      bmt = zero
	      bmb = zero
            ENDIF 
!
            IF (ROP_S(IJK,M)>ZERO .AND. SUM_R_S(IJK,M)<ZERO) THEN 
               SRC = VOL(IJK)*(-SUM_R_S(IJK,M))/ROP_S(IJK,M) 
            ELSE 
               SRC = ZERO 
            ENDIF 
!
            A_M(IJK,0,0) = -(A_M(IJK,0,0)+VOL(IJK)*ODT*RO_S(M)+SRC*RO_S(M)) 
!
            bma = (ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT
	    bmr = SUM_R_S(IJK,M)*VOL(IJK)
            B_M(IJK,0) = -(B_M(IJK,0) - bma + bmr) 
            B_MMAX(IJK,0) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb), abs(bmr) ) 
!
            IF ((-A_M(IJK,0,0)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,0) = -ONE            ! Equation is undefined. 
                  B_M(IJK,0) = ZERO              ! Use existing value 
               ELSE 
!$omp             critical
                  WRITE (LINE, '(A,I6,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' A = 0 and b = ', B_M(IJK,0) 
                  CALL WRITE_ERROR ('CONV_SOURCE_EPp0', LINE, 1) 
!$omp             end critical
               ENDIF 
            ENDIF 
         ELSE 
            A_M(IJK,0,0) = -ONE 
            B_M(IJK,0) = ZERO 
         ENDIF 
      END DO 
!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_S = SUM_R_S_temp
      END IF
!     CHEM & ISAT end (nan xie)
!      
      RETURN  
      END SUBROUTINE CONV_SOURCE_EPP0 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_SOURCE_EPp1(A_m, B_m, B_MMAX, IER)                C
!  Purpose: Determine convection terms for solids volume fraction      C
!           correction equation.  Higher order scheme.                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-SEP-96  C
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
      SUBROUTINE CONV_SOURCE_EPP1(A_M, B_M, B_MMAX, IER) 
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
      USE fldvar
      USE run
      USE parallel 
      USE matrix 
      USE constant
      USE physprop
      USE rxns
      USE geometry
      USE indices
      USE pgcor
      USE pscor
      USE xsi_array
      USE vshear
      USE compar
      USE sendrecv
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
! 
!                      maximum term in b_m expression
      DOUBLE PRECISION B_mmax(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      phase index 
      INTEGER          M 
! 
!                      Indices 
      INTEGER          I, J, K, IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS 
      INTEGER          IJKB, IJKT 
!
! loezos

	INTEGER  incr
! loezos      
 
!                      dPodEP_s(EP_s(IJK, M)) 
      DOUBLE PRECISION K_P 
! 
!                      Mass source 
      DOUBLE PRECISION Src 
! 
!                      face value of ROP_s 
      DOUBLE PRECISION ROP_sf 
! 
!                      error message 
      CHARACTER*80     LINE(1) 
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
! 
!                      terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb, bmr
!
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'


      call lock_xsi_array

!
      M = MCP 
!
!  Calculate convection factors
!

! loezos
	incr=0		
! loezos

      CALL CALC_XSI (DISCRETIZE(2), ROP_S(1,M), U_S(1,M), V_S(1,M), W_S(1,M), &
         XSI_E, XSI_N, XSI_T,incr) 


! loezos
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)  
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V_S(IJK,m)=V_s(IJK,m)+VSH(IJK)	
         END IF
       END DO 
 
      END IF

! loezos

!
!
!     CHEM & ISAT begin (nan xie)
!
! Set the source terms zero
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_S_temp = SUM_R_S
         SUM_R_S = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
!
!     Calculate convection-diffusion fluxes through each of the faces
!
!$omp parallel do                                                      &   
!$omp&   private(I, J, K, IJK, IPJK, IJPK, IJKP,                &
!$omp&           IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB, &
!$omp&           K_P,ROP_SF,SRC, bma, bme, bmw, bmn, bms, bmt, bmb, bmr )
      DO IJK = ijkstart3, ijkend3
!// Determine if IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
            IJKE = EAST_OF(IJK) 
            IJKW = WEST_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IJKS = SOUTH_OF(IJK) 
            IJKT = TOP_OF(IJK)
	    IJKB = BOTTOM_OF(IJK) 
!
            A_M(IJK,0,0) = ZERO 
            B_M(IJK,0) = ZERO 
            K_P = K_CP(IJK) 
!
!
!         East face (i+1/2, j, k)
            ROP_SF = ROP_S(IJKE,M)*XSI_E(IJK) + ROP_S(IJK,M)*(ONE - XSI_E(IJK)) 
!
            A_M(IJK,E,0) = (ROP_SF*E_E(IJK)*K_CP(IJKE)-RO_S(M)*U_S(IJK,M)*XSI_E&
               (IJK))*AYZ(IJK) 
!
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_E(IJK)*K_P+RO_S(M)*U_S(IJK,&
               M)*(ONE-XSI_E(IJK)))*AYZ(IJK) 
!
            bme = (-ROP_SF*U_S(IJK,M))*AYZ(IJK)
            B_M(IJK,0) = B_M(IJK,0) + bme 
!
!
!
!         West face (i-1/2, j, k)
            ROP_SF=ROP_S(IJK,M)*XSI_E(IMJK)+ROP_S(IJKW,M)*(ONE-XSI_E(IMJK)) 
!
            A_M(IJK,W,0) = (ROP_SF*E_E(IMJK)*K_CP(IJKW)+RO_S(M)*U_S(IMJK,M)*(&
               ONE-XSI_E(IMJK)))*AYZ(IMJK) 
!
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_E(IMJK)*K_P-RO_S(M)*U_S(&
               IMJK,M)*XSI_E(IMJK))*AYZ(IMJK) 
!
            bmw = (ROP_SF*U_S(IMJK,M))*AYZ(IMJK)
            B_M(IJK,0) = B_M(IJK,0) +  bmw
!
!
!         North face (i, j+1/2, k)
            ROP_SF = ROP_S(IJKN,M)*XSI_N(IJK) + ROP_S(IJK,M)*(ONE - XSI_N(IJK)) 
!
            A_M(IJK,N,0) = (ROP_SF*E_N(IJK)*K_CP(IJKN)-RO_S(M)*V_S(IJK,M)*XSI_N&
               (IJK))*AXZ(IJK) 
!
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_N(IJK)*K_P+RO_S(M)*V_S(IJK,&
               M)*(ONE-XSI_N(IJK)))*AXZ(IJK) 
!
            bmn = (-ROP_SF*V_S(IJK,M))*AXZ(IJK) 
            B_M(IJK,0) = B_M(IJK,0) + bmn
!
!
!         South face (i, j-1/2, k)
            ROP_SF=ROP_S(IJK,M)*XSI_N(IJMK)+ROP_S(IJKS,M)*(ONE-XSI_N(IJMK)) 
!
            A_M(IJK,S,0) = (ROP_SF*E_N(IJMK)*K_CP(IJKS)+RO_S(M)*V_S(IJMK,M)*(&
               ONE-XSI_N(IJMK)))*AXZ(IJMK) 
!
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_N(IJMK)*K_P-RO_S(M)*V_S(&
               IJMK,M)*XSI_N(IJMK))*AXZ(IJMK) 
!
            bms = (ROP_SF*V_S(IJMK,M))*AXZ(IJMK)
            B_M(IJK,0) = B_M(IJK,0) + bms
!
            IF (DO_K) THEN 
!
!           Top face (i, j, k+1/2)
               ROP_SF=ROP_S(IJKT,M)*XSI_T(IJK)+ROP_S(IJK,M)*(ONE-XSI_T(IJK)) 
!
               A_M(IJK,T,0) = (ROP_SF*E_T(IJK)*K_CP(IJKT)-RO_S(M)*W_S(IJK,M)*&
                  XSI_T(IJK))*AXY(IJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_T(IJK)*K_P+RO_S(M)*W_S(&
                  IJK,M)*(ONE-XSI_T(IJK)))*AXY(IJK) 
!
               bmt = (-ROP_SF*W_S(IJK,M))*AXY(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmt
!
!
!           Bottom face (i, j, k-1/2)
               ROP_SF = ROP_S(IJK,M)*XSI_T(IJKM) + ROP_S(IJKB,M)*(ONE - XSI_T(&
                  IJKM)) 
!
               A_M(IJK,B,0) = (ROP_SF*E_T(IJKM)*K_CP(IJKB)+RO_S(M)*W_S(IJKM,M)*&
                  (ONE-XSI_T(IJKM)))*AXY(IJKM) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_T(IJKM)*K_P-RO_S(M)*W_S(&
                  IJKM,M)*XSI_T(IJKM))*AXY(IJKM) 
!
               bmb = (ROP_SF*W_S(IJKM,M))*AXY(IJKM)
               B_M(IJK,0) = B_M(IJK,0) + bmb
!
            ELSE
	      bmt = zero
	      bmb = zero
            ENDIF 
!
            IF (ROP_S(IJK,M)>ZERO .AND. SUM_R_S(IJK,M)<ZERO) THEN 
               SRC = VOL(IJK)*(-SUM_R_S(IJK,M))/ROP_S(IJK,M) 
            ELSE 
               SRC = ZERO 
            ENDIF 
!
            A_M(IJK,0,0) = -(A_M(IJK,0,0)+VOL(IJK)*ODT*RO_S(M)+SRC*RO_S(M)) 
!
            bma = (ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT
	    bmr = SUM_R_S(IJK,M)*VOL(IJK)
            B_M(IJK,0) = -(B_M(IJK,0)- bma + bmr) 
            B_MMAX(IJK,0) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb), abs(bmr) ) !
            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,0) = -ONE            ! Equation is undefined. 
                  B_M(IJK,0) = ZERO              ! Use existing value 
               ELSE 
!$omp             critical
                  WRITE (LINE(1), '(A,I6,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' A = 0 and b = ', B_M(IJK,0) 
!//SP Having problem to compile this statement on SGI
                  CALL WRITE_ERROR ('CONV_SOURCE_EPp1', LINE, 1) 
!$omp             end critical
               ENDIF 
            ENDIF 
         ELSE 
            A_M(IJK,0,0) = -ONE 
            B_M(IJK,0) = ZERO 
         ENDIF 
      END DO
!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_S = SUM_R_S_temp
      END IF
!     CHEM & ISAT end (nan xie)
!
! loezos
      IF (SHEAR) THEN
!$omp parallel do private(IJK)  
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V_S(IJK,m)=V_s(IJK,m)-VSH(IJK)	
         END IF
       END DO 
      END IF
! loezos       
      call unlock_xsi_array
      
      
      RETURN  
      END SUBROUTINE CONV_SOURCE_EPP1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 360 Check if ijk falls within 1 ghost layer

