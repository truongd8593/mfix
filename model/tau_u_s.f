!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_U_s(TAU_U_s, IER)                             C
!  Purpose: Cross terms in the gradient of stress in U_s momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
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
      SUBROUTINE CALC_TAU_U_S(TAU_U_S, IER) 
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
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE vshear
      USE compar        !//d
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
!                      TAU_U_s 
      DOUBLE PRECISION TAU_U_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Indices 
      INTEGER          I, J, JM, K, KM, IJK, IJKE, IPJK, IP, IMJK, IJKN,& 
                       IJKNE, IJKS, IJKSE, IPJMK, IJMK, IJKT, IJKTE,& 
                       IJKB, IJKBE, IJKM, IPJKM 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPGA 
! 
!                      Average density 
      DOUBLE PRECISION ROPGA 
! 
!                      Average viscosity 
      DOUBLE PRECISION MUGA 
! 
!                      Average EP_s*viscosity 
      DOUBLE PRECISION EPMU_ste, EPMU_sbe, EPMUGA 
! 
!                      Average dW/Xdz 
      DOUBLE PRECISION dWoXdz 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION Vtzb 
! 
!                      error message 
      CHARACTER*80     LINE 
! 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
! loezos
       INTEGER I1,J1    
! loezos
!


      DO M = 1, MMAX 

! loezos
! update to true velocity
      IF (SHEAR) THEN        
!$omp  parallel do private(IJK)
!//SP
	 DO IJK = IJKSTART3, IJKEND3
          IF (FLUID_AT(IJK)) THEN 
	   V_S(IJK,m)=V_S(IJK,m)+VSH(IJK)
          END IF
        END DO     
      END IF
! loezos

	
!!$omp  parallel do private( IJK, I, IJKE, EPGA, IP, J, JM, K, KM,  &
!!$omp&  IPJK,IMJK,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJMK,IJKT,IJKTE,  &
!!$omp&  IJKB,IJKBE,IJKM,IPJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ, EPMU_STE,EPMU_SBE, &
!!$omp&  EPMUGA,DWOXDZ,VTZB ) &
!!$omp&  schedule(static)

!//SP
      DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            EPGA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
            IF ( .NOT.SIP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN 
               IP = IP1(I) 
               J = J_OF(IJK) 
               JM = JM1(J) 
               K = K_OF(IJK) 
               KM = KM1(K) 
               IPJK = IP_OF(IJK) 
               IMJK = IM_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               IJKNE = EAST_OF(IJKN) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
               IPJMK = JM_OF(IPJK) 
               IJMK = JM_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB) 
               IJKM = KM_OF(IJK) 
               IPJKM = IP_OF(IJKM) 
!
!       Surface forces
!
!         bulk viscosity term
               SBV = (LAMBDA_S(IJKE,M)*TRD_S(IJKE,M)-LAMBDA_S(IJK,M)*TRD_S(IJK,&
                  M))*AYZ(IJK) 
!
!         shear stress terms
               SSX = MU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*ODX(IP)*AYZ_U(IJK)&
                   - MU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*ODX(I)*AYZ_U(IMJK) 
               SSY = AVG_X_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
                  IJKE,M),MU_S(IJKNE,M),J),I)*(V_S(IPJK,M)-V_S(IJK,M))*ODX_E(I)&
                  *AXZ_U(IJK) - AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
                  AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*(V_S(IPJMK,M)-V_S(&
                  IJMK,M))*ODX_E(I)*AXZ_U(IJMK) 
               EPMU_STE = AVG_X_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(&
                  MU_S(IJKE,M),MU_S(IJKTE,M),K),I) 
               EPMU_SBE = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),AVG_Z_H(&
                  MU_S(IJKBE,M),MU_S(IJKE,M),KM),I) 
               SSZ = EPMU_STE*(W_S(IPJK,M)-W_S(IJK,M))*ODX_E(I)*AXY_U(IJK) - &
                  EPMU_SBE*(W_S(IPJKM,M)-W_S(IJKM,M))*ODX_E(I)*AXY_U(IJKM) 
!
!
!         Special terms for cylindrical coordinates
               IF (CYLINDRICAL) THEN 
!
!           modify Ssz
                  SSZ = SSZ - (EPMU_STE*(HALF*(W_S(IPJK,M)+W_S(IJK,M))*OX_E(I))&
                     *AXY_U(IJK)-EPMU_SBE*(HALF*(W_S(IPJKM,M)+W_S(IJKM,M))*OX_E&
                     (I))*AXY_U(IJKM)) 
!
!           Tau_zz/X
                  EPMUGA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                  DWOXDZ = HALF*((W_S(IJK,M)-W_S(IJKM,M))*OX(I)*ODZ(K)+(W_S(&
                     IPJK,M)-W_S(IPJKM,M))*OX(IP)*ODZ(K)) 
                  VTZB = -2.*EPMUGA*OX_E(I)*DWOXDZ 
               ELSE 
                  VTZB = ZERO 
               ENDIF 
!
!         Add the terms
               TAU_U_S(IJK,M) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK) 
            ELSE 
               TAU_U_S(IJK,M) = ZERO 
            ENDIF 
         END DO 
! loezos 
       IF (SHEAR) THEN
!$omp  parallel do private(IJK) 
!//SP
	 DO IJK = IJKSTART3, IJKEND3
          IF (FLUID_AT(IJK)) THEN  	 
	   V_S(IJK,m)=V_S(IJK,m)-VSH(IJK)	
	  END IF
         END DO	
        END IF
! loezos      

      END DO 
      RETURN  
      END SUBROUTINE CALC_TAU_U_S 
