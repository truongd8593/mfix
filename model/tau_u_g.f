!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_U_g(TAU_U_g, IER)                             C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
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
      SUBROUTINE CALC_TAU_U_G(TAU_U_G, IER) 
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
      USE visc_g
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
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
!                      TAU_U_g 
      DOUBLE PRECISION TAU_U_g(DIMENSION_3) 
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
!                      Average viscosity 
      DOUBLE PRECISION EPMU_gte, EPMU_gbe, EPMUGA 
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
!
!
!!$omp  parallel do private( IJK, I, IJKE, EPGA, IP, J, JM, K, KM,  &
!!$omp&  IPJK,IMJK,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJMK,IJKT,IJKTE,  &
!!$omp&  IJKB,IJKBE,IJKM,IPJKM, &
!!$omp&  SBV,  SSX,SSY,  EPMU_GTE,EPMU_GBE, SSZ, &
!!$omp&  EPMUGA,DWOXDZ,VTZB )
!!$omp&  schedule(static)
      DO IJK = 1, IJKMAX2 
         I = I_OF(IJK) 
         IJKE = EAST_OF(IJK) 
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
         IF ( .NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN 
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
            SBV = (LAMBDA_GT(IJKE)*TRD_G(IJKE)-LAMBDA_GT(IJK)*TRD_G(IJK))*AYZ(&
               IJK) 
!
!         shear stress terms
            SSX = MU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ODX(IP)*AYZ_U(IJK) - MU_GT(&
               IJK)*(U_G(IJK)-U_G(IMJK))*ODX(I)*AYZ_U(IMJK) 
            SSY = AVG_X_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE)&
               ,MU_GT(IJKNE),J),I)*(V_G(IPJK)-V_G(IJK))*ODX_E(I)*AXZ_U(IJK) - &
               AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(IJKSE),&
               MU_GT(IJKE),JM),I)*(V_G(IPJMK)-V_G(IJMK))*ODX_E(I)*AXZ_U(IJMK) 
            EPMU_GTE = AVG_X_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
               IJKE),MU_GT(IJKTE),K),I) 
            EPMU_GBE = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT&
               (IJKBE),MU_GT(IJKE),KM),I) 
            SSZ = EPMU_GTE*(W_G(IPJK)-W_G(IJK))*ODX_E(I)*AXY_U(IJK) - EPMU_GBE*&
               (W_G(IPJKM)-W_G(IJKM))*ODX_E(I)*AXY_U(IJKM) 
!
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           modify Ssz
               SSZ = SSZ - (EPMU_GTE*(HALF*(W_G(IPJK)+W_G(IJK))*OX_E(I))*AXY_U(&
                  IJK)-EPMU_GBE*(HALF*(W_G(IPJKM)+W_G(IJKM))*OX_E(I))*AXY_U(&
                  IJKM)) 
!
!           Tau_zz/X
               EPMUGA = AVG_X(MU_G(IJK),MU_G(IJKE),I) 
               DWOXDZ = HALF*((W_G(IJK)-W_G(IJKM))*OX(I)*ODZ(K)+(W_G(IPJK)-W_G(&
                  IPJKM))*OX(IP)*ODZ(K)) 
               VTZB = -2.*EPMUGA*OX_E(I)*DWOXDZ 
            ELSE 
               VTZB = ZERO 
            ENDIF 
!
!         Add the terms
            TAU_U_G(IJK) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK) 
         ELSE 
            TAU_U_G(IJK) = ZERO 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_TAU_U_G 
