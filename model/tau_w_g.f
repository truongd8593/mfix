!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_W_g(TAU_W_g, IER)                             C
!  Purpose: Cross terms in the gradient of stress in W_g momentum      c
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
      SUBROUTINE CALC_TAU_W_G(TAU_W_G, IER) 
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
!                      TAU_W_g 
      DOUBLE PRECISION TAU_W_g(DIMENSION_3) 
! 
!                      Indices 
      INTEGER          IJK, J, I, IM, IJKP, IMJK, IJKN, IJKNT, IJKS,& 
                       IJKST, IJMKP, IJMK, IJKE, IJKTE, IJKW, IJKTW,& 
                       IMJKP, K, IJKT, JM, KP, IJKM 
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
!                      Average gradients 
      DOUBLE PRECISION dWoXdz, duodz 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION Vxz 
! 
!                      error message 
      CHARACTER*80     LINE 
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!
!!$omp  parallel do private( IJK, I, IJKE, EPGA,  J,  K,   &
!!$omp&  JM,IJKP,IJKNT,IJKS,IJKST,IJMKP,IJKTE,IJKTW,IMJKP,KP, &
!!$omp&  DUODZ,VXZ, &
!!$omp&  IM,IJKW, &
!!$omp&  IMJK,IJKN,IJMK,IJKT,  &
!!$omp&  IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ) &
!!$omp&  schedule(static)

!//SP
      DO IJK = IJKSTART3, IJKEND3
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K) 
         IF ( .NOT.IP_AT_T(IJK) .AND. EPGA>DIL_EP_S) THEN 
            J = J_OF(IJK) 
            I = I_OF(IJK) 
            IM = IM1(I) 
            JM = JM1(J) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKNT = TOP_OF(IJKN) 
            IJKS = SOUTH_OF(IJK) 
            IJKST = TOP_OF(IJKS) 
            IJMKP = JM_OF(IJKP) 
            IJMK = JM_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKTE = EAST_OF(IJKT) 
            IJKW = WEST_OF(IJK) 
            IJKTW = WEST_OF(IJKT) 
            IMJKP = KP_OF(IMJK) 
            KP = KP1(K) 
            IJKM = KM_OF(IJK) 
!
!       Surface forces
!
!         bulk viscosity term
            SBV = (LAMBDA_GT(IJKT)*TRD_G(IJKT)-LAMBDA_GT(IJK)*TRD_G(IJK))*AXY(&
               IJK) 
!
!         shear stress terms
            SSX = AVG_Z_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT)&
               ,MU_GT(IJKTE),I),K)*(U_G(IJKP)-U_G(IJK))*OX_E(I)*ODZ_T(K)*AYZ_W(&
               IJK) - AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(&
               IJKTW),MU_GT(IJKT),IM),K)*(U_G(IMJKP)-U_G(IMJK))*ODZ_T(K)*DY(J)*&
               (HALF*(DZ(K)+DZ(KP))) 
            SSY = AVG_Z_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT)&
               ,MU_GT(IJKNT),J),K)*(V_G(IJKP)-V_G(IJK))*OX(I)*ODZ_T(K)*AXZ_W(&
               IJK) - AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(&
               IJKST),MU_GT(IJKT),JM),K)*(V_G(IJMKP)-V_G(IJMK))*OX(I)*ODZ_T(K)*&
               AXZ_W(IJMK) 
!
            SSZ = MU_GT(IJKT)*(W_G(IJKP)-W_G(IJK))*OX(I)*ODZ(KP)*AXY_W(IJK) - &
               MU_GT(IJK)*(W_G(IJK)-W_G(IJKM))*OX(I)*ODZ(K)*AXY_W(IJKM) 
!
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           modify Szz
               SSZ = SSZ + MU_GT(IJKT)*(U_G(IJKP)+U_G(IMJKP))*OX(I)*AXY_W(IJK)&
                   - MU_GT(IJK)*(U_G(IJK)+U_G(IMJK))*OX(I)*AXY_W(IJKM) 
!
!           part of tau_xz/X
               IF (OX_E(IM) /= UNDEFINED) THEN 
                  DUODZ = (U_G(IMJKP)-U_G(IMJK))*OX_E(IM)*ODZ_T(K) 
               ELSE 
                  DUODZ = ZERO 
               ENDIF 
               VXZ = AVG_Z(MU_GT(IJK),MU_GT(IJKT),K)*OX(I)*HALF*((U_G(IJKP)-U_G&
                  (IJK))*OX_E(I)*ODZ_T(K)+DUODZ) 
!
            ELSE 
               VXZ = ZERO 
            ENDIF 
!
!         Add the terms
            TAU_W_G(IJK) = SBV + SSX + SSY + SSZ + VXZ*VOL_W(IJK) 
         ELSE 
            TAU_W_G(IJK) = ZERO 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_TAU_W_G 
