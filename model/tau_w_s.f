!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_W_s(TAU_W_s, IER)                             C
!  Purpose: Cross terms in the gradient of stress in W_s momentum      c
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
      SUBROUTINE CALC_TAU_W_S(TAU_W_S, IER) 
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
!                      TAU_W_s 
      DOUBLE PRECISION TAU_W_s(DIMENSION_3, DIMENSION_M) 
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
!                      Average velocity gradients 
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
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
!
      DO M = 1, MMAX 

!!$omp  parallel do private( IJK, I, IJKE, EPGA,  J,  K,   &
!!$omp&  JM,IJKP,IJKNT,IJKS,IJKST,IJMKP,IJKTE,IJKTW,IMJKP,KP, &
!!$omp&  DUODZ,VXZ, &
!!$omp&  IM,IJKW, &
!!$omp&  IMJK,IJKN,IJMK,IJKT,  &
!!$omp&  IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ) &
!!$omp&  schedule(static)

         DO IJK = 1, IJKMAX2 
            K = K_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            EPGA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K) 
            IF ( .NOT.SIP_AT_T(IJK) .AND. EPGA>DIL_EP_S) THEN 
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
               SBV = (LAMBDA_S(IJKT,M)*TRD_S(IJKT,M)-LAMBDA_S(IJK,M)*TRD_S(IJK,&
                  M))*AXY(IJK) 
!
!         shear stress terms
               SSX = AVG_Z_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
                  IJKT,M),MU_S(IJKTE,M),I),K)*(U_S(IJKP,M)-U_S(IJK,M))*OX_E(I)*&
                  ODZ_T(K)*AYZ_W(IJK) - AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M&
                  ),IM),AVG_X_H(MU_S(IJKTW,M),MU_S(IJKT,M),IM),K)*(U_S(IMJKP,M)&
                  -U_S(IMJK,M))*ODZ_T(K)*DY(J)*(HALF*(DZ(K)+DZ(KP))) 
!
               SSY = AVG_Z_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
                  IJKT,M),MU_S(IJKNT,M),J),K)*(V_S(IJKP,M)-V_S(IJK,M))*OX(I)*&
                  ODZ_T(K)*AXZ_W(IJK) - AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M&
                  ),JM),AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*(V_S(IJMKP,M)&
                  -V_S(IJMK,M))*OX(I)*ODZ_T(K)*AXZ_W(IJMK) 
!
               SSZ = MU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*OX(I)*ODZ(KP)*AXY_W(&
                  IJK) - MU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*OX(I)*ODZ(K)*&
                  AXY_W(IJKM) 
!
!
!         Special terms for cylindrical coordinates
               IF (CYLINDRICAL) THEN 
!
!           modify Szz
                  SSZ = SSZ + MU_S(IJKT,M)*(U_S(IJKP,M)+U_S(IMJKP,M))*OX(I)*&
                     AXY_W(IJK) - MU_S(IJK,M)*(U_S(IJK,M)+U_S(IMJK,M))*OX(I)*&
                     AXY_W(IJKM) 
!
!           part of tau_xz/X
                  IF (OX_E(IM) /= UNDEFINED) THEN 
                     DUODZ = (U_S(IMJKP,M)-U_S(IMJK,M))*OX_E(IM)*ODZ_T(K) 
                  ELSE 
                     DUODZ = ZERO 
                  ENDIF 
                  VXZ = AVG_Z(MU_S(IJK,M),MU_S(IJKT,M),K)*OX(I)*HALF*((U_S(IJKP&
                     ,M)-U_S(IJK,M))*OX_E(I)*ODZ_T(K)+DUODZ) 
!
               ELSE 
                  VXZ = ZERO 
               ENDIF 
!
!         Add the terms
               TAU_W_S(IJK,M) = SBV + SSX + SSY + SSZ + VXZ*VOL_W(IJK) 
            ELSE 
               TAU_W_S(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE CALC_TAU_W_S 
