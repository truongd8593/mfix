!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_V_g(TAU_V_g, IER)                             C
!  Purpose: Cross terms in the gradient of stress in V_g momentum      c
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
      SUBROUTINE CALC_TAU_V_G(TAU_V_G, IER) 
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
      USE sendrecv    
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
!                      Error index 
      INTEGER          IER 
! 
!                      TAU_V_g 
      DOUBLE PRECISION TAU_V_g(DIMENSION_3) 
! 
!                      Indices 
      INTEGER          I, J, K, IJK, IJKN, JP, IM,  KM, IJPK, IJMK,& 
                       IJKE, IJKNE, IJKW, IJKNW, IMJPK, IMJK, IJKT,& 
                       IJKTN, IJKB, IJKBN, IJKM, IJPKM 
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
!                      Source terms (Surface) 
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!

!!$omp  parallel do private( IJK, I, IJKE, EPGA,  J,  K, KM,  &
!!$omp&  IMJK,IJKN,IJKNE,IJMK,IJKT,  &
!!$omp&  JP,IM,IJPK,IJKW,IJKNW,IMJPK,IJKTN,IJKBN,IJPKM,IJKB,IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ)  &
!!$omp&  schedule(static)
      DO IJK = IJKSTART3, IJKEND3
         J = J_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
         IF ( .NOT.IP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN 
            JP = JP1(J) 
            I = I_OF(IJK) 
            IM = IM1(I) 
            K = K_OF(IJK) 
            KM = KM1(K) 
            IJPK = JP_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
            IJKW = WEST_OF(IJK) 
            IJKNW = NORTH_OF(IJKW) 
            IMJPK = IM_OF(IJPK) 
            IMJK = IM_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKB = BOTTOM_OF(IJK) 
            IJKBN = NORTH_OF(IJKB) 
            IJKM = KM_OF(IJK) 
            IJPKM = JP_OF(IJKM) 
!
!       Surface forces
!
!         bulk viscosity term
            SBV = (LAMBDA_GT(IJKN)*TRD_G(IJKN)-LAMBDA_GT(IJK)*TRD_G(IJK))*AXZ(&
               IJK) 
!
!         shear stress terms
            SSX = AVG_Y_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN)&
               ,MU_GT(IJKNE),I),J)*(U_G(IJPK)-U_G(IJK))*ODY_N(J)*AYZ_V(IJK) - &
               AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(IJKNW),&
               MU_GT(IJKN),IM),J)*(U_G(IMJPK)-U_G(IMJK))*ODY_N(J)*AYZ_V(IMJK) 
            SSY = MU_GT(IJKN)*(V_G(IJPK)-V_G(IJK))*ODY(JP)*AXZ_V(IJK) - MU_GT(&
               IJK)*(V_G(IJK)-V_G(IJMK))*ODY(J)*AXZ_V(IJMK) 
            SSZ = AVG_Y_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(IJKN)&
               ,MU_GT(IJKTN),K),J)*(W_G(IJPK)-W_G(IJK))*ODY_N(J)*AXY_V(IJK) - &
               AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT(IJKBN),&
               MU_GT(IJKN),KM),J)*(W_G(IJPKM)-W_G(IJKM))*ODY_N(J)*AXY_V(IJKM) 
!
!
!         Add the terms
            TAU_V_G(IJK) = SBV + SSX + SSY + SSZ 
         ELSE 
            TAU_V_G(IJK) = ZERO 
         ENDIF 
      END DO 
      call send_recv(tau_v_g,2)
      RETURN  
      END SUBROUTINE CALC_TAU_V_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
