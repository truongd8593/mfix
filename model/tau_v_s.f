!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_V_s(TAU_V_s, IER)                             C
!  Purpose: Cross terms in the gradient of stress in V_s momentum      c
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
      SUBROUTINE CALC_TAU_V_S(TAU_V_S, IER) 
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
!                      TAU_V_s 
      DOUBLE PRECISION TAU_V_s(DIMENSION_3, DIMENSION_M) 
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


!!$omp  parallel do private( IJK, I, IJKE, EPGA,  J,  K, KM,  &
!!$omp& JP,IM,IJPK,IJKW,IJKNW,IMJPK,IJKTN,IJKBN,IJPKM, &
!!$omp&  IMJK,IJKN,IJKNE,IJMK,IJKT,  &
!!$omp&  IJKB,IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ) &
!!$omp&  schedule(static)

         DO IJK = 1, IJKMAX2 
            J = J_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            EPGA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J) 
            IF ( .NOT.SIP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN 
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
               SBV = (LAMBDA_S(IJKN,M)*TRD_S(IJKN,M)-LAMBDA_S(IJK,M)*TRD_S(IJK,&
                  M))*AXZ(IJK) 
!
!         shear stress terms
               SSX = AVG_Y_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
                  IJKN,M),MU_S(IJKNE,M),I),J)*(U_S(IJPK,M)-U_S(IJK,M))*ODY_N(J)&
                  *AYZ_V(IJK) - AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                  AVG_X_H(MU_S(IJKNW,M),MU_S(IJKN,M),IM),J)*(U_S(IMJPK,M)-U_S(&
                  IMJK,M))*ODY_N(J)*AYZ_V(IMJK) 
               SSY = MU_S(IJKN,M)*(V_S(IJPK,M)-V_S(IJK,M))*ODY(JP)*AXZ_V(IJK)&
                   - MU_S(IJK,M)*(V_S(IJK,M)-V_S(IJMK,M))*ODY(J)*AXZ_V(IJMK) 
               SSZ = AVG_Y_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKN,M),MU_S(IJKTN,M),K),J)*(W_S(IJPK,M)-W_S(IJK,M))*ODY_N(J)&
                  *AXY_V(IJK) - AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),&
                  AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*(W_S(IJPKM,M)-W_S(&
                  IJKM,M))*ODY_N(J)*AXY_V(IJKM) 
!
!
!         Add the terms
               TAU_V_S(IJK,M) = SBV + SSX + SSY + SSZ 
            ELSE 
               TAU_V_S(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE CALC_TAU_V_S 
