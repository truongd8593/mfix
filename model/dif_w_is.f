!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DIF_W_IS(Dif, A_m, B_m, M, IER)                        C
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
      SUBROUTINE DIF_W_IS(DIF, A_M, B_M, M, IER) 
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
      USE toleranc 
      USE run
      USE physprop
      USE fldvar
      USE geometry
      USE output
      USE indices
      USE is
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
!                      Internal surface 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,& 
                       IJKE, IJKN, IJKT, IJPK, IPJK, IJKTN, IJKTE 
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
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
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
!
            IF (IS_PLANE(L) == 'N') THEN 
!
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
                        IJK = FUNIJK(I,J,K) 
                        IJKT = TOP_OF(IJK) 
                        IJKN = NORTH_OF(IJK) 
                        IJKTN = TOP_OF(IJKN) 
                        IJPK = JP_OF(IJK) 
!
                        D_F = AVG_Z_H(AVG_Y_H(DIF(IJK),DIF(IJKN),J),AVG_Y_H(DIF&
                           (IJKT),DIF(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK) 
!
                        A_M(IJK,N,M) = A_M(IJK,N,M) - D_F 
                        A_M(IJPK,S,M) = A_M(IJPK,S,M) - D_F 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (IS_PLANE(L) == 'E') THEN 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
                        IJK = FUNIJK(I,J,K) 
                        IJKE = EAST_OF(IJK) 
                        IJKT = TOP_OF(IJK) 
                        IJKTE = EAST_OF(IJKT) 
                        IPJK = KP_OF(IJK) 
!
                        D_F = AVG_Z_H(AVG_X_H(DIF(IJK),DIF(IJKE),I),AVG_X_H(DIF&
                           (IJKT),DIF(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK) 
!
                        A_M(IJK,E,M) = A_M(IJK,E,M) - D_F 
                        A_M(IPJK,W,M) = A_M(IPJK,W,M) - D_F 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE DIF_W_IS 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
