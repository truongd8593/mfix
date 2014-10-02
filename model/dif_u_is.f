!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DIF_U_IS(Dif, A_m, B_m, M, IER)                        C
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
      SUBROUTINE DIF_U_IS(DIF, A_M, B_M, M, IER)
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
      USE fun_avg
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
                       IJKE, IJKN, IJKT, IJPK, IJKP, IJKNE, IJKTE
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
      INCLUDE 'function.inc'
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
            IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

            IF (IS_PLANE(L) == 'N') THEN
!
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IJKE = EAST_OF(IJK)
                        IJKN = NORTH_OF(IJK)
                        IJKNE = EAST_OF(IJKN)
                        IJPK = JP_OF(IJK)
!
                        D_F = AVG_X_H(AVG_Y_H(DIF(IJK),DIF(IJKN),J),AVG_Y_H(DIF&
                           (IJKE),DIF(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK)
!
                        A_M(IJK,N,M) = A_M(IJK,N,M) - D_F
                        A_M(IJPK,S,M) = A_M(IJPK,S,M) - D_F
                     END DO
                  END DO
               END DO
            ELSE IF (IS_PLANE(L) == 'T') THEN
               IF (DO_K) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IJKE = EAST_OF(IJK)
                           IJKT = TOP_OF(IJK)
                           IJKTE = EAST_OF(IJKT)
                           IJKP = KP_OF(IJK)
!
                           D_F = AVG_X_H(AVG_Z_H(DIF(IJK),DIF(IJKT),K),AVG_Z_H(&
                              DIF(IJKE),DIF(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U&
                              (IJK)
!
                           A_M(IJK,T,M) = A_M(IJK,T,M) - D_F
                           A_M(IJKP,B,M) = A_M(IJKP,B,M) - D_F
                        END DO
                     END DO
                  END DO
               ENDIF
            ENDIF
         ENDIF
      END DO
!

      RETURN
      END SUBROUTINE DIF_U_IS

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 300 Limit I1, I2 and all to local processor first ghost layer
