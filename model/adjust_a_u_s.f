!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_A_U_s(A_m, B_m, IER)                            C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_s momentum eq. becoming zero.                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
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
      SUBROUTINE ADJUST_A_U_S(A_M, B_M, IER) 
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
      USE fldvar
      USE physprop
      USE geometry
      USE run
      USE indices
      USE compar
      USE sendrecv  !// 400      
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
!                      Indices
      INTEGER          I, IP, IJK, IJKE, IMJK
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
!
!// 350 1229 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

!!$omp     parallel do private(I, IP, IJK, IJKE, IMJK )
            DO IJK = ijkstart3, ijkend3
               IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  IF (B_M(IJK,M) < ZERO) THEN 
                     IJKE = EAST_OF(IJK) 
                     IP = IP1(I_OF(IJK)) 
                     IF (ROP_S(IJKE,M) > SMALL_NUMBER) THEN 
                        B_M(IJK,M) = SQRT((-B_M(IJK,M)/(ROP_S(IJKE,M)*AVG_X_E(&
                           ONE,ZERO,IP)*AYZ_U(IJK)))) 
                     ELSE 
                        B_M(IJK,M) = ZERO 
                     ENDIF 
                  ELSE IF (B_M(IJK,M) > ZERO) THEN 
                     I = I_OF(IJK) 
                     IMJK = IM_OF(IJK) 
                     IF (ROP_S(IJK,M) > SMALL_NUMBER) THEN 
                        B_M(IJK,M) = SQRT(B_M(IJK,M)/(ROP_S(IJK,M)*AVG_X_E(ZERO&
                           ,ONE,I)*AYZ_U(IMJK))) 
                     ELSE 
                        B_M(IJK,M) = ZERO 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

!// 400 1229 COMM A_M & B_M
!!!      call send_recv(A_M,2)
!!!      call send_recv(B_M,2)            
      
      RETURN  
      END SUBROUTINE ADJUST_A_U_S 
