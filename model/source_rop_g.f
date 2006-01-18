!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_ROP_g(A_m, B_m, IER)                            C
!  Purpose: Determine source terms for continuity equation.            C
!  The off-diagonal coefficients are                                   C
!  positive. The center coefficient and the source vector are          C
!  negative.                                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2 -JUL-96  C
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
      SUBROUTINE SOURCE_ROP_G(A_M, B_M, IER) 
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
      USE rxns
      USE run
      USE geometry
      USE indices
      USE pgcor
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
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      DEL dot V 
      DOUBLE PRECISION DEL_V 
! 
!                      Mass source 
      DOUBLE PRECISION Src 
! 
!                      Indices 
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM 
! 
!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_G_temp = SUM_R_G
         SUM_R_G = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
!
!!$omp  parallel do private( I, J, K, IJK, IMJK, IJMK, IJKM,  DEL_V, &
!!$omp&  Src, LINE) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK) .AND. PHASE_4_P_G(IJK)/=0) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
            DEL_V = U_G(IJK)*AYZ(IJK) - U_G(IMJK)*AYZ(IMJK) + V_G(IJK)*AXZ(IJK)&
                - V_G(IJMK)*AXZ(IJMK) + W_G(IJK)*AXY(IJK) - W_G(IJKM)*AXY(IJKM) 
!
            IF (ROP_G(IJK) > ZERO) THEN 
               SRC = VOL(IJK)*ZMAX((-SUM_R_G(IJK)))/ROP_G(IJK) 
            ELSE 
               SRC = ZERO 
            ENDIF 
!
            A_M(IJK,0,0) = -(A_M(IJK,E,0)+A_M(IJK,W,0)+A_M(IJK,N,0)+A_M(IJK,S,0&
               )+A_M(IJK,T,0)+A_M(IJK,B,0)+VOL(IJK)*ODT+ZMAX(DEL_V)+SRC) 
            B_M(IJK,0) = -(ROP_GO(IJK)*VOL(IJK)*ODT+ZMAX((-DEL_V))*ROP_G(IJK)+&
               ZMAX(SUM_R_G(IJK))*VOL(IJK)) 
            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,0) = -ONE 
                  B_M(IJK,0) = ZERO 
               ELSE 
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0) 
                  CALL WRITE_ERROR ('SOURCE_ROP_g', LINE, 1) 
!!$omp             end critical
               ENDIF 
            ENDIF 
         ELSE 
            A_M(IJK,E,0) = ZERO 
            A_M(IJK,W,0) = ZERO 
            A_M(IJK,N,0) = ZERO 
            A_M(IJK,S,0) = ZERO 
            A_M(IJK,T,0) = ZERO 
            A_M(IJK,B,0) = ZERO 
            A_M(IJK,0,0) = -ONE 
            B_M(IJK,0) = -ROP_G(IJK) 
         ENDIF 
      END DO 
!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_G = SUM_R_G_temp
      END IF 
!     CHEM & ISAT end (nan xie) 
!      
      RETURN  
      END SUBROUTINE SOURCE_ROP_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
