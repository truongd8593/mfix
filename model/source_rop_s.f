!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_ROP_s(A_m, B_m, M, IER)                         C
!  Purpose: Determine source terms for solids continuity equation.     C
!  The off-diagonal coefficients are                                   C
!  positive. The center coefficient and the source vector are          C
!  negative.                                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3 -JUL-96  C
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
      SUBROUTINE SOURCE_ROP_S(A_M, B_M, M, IER) 
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
      USE pscor
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
!                      phase index 
      INTEGER          M 
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
!                      FOR CALL_CHEM and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!!$omp  parallel do private( I, J, K, IJK, IMJK, IJMK, IJKM,  DEL_V, &
!!$omp&  Src, LINE) &
!!$omp&  schedule(static)
!
!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
      IF (CALL_CHEM .or. CALL_ISAT) THEN
         SUM_R_S_temp = SUM_R_S
         SUM_R_S = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
!
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK) .AND. PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M&
            ) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
            DEL_V = U_S(IJK,M)*AYZ(IJK) - U_S(IMJK,M)*AYZ(IMJK) + V_S(IJK,M)*&
               AXZ(IJK) - V_S(IJMK,M)*AXZ(IJMK) + W_S(IJK,M)*AXY(IJK) - W_S(&
               IJKM,M)*AXY(IJKM) 
!
            IF (ROP_S(IJK,M) > ZERO) THEN 
               SRC = VOL(IJK)*ZMAX((-SUM_R_S(IJK,M)))/ROP_S(IJK,M) 
            ELSE 
               SRC = ZERO 
            ENDIF 
!
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+VOL(IJK)*ODT+ZMAX(DEL_V)+SRC) 
            B_M(IJK,M) = -(ROP_SO(IJK,M)*VOL(IJK)*ODT+ZMAX((-DEL_V))*ROP_S(IJK,&
               M)+ZMAX(SUM_R_S(IJK,M))*VOL(IJK)) 
            IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,M)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,M) = -ONE            ! Equation is undefined. 
                  B_M(IJK,M) = -ROP_S(IJK,M)     ! Use existing value 
               ELSE 
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', M, ' A = 0 and b = ', B_M(IJK,M) 
                  CALL WRITE_ERROR ('SOURCE_ROP_s', LINE, 1) 
!!$omp             end critical
               ENDIF 
            ENDIF 
         ELSE 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = -ROP_S(IJK,M) 
         ENDIF 
      END DO 
!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_CHEM .or. CALL_ISAT) THEN
         SUM_R_S = SUM_R_S_temp
      END IF
!     CHEM & ISAT end (nan xie)
! 
      
      RETURN  
      END SUBROUTINE SOURCE_ROP_S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
