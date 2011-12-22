!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_Pp_g(A_m, B_m, IER)
!  Purpose: Determine convection terms for Pressure                    C
!  correction equation.  The off-diagonal coefficients calculated here C
!  must be positive. The center coefficient and the source vector are  C
!  negative. Multiplication with factors d_e, d_n, and d_t are carried C
!  out in source_pp_g.  Constant pressure boundaries are handled by
!  holding the Pp_g at the boundaries zero.  For specified mass flow   C
!  boundaries (part of) a's are calculated here since b is calculated  C
!  from a's in source_pp_g.  After calculating b, a's are multiplied by
!  d and at the flow boundaries get are set to zero.                   C                                !  Author: M. Syamlal                                 Date: 20-JUN-96  C
!  Reviewer:                                          Date:            C
!   Revision to use face densities calculated in CONV_ROP              C
!  Author: M. Syamlal                                 Date: 1-JUN-05  C
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
      SUBROUTINE CONV_PP_G(A_M, B_M, IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE run 
      USE parallel 
      USE matrix 
      USE physprop
      USE geometry
      USE indices
      USE pgcor
      USE compar
      USE mflux   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
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
!                      Indices 
      INTEGER          IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM 
      INTEGER          M 
! 
!                      local value of A_m 
      DOUBLE PRECISION am 
!-----------------------------------------------
      INCLUDE 'function.inc'
      

!
!  Calculate convection fluxes through each of the faces

!!!$omp  parallel do private( IJK, IPJK, IJPK, M, AM, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            AM = ROP_GE(IJK)*AYZ(IJK) 
            A_M(IJK,E,0) = AM 
            A_M(IPJK,W,0) = AM 
!
!         North face (i, j+1/2, k)
            AM = ROP_GN(IJK)*AXZ(IJK) 
            A_M(IJK,N,0) = AM 
            A_M(IJPK,S,0) = AM 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               AM = ROP_GT(IJK)*AXY(IJK) 
               A_M(IJK,T,0) = AM 
               A_M(IJKP,B,0) = AM 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               AM = ROP_GE(IMJK)*AYZ(IMJK) 
               A_M(IJK,W,0) = AM 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               AM = ROP_GN(IJMK)*AXZ(IJMK) 
               A_M(IJK,S,0) = AM 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  AM = ROP_GT(IJKM)*AXY(IJKM) 
                  A_M(IJK,B,0) = AM 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      DO M = 1, MMAX 
         IF (.NOT.CLOSE_PACKED(M)) THEN
            DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN 
                  IPJK = IP_OF(IJK) 
                  IJPK = JP_OF(IJK) 
                  IJKP = KP_OF(IJK) 
!
!             East face (i+1/2, j, k)
                  AM = ROP_SE(IJK,M)*AYZ(IJK) 
                  A_M(IJK,E,M) = AM 
                  A_M(IPJK,W,M) = AM 
!
!             North face (i, j+1/2, k)
                  AM = ROP_SN(IJK,M)*AXZ(IJK) 
                  A_M(IJK,N,M) = AM 
                  A_M(IJPK,S,M) = AM 
!
!             Top face (i, j, k+1/2)
                  IF (DO_K) THEN 
                     AM = ROP_ST(IJK,M)*AXY(IJK) 
                     A_M(IJK,T,M) = AM 
                     A_M(IJKP,B,M) = AM 
                  ENDIF 
!
!             West face (i-1/2, j, k)
                  IMJK = IM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IMJK)) THEN 
                     AM = ROP_SE(IMJK,M)*AYZ(IMJK) 
                     A_M(IJK,W,M) = AM 
                  ENDIF 
!
!             South face (i, j-1/2, k)
                  IJMK = JM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IJMK)) THEN 
                     AM = ROP_SN(IJMK,M)*AXZ(IJMK) 
                     A_M(IJK,S,M) = AM 
                  ENDIF 
!
!             Bottom face (i, j, k-1/2)
                  IF (DO_K) THEN 
                     IJKM = KM_OF(IJK) 
                     IF (.NOT.FLUID_AT(IJKM)) THEN 
                        AM = ROP_ST(IJKM,M)*AXY(IJKM) 
                        A_M(IJK,B,M) = AM 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CONV_PP_G 
