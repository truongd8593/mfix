!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_g(A_m, B_m, IER)                              C
!  Purpose: Determine convection terms for gas continuity              C
!            equation.  Master routine.                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-MAR-97  C
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
      SUBROUTINE CONV_ROP_G(A_M, B_M, IER) 
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
      USE fldvar
      USE run
      USE compar    
  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!! 
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
!
!
      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_ROP_G0 (A_M, B_M, IER) 
      ELSE 
         CALL CONV_ROP_G1 (A_M, B_M, IER) 
      ENDIF 
      
            
      RETURN  
      END SUBROUTINE CONV_ROP_G 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_g0(A_m, B_m, IER)                             C
!  Purpose: Determine convection terms for gas continuity              C
!           equation.  The off-diagonal coefficients calculated here   C
!  must be positive. The center coefficient and the source vector are  C
!  negative.  -- First order upwinding                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2-JUL-96   C
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
      SUBROUTINE CONV_ROP_G0(A_M, B_M, IER) 
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
      USE fldvar
      USE run 
      USE parallel 
      USE matrix 
      USE physprop
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
      INTEGER          I, J, K, IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'


!
!  Calculate convection fluxes through each of the faces
!
!!$omp  parallel do private( J, K, IJK, IPJK, IJPK, IJKP, &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3 
         IF (PHASE_4_P_G(IJK) /= 0) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
!         East face (i+1/2, j, k)
            A_M(IJK,E,0) = ZMAX((-U_G(IJK)))*AYZ(IJK) 
            A_M(IPJK,W,0) = ZMAX(U_G(IJK))*AYZ(IJK) 
!
!         North face (i, j+1/2, k)
            A_M(IJK,N,0) = ZMAX((-V_G(IJK)))*AXZ(IJK) 
            A_M(IJPK,S,0) = ZMAX(V_G(IJK))*AXZ(IJK) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               A_M(IJK,T,0) = ZMAX((-W_G(IJK)))*AXY(IJK) 
               A_M(IJKP,B,0) = ZMAX(W_G(IJK))*AXY(IJK) 
            ENDIF 
            IF(PHASE_4_P_G(IMJK)==0)A_M(IJK,W,0)=ZMAX(U_G(IMJK))*AYZ(IMJK) 
            IF(PHASE_4_P_G(IJMK)==0)A_M(IJK,S,0)=ZMAX(V_G(IJMK))*AXZ(IJMK) 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IF(PHASE_4_P_G(IJKM)==0)A_M(IJK,B,0)=ZMAX(W_G(IJKM))*AXY(IJKM) 
            ENDIF 
         ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE CONV_ROP_G0 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_g1(A_m, B_m, IER)                             C
!  Purpose: Determine convection terms for gas continuity              C
!           equation.  The off-diagonal coefficients calculated here   C
!  must be positive. The center coefficient and the source vector are  C
!  negative.  -- Higher order schemes                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date:19-MAR-97   C
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
      SUBROUTINE CONV_ROP_G1(A_M, B_M, IER) 
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
      USE fldvar
      USE run
      USE parallel 
      USE matrix 
      USE physprop
      USE geometry
      USE indices
      USE pgcor
      Use xsi_array
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
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
! 
!                      Indices 
      INTEGER          I, J, K, IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'

! loezos
	INTEGER incr
! loezos     
      call lock_xsi_array
!
!  Calculate convection factors
!

! loezos
	incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(1), ROP_G, U_G, V_G, W_G, XSI_E, XSI_N,& 
	XSI_T,incr) 
!
!  Calculate convection fluxes through each of the faces
!
!!$omp  parallel do private( J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK) /= 0) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
!         East face (i+1/2, j, k)
            A_M(IJK,E,0) = -XSI_E(IJK)*U_G(IJK)*AYZ(IJK) 
            A_M(IPJK,W,0) = (ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK) 
!
!         North face (i, j+1/2, k)
            A_M(IJK,N,0) = -XSI_N(IJK)*V_G(IJK)*AXZ(IJK) 
            A_M(IJPK,S,0) = (ONE - XSI_N(IJK))*V_G(IJK)*AXZ(IJK) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               A_M(IJK,T,0) = -XSI_T(IJK)*W_G(IJK)*AXY(IJK) 
               A_M(IJKP,B,0) = (ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK) 
            ENDIF 
            IF (PHASE_4_P_G(IMJK) == 0) A_M(IJK,W,0) = (ONE - XSI_E(IMJK))*U_G(&
               IMJK)*AYZ(IMJK) 
            IF (PHASE_4_P_G(IJMK) == 0) A_M(IJK,S,0) = (ONE - XSI_N(IJMK))*V_G(&
               IJMK)*AXZ(IJMK) 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IF (PHASE_4_P_G(IJKM) == 0) A_M(IJK,B,0) = (ONE - XSI_T(IJKM))*&
                  W_G(IJKM)*AXY(IJKM) 
            ENDIF 
         ENDIF 
      END DO 
      
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE CONV_ROP_G1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
