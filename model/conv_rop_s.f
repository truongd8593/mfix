!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_s(A_m, B_m, M, IER)                           C
!  Purpose: Determine convection terms for solids continuity           C
!            equation.  Master routine.                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
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
      SUBROUTINE CONV_ROP_S(A_M, B_M, M, IER) 
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
      USE compar     !//
      USE sendrecv   !// 400
      
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
!                      phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!
      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_ROP_S0 (A_M, B_M, M, IER) 
      ELSE 
         CALL CONV_ROP_S1 (A_M, B_M, M, IER) 
      ENDIF 

!// 400 1218 Communicate boundaries   
      CALL SEND_RECV(A_M, 2)
      CALL SEND_RECV(B_M, 2)
      
      RETURN  
      END SUBROUTINE CONV_ROP_S 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_s0(A_m, B_m, M, IER)                          C
!  Purpose: Determine convection terms for solids continuity           C
!           equation.  The off-diagonal coefficients calculated here   C
!  must be positive. The center coefficient and the source vector are  C
!  negative. -- First order upwinding                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-JUL-96   C
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
      SUBROUTINE CONV_ROP_S0(A_M, B_M, M, IER) 
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
      USE pscor
      USE compar        !//d
      USE sendrecv      !// 400
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
!                      Indices 
      INTEGER          I, J, K, IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!     Calculate convection-diffusion fluxes through each of the faces
!
!//? make sure all vars used in following section to calc. A_M are uptodate
!//  at the boundaries, if necessary insert redundant COMM of them here.

!// 350 1218 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!!$omp  parallel do private( I, J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3 
         IF (PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
!// 360 1218 Check if current i,j,k resides on this PE	    
            IF(.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	    	    
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
!         East face (i+1/2, j, k)
            A_M(IJK,E,M) = ZMAX((-U_S(IJK,M)))*AYZ(IJK) 
            A_M(IPJK,W,M) = ZMAX(U_S(IJK,M))*AYZ(IJK) 
!
!         North face (i, j+1/2, k)
            A_M(IJK,N,M) = ZMAX((-V_S(IJK,M)))*AXZ(IJK) 
            A_M(IJPK,S,M) = ZMAX(V_S(IJK,M))*AXZ(IJK) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               A_M(IJK,T,M) = ZMAX((-W_S(IJK,M)))*AXY(IJK) 
               A_M(IJKP,B,M) = ZMAX(W_S(IJK,M))*AXY(IJK) 
            ENDIF 
            IF (PHASE_4_P_G(IMJK)==M .OR. PHASE_4_P_S(IMJK)==M) A_M(IJK,W,M) = &
               ZMAX(U_S(IMJK,M))*AYZ(IMJK) 
            IF (PHASE_4_P_G(IJMK)==M .OR. PHASE_4_P_S(IJMK)==M) A_M(IJK,S,M) = &
               ZMAX(V_S(IJMK,M))*AXZ(IJMK) 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IF (PHASE_4_P_G(IJKM)==M .OR. PHASE_4_P_S(IJKM)==M) A_M(IJK,B,M)&
                   = ZMAX(W_S(IJKM,M))*AXY(IJKM) 
            ENDIF 
         ENDIF 
      END DO 

!//S This COMM could be removed as it is redundant but inserted due fool-proof      
!// 400 1218 Communicate boundaries         
      call send_recv(A_M,2)      
      
      RETURN  
      END SUBROUTINE CONV_ROP_S0 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP_s1(A_m, B_m, M, IER)                          C
!  Purpose: Determine convection terms for solids continuity           C
!           equation.  The off-diagonal coefficients calculated here   C
!  must be positive. The center coefficient and the source vector are  C
!  negative. -- Higher order methods                                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
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
      SUBROUTINE CONV_ROP_S1(A_M, B_M, M, IER) 
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
      USE pscor
      Use xsi_array
      USE compar        !//d
      USE sendrecv   !// 400      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
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
!
!  Calculate convection factors
!

! loezos
	 incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(2), ROP_S(1,M), U_S(1,M), V_S(1,M), W_S(1,M), &
         XSI_E, XSI_N, XSI_T,incr) 
!
!     Calculate convection-diffusion fluxes through each of the faces
!
!// 350 1218 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!!$omp  parallel do private( I, J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
!// 360 1218 Check if current i,j,k resides on this PE	    
            IF(.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	    	    
	    
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
!         East face (i+1/2, j, k)
            A_M(IJK,E,M) = -XSI_E(IJK)*U_S(IJK,M)*AYZ(IJK) 
            A_M(IPJK,W,M) = (ONE - XSI_E(IJK))*U_S(IJK,M)*AYZ(IJK) 
!
!         North face (i, j+1/2, k)
            A_M(IJK,N,M) = -XSI_N(IJK)*V_S(IJK,M)*AXZ(IJK) 
            A_M(IJPK,S,M) = (ONE - XSI_N(IJK))*V_S(IJK,M)*AXZ(IJK) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               A_M(IJK,T,M) = -XSI_T(IJK)*W_S(IJK,M)*AXY(IJK) 
               A_M(IJKP,B,M) = (ONE - XSI_T(IJK))*W_S(IJK,M)*AXY(IJK) 
            ENDIF 
            IF (PHASE_4_P_G(IMJK)==M .OR. PHASE_4_P_S(IMJK)==M) A_M(IJK,W,M) = &
               (ONE - XSI_E(IMJK))*U_S(IMJK,M)*AYZ(IMJK) 
            IF (PHASE_4_P_G(IJMK)==M .OR. PHASE_4_P_S(IJMK)==M) A_M(IJK,S,M) = &
               (ONE - XSI_N(IJMK))*V_S(IJMK,M)*AXZ(IJMK) 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IF (PHASE_4_P_G(IJKM)==M .OR. PHASE_4_P_S(IJKM)==M) A_M(IJK,B,M)&
                   = (ONE - XSI_T(IJKM))*W_S(IJKM,M)*AXY(IJKM) 
            ENDIF 
         ENDIF 
      END DO 

!//S This COMM could be removed as it is redundant but inserted due fool-proof      
!// 400 1218 Communicate boundaries         
      call send_recv(A_M,2)      
      
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE CONV_ROP_S1 
