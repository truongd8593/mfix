!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_Pp_g(A_m, B_m, IER)                               C
!  Purpose: Determine convection terms for pressure correction         C
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
      SUBROUTINE CONV_PP_G(A_M, B_M, IER) 
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
!
!
      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_PP_G0 (A_M, B_M, IER) 
      ELSE 
         CALL CONV_PP_G1 (A_M, B_M, IER) 
      ENDIF 
      RETURN  
      END SUBROUTINE CONV_PP_G 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_Pp_g0(A_m, B_m, IER)
!  Purpose: Determine convection terms for Pressure                    C
!  correction equation.  The off-diagonal coefficients calculated here C
!  must be positive. The center coefficient and the source vector are  C
!  negative. Multiplication with factors d_e, d_n, and d_t are carried C
!  out in source_pp_g.  Constant pressure boundaries are handled by
!  holding the Pp_g at the boundaries zero.  For specified mass flow   C
!  boundaries (part of) a's are calculated here since b is calculated  C
!  from a's in source_pp_g.  After calculating b, a's are multiplied by
!  d and at the flow boundaries get are set to zero.                   C                                                    C
!   FIrst order upwinding                                              C
!  Author: M. Syamlal                                 Date: 20-JUN-96  C
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
      SUBROUTINE CONV_PP_G0(A_M, B_M, IER) 
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
      INTEGER          I, J, K, IJK, IJKE, IJKN, IJKT, IPJK, IJPK, IJKP 
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
      INTEGER          M 
! 
!                      local value of A_m 
      DOUBLE PRECISION am 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      

!
!  Calculate convection fluxes through each of the faces

!$omp  parallel do private( I, J, K, IJK, IJKE, IJKN, IJKT, IPJK, IJPK, &
!$omp&  IJKP, M, AM, &
!$omp&  IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            IF (U_G(IJK) >= ZERO) THEN 
               AM = ROP_G(IJK)*AYZ(IJK) 
            ELSE 
               AM = ROP_G(IJKE)*AYZ(IJK) 
            ENDIF 
            A_M(IJK,E,0) = AM 
            A_M(IPJK,W,0) = AM 
!
!         North face (i, j+1/2, k)
            IF (V_G(IJK) >= ZERO) THEN 
               AM = ROP_G(IJK)*AXZ(IJK) 
            ELSE 
               AM = ROP_G(IJKN)*AXZ(IJK) 
            ENDIF 
            A_M(IJK,N,0) = AM 
            A_M(IJPK,S,0) = AM 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IF (W_G(IJK) >= ZERO) THEN 
                  AM = ROP_G(IJK)*AXY(IJK) 
               ELSE 
                  AM = ROP_G(IJKT)*AXY(IJK) 
               ENDIF 
               A_M(IJK,T,0) = AM 
               A_M(IJKP,B,0) = AM 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               IF (U_G(IMJK) >= ZERO) THEN 
                  AM = ROP_G(IJKW)*AYZ(IMJK) 
               ELSE 
                  AM = ROP_G(IJK)*AYZ(IMJK) 
               ENDIF 
               A_M(IJK,W,0) = AM 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               IF (V_G(IJMK) >= ZERO) THEN 
                  AM = ROP_G(IJKS)*AXZ(IJMK) 
               ELSE 
                  AM = ROP_G(IJK)*AXZ(IJMK) 
               ENDIF 
               A_M(IJK,S,0) = AM 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  IF (W_G(IJKM) >= ZERO) THEN 
                     AM = ROP_G(IJKB)*AXY(IJKM) 
                  ELSE 
                     AM = ROP_G(IJK)*AXY(IJKM) 
                  ENDIF 
                  A_M(IJK,B,0) = AM 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      DO M = 1, MMAX 
         IF (.NOT.CLOSE_PACKED(M)) THEN
            DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN 
                  I = I_OF(IJK) 
                  J = J_OF(IJK) 
                  K = K_OF(IJK) 
                  IJKE = EAST_OF(IJK) 
                  IJKN = NORTH_OF(IJK) 
                  IJKT = TOP_OF(IJK) 
                  IPJK = IP_OF(IJK) 
                  IJPK = JP_OF(IJK) 
                  IJKP = KP_OF(IJK) 
!
!             East face (i+1/2, j, k)
                  IF (U_S(IJK,M) >= ZERO) THEN 
                     AM = ROP_S(IJK,M)*AYZ(IJK) 
                  ELSE 
                     AM = ROP_S(IJKE,M)*AYZ(IJK) 
                  ENDIF 
                  A_M(IJK,E,M) = AM 
                  A_M(IPJK,W,M) = AM 
!
!             North face (i, j+1/2, k)
                  IF (V_S(IJK,M) >= ZERO) THEN 
                     AM = ROP_S(IJK,M)*AXZ(IJK) 
                  ELSE 
                     AM = ROP_S(IJKN,M)*AXZ(IJK) 
                  ENDIF 
                  A_M(IJK,N,M) = AM 
                  A_M(IJPK,S,M) = AM 
!
!             Top face (i, j, k+1/2)
                  IF (DO_K) THEN 
                     IF (W_S(IJK,M) >= ZERO) THEN 
                        AM = ROP_S(IJK,M)*AXY(IJK) 
                     ELSE 
                        AM = ROP_S(IJKT,M)*AXY(IJK) 
                     ENDIF 
                     A_M(IJK,T,M) = AM 
                     A_M(IJKP,B,M) = AM 
                  ENDIF 
!
!             West face (i-1/2, j, k)
                  IMJK = IM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IMJK)) THEN 
                     IJKW = WEST_OF(IJK) 
                     IF (U_S(IMJK,M) >= ZERO) THEN 
                        AM = ROP_S(IJKW,M)*AYZ(IMJK) 
                     ELSE 
                        AM = ROP_S(IJK,M)*AYZ(IMJK) 
                     ENDIF 
                     A_M(IJK,W,M) = AM 
                  ENDIF 
!
!             South face (i, j-1/2, k)
                  IJMK = JM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IJMK)) THEN 
                     IJKS = SOUTH_OF(IJK) 
                     IF (V_S(IJMK,M) >= ZERO) THEN 
                        AM = ROP_S(IJKS,M)*AXZ(IJMK) 
                     ELSE 
                        AM = ROP_S(IJK,M)*AXZ(IJMK) 
                     ENDIF 
                     A_M(IJK,S,M) = AM 
                  ENDIF 
!
!             Bottom face (i, j, k-1/2)
                  IF (DO_K) THEN 
                     IJKM = KM_OF(IJK) 
                     IF (.NOT.FLUID_AT(IJKM)) THEN 
                        IJKB = BOTTOM_OF(IJK) 
                        IF (W_S(IJKM,M) >= ZERO) THEN 
                           AM = ROP_S(IJKB,M)*AXY(IJKM) 
                        ELSE 
                           AM = ROP_S(IJK,M)*AXY(IJKM) 
                        ENDIF 
                        A_M(IJK,B,M) = AM 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CONV_PP_G0 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_Pp_g1(A_m, B_m, IER)
!  Purpose: Determine convection terms for Pressure                    C
!  correction equation.  The off-diagonal coefficients calculated here C
!  must be positive. The center coefficient and the source vector are  C
!  negative. Multiplication with factors d_e, d_n, and d_t are carried C
!  out in source_pp_g.  Constant pressure boundaries are handled by
!  holding the Pp_g at the boundaries zero.  For specified mass flow   C
!  boundaries (part of) a's are calculated here since b is calculated  C
!  from a's in source_pp_g.  After calculating b, a's are multiplied by
!  d and at the flow boundaries get are set to zero.                   C                                                    C
!   Higher order scheme                                                C
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
      SUBROUTINE CONV_PP_G1(A_M, B_M, IER) 
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
      INTEGER          I, J, K, IJK, IJKE, IJKN, IJKT, IPJK, IJPK, IJKP 
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
      INTEGER          M 
! 
!                      local value of A_m 
      DOUBLE PRECISION am 
      
      Integer          incr
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      

      call lock_xsi_array
!
!  Calculate convection factors
!
   
       incr=0	
       CALL CALC_XSI (DISCRETIZE(1), ROP_G, U_G, V_G, W_G, XSI_E, XSI_N,& 
	XSI_T,incr) 
!
!
!  Calculate convection fluxes through each of the faces
!
!$omp  parallel do private( I, J, K, IJK, IJKE, IJKN, IJKT, IPJK, IJPK, &
!$omp&  IJKP, M, am, &
!$omp&  IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLUID_AT(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            AM=((ONE-XSI_E(IJK))*ROP_G(IJK)+XSI_E(IJK)*ROP_G(IJKE))*AYZ(IJK) 
            A_M(IJK,E,0) = AM 
            A_M(IPJK,W,0) = AM 
!
!         North face (i, j+1/2, k)
            AM=((ONE-XSI_N(IJK))*ROP_G(IJK)+XSI_N(IJK)*ROP_G(IJKN))*AXZ(IJK) 
            A_M(IJK,N,0) = AM 
            A_M(IJPK,S,0) = AM 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               AM = ((ONE - XSI_T(IJK))*ROP_G(IJK)+XSI_T(IJK)*ROP_G(IJKT))*AXY(&
                  IJK) 
               A_M(IJK,T,0) = AM 
               A_M(IJKP,B,0) = AM 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               AM = ((ONE - XSI_E(IMJK))*ROP_G(IJKW)+XSI_E(IMJK)*ROP_G(IJK))*&
                  AYZ(IMJK) 
               A_M(IJK,W,0) = AM 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               AM = ((ONE - XSI_N(IJMK))*ROP_G(IJKS)+XSI_N(IJMK)*ROP_G(IJK))*&
                  AXZ(IJMK) 
               A_M(IJK,S,0) = AM 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  AM = ((ONE - XSI_T(IJKM))*ROP_G(IJKB)+XSI_T(IJKM)*ROP_G(IJK))&
                     *AXY(IJKM) 
                  A_M(IJK,B,0) = AM 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      DO M = 1, MMAX 
         IF (.NOT.CLOSE_PACKED(M)) THEN 
!
!         Calculate convection factors
!
	    incr=0
            CALL CALC_XSI (DISCRETIZE(2), ROP_S(1,M), U_S(1,M), V_S(1,M), W_S(1&
               ,M), XSI_E, XSI_N, XSI_T,incr) 

            DO IJK = ijkstart3, ijkend3
!
               IF (FLUID_AT(IJK)) THEN 
                  I = I_OF(IJK) 
                  J = J_OF(IJK) 
                  K = K_OF(IJK) 
                  IJKE = EAST_OF(IJK) 
                  IJKN = NORTH_OF(IJK) 
                  IJKT = TOP_OF(IJK) 
                  IPJK = IP_OF(IJK) 
                  IJPK = JP_OF(IJK) 
                  IJKP = KP_OF(IJK) 
!
!             East face (i+1/2, j, k)
                  AM = ((ONE - XSI_E(IJK))*ROP_S(IJK,M)+XSI_E(IJK)*ROP_S(IJKE,M&
                     ))*AYZ(IJK) 
                  A_M(IJK,E,M) = AM 
                  A_M(IPJK,W,M) = AM 
!
!             North face (i, j+1/2, k)
                  AM = ((ONE - XSI_N(IJK))*ROP_S(IJK,M)+XSI_N(IJK)*ROP_S(IJKN,M&
                     ))*AXZ(IJK) 
                  A_M(IJK,N,M) = AM 
                  A_M(IJPK,S,M) = AM 
!
!             Top face (i, j, k+1/2)
                  IF (DO_K) THEN 
                     AM = ((ONE - XSI_T(IJK))*ROP_S(IJK,M)+XSI_T(IJK)*ROP_S(&
                        IJKT,M))*AXY(IJK) 
                     A_M(IJK,T,M) = AM 
                     A_M(IJKP,B,M) = AM 
                  ENDIF 
!
!             West face (i-1/2, j, k)
                  IMJK = IM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IMJK)) THEN 
                     IJKW = WEST_OF(IJK) 
                     AM = ((ONE - XSI_E(IMJK))*ROP_S(IJKW,M)+XSI_E(IMJK)*ROP_S(&
                        IJK,M))*AYZ(IMJK) 
                     A_M(IJK,W,M) = AM 
                  ENDIF 
!
!             South face (i, j-1/2, k)
                  IJMK = JM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IJMK)) THEN 
                     IJKS = SOUTH_OF(IJK) 
                     AM = ((ONE - XSI_N(IJMK))*ROP_S(IJKS,M)+XSI_N(IJMK)*ROP_S(&
                        IJK,M))*AXZ(IJMK) 
                     A_M(IJK,S,M) = AM 
                  ENDIF 
!
!             Bottom face (i, j, k-1/2)
                  IF (DO_K) THEN 
                     IJKM = KM_OF(IJK) 
                     IF (.NOT.FLUID_AT(IJKM)) THEN 
                        IJKB = BOTTOM_OF(IJK) 
                        AM = ((ONE - XSI_T(IJKM))*ROP_S(IJKB,M)+XSI_T(IJKM)*&
                           ROP_S(IJK,M))*AXY(IJKM) 
                        A_M(IJK,B,M) = AM 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 
      
      call unlock_xsi_array
 
           
      RETURN  
      END SUBROUTINE CONV_PP_G1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
