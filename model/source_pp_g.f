!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_Pp_g(A_m, B_m, IER)
!  Purpose: Determine source terms for Pressure                        C
!  correction equation.  The off-diagonal coefficients are             C
!   positive. The center coefficient and the source vector are         C
!  negative.                                                           C
!  See conv_Pp_g
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
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
      SUBROUTINE SOURCE_PP_G(A_M, B_M, IER) 
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
      USE physprop
      USE fldvar
      USE rxns
      USE run
      USE geometry
      USE indices
      USE pgcor
      USE bc
      USE vshear
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
      INTEGER          IJK, IMJK, IJMK, IJKM, I, J, K 
      INTEGER          M, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB 
!
! loezos
      integer incr
 
!                      error message 
      CHARACTER*80     LINE(1) 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DROODP_G 
!-----------------------------------------------
 
!// Need to extract i, j, k from ijk_p_g to determine the processor which
!     acts on ijk_p_g to fix the value of pressure
!       ----------------------------------------------------------
!       inline functions for determining i, j, k for global ijk_p_g
!       -----------------------------------------------------------
        integer i_of_g,j_of_g,k_of_g
        integer ijk_p_g_local

      INCLUDE 'function.inc'

        k_of_g(ijk) = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
        j_of_g(ijk) = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                      - 1)/(imax3-imin3+1)) + jmin3
        i_of_g(ijk) = ijk - (j_of_g(ijk)-jmin3)*(imax3-imin3+1) - &
                      (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + imin3

! loezos
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK) 
	 DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN  
	   V_G(IJK)=V_G(IJK)+VSH(IJK)	
         END IF
        END DO 
 
      END IF

      call lock_xsi_array

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp    parallel do private(IJK, IMJK, IJMK, IJKM, M)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
            B_M(IJK,0) = -((-((ROP_G(IJK)-ROP_GO(IJK))*VOL(IJK)*ODT+A_M(IJK,E,0&
               )*U_G(IJK)-A_M(IJK,W,0)*U_G(IMJK)+A_M(IJK,N,0)*V_G(IJK)-A_M(IJK,&
               S,0)*V_G(IJMK)+A_M(IJK,T,0)*W_G(IJK)-A_M(IJK,B,0)*W_G(IJKM)))+&
               SUM_R_G(IJK)*VOL(IJK)) 
            A_M(IJK,E,0) = A_M(IJK,E,0)*D_E(IJK,0) 
            A_M(IJK,W,0) = A_M(IJK,W,0)*D_E(IMJK,0) 
            A_M(IJK,N,0) = A_M(IJK,N,0)*D_N(IJK,0) 
            A_M(IJK,S,0) = A_M(IJK,S,0)*D_N(IJMK,0) 
            A_M(IJK,T,0) = A_M(IJK,T,0)*D_T(IJK,0) 
            A_M(IJK,B,0) = A_M(IJK,B,0)*D_T(IJKM,0) 
            DO M = 1, MMAX 
               IF (.NOT.CLOSE_PACKED(M)) THEN 
                  B_M(IJK,0) = B_M(IJK,0) - ((-((ROP_S(IJK,M)-ROP_SO(IJK,M))*&
                     VOL(IJK)*ODT+A_M(IJK,E,M)*U_S(IJK,M)-A_M(IJK,W,M)*U_S(IMJK&
                     ,M)+A_M(IJK,N,M)*V_S(IJK,M)-A_M(IJK,S,M)*V_S(IJMK,M)+A_M(&
                     IJK,T,M)*W_S(IJK,M)-A_M(IJK,B,M)*W_S(IJKM,M)))+SUM_R_S(IJK&
                     ,M)*VOL(IJK)) 
                  A_M(IJK,E,0) = A_M(IJK,E,0) + A_M(IJK,E,M)*D_E(IJK,M) 
                  A_M(IJK,W,0) = A_M(IJK,W,0) + A_M(IJK,W,M)*D_E(IMJK,M) 
                  A_M(IJK,N,0) = A_M(IJK,N,0) + A_M(IJK,N,M)*D_N(IJK,M) 
                  A_M(IJK,S,0) = A_M(IJK,S,0) + A_M(IJK,S,M)*D_N(IJMK,M) 
                  A_M(IJK,T,0) = A_M(IJK,T,0) + A_M(IJK,T,M)*D_T(IJK,M) 
                  A_M(IJK,B,0) = A_M(IJK,B,0) + A_M(IJK,B,M)*D_T(IJKM,M) 
               ENDIF 
            END DO 
            A_M(IJK,0,0) = -(A_M(IJK,E,0)+A_M(IJK,W,0)+A_M(IJK,N,0)+A_M(IJK,S,0&
               )+A_M(IJK,T,0)+A_M(IJK,B,0)) 
!
            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,0) = -ONE 
                  B_M(IJK,0) = ZERO 
               ELSE 
!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0) 
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1) 
!$omp             end critical
               ENDIF 
            ENDIF 
!
         ELSE 
            A_M(IJK,E,0) = ZERO 
            A_M(IJK,W,0) = ZERO 
            A_M(IJK,N,0) = ZERO 
            A_M(IJK,S,0) = ZERO 
            A_M(IJK,T,0) = ZERO 
            A_M(IJK,B,0) = ZERO 
            A_M(IJK,0,0) = -ONE 
            B_M(IJK,0) = ZERO 
         ENDIF 
      END DO 

! loezos
      IF (SHEAR) THEN
!$omp parallel do private(IJK) 
	 DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN  
	   V_G(IJK)=V_G(IJK)-VSH(IJK)	
         END IF
       END DO 
      END IF
! loezos

      IF (RO_G0 == UNDEFINED) THEN 

! loezos
	incr=0		
! loezos
         CALL CALC_XSI(DISCRETIZE(1),ROP_G,U_G,V_G,W_G,XSI_E,XSI_N,XSI_T,incr) 
	
!$omp    parallel do                                                     &
!$omp&   private(IJK,I,J,K,                                       &
!$omp&            IMJK,IJMK,IJKM,IJKE,IJKW,IJKN,IJKS,IJKT,IJKB)
         DO IJK = ijkstart3, ijkend3 
            IF (FLUID_AT(IJK)) THEN 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IMJK = IM_OF(IJK) 
               IJMK = JM_OF(IJK) 
               IJKM = KM_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               IJKW = WEST_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               IJKS = SOUTH_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKB = BOTTOM_OF(IJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) - DROODP_G(RO_G(IJK),P_G(IJK))*EP_G(&
                  IJK)*((ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK)-XSI_E(IMJK)*U_G(&
                  IMJK)*AYZ(IMJK)+(ONE-XSI_N(IJK))*V_G(IJK)*AXZ(IJK)-XSI_N(IJMK&
                  )*V_G(IJMK)*AXZ(IJMK)+VOL(IJK)*ODT) 
!
               A_M(IJK,E,0) = A_M(IJK,E,0) - EP_G(IJKE)*DROODP_G(RO_G(IJKE),P_G&
                  (IJKE))*XSI_E(IJK)*U_G(IJK)*AYZ(IJK) 
               A_M(IJK,W,0) = A_M(IJK,W,0) + EP_G(IJKW)*DROODP_G(RO_G(IJKW),P_G&
                  (IJKW))*(ONE - XSI_E(IMJK))*U_G(IMJK)*AYZ(IMJK) 
               A_M(IJK,N,0) = A_M(IJK,N,0) - EP_G(IJKN)*DROODP_G(RO_G(IJKN),P_G&
                  (IJKN))*XSI_N(IJK)*V_G(IJK)*AXZ(IJK) 
               A_M(IJK,S,0) = A_M(IJK,S,0) + EP_G(IJKS)*DROODP_G(RO_G(IJKS),P_G&
                  (IJKS))*(ONE - XSI_N(IJMK))*V_G(IJMK)*AXZ(IJMK) 
               IF (DO_K) THEN 
                  A_M(IJK,0,0) = A_M(IJK,0,0) - DROODP_G(RO_G(IJK),P_G(IJK))*&
                     EP_G(IJK)*((ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK)-XSI_T(IJKM&
                     )*W_G(IJKM)*AXY(IJKM)) 
                  A_M(IJK,T,0) = A_M(IJK,T,0) - EP_G(IJKT)*DROODP_G(RO_G(IJKT),&
                     P_G(IJKT))*XSI_T(IJK)*W_G(IJK)*AXY(IJK) 
                  A_M(IJK,B,0) = A_M(IJK,B,0) + EP_G(IJKB)*DROODP_G(RO_G(IJKB),&
                     P_G(IJKB))*(ONE - XSI_T(IJKM))*W_G(IJKM)*AXY(IJKM) 
               ENDIF 
!
            ENDIF 
         END DO 
      ENDIF 
!
!
!  Specify P' to zero at a certain location for incompressible flows and
!  cyclic boundary conditions.
!
!// Parallel implementation of fixing a pressure at a point
   I = I_OF_G(IJK_P_G)
   J = J_OF_G(IJK_P_G)
   K = K_OF_G(IJK_P_G)
   IF(IS_ON_myPE_OWNS(I,J,K)) THEN
      IF (IJK_P_G /= UNDEFINED_I) THEN 
         IJK_P_G_LOCAL = FUNIJK(I,J,K)
         A_M(IJK_P_G_LOCAL,E,0) = ZERO 
         A_M(IJK_P_G_LOCAL,W,0) = ZERO 
         A_M(IJK_P_G_LOCAL,N,0) = ZERO 
         A_M(IJK_P_G_LOCAL,S,0) = ZERO 
         A_M(IJK_P_G_LOCAL,T,0) = ZERO 
         A_M(IJK_P_G_LOCAL,B,0) = ZERO 
         A_M(IJK_P_G_LOCAL,0,0) = -ONE 
         B_M(IJK_P_G_LOCAL,0) = ZERO 
      ENDIF 
   ENDIF
!
      call unlock_xsi_array

      RETURN  
      END SUBROUTINE SOURCE_PP_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

