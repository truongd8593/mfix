!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP(IER)                                          C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. Master routine.                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
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
      SUBROUTINE CONV_ROP(IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE mflux
      USE physprop
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
!                      solids phase index 
      INTEGER          M 
! 
! 
!
!
      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_ROP0 (ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, IER) 
      ELSE 
         CALL CONV_ROP1 (DISCRETIZE(1), ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, IER) 
      ENDIF 
      
      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding 
        DO M = 1, MMAX
          CALL CONV_ROP0 (ROP_s(1, M), U_s(1, M), V_s(1, M), W_s(1, M), &
	                  ROP_sE(1, M), ROP_sN(1, M), ROP_sT(1, M), IER) 
        ENDDO
      ELSE 
        DO M = 1, MMAX
          CALL CONV_ROP1 (DISCRETIZE(2), ROP_s(1, M), U_s(1, M), V_s(1, M), W_s(1, M), &
	                  ROP_sE(1, M), ROP_sN(1, M), ROP_sT(1, M), IER) 
        ENDDO
      ENDIF 
      
      RETURN  
      END SUBROUTINE CONV_ROP 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER)      C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. FOU routine.                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
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
      SUBROUTINE CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run 
      USE parallel 
      USE physprop
      USE geometry
      USE indices
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
!                      macroscopic density (rho_prime)
      DOUBLE PRECISION ROP(DIMENSION_3) 
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3) 
!
!                      Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION ROP_E(DIMENSION_3), ROP_N(DIMENSION_3), ROP_T(DIMENSION_3) 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'
      

!
!  Interpolate the face value of density for calculating the convection fluxes 
!!!$omp  parallel do private( IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            IF (U(IJK) >= ZERO) THEN 
               ROP_E(IJK) = ROP(IJK) 
            ELSE 
               ROP_E(IJK) = ROP(IJKE) 
            ENDIF 
!
!         North face (i, j+1/2, k)
            IF (V(IJK) >= ZERO) THEN 
               ROP_N(IJK) = ROP(IJK) 
            ELSE 
               ROP_N(IJK) = ROP(IJKN) 
            ENDIF 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  ROP_T(IJK) = ROP(IJK) 
               ELSE 
                  ROP_T(IJK) = ROP(IJKT) 
               ENDIF 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               IF (U(IMJK) >= ZERO) THEN 
                  ROP_E(IMJK) = ROP(IJKW) 
               ELSE 
                  ROP_E(IMJK) = ROP(IJK) 
               ENDIF 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               IF (V(IJMK) >= ZERO) THEN 
                 ROP_N(IJMK) = ROP(IJKS) 
               ELSE 
                 ROP_N(IJMK) = ROP(IJK) 
               ENDIF 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  IF (W(IJKM) >= ZERO) THEN 
                     ROP_T(IJKM) = ROP(IJKB) 
                  ELSE 
                     ROP_T(IJKM) = ROP(IJK) 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CONV_ROP0 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER)C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. HR routine.                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER) 
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE parallel 
      USE physprop
      USE geometry
      USE indices
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
!
!                      Discretization scheme 
      INTEGER          DISC
!
!                      macroscopic density (rho_prime)
      DOUBLE PRECISION ROP(DIMENSION_3) 
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3) 
!
!                      Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION ROP_E(DIMENSION_3), ROP_N(DIMENSION_3), ROP_T(DIMENSION_3) 
!
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT 
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
      
      Integer          incr
!-----------------------------------------------
      INCLUDE 'function.inc'
      

      call lock_xsi_array
!
!  Calculate convection factors
!
   
       incr=0	
       CALL CALC_XSI (DISC, ROP, U, V, W, XSI_E, XSI_N,	XSI_T,incr) 
!
!
!  Calculate convection fluxes through each of the faces
!
!!!$omp  parallel do private(IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLUID_AT(IJK)) THEN 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            ROP_E(IJK) = ((ONE-XSI_E(IJK))*ROP(IJK)+XSI_E(IJK)*ROP(IJKE)) 
!
!         North face (i, j+1/2, k)
            ROP_N(IJK) = ((ONE-XSI_N(IJK))*ROP(IJK)+XSI_N(IJK)*ROP(IJKN)) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               ROP_T(IJK) = ((ONE - XSI_T(IJK))*ROP(IJK)+XSI_T(IJK)*ROP(IJKT))
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               ROP_E(IMJK) = ((ONE - XSI_E(IMJK))*ROP(IJKW)+XSI_E(IMJK)*ROP(IJK))
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               ROP_N(IJMK) = ((ONE - XSI_N(IJMK))*ROP(IJKS)+XSI_N(IJMK)*ROP(IJK))
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  ROP_T(IJKM) = ((ONE - XSI_T(IJKM))*ROP(IJKB)+XSI_T(IJKM)*ROP(IJK))
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      
      call unlock_xsi_array
 
           
      RETURN  
      END SUBROUTINE CONV_ROP1 
