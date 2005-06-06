!
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0(IER)                                         C
!  Purpose: Correct the fluid pressure and gas and solids velocities   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
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
      SUBROUTINE CORRECT_0(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE pgcor
      USE ur_facs 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                      error index
!
      CALL CORRECT_0G (PP_G, UR_FAC(1), D_E, D_N, D_T, P_G, U_G, V_G, W_G, IER) 
!      CALL CORRECT_0S (PP_G, D_E, D_N, D_T, U_S, V_S, W_S, IER) 
      RETURN  
      END SUBROUTINE CORRECT_0 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0g(Pp_g, UR_fac, d_e, d_n, d_t,                C
!     &                P_g, U_g, V_g, W_g, IER)                        C
!  Purpose: Correct the fluid pressure and velocities.                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
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
      SUBROUTINE CORRECT_0G(PP_G,UR_FAC,D_E,D_N,D_T,P_G,U_G,V_G,W_G,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Pressure correction 
      DOUBLE PRECISION Pp_g(DIMENSION_3) 
! 
!                      Under relaxation factor for Pressure correction 
      DOUBLE PRECISION UR_fac 
! 
!                      Pressure correction coefficient -- East 
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- North 
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- Top 
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure 
      DOUBLE PRECISION P_g(DIMENSION_3) 
! 
!                      Velocity components 
      DOUBLE PRECISION U_g(DIMENSION_3), V_g(DIMENSION_3),& 
                       W_g(DIMENSION_3) 
! 
!                      error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT, M 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Underrelax pressure correction.  Velocity corrections should not be
!  underrelaxed, so that the continuity eq. is satisfied.
!
!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
      DO IJK = ijkstart3, ijkend3 
         IF (FLUIDORP_FLOW_AT(IJK)) THEN 
            P_G(IJK) = P_G(IJK) + UR_FAC*PP_G(IJK) 
!
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            U_G(IJK) = U_G(IJK) - D_E(IJK,0)*(PP_G(IJKE)-PP_G(IJK)) 
            V_G(IJK) = V_G(IJK) - D_N(IJK,0)*(PP_G(IJKN)-PP_G(IJK)) 
            IF (DO_K) THEN 
               IJKT = TOP_OF(IJK) 
               W_G(IJK) = W_G(IJK) - D_T(IJK,0)*(PP_G(IJKT)-PP_G(IJK)) 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CORRECT_0G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0s(Pp_g, d_e, d_n, d_t, U_s, V_s, W_s, IER)    C
!  Purpose: Correct the solids velocities.                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
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
      SUBROUTINE CORRECT_0S(PP_G, D_E, D_N, D_T, U_S, V_S, W_S, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Pressure correction 
      DOUBLE PRECISION Pp_g(DIMENSION_3) 
! 
!                      Pressure correction coefficient -- East 
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- North 
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- Top 
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Velocity components 
      DOUBLE PRECISION U_s(DIMENSION_3, DIMENSION_M),& 
                       V_s(DIMENSION_3, DIMENSION_M),& 
                       W_s(DIMENSION_3, DIMENSION_M) 
! 
!                      error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT, M 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Velocity corrections should not be
!  underrelaxed, so that the continuity eq. is satisfied.
!
      DO M = 1, MMAX 
!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
         DO IJK = ijkstart3, ijkend3
            IF (FLUIDORP_FLOW_AT(IJK)) THEN 
!
               IJKE = EAST_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               U_S(IJK,M) = U_S(IJK,M) - D_E(IJK,M)*(PP_G(IJKE)-PP_G(IJK)) 
               V_S(IJK,M) = V_S(IJK,M) - D_N(IJK,M)*(PP_G(IJKN)-PP_G(IJK)) 
               IF (DO_K) THEN 
                  IJKT = TOP_OF(IJK) 
                  W_S(IJK,M) = W_S(IJK,M) - D_T(IJK,M)*(PP_G(IJKT)-PP_G(IJK)) 
               ENDIF 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE CORRECT_0S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
