!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_e_e(A_m,  e_e, IER)                               C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- East                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
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
      SUBROUTINE CALC_E_E(A_M, MCP, E_E, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar        !//d
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
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

!                      Index of close packed solids phase
      INTEGER          Mcp
!                      
      DOUBLE PRECISION e_e(DIMENSION_3)
!
!                      Indices
      INTEGER          IJK
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (.NOT.MOMENTUM_X_EQ(MCP)) RETURN  
!
!
!!$omp parallel do private(ijk)
      DO IJK = 1, IJKMAX2 
         IF (SIP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN 
            E_E(IJK) = ZERO 
         ELSE 
!
            IF (A_M(IJK,0,MCP) /= ZERO) THEN 
               E_E(IJK) = AYZ(IJK)/(-A_M(IJK,0,MCP)) 
            ELSE 
               E_E(IJK) = LARGE_NUMBER 
            ENDIF 
!
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_E_E 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_e_n(A_m,  e_n, IER)                               C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- North                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
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
      SUBROUTINE CALC_E_N(A_M, MCP, E_N, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar    !//d
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
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

!                      Index of close packed solids phase
      INTEGER          Mcp
!                      
      DOUBLE PRECISION e_n(DIMENSION_3)
!
!                      Indices
      INTEGER          I, K, IJK
!
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (.NOT.MOMENTUM_Y_EQ(MCP)) RETURN  
!
!
!!$omp parallel do private(IJK,I,K)
      DO IJK = 1, IJKMAX2 
         IF (SIP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
            E_N(IJK) = ZERO 
         ELSE 
            I = I_OF(IJK) 
            K = K_OF(IJK) 
!
            IF ((-A_M(IJK,0,MCP)) /= ZERO) THEN 
               E_N(IJK) = AXZ(IJK)/(-A_M(IJK,0,MCP)) 
            ELSE 
               E_N(IJK) = LARGE_NUMBER 
            ENDIF 
!
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_E_N 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_e_t(A_m,  e_t, IER)                               C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- Top                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
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
      SUBROUTINE CALC_E_T(A_M, MCP, E_T, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar   !//d
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
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

!                      Index of close packed solids phase
      INTEGER          Mcp
!                      
      DOUBLE PRECISION e_t(DIMENSION_3)
!
!                      Indices
      INTEGER          I, J, IJK
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (.NOT.MOMENTUM_Z_EQ(MCP)) RETURN  
!
!
!!$omp parallel do private(I,J,IJK)
      DO IJK = 1, IJKMAX2 
         IF (SIP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN 
            E_T(IJK) = ZERO 
         ELSE 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            IF ((-A_M(IJK,0,MCP)) /= ZERO) THEN 
               E_T(IJK) = AXY(IJK)/(-A_M(IJK,0,MCP)) 
            ELSE 
               E_T(IJK) = LARGE_NUMBER 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CALC_E_T 
