!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_OUTFLOW(L)                                        C
!  Purpose: Calculate mass and volumetric out flow rates from an       C
!           outflow boundary                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-OCT-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE CALC_OUTFLOW(L) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE bc
      USE fldvar
      USE indices
      USE physprop
      USE geometry
      USE compar    !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK
!
!                      Solids phase
      INTEGER          M
!
      INTEGER IJK2
!
!  Passed variables
!                      Boundary condition number
      INTEGER          L
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      BC_OUT_N(L) = BC_OUT_N(L) + 1 
      DO K = BC_K_B(L), BC_K_T(L) 
         DO J = BC_J_S(L), BC_J_N(L) 
            DO I = BC_I_W(L), BC_I_E(L) 
               IJK = FUNIJK(I,J,K) 
               SELECT CASE (BC_PLANE(L))  
               CASE ('W')  
                  IJK2 = IM_OF(IJK) 
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY(J)*X_E(I-1)*DZ(K)*U_G(IJK2)*&
                     ROP_G(IJK2) 
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY(J)*X_E(I-1)*DZ(K)*U_G(IJK2)*&
                     EP_G(IJK2) 
               CASE ('E')  
                  IJK2 = IP_OF(IJK) 
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DY(J)*X_E(I)*DZ(K)*U_G(IJK)*&
                     ROP_G(IJK2) 
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DY(J)*X_E(I)*DZ(K)*U_G(IJK)*&
                     EP_G(IJK2) 
               CASE ('S')  
                  IJK2 = JM_OF(IJK) 
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX(I)*X(I)*DZ(K)*V_G(IJK2)*&
                     ROP_G(IJK2) 
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX(I)*X(I)*DZ(K)*V_G(IJK2)*EP_G&
                     (IJK2) 
               CASE ('N')  
                  IJK2 = JP_OF(IJK) 
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX(I)*X(I)*DZ(K)*V_G(IJK)*ROP_G&
                     (IJK2) 
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX(I)*X(I)*DZ(K)*V_G(IJK)*EP_G(&
                     IJK2) 
               CASE ('B')  
                  IJK2 = KM_OF(IJK) 
                  BC_MOUT_G(L) = BC_MOUT_G(L) + DX(I)*DY(J)*W_G(IJK2)*ROP_G(&
                     IJK2) 
                  BC_VOUT_G(L)=BC_VOUT_G(L)+DX(I)*DY(J)*W_G(IJK2)*EP_G(IJK2) 
               CASE ('T')  
                  IJK2 = KP_OF(IJK) 
                  BC_MOUT_G(L)=BC_MOUT_G(L)+DX(I)*DY(J)*W_G(IJK)*ROP_G(IJK2) 
                  BC_VOUT_G(L) = BC_VOUT_G(L) + DX(I)*DY(J)*W_G(IJK)*EP_G(IJK2) 
               END SELECT 
               DO M = 1, MMAX 
                  SELECT CASE (BC_PLANE(L))  
                  CASE ('W')  
                     IJK2 = IM_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DY(J)*X_E(I-1)*DZ(K)*U_S&
                        (IJK2,M)*ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DY(J)*X_E(I-1)*DZ(K)*U_S&
                        (IJK2,M)*EP_S(IJK2,M) 
                  CASE ('E')  
                     IJK2 = IP_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DY(J)*X_E(I)*DZ(K)*U_S(&
                        IJK,M)*ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DY(J)*X_E(I)*DZ(K)*U_S(&
                        IJK,M)*EP_S(IJK2,M) 
                  CASE ('S')  
                     IJK2 = JM_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DX(I)*X(I)*DZ(K)*V_S(&
                        IJK2,M)*ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DX(I)*X(I)*DZ(K)*V_S(&
                        IJK2,M)*EP_S(IJK2,M) 
                  CASE ('N')  
                     IJK2 = JP_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DX(I)*X(I)*DZ(K)*V_S(IJK&
                        ,M)*ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DX(I)*X(I)*DZ(K)*V_S(IJK&
                        ,M)*EP_S(IJK2,M) 
                  CASE ('B')  
                     IJK2 = KM_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DX(I)*DY(J)*W_S(IJK2,M)*&
                        ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DX(I)*DY(J)*W_S(IJK2,M)*&
                        EP_S(IJK2,M) 
                  CASE ('T')  
                     IJK2 = KP_OF(IJK) 
                     BC_MOUT_S(L,M) = BC_MOUT_S(L,M) + DX(I)*DY(J)*W_S(IJK,M)*&
                        ROP_S(IJK2,M) 
                     BC_VOUT_S(L,M) = BC_VOUT_S(L,M) + DX(I)*DY(J)*W_S(IJK,M)*&
                        EP_S(IJK2,M) 
                  END SELECT 
               END DO 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE CALC_OUTFLOW 
