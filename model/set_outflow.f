!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2)                C
!  Purpose: Set specified pressure outflow bc for a specified range of C
!           cells                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX                                          C
!                                                                      C
!  Variables modified: I, J, K, RO_g, ROP_g,                           C
!                      EP_g, C
!                      T_g, T_s,  M, ROP_s, U_g, U_s, V_g, V_s,  C
!                      W_g, W_s,
!                                                                      C
!  Local variables: IJK, LFLUID                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2) 
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
      USE geometry
      USE indices
      USE physprop
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Starting I index 
      INTEGER          I1 
! 
!                      Ending I index 
      INTEGER          I2 
! 
!                      Starting J index 
      INTEGER          J1 
! 
!                      Ending J index 
      INTEGER          J2 
! 
!                      Starting K index 
      INTEGER          K1 
! 
!                      Ending K index 
      INTEGER          K2 
! 
!                      indices 
      INTEGER          I, J, K, M, N 
! 
! 
!                      Local index for boundary cell 
      INTEGER          IJK 
! 
!                      Boundary condition number 
      INTEGER          BCV 
! 
!                      Locall index for a fluid cell near the boundary cell 
      INTEGER          LFLUID 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
!//SP Check if current i,j,k resides on this PE
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K) 
!
! Fluid cell at West
!
               IF (FLUID_AT(IM_OF(IJK))) THEN 
                  LFLUID = IM_OF(IJK) 
                  IF (U_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (U_S(LFLUID,M) >= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
!
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     U_G(IJK) = ZERO 
                  ENDIF 
                  V_G(IJK) = V_G(LFLUID) 
                  W_G(IJK) = W_G(LFLUID) 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     WHERE (ROP_S(IJK,:MMAX) > ZERO)  
                        U_S(IJK,:MMAX) = ROP_S(LFLUID,:MMAX)*U_S(LFLUID,:MMAX)/&
                           ROP_S(IJK,:MMAX) 
                     ELSEWHERE 
                        U_S(IJK,:MMAX) = ZERO 
                     END WHERE 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     M = MMAX + 1 
                  ENDIF 
               ENDIF 
!
! Fluid cell at East
!
               IF (FLUID_AT(IP_OF(IJK))) THEN 
                  LFLUID = IP_OF(IJK) 
                  IF (U_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (U_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  IF (U_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        U_G(IJK) = ROP_G(LFLUID)*U_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        U_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  V_G(IJK) = V_G(LFLUID) 
                  W_G(IJK) = W_G(LFLUID) 
                  DO M = 1, MMAX 
                     IF (U_S(IJK,M) == UNDEFINED) THEN 
                        IF (ROP_S(IJK,M) > ZERO) THEN 
                           U_S(IJK,M) = ROP_S(LFLUID,M)*U_S(LFLUID,M)/ROP_S(IJK&
                              ,M) 
                        ELSE 
                           U_S(IJK,M) = ZERO 
                        ENDIF 
                     ENDIF 
                     V_S(IJK,M) = V_S(LFLUID,M) 
                     W_S(IJK,M) = W_S(LFLUID,M) 
                  END DO 
               ENDIF 
!
! Fluid cell at South
!
               IF (FLUID_AT(JM_OF(IJK))) THEN 
                  LFLUID = JM_OF(IJK) 
                  IF (V_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (V_S(LFLUID,M) >= 0.) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     V_G(IJK) = ZERO 
                  ENDIF 
                  W_G(IJK) = W_G(LFLUID) 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     M = MMAX + 1 
                  ENDIF 
               ENDIF 
!
! Fluid cell at North
!
               IF (FLUID_AT(JP_OF(IJK))) THEN 
                  LFLUID = JP_OF(IJK) 
                  IF (V_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (V_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  IF (V_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        V_G(IJK) = ROP_G(LFLUID)*V_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        V_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  W_G(IJK) = W_G(LFLUID) 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     WHERE (V_S(IJK,:MMAX) == UNDEFINED) V_S(IJK,:MMAX) = V_S(&
                        LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     M = MMAX + 1 
                  ENDIF 
               ENDIF 
!
! Fluid cell at Bottom
!
               IF (FLUID_AT(KM_OF(IJK))) THEN 
                  LFLUID = KM_OF(IJK) 
                  IF (W_G(LFLUID)>=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (W_S(LFLUID,M) >= 0.) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  V_G(IJK) = V_G(LFLUID) 
                  IF (ROP_G(IJK) > ZERO) THEN 
                     W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK) 
                  ELSE 
                     W_G(IJK) = ZERO 
                  ENDIF 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     W_S(IJK,:MMAX) = W_S(LFLUID,:MMAX) 
                     M = MMAX + 1 
                  ENDIF 
               ENDIF 
!
! Fluid cell at Top
!
               IF (FLUID_AT(KP_OF(IJK))) THEN 
                  LFLUID = KP_OF(IJK) 
                  IF (W_G(IJK)<=ZERO .OR. EP_G(IJK)==UNDEFINED) THEN 
                     IF (BC_TYPE(BCV) /= 'P_OUTFLOW') P_G(IJK) = P_G(LFLUID) 
                     T_G(IJK) = T_G(LFLUID) 
                     N = 1 
                     IF (NMAX(0) > 0) THEN 
                        X_G(IJK,:NMAX(0)) = X_G(LFLUID,:NMAX(0)) 
                        N = NMAX(0) + 1 
                     ENDIF 
                     MW_MIX_G(IJK) = MW_MIX_G(LFLUID) 
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G(IJK),P_G&
                        (IJK),T_G(IJK)) 
                  ENDIF 
                  P_STAR(IJK) = P_STAR(LFLUID) 
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE 
                  DO M = 1, MMAX 
                     P_S(IJK,M) = P_S(LFLUID,M) 
                     IF (W_S(IJK,M) <= ZERO) THEN 
                        ROP_S(IJK,M) = ROP_S(LFLUID,M) 
                        T_S(IJK,M) = T_S(LFLUID,M) 
                     ELSE 
                        ROP_S(IJK,M) = ZERO 
                     ENDIF 
!
                     IF(BC_ROP_S(BCV,M)/=UNDEFINED)ROP_S(IJK,M)=BC_ROP_S(BCV,M) 
!
                     IF(BC_EP_G(BCV)==UNDEFINED)EP_G(IJK)=EP_G(IJK)-EP_S(IJK,M) 
                     N = 1 
                     IF (NMAX(M) > 0) THEN 
                        X_S(IJK,M,:NMAX(M)) = X_S(LFLUID,M,:NMAX(M)) 
                        N = NMAX(M) + 1 
                     ENDIF 
                  END DO 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
                  U_G(IJK) = U_G(LFLUID) 
                  V_G(IJK) = V_G(LFLUID) 
                  IF (W_G(IJK) == UNDEFINED) THEN 
                     IF (ROP_G(IJK) > ZERO) THEN 
                        W_G(IJK) = ROP_G(LFLUID)*W_G(LFLUID)/ROP_G(IJK) 
                     ELSE 
                        W_G(IJK) = ZERO 
                     ENDIF 
                  ENDIF 
                  M = 1 
                  IF (MMAX > 0) THEN 
                     U_S(IJK,:MMAX) = U_S(LFLUID,:MMAX) 
                     V_S(IJK,:MMAX) = V_S(LFLUID,:MMAX) 
                     WHERE (W_S(IJK,:MMAX) == UNDEFINED) W_S(IJK,:MMAX) = W_S(&
                        LFLUID,:MMAX) 
                     M = MMAX + 1 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE SET_OUTFLOW 
