CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_MU_s(M, IER)                                      C
C  Purpose: Calculate granular stress terms: THETA, P_s, LAMBDA_s, MU_sC
C                                                                      C
C  Author: W. Rogers                                  Date: 04-mar-92  C
C  Reviewer: M. Syamlal                               Date: 16-MAR-92  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Modifications for cylindrical geometry                     C
C  Author: M. Syamlal                                 Date: 15-MAY-92  C
C  Revision Number: 2                                                  C
C  Purpose: Add volume-weighted averaging statement functions for      C
C           variable grid capability                                   C
C  Author:  W. Rogers                                 Date: 21-JUL-92  C
C  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
C  Revision Number: 3                                                  C
C  Purpose: Add plastic-flow stress terms                              C
C  Author: M. Syamlal                                 Date: 10-FEB-93  C
C  Revision Number: 4                                                  C
C  Purpose: Add Boyle-Massoudi stress terms                            C
C  Author: M. Syamlal                                 Date: 2-NOV-95   C
C  Revision Number: 5                                                  C
C  Purpose: MFIX 2.0 mods  (old name CALC_THETA)                       C
C  Author: M. Syamlal                                 Date: 24-APR_96  C
C  Author: Kapil Agrawal, Princeton University        Date: 6-FEB-98   C
C  Revision Number: 6                                                  C
C  Purpose: Add calculation of viscosities and conductivities for use  C
C           with granular temperature PDE. New common block contained  C
C           in 'trace.inc' contains trD_s_C(DIMENSION_3, DIMENSION_M)  C
C           and trD_s2(DIMENSION_3, DIMENSION_M)                       C
C                                                                      C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: U_s, V_s, W_s, IMAX2, JMAX2, KMAX2, DX, DY,   C
C                        DZ, IMJPK, IMJK, IPJMK, IPJK, IJMK, IJKP,     C
C                        IMJKP, IPJKM, IJKM, IJMKP, IJPK, IJPKM, IJMK, C
C                        M,  RO_s, C_e, D_p, Pi, G_0, X                C
C                                                                      C
C  Variables modified: I, J, K, IJK, MU_s, LAMBDA_s, P_s               C
C                                                                      C
C  Local variables: K_1m, K_2m, K_3m, K_4m, D_s, U_s_N, U_s_S, V_s_E,  C
C                   V_s_W, U_s_T, U_s_B, W_s_E, W_s_W, V_s_T, V_s_B,   C
C                   W_s_N, W_s_S, trD_s_C, W_s_C                       C
C                   trD_s2, EP_s2xTHETA, EP_sxSQRTHETA, I1, I2, U_s_C, C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CALC_MU_s(M, IER)
C
      IMPLICIT NONE
C
C                      Maximum value of solids viscosity in poise
      DOUBLE PRECISION MAX_MU_s
      PARAMETER (MAX_MU_s = 1000.)
C
C     Include param.inc file to specify parameter values
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'parallel.inc'
C
C     Include physical and numerical parameters section
C
      INCLUDE 'physprop.inc'
      INCLUDE 'drag.inc'
      INCLUDE 'run.inc'
C
C     Include geometry and discretization section
C
      INCLUDE 'geometry.inc'
C
C     Include Field Variables
C
      INCLUDE 'fldvar.inc'
      INCLUDE 'visc_g.inc'
      INCLUDE 'visc_s.inc'
      INCLUDE 'trace.inc'
C
C     Include indices
C
      INCLUDE 'indices.inc'
C
C     Include constants
C
      INCLUDE 'constant.inc'
C
C  Function subroutines
C
      DOUBLE PRECISION G_0
C
C     Local Variables
C
C                      Error index
      INTEGER          IER
C
C                      Denotes cell class, column of STORE_INCREMENTS
      INTEGER          ICLASS
C
C                      Constant in equation for mth solids phase pressure
      DOUBLE PRECISION K_1m
C
C                      Constant in equation for mth solids phase bulk viscosity
      DOUBLE PRECISION K_2m
C
C                      Constant in equation for mth solids phase viscosity
      DOUBLE PRECISION K_3m
C
C                      Constant in equation for mth solids phase dissipation
      DOUBLE PRECISION K_4m
C
C                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_s(3,3)
C
C                      U_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_s_N
C
C                      U_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_s_S
C
C                      U_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_s_T
C
C                      U_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_s_B
C
C                      U_s at the center of the THETA cell-(i, j, k)
C                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_s_C
C
C                      V_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_s_E
C
C                      V_s at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_s_W
C
C                      V_s at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_s_T
C
C                      V_s at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_s_B
C
C                      W_s at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_s_E
C
C                      W_s at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_s_W
C
C                      W_s at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_s_N
C
C                      W_s at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_s_S
C
C                      W_s at the center of the THETA cell-(i, j, k).
C                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_s_C
C
C                      Value of EP_s * SQRT( THETA )for Mth solids phase
C                      continuum
      DOUBLE PRECISION EP_sxSQRTHETA
C
C                      Value of EP_s * EP_s * THETA for Mth solids phase
C                      continuum
      DOUBLE PRECISION EP_s2xTHETA
C
C                      Local DO-LOOP counters
      INTEGER          I1, I2
C
C                      Second invariant of the deviator of D_s
      DOUBLE PRECISION I2_devD_s
C
C                      Factor in plastic-flow stress terms
      DOUBLE PRECISION qxP_s
C
C                      Coefficients of quadratic equation
      DOUBLE PRECISION aq, bq, cq
C
C                      Constant in Boyle-Massoudi stress term
      DOUBLE PRECISION K_5m
C
C                      d(EP_sm)/dX
      DOUBLE PRECISION DEP_soDX
C
C                      d(EP_sm)/dY
      DOUBLE PRECISION DEP_soDY
C
C                      d(EP_sm)/XdZ
      DOUBLE PRECISION DEP_soXDZ
C
C                      Solids volume fraction gradient tensor
      DOUBLE PRECISION M_s(3,3)
C
C                      Trace of M_s
      DOUBLE PRECISION trM_s
C
C                      Trace of (D_s)(M_s)
      DOUBLE PRECISION trDM_s
C
C                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,
     &                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,
     &                 IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP,
     &                 IJMKM, IJPKM
C
C                      Solids phase
      INTEGER          M
C
C                      Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION Mu, Mu_b, Mu_star,Kth,Kth_star

C     SWITCH enables us to turn on/off the modification to the
C     particulate phase viscosity. If we want to simulate gas-particle
C     flow then SWITCH=1 to incorporate the effect of drag on the
C     particle viscosity. If we want to simulate granular flow
C     without the effects of an interstitial gas, SWITCH=0.
C     (Same for conductivity)

C  Function subroutines
C
C
C                      dg0/dep
      DOUBLE PRECISION DG_0DNU
C
C     Include statement functions
C
C
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
C
C
C
C$DOACROSS LOCAL(I, J, K, IJK, ICLASS, IMJK, IPJK, IJMK, IJPK, IJKM,
C$&  IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, IM, JM, KM, 
C$&  U_s_N, U_s_S, U_s_T, U_s_B, V_s_E, V_s_W, V_s_T, V_s_B, W_s_N,
C$&  W_s_S, W_s_E, W_s_W, U_s_C, W_s_C, D_s, I2_devD_s, trD_s_C,
C$&  qxP_s, trD_s2, K_1m, K_2m, K_3m, K_4m, K_5m, aq, bq, cq, 
C$&  DEP_soDX, DEP_soDY, DEP_soXDZ, M_s, trM_s, trDM_s, I1, I2,
C$&  EP_sxSQRTHETA, EP_s2xTHETA 
C$&  ),  
C & mp_SCHEDTYPE = INTERLEAVE, CHUNK = chunk_size
C$& mp_SCHEDTYPE = GSS
      DO 200 IJK = IJKMIN1, IJKMAX1
C
        IF ( FLUID_AT(IJK) ) THEN
C
C------------------------------------------------------------------------
c          CALL SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,
c     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
c     &                       IJKB, IJKT, IM, JM, KM)
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)
          IM = Im1(I)
          JM = Jm1(J)
          KM = Km1(K)
          ICLASS = CELL_CLASS(IJK) !Determine the class of cell IJK
          IJKW  = IJK + INCREMENT_FOR_w (ICLASS)
          IJKE  = IJK + INCREMENT_FOR_e (ICLASS)
          IJKS  = IJK + INCREMENT_FOR_s (ICLASS)
          IJKN  = IJK + INCREMENT_FOR_n (ICLASS)
          IJKB  = IJK + INCREMENT_FOR_b (ICLASS)
          IJKT  = IJK + INCREMENT_FOR_t (ICLASS)
          IMJK  = IJK + INCREMENT_FOR_im(ICLASS)
          IPJK  = IJK + INCREMENT_FOR_ip(ICLASS)
          IJMK  = IJK + INCREMENT_FOR_jm(ICLASS)
          IJPK  = IJK + INCREMENT_FOR_jp(ICLASS)
          IJKM  = IJK + INCREMENT_FOR_km(ICLASS)
          IJKP  = IJK + INCREMENT_FOR_kp(ICLASS)
          IMJPK = IM_OF(IJPK)
          IMJMK = IM_OF(IJMK)
          IMJKP = IM_OF(IJKP)
          IMJKM = IM_OF(IJKM)
          IPJKM = IP_OF(IJKM)
          IPJMK = IP_OF(IJMK)
          IJMKP = JM_OF(IJKP)
          IJMKM = JM_OF(IJKM)
          IJPKM = JP_OF(IJKM)

          U_s_N = AVG_Y(                                 !i, j+1/2, k
     &             AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),
     &             AVG_X_E(U_s(IMJPK, M), U_s(IJPK, M), I), J
     &           )
          U_s_S = AVG_Y(                                 !i, j-1/2, k
     &             AVG_X_E(U_s(IMJMK, M), U_s(IJMK, M), I),
     &             AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), JM
     &           )
          U_s_T = AVG_Z(                                 !i, j, k+1/2
     &             AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),
     &             AVG_X_E(U_s(IMJKP, M), U_s(IJKP, M), I), K
     &           )
          U_s_B = AVG_Z(                                 !i, j, k-1/2
     &             AVG_X_E(U_s(IMJKM, M), U_s(IJKM, M), I),
     &             AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), KM
     &           )
          V_s_E = AVG_X(                                 !i+1/2, j, k
     &             AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),
     &             AVG_Y_N(V_s(IPJMK, M), V_s(IPJK, M)), I
     &           )
          V_s_W = AVG_X(                                 !i-1/2, j, k
     &             AVG_Y_N(V_s(IMJMK, M), V_s(IMJK, M)),
     &             AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM
     &           )
          V_s_T = AVG_Z(                                 !i, j, k+1/2
     &             AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),
     &             AVG_Y_N(V_s(IJMKP, M), V_s(IJKP, M)), K
     &           )
          V_s_B = AVG_Z(                                 !i, j, k-1/2
     &             AVG_Y_N(V_s(IJMKM, M), V_s(IJKM, M)),
     &             AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), KM
     &           )
          W_s_N = AVG_Y(                                 !i, j+1/2, k
     &             AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),
     &             AVG_Z_T(W_s(IJPKM, M), W_s(IJPK, M)), J
     &           )
          W_s_S = AVG_Y(                                 !i, j-1/2, k
     &             AVG_Z_T(W_s(IJMKM, M), W_s(IJMK, M)),
     &             AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), JM
     &           )
          W_s_E = AVG_X(                                 !i+1/2, j, k
     &             AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),
     &             AVG_Z_T(W_s(IPJKM, M), W_s(IPJK, M)), I
     &           )
          W_s_W = AVG_X(                                 !i-1/2, j, k
     &             AVG_Z_T(W_s(IMJKM, M), W_s(IMJK, M)),
     &             AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), IM
     &           )
C
          IF(CYLINDRICAL) THEN
            U_s_C = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)    !i, j, k
            W_s_C = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))    !i, j, k
          ELSE
            U_s_C = ZERO
            W_s_C = ZERO
          ENDIF
C
C         Find components of Mth solids phase continuum strain rate
C         tensor, D_s, at center of THETA cell-(i, j, k)
          D_s(1,1) = ( U_s(IJK,M) - U_s(IMJK,M) ) * oDX(I)
          D_s(1,2) = HALF * ( (U_s_N - U_s_S) * oDY(J) +
     &                         (V_s_E - V_s_W) * oDX(I) )
          D_s(1,3) = HALF * ( (W_s_E - W_s_W) * oDX(I) +
     &                         (U_s_T - U_s_B) * (oX(I)*oDZ(K)) -
     &                          W_s_C * oX(I) )
          D_s(2,1) = D_s(1,2)
          D_s(2,2) = ( V_s(IJK,M) - V_s(IJMK,M) ) * oDY(J)
          D_s(2,3) = HALF * ( (V_s_T - V_s_B) * (oX(I)*oDZ(K)) +
     &                         (W_s_N - W_s_S) * oDY(J) )
          D_s(3,1) = D_s(1,3)
          D_s(3,2) = D_s(2,3)
          D_s(3,3) = ( W_s(IJK,M) - W_s(IJKM,M) ) * (oX(I)*oDZ(K)) +
     &                U_s_C * oX(I)
C
C         Calculate the trace of D_s
          trD_s_C(IJK,M) = D_s(1,1) + D_s(2,2) + D_s(3,3)
C
          IF(EP_g(IJK) .LT. EP_star) THEN
            P_star(IJK) = Neg_H(EP_g(IJK))
C
C  Plastic-flow stress tensor
C
C           Calculate the second invariant of the deviator of D_s
            I2_devD_s = ( (D_s(1,1)-D_s(2,2))**2 
     &                    +(D_s(2,2)-D_s(3,3))**2
     &                    +(D_s(3,3)-D_s(1,1))**2 )/6. 
     &                    + D_s(1,2)**2 + D_s(2,3)**2 + D_s(3,1)**2
C
C-----------------------------------------------------------------------
C            Gray and Stiles (1988)
C            IF(Sin2_Phi .GT. SMALL_NUMBER) THEN
C              qxP_s = SQRT( (4. * Sin2_Phi) * I2_devD_s 
C     &                       + trD_s_C(IJK,M) * trD_s_C(IJK,M))
C              MU_s(IJK, M)     = P_star(IJK) * Sin2_Phi
c     &                         / (qxP_s + SMALL_NUMBER)
c              MU_s(IJK, M)     = MIN(MU_s(IJK, M), MAX_MU_s)
C              LAMBDA_s(IJK, M) = P_star(IJK) * F_Phi
c     &                         / (qxP_s + SMALL_NUMBER)
c              LAMBDA_s(IJK, M) = MIN(LAMBDA_s(IJK, M), MAX_MU_s)
C            ELSE
C              MU_s(IJK, M)     = ZERO
C              LAMBDA_s(IJK, M) = ZERO
C            ENDIF
C-----------------------------------------------------------------------
C           Schaeffer (1987)
C
            qxP_s            = SQRT( (4. * Sin2_Phi) * I2_devD_s)
            MU_s(IJK, M)     = P_star(IJK) * Sin2_Phi
     &                         / (qxP_s + SMALL_NUMBER)
            MU_s(IJK, M)     = MIN(MU_s(IJK, M), MAX_MU_s)

            LAMBDA_s(IJK, M) = ZERO
            ALPHA_s(IJK, M)  = ZERO
          ELSE

C
C  Viscous-flow stress tensor
C
C
C           Calculate trace of the square of D_s
            trD_s2(IJK,M) = 0.D0  !Initialize the totalizer
            DO 20 I1 = 1,3
              DO 10 I2 = 1,3
                trD_s2(IJK,M) = trD_s2(IJK,M) + D_s(I1,I2)*D_s(I1,I2)
   10         CONTINUE
   20       CONTINUE

            IF(.NOT.GRANULAR_ENERGY) THEN  !algebraic granular energy equation
C
C
C             Calculate K_1m, K_2m, K_3m, K_4m
              K_1m = 2.D0 * (ONE + C_e) * RO_s(M) * G_0(IJK, M,M)
              K_3m = HALF * D_p(M) * RO_s(M) * (
     &            ( (SQRT_PI / (3.D0*(3.D0 - C_e))) *
     &            (ONE + 0.4D0*(ONE + C_e)*(3.D0*C_e - ONE)*
     &            EP_s(IJK,M)*G_0(IJK, M,M)) ) +
     &            8.D0*EP_s(IJK,M)*G_0(IJK, M,M)*(ONE + C_e) /
     &            (5.D0*SQRT_PI) )
              K_2m = 4.D0 * D_p(M) * RO_s(M) * (ONE + C_e) *
     &            EP_s(IJK,M) * G_0(IJK, M,M) / (3.D0 * SQRT_PI) -
     &            2.D0/3.D0 * K_3m
              K_4m = 12.D0 * (ONE - C_e*C_e) *
     &            RO_s(M) * G_0(IJK, M,M) / (D_p(M) * SQRT_PI)
              aq   = K_4m*EP_s(IJK,M)
              bq   = K_1m*EP_s(IJK,M)*trD_s_C(IJK,M)
              cq   = -(K_2m*trD_s_C(IJK,M)*trD_s_C(IJK,M) 
     &               + 2.D0*K_3m*trD_s2(IJK,M))
C
C             Boyle-Massoudi Stress term
C
              IF(V_ex .NE. ZERO) THEN        
                K_5m = 0.4 * (ONE + C_e) * G_0(IJK, M,M) * RO_s(M) *
     &            ( (V_ex * D_p(M)) / (ONE - EP_s(IJK,M) * V_ex) )**2
                DEP_soDX  = ( EP_s(IJKE, M) - EP_s(IJK, M) ) * oDX_E(I)
     &           * ( ONE / ( (oDX_E(IM)/oDX_E(I)) + ONE ) ) +
     &           ( EP_s(IJK, M) - EP_s(IJKW, M) ) * oDX_E(IM)
     &           * ( ONE / ( (oDX_E(I)/oDX_E(IM)) + ONE ) ) 
                DEP_soDY  = ( EP_s(IJKN, M) - EP_s(IJK, M) ) * oDY_N(J)
     &           * ( ONE / ( (oDY_N(JM)/oDY_N(J)) + ONE ) ) +
     &           ( EP_s(IJK, M) - EP_s(IJKS, M) ) * oDY_N(JM)
     &           * ( ONE / ( (oDY_N(J)/oDY_N(JM)) + ONE ) ) 
                DEP_soXDZ  = (( EP_s(IJKT, M) - EP_s(IJK, M) ) 
     &           * oX(I)*oDZ_T(K)
     &           * ( ONE / ( (oDZ_T(KM)/oDZ_T(K)) + ONE ) ) +
     &           ( EP_s(IJK, M) - EP_s(IJKB, M) ) * oX(I)*oDZ_T(KM)
     &           * ( ONE / ( (oDZ_T(K)/oDZ_T(KM)) + ONE ) ) ) /
     &           X(I)
                M_s(1,1) = DEP_soDX * DEP_soDX
                M_s(1,2) = DEP_soDX * DEP_soDY
                M_s(1,3) = DEP_soDX * DEP_soXDZ
                M_s(2,1) = DEP_soDX * DEP_soDY
                M_s(2,2) = DEP_soDY * DEP_soDY
                M_s(2,3) = DEP_soDY * DEP_soXDZ
                M_s(3,1) = DEP_soDX * DEP_soXDZ
                M_s(3,2) = DEP_soDY * DEP_soXDZ
                M_s(3,3) = DEP_soXDZ * DEP_soXDZ
                trM_s    = M_s(1,1) + M_s(2,2) + M_s(3,3)
                trDM_s = ZERO
                DO 40 I1 = 1,3
                  DO 30 I2 = 1,3
                    trDM_s = trDM_s + D_s(I1,I2)*M_s(I1,I2)
   30             CONTINUE
   40           CONTINUE
                bq   = bq + EP_s(IJK,M) * K_5m * (trM_s + 2. * trDM_s)
              ELSE
                K_5m = ZERO
              ENDIF
C
C             Calculate EP_sxSQRTHETA and EP_s2xTHETA
              EP_sxSQRTHETA = (-bq + SQRT(bq**2 - 4. * aq * cq ))
     &                       / ( 2. * K_4m )
              EP_s2xTHETA = EP_sxSQRTHETA * EP_sxSQRTHETA

              IF(EP_s(IJK,M) .NE. ZERO)THEN
cstart      kapil&anuj 01/19/98
c               Find pseudo-thermal temperature in the Mth solids phase
                THETA_m(IJK,M) = EP_s2xTHETA/(EP_s(IJK,M)*EP_s(IJK,M))
cend      kapil&anuj 01/19/98
              ELSE
                THETA_m(IJK,M) = ZERO
              ENDIF
C
C             Find pressure in the Mth solids phase
              P_s(IJK,M) = K_1m * EP_s2xTHETA
C
C             bulk viscosity in Mth solids phase 
              LAMBDA_s(IJK,M) = K_2m * EP_sxSQRTHETA 
C
C             shear viscosity in Mth solids phase
              MU_s(IJK,M) = K_3m * EP_sxSQRTHETA
C
C             Boyle-Massoudi stress coefficient
              ALPHA_s(IJK, M) = -K_5m * EP_s2xTHETA

            ELSE	!granular energy transport equation

C             Find pressure in the Mth solids phase
              P_s(IJK,M) = RO_s(M)*EP_S(IJK,M)*(1d0+2d0*(1+C_e)*
     &                     EP_s(IJK,M)*G_0(IJK,M,M))*Theta_m(IJK,M)
C              
              Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(M)*RO_s(M))/96d0

              Mu_b = (256d0*Mu*EP_s(IJK,M)*EP_s(IJK,M)*G_0(IJK,M,M))
     &               /(5d0*Pi)

              Mu_star = Mu/(1.+(2d0*SWITCH*F_gs(IJK,M)*Mu/
     &                  (RO_s(M)*RO_s(M)*EP_s(IJK,M)*EP_s(IJK,M)
     &                  *G_0(IJK,M,M)*Theta_m(IJK,M))))
C
C             shear viscosity in Mth solids phase     
              Mu_s(IJK,M) = 
     &             ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*
     &             G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*
     &             G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*
     &             EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b*Eta))
C
C             bulk viscosity in Mth solids phase 
              LAMBDA_s(IJK,M) = Eta*Mu_b - (2d0*MU_s(IJK,M)/3d0)
              
              Kth=75*RO_s(M)*D_p(M)*DSQRT(Pi*Theta_m(IJK,M))/
     &            (48*Eta*(41d0-33*Eta))

              Kth_star=Kth/(1.+(1.2d0*SWITCH*F_gs(IJK,M)*Kth/
     &                 (RO_s(M)*RO_s(M)*EP_s(IJK,M)*EP_s(IJK,M)
     &                 *G_0(IJK,M,M)*Theta_m(IJK,M))))
C
C             granular conductivity in Mth solids phase
              Kth_s(IJK,M) = Kth_star*(
     &              ( (1d0/G_0(IJK,M,M)) + (12d0/5.)*Eta*EP_s(IJK,M) )
     &            * ( (1d0) + (12d0/5.)*Eta*Eta*(4d0*Eta-3d0)
     &                *EP_s(IJK,M)*G_0(IJK,M,M) )
     &            + (64d0/(25d0*Pi)) * (41d0-33d0*Eta) *
     &               (Eta*EP_s(IJK,M))**2 * G_0(IJK,M,M)
     &        )
C
C     granular 'conductivity' in the Mth solids phase associated
C     with gradient in volume fraction
              Kphi_s(IJK,M) = ZERO
c     &            (Kth_star/(G_0(IJK,M,M)))*(12d0/5.)*Eta*(Eta-1.)*
c     &            (2.*Eta-1.)*(1.+(12d0/5.)*Eta*EP_s(IJK,M)*
c     &            G_0(IJK,M,M))*(EP_s(IJK,M)*
c     &            DG_0DNU(EP_s(IJK,M)) 
c     &            + 2*G_0(IJK,M,M))*Theta_m(IJK,M)
C
C             Boyle-Massoudi stress coefficient
              ALPHA_s(IJK, M) = ZERO
            ENDIF
          ENDIF
C
        END IF
C
  200 CONTINUE
C
      RETURN
      END
