MODULE fun_avg

  USE param1
  USE geometry

CONTAINS

  !   average at i+1/2 location
  ! Arithmetic averages
  !     i+1/2, j, k. Use IJKE for scalar
  !          F(IJK) F(IPJK)  I
  DOUBLE PRECISION FUNCTION AVG_X(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_X = FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp
  END FUNCTION AVG_X

  !   average (of U) at i+1 location
  !
  !     i+1, j, k
  !          F(IJK) F(IPJK)  IP
  DOUBLE PRECISION FUNCTION AVG_X_E(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_X_E = FX(xL) * XXXm + FX_bar(xL) * XXXp
  END FUNCTION AVG_X_E

  !   average at j+1/2 location
  !
  !     i, j+1/2, k. Use IJKN for scalar
  !           F(IJK) F(IJPK) J
  DOUBLE PRECISION FUNCTION AVG_Y(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Y = FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp
  END FUNCTION AVG_Y

  !   average (of V) at j+1 location
  !
  !     i, j+1, k
  !           F(IJK) F(IJPK)
  DOUBLE PRECISION FUNCTION AVG_Y_N(XXXm,  XXXp)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm

    AVG_Y_N = HALF *( XXXm + XXXp )
  END FUNCTION AVG_Y_N

  !   average at k+1/2 location
  !
  !     i, j, k+1/2. Use IJKT for scalar
  !           F(IJK) F(IJKP) K
  DOUBLE PRECISION FUNCTION AVG_Z(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Z = FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp
  END FUNCTION AVG_Z

  !   average (of W) at k+1 location
  !
  !     i, j, k+1
  !           F(IJK) F(IJKP)
  DOUBLE PRECISION FUNCTION AVG_Z_T(XXXm,  XXXp)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm

    AVG_Z_T = HALF *( XXXm + XXXp )
  END FUNCTION AVG_Z_T
  !
  !   Harmonic average at i+1/2 location
  !
  !     i+1/2, j, k. Use IJKE for scalar
  !             F(IJK) F(IPJK)  I
  DOUBLE PRECISION FUNCTION AVG_X_h(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_X_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
  END FUNCTION AVG_X_h
  !
  !   Harmonic averageat j+1/2 location
  !
  !     i, j+1/2, k. Use IJKN for scalar
  !             F(IJK) F(IJPK) J
  DOUBLE PRECISION FUNCTION AVG_Y_h(XXXm,  XXXp,   xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Y_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
  END FUNCTION AVG_Y_h
  !
  !   Harmonic average at k+1/2 location
  !
  !     i, j, k+1/2. Use IJKT for scalar
  !             F(IJK) F(IJKP) K
  DOUBLE PRECISION FUNCTION AVG_Z_h(XXXm,  XXXp,   xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Z_h = XXXm * XXXp / &
         MAX(SMALL_NUMBER, ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )
  END FUNCTION AVG_Z_h
  !
  ! Harmonic averages for possibly negative scalars (sof, Aug 22 2006)
  !   Harmonic average at i+1/2 location
  !
  !     i+1/2, j, k. Use IJKE for scalar
  !             F(IJK) F(IPJK)  I
  DOUBLE PRECISION FUNCTION AVG_X_S(XXXm,  XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_X_S = XXXm * XXXp / &
         (1D-30 + ( FX_E(xL) * XXXm + FX_E_bar(xL) * XXXp) )
  END FUNCTION AVG_X_S

  !   Harmonic averageat j+1/2 location
  !
  !     i, j+1/2, k. Use IJKN for scalar
  !             F(IJK) F(IJPK) J
  DOUBLE PRECISION FUNCTION AVG_Y_S(XXXm, XXXp, xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Y_S = XXXm * XXXp / &
         (1D-30 + ( FY_N(xL) * XXXm + FY_N_bar(xL) * XXXp) )
  END FUNCTION AVG_Y_S

  !   Harmonic average at k+1/2 location
  !
  !     i, j, k+1/2. Use IJKT for scalar
  !             F(IJK) F(IJKP) K
  DOUBLE PRECISION FUNCTION AVG_Z_S(XXXm,  XXXp,   xL)
    IMPLICIT NONE
    DOUBLE PRECISION XXXp, XXXm
    INTEGER          xL

    AVG_Z_S = XXXm * XXXp / &
         (1D-30 + ( FZ_T(xL) * XXXm + FZ_T_bar(xL) * XXXp) )
  END FUNCTION AVG_Z_S

END MODULE fun_avg
