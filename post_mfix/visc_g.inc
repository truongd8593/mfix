C
C                      trace of D_g at i, j, k
      DOUBLE PRECISION trD_g(DIMENSION_3)
C
C                      Turbulent viscosity of fluid phase: sum of molecular
C                      and eddy viscosities
      DOUBLE PRECISION MU_gt (DIMENSION_3)
C
C                      Granular second coefficient of viscosity
      DOUBLE PRECISION LAMBDA_gt (DIMENSION_3)
C
C                      Characteristic length for turbulence
      DOUBLE PRECISION L_scale (DIMENSION_3)

      COMMON /visc_g_DP/
     &  trD_g, MU_gt, LAMBDA_gt, L_scale

