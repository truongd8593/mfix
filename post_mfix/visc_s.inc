C
C                      trace of D_s at i, j, k
      DOUBLE PRECISION trD_s(DIMENSION_3, DIMENSION_M)
C
C                      Dissipation term for granular temperature
      DOUBLE PRECISION GAMMA_m
C
C                      Granular first coefficient of (shear) viscosity(i,j,k)
      DOUBLE PRECISION MU_s (DIMENSION_3, DIMENSION_M)
C
C                      Granular second coefficient of viscosity
      DOUBLE PRECISION LAMBDA_s (DIMENSION_3, DIMENSION_M)
C
C                      Boyle-Massoudi stress tensor (i,j,k)
      DOUBLE PRECISION ALPHA_s (DIMENSION_3, DIMENSION_M)

      COMMON /visc_s_DP/
     &  trD_s, GAMMA_m,
     &  MU_s, LAMBDA_s, 
     &  ALPHA_s
