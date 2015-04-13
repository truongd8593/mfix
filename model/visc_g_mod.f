      MODULE visc_g

! trace of D_g at i, j, k
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  trD_g

! Turbulent viscosity of fluid phase: sum of molecular
! and eddy viscosities
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_gt

! Granular second coefficient of viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  LAMBDA_gt

! Characteristic length for turbulence
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  L_scale

! gas viscosity needs to be updated during iterations?
      LOGICAL :: Recalc_visc_g

      END MODULE visc_g
