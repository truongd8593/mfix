      MODULE visc_g

! trace of D_g at i, j, k
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  trD_g

! Turbulent viscosity of fluid phase: sum of molecular
! and eddy viscosities
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_gt

! Turbulent viscosity of fluid phase: sum of molecular
! and eddy viscosities x the void fraction of fluid phase
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  EPMU_gt


! Granular second coefficient of viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  LAMBDA_gt

! Characteristic length for turbulence
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  L_scale

! gas viscosity needs to be updated during iterations?
      LOGICAL :: Recalc_visc_g

! diffusive component of conv-dif 
! stores diffusive term in matrix for x, y, z momentum cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DF_gu
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DF_gv
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DF_gw

      END MODULE visc_g
