      MODULE visc_g


      Use param
      Use param1


!
!                      trace of D_g at i, j, k
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  trD_g
!
!                      Turbulent viscosity of fluid phase: sum of molecular
!                      and eddy viscosities
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_gt
!
!                      Granular second coefficient of viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  LAMBDA_gt
!
!                      Characteristic length for turbulence
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  L_scale

!
!                      gas viscosity needs to be updated during iterations?
      LOGICAL          Recalc_visc_g



!!!HPF$ align trD_g(:) with TT(:)
!!!HPF$ align MU_gt(:) with TT(:)
!!!HPF$ align LAMBDA_gt(:) with TT(:)
!!!HPF$ align L_scale(:) with TT(:)

      END MODULE visc_g
