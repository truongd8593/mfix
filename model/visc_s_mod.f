      MODULE visc_s
 
 
      Use param
      Use param1
 
 
!
!                      trace of D_s at i, j, k
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s
!
!                      Dissipation term for granular temperature
      DOUBLE PRECISION GAMMA_m
!
!                      Granular first coefficient of (shear) viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s 
!
!                      Granular second coefficient of viscosity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s 
!
!                      Boyle-Massoudi stress tensor (i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ALPHA_s 
 
!                      Granular first coefficient of (shear) viscosity
!                      (i,j,k) - Collisional contribution
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s_c 
!
!                      Granular second coefficient of viscosity
!                      Collisional contribution
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s_c 
 
 
 
!!!HPF$ align trD_s(:, *) with TT(:)
!!!HPF$ align MU_s(:, *) with TT(:)
!!!HPF$ align LAMBDA_s(:, *) with TT(:)
!!!HPF$ align ALPHA_s(:, *) with TT(:)
!!!HPF$ align MU_s_c(:, *) with TT(:)
!!!HPF$ align LAMBDA_s_c(:, *) with TT(:)

      END MODULE visc_s
