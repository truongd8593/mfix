      MODULE energy
 
 
      Use param
      Use param1
 
 
!
!                      Gas-phase heat of reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  HOR_g 
!
!                      Solids-phase heat of reaction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  HOR_s 
!
!                      Gas-solids heat transfer coefficient
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_gs 
!
!                      Gas-phase radiation coefficient
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  GAMA_Rg 
!
!                      Solids-phase radiation coefficient
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_Rs 
!
!                      Gas-phase radiation temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  T_Rg 
!
!                      Solids-phase radiation temperature
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  T_Rs 
!
 
 
!!!HPF$ align HOR_g(:) with TT(:)
!!!HPF$ align HOR_s(:, *) with TT(:)
!!!HPF$ align GAMA_gs(:, *) with TT(:)
!!!HPF$ align GAMA_Rg(:) with TT(:)
!!!HPF$ align GAMA_Rs(:, *) with TT(:)
!!!HPF$ align T_Rg(:) with TT(:)
!!!HPF$ align T_Rs(:, *) with TT(:)

      END MODULE energy
