      MODULE visc_s
      
      
      Use param
      Use param1
      
!     
!     trace of D_s at i, j, k
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s
!     
!     Dissipation term for granular temperature
      DOUBLE PRECISION GAMMA_m
!     
!     Granular first coefficient of (shear) viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s 
!     
!     Granular second coefficient of viscosity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s 
!     
!     Boyle-Massoudi stress tensor (i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ALPHA_s 
      
!     Granular first coefficient of (shear) viscosity
!     (i,j,k) - Collisional contribution
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s_c 
!     
!     Granular second coefficient of viscosity
!     Collisional contribution
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s_c 

!     Granular first coefficient of (shear) viscosity(i,j,k) - Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_v 
!     
!     Granular first coefficient of (shear) viscosity(i,j,k) - Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_p 
!     
!     Granular first coefficient of (shear) viscosity(i,j,k) - Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_f 
!     
!     Granular second coefficient of (shear) viscosity(i,j,k) - Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_v 
!     
!     Granular second coefficient of (shear) viscosity(i,j,k) - Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_p 
!     
!     Granular second coefficient of (shear) viscosity(i,j,k) - Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_f 
!     
!     Boyle-Massoudi stress tensor (i,j,k) - Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Alpha_s_v 
!     
!     Boyle-Massoudi stress tensor (i,j,k)- Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Alpha_s_p 
!     
!     Boyle-Massoudi stress tensor (i,j,k)- Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Alpha_s_f 
!     
!              Bulk viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mu_b_v

      
!     Packed bed (close packed) void fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_star_array

!     Second invariant of the deviator of D_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: I2_devD_s
!     
!     Trace of M_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trM_s
!     
!     Trace of (D_s)(M_s)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trDM_s
!     
!     Relative velocity array
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VREL_array

      
      
      
!!!   HPF$ align trD_s(:, *) with TT(:)
!!!   HPF$ align MU_s(:, *) with TT(:)
!!!   HPF$ align LAMBDA_s(:, *) with TT(:)
!!!   HPF$ align ALPHA_s(:, *) with TT(:)
!!!   HPF$ align MU_s_c(:, *) with TT(:)
!!!   HPF$ align LAMBDA_s_c(:, *) with TT(:)

      END MODULE visc_s
