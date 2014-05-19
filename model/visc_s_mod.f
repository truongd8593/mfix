      MODULE visc_s
      
          
! Granular first coefficient of (shear) viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s 

! Granular first coefficient of (shear) viscosity(i,j,k) 
! Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_v 
           
! Granular first coefficient of (shear) viscosity(i,j,k) 
! Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_p 
     
! Granular first coefficient of (shear) viscosity(i,j,k) 
! Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_s_f 

! Bulk viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mu_b_v

! Granular first coefficient of (shear) viscosity(i,j,k)
! Collisional contribution of viscous component
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  MU_s_c

! Granular second coefficient of viscosity(i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s 
      
! Granular second coefficient of (shear) viscosity(i,j,k)
! Viscous contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_v 
     
! Granular second coefficient of (shear) viscosity(i,j,k) 
! Plastic contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_p 
     
! Granular second coefficient of (shear) viscosity(i,j,k) 
! Frictional contribution
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Lambda_s_f 
     
! Granular second coefficient of viscosity
! Collisional contribution of viscous component
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  LAMBDA_s_c 

! Packed bed (close packed) void fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_star_array

! Start of Blending
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_g_blend_start

! End of Blending
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EP_g_blend_end

! trace of D_s (rate of strain tensor) at i, j, k
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s

! Second invariant of the deviator of D_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: I2_devD_s

! Boyle-Massoudi stress tensor coefficient (i,j,k)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ALPHA_s 

! For Boyle-Massoudi stress tensor
! Trace of M_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trM_s
     
! Trace of (D_s)(M_s)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trDM_s
     
! Relative velocity array
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VREL_array

      
            END MODULE visc_s
