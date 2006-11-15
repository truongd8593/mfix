  MODULE kintheory
 
 
      Use param
      Use param1
 
!
!     coefficient terms needed for stress (in addition to: 
!     mu_s, lambda_s which are defined in visc_s_mod, allocated in
!          allocate_arrays, and initialized in set_constprop and
!     P_s which is defined in fldvar_mod, allocated in 
!          allocate_arrays, and initialized in init_fvars)
!
!     stress term with gradient in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sM_ip
!     stress term with gradient in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sL_ip
!     stress term with trace in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sM_ip
!     stress term with trace in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sL_ip
!
!
!
!     coefficient terms needed for momentum source (in addition to:
!     F_SS which is defined in drag_mod, allocated in allocate_arrays
!          and initialized in set_constprop)
!
!     momentum source term with gradient in number density
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Fnu_s_ip
!     momentum source term with gradient in mixture temperature or
!     with the gradient in temperature of species M
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sM_ip
!     momentum source term with gradient in temperature of species L
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sL_ip
!
!
!
!     coefficient terms needed for heat flux (in addition to: 
!     kth_s which is defined in set_constprop, allocated 
!          allocate_arrays and initialized in init_fvars)
!
!     heat flux term with gradient in granular temperature of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kth_sL_ip
!     heat flux term with gradient in number density of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sM_ip
!     heat flux term with gradient in number density of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sL_ip
!     heat flux term with velocity difference
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kvel_s_ip
!
!
!
!     coefficient terms needed for energy dissipation 
!
!     energy dissipation with difference in species granular
!     temperature: transfer between solid solid phases
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ED_ss_ip
!
!     energy dissipation term
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDT_s_ip
!     energy dissipation with divergence of velocity of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sM_ip
!     energy dissipation with divergence of velocity of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sL_ip


      END MODULE kintheory
