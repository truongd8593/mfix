!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: physprop                                               C
!  Purpose: Common block containing physical property data             C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE physprop
 
!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use param, only: dim_m, dim_n, dim_n_g, dim_n_s
!-----------------------------------------------


! Scale factor for gas turbulence length scale
      DOUBLE PRECISION :: K_scale
 
! Number of solids phases
      INTEGER :: MMAX

! Real number of solids phases for GHD theory
      INTEGER :: SMAX

! Particle diameters
      DOUBLE PRECISION :: D_p0(DIM_M)

! Constant solids phase densities.
      DOUBLE PRECISION :: RO_s0(DIM_M)

! Baseline/Unreacted solids phase phase density. This value is only used
! for variable solids density simulations in EOSS.
      DOUBLE PRECISION :: BASE_ROs(DIM_M)

! Constant solids phase species mass fractions. These values delinate
! the baseline/initial solids phase composition for ariable density
      DOUBLE PRECISION :: X_S0(DIM_M, DIM_N_s)

! Density of solid species (constant)
      DOUBLE PRECISION :: RO_Xs0(DIM_M, DIM_N_s)

! The index of an inert solids phase species. This is needed for 
! calculating the variable solids phase density.
      INTEGER :: INERT_SPECIES(DIM_M)

! Particle shape factor
      DOUBLE PRECISION :: PSI_s(DIM_M)

! Specified constant solids viscosity
      DOUBLE PRECISION MU_s0

! Flag indicates whether the phase becomes close-packed at ep_star
      LOGICAL :: CLOSE_PACKED (DIM_M)

! Specified constant gas density
      DOUBLE PRECISION RO_g0

! Specified constant gas viscosity
      DOUBLE PRECISION MU_g0

! Virtual (added) mass coefficient Cv
      DOUBLE PRECISION Cv

! Gas viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_g 

! Average molecular weight of gas
      DOUBLE PRECISION MW_AVG

! Constant constant-pressure specific heat of gas
      DOUBLE PRECISION C_pg0

! Reference temperature for enthalpy calculations (K)
      DOUBLE PRECISION, PARAMETER :: T_ref = 298

! Constant pressure specific heat of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  C_pg 

! Constant constant-pressure specific heat of solids
      DOUBLE PRECISION C_ps0

! Constant pressure specific heat of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  C_ps 

! Specified constant gas conductivity
      DOUBLE PRECISION :: K_g0

! Conductivity of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_g 

! Specified constant solids conductivity
      DOUBLE PRECISION K_s0

! Conductivity of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  K_s 

! Granular Temperature Conductivity (associated with temperature grad)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kth_s 

! Granular Temperature Conductivity (associated with volume fraction grad)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kphi_s 

! Specified constant gas diffusivity
      DOUBLE PRECISION DIF_g0

! Diffusivity of gas species N
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  DIF_g 

! Specified constant solids diffusivity
      DOUBLE PRECISION DIF_s0

! Diffusivity of solids species N
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  DIF_s 

! Total number of gas or solids species
      INTEGER :: NMAX(0:DIM_M) ! Runtime (all phases)
      INTEGER :: NMAX_g        ! Number of gas phase species
      INTEGER :: NMAX_s(DIM_M) ! Number of solids phase species

! Molecular weight of gas species
      DOUBLE PRECISION :: MW_g (DIM_N_g)

! Molecular weight of solids species
      DOUBLE PRECISION :: MW_s (DIM_M, DIM_N_s)

! Molecular weight of gas mixture
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MW_MIX_g

! Logical for reading thermochemical database.
      LOGICAL :: DATABASE_READ = .TRUE.

! Polynomical coefficients for calculating specific heat.
      DOUBLE PRECISION Alow (7,0:DIM_M, DIM_N) ! Tlow --> Tcom
      DOUBLE PRECISION Ahigh(7,0:DIM_M, DIM_N) ! Tcom --> Thigh

! Range where the polynomials are valid.
      DOUBLE PRECISION Thigh(0:DIM_M, DIM_N) ! Upper bound
      DOUBLE PRECISION Tlow (0:DIM_M, DIM_N) ! Lower bound
      DOUBLE PRECISION Tcom (0:DIM_M, DIM_N) ! Switch from low to high

! Heat of formation at Tref divided by the gas constant.
      DOUBLE PRECISION HfrefoR(0:DIM_M, DIM_N)

! Reference values.
      DOUBLE PRECISION ICpoR_l(0:DIM_M, DIM_N)
      DOUBLE PRECISION ICpoR_h(0:DIM_M, DIM_N)

      END MODULE physprop
