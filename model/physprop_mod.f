!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: physical_prop.inc                                      C
!  Purpose: Common block containing physical property data             C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add EP_zero, D_p3, oD_p3, MASS_s                     C
!  Author: W. Sams                                    Date: 22-JUL-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Numner: 2                                                  C
!  Purpose: Add K_scale                                                C
!  Author: W. Sams                                    Date: 26-APR-94  C
!  Reviewer:                                                           C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
 
 
      MODULE physprop
 
 
      Use param
      Use param1
 
 
!
!                      Number of solids phases
      INTEGER          MMAX
!
!                      Scale factor for gas turbulence length scale
      DOUBLE PRECISION K_scale
!
!                      Particle diameters
      DOUBLE PRECISION D_p0 (DIM_M)
!
!                      index to rearrange particles from coarsest to finest
!                      for use in function CALC_ep_star(IJK,IER)
      INTEGER          M_MAX (DIM_M)
!
!                      Particle densities
      DOUBLE PRECISION RO_s  (DIM_M)
!
!                      Specified constant solids viscosity
      DOUBLE PRECISION MU_s0
!
!                      Flag indicates whether the phase becomes close-packed
!                      at ep_star
      LOGICAL          CLOSE_PACKED (DIM_M)
!
!                      Specified constant gas density
      DOUBLE PRECISION RO_g0
!
!                      Specified constant gas viscosity
      DOUBLE PRECISION MU_g0
!
!                      gas viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_g 
!
!                      average molecular weight of gas
      DOUBLE PRECISION MW_AVG
!
!                      Constant constant-pressure specific heat of gas
      DOUBLE PRECISION C_pg0
!
!                      Reference temperature for enthalpy calculations (K)
      DOUBLE PRECISION, PARAMETER :: T_ref = 298
!
!                      Constant pressure specific heat of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  C_pg 
!
!                      Constant constant-pressure specific heat of solids
      DOUBLE PRECISION C_ps0
!
!                      Constant pressure specific heat of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  C_ps 
!
!                      Specified constant gas conductivity
      DOUBLE PRECISION K_g0
!
!                      Conductivity of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_g 
!
!                      Specified constant solids conductivity
      DOUBLE PRECISION K_s0
!
!                      Conductivity of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  K_s 
!
!		       Granular Temperature Conductivity (associated
!                      with temperature gradient)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kth_s 
!
!		       Granular Temperature Conductivity (associated
!                      with volume fraction gradient)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kphi_s 
!
!                      Specified constant gas diffusivity
      DOUBLE PRECISION DIF_g0
!
!                      Diffusivity of gas species N
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  DIF_g 
!
!                      Specified constant solids diffusivity
      DOUBLE PRECISION DIF_s0
!
!                      Diffusivity of solids species N
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  DIF_s 
!
!                      Total number of gas or solids species
      INTEGER          NMAX(0:DIM_M)
!
!                      Molecular weight of gas species
      DOUBLE PRECISION MW_g (DIM_N_g)
!
!                      Molecular weight of solids species
      DOUBLE PRECISION MW_s (DIM_M, DIM_N_s)
!
!                      Molecular weight of gas mixture
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MW_MIX_g 
!
!                      Coefficients for high and low temperature ranges, heat of formation at T_ref K
!                      and the temperature ranges for calculating thermochemical properties
      LOGICAL :: DATABASE_READ = .FALSE.
      DOUBLE PRECISION Ahigh_g(7, DIM_N_g), Alow_g(7, DIM_N_g), HfrefoR_g(DIM_N_g)
      DOUBLE PRECISION Thigh_g(DIM_N_g), Tlow_g(DIM_N_g), Tcom_g(DIM_N_g), IC_PGrefoR(DIM_N_g)
      DOUBLE PRECISION Ahigh_s(7, DIM_M, DIM_N_s), Alow_s(7, DIM_M, DIM_N_s), HfrefoR_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION Thigh_s(DIM_M, DIM_N_s), Tlow_s(DIM_M, DIM_N_s), Tcom_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION IC_PsrefoR(DIM_M, DIM_N_s)
 
 
!!!HPF$ align MU_g(:) with TT(:)
!!!HPF$ align C_pg(:) with TT(:)
!!!HPF$ align C_ps(:, *) with TT(:)
!!!HPF$ align K_g(:) with TT(:)
!!!HPF$ align K_s(:, *) with TT(:)
!!!HPF$ align Kth_s(:, *) with TT(:)
!!!HPF$ align Kphi_s(:, *) with TT(:)
!!!HPF$ align DIF_g(:, *) with TT(:)
!!!HPF$ align DIF_s(:, *, *) with TT(:)
!!!HPF$ align MW_MIX_g(:) with TT(:)

      END MODULE physprop
