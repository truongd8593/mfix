CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: physical_prop.inc                                      C
C  Purpose: Common block containing physical property data             C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Add PHIP, EP_zero, D_p3, oD_p3, MASS_s                     C
C  Author: W. Sams                                    Date: 22-JUL-93  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Numner: 2                                                  C
C  Purpose: Add K_scale                                                C
C  Author: W. Sams                                    Date: 26-APR-94  C
C  Reviewer:                                                           C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      Number of solids phases
      INTEGER          MMAX
C
C                      Specularity factor
      DOUBLE PRECISION PHIP
C
C                      Scale factor for gas turbulence length scale
      DOUBLE PRECISION K_scale
C
C                      Particle diameters
      DOUBLE PRECISION D_p (DIMENSION_M)
C
C                      Particle densities
      DOUBLE PRECISION RO_s (DIMENSION_M)
C
C                      Specified constant solids viscosity
      DOUBLE PRECISION MU_s0
C
C                      Flag indicates whether the phase becomes close-packed
C                      at ep_star 
      LOGICAL          CLOSE_PACKED (DIMENSION_M)
C
C                      Specified constant gas density
      DOUBLE PRECISION RO_g0
C
C                      Specified constant gas viscosity
      DOUBLE PRECISION MU_g0
C
C                      gas viscosity
      DOUBLE PRECISION MU_g (DIMENSION_3)
C
C                      average molecular weight of gas
      DOUBLE PRECISION MW_AVG
C
C                      Constant constant-pressure specific heat of gas
      DOUBLE PRECISION C_pg0
C
C                      Constant pressure specific heat of gas
      DOUBLE PRECISION C_pg (DIMENSION_3)
C
C                      Constant constant-pressure specific heat of solids
      DOUBLE PRECISION C_ps0
C
C                      Constant pressure specific heat of solids
      DOUBLE PRECISION C_ps (DIMENSION_3, DIMENSION_M)
C
C                      Specified constant gas conductivity
      DOUBLE PRECISION K_g0
C
C                      Conductivity of gas
      DOUBLE PRECISION K_g (DIMENSION_3)
C
C                      Specified constant solids conductivity
      DOUBLE PRECISION K_s0
C
C                      Conductivity of solids
      DOUBLE PRECISION K_s (DIMENSION_3, DIMENSION_M)
C
C		       Granular Temperature Conductivity (associated
C                      with temperature gradient)
      DOUBLE PRECISION Kth_s (DIMENSION_3, DIMENSION_M)
C
C		       Granular Temperature Conductivity (associated
C                      with volume fraction gradient)
      DOUBLE PRECISION Kphi_s (DIMENSION_3, DIMENSION_M)
C
C                      Specified constant gas diffusivity
      DOUBLE PRECISION DIF_g0
C
C                      Diffusivity of gas species N
      DOUBLE PRECISION DIF_g (DIMENSION_3, DIMENSION_N_g)
C
C                      Specified constant solids diffusivity
      DOUBLE PRECISION DIF_s0
C
C                      Diffusivity of solids species N
      DOUBLE PRECISION DIF_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s)
C
C                      Total number of gas or solids species
      INTEGER          NMAX(0:DIMENSION_M)
C
C                      Molecular weight of gas species
      DOUBLE PRECISION MW_g (DIMENSION_N_g)
C
C                      Molecular weight of solids species
      DOUBLE PRECISION MW_s (DIMENSION_M, DIMENSION_N_s)
C
C                      Molecular weight of gas mixture
      DOUBLE PRECISION MW_MIX_g (DIMENSION_3)
C
      COMMON / PHYSICAL_PROP_DP /
     &       PHIP, K_scale,
     &       D_p, RO_s, MU_s0,
     &       RO_g0, MU_g0, MU_g, MW_AVG,
     &       C_pg, C_ps, C_pg0, C_ps0,
     &       K_g, K_s, Kth_s, Kphi_s, MW_g, MW_s, MW_MIX_g, K_g0, K_s0,
     &       DIF_g, DIF_s, DIF_g0, DIF_s0
      COMMON / PHYSICAL_PROP_I /
     &       MMAX, NMAX

      COMMON / PHYSICAL_PROP_L /
     &       CLOSE_PACKED
