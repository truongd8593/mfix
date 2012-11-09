!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_SPECIES                                            !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM species transfer.             !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_RXNS

      USE param
      USE rxn_com


! Run time options:
!---------------------------------------------------------------------//
! DES - Species Equation
      LOGICAL DES_SPECIES_EQ(DIM_M)
      LOGICAL ANY_DES_SPECIES_EQ

! This indicates how a particle will respond to a chemical reaction.
!  1) [VARIABLE_DENSITY] : The size of the particle remains constant
!        but the density of the particle is changed to reflect the 
!        consumption of the components of the particle.
!  2) SHRINKING_PARTICLE : As components of the particle are consumed,
!        the size of the particle reduces to maintain a constant
!        density.
      CHARACTER*64 REACTION_MODEL
      LOGICAL RM_VARIABLE_DENSITY
      LOGICAL RM_SHRINKING_PARTICLE

! Logical flags indicating that the database was read for a particular
! species. DIM_N_g is used (as opposed to DIM_N_s) to ensure the array
! sizes for DEM/TFM are the same. Note that both are set to 100 in
! param_mod.f.
      LOGICAL DES_rDatabase(DIM_M, DIM_N_g)

! Data Storate:
!---------------------------------------------------------------------//
! maximum number of species over all solids phases
      INTEGER MAX_DES_NMAX

! molecular weight of discrete solids species
      DOUBLE PRECISION DES_MW_s(DIM_M, DIM_N_s)

! total number of discrete solids species for each phase
      INTEGER DES_NMAX_s(DIM_M)

! Solids phase species names (database) and aliases
      CHARACTER(len=18) DES_SPECIES_s(DIM_M, DIM_N_s) ! database name
      CHARACTER(len=32)  DES_SPECIES_ALIAS_s(DIM_M, DIM_N_s) ! alias

! discrete solids species mass fractions (PARTICLES, 0:MAX_DES_NMAX))
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_X_s

! Rate of production of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_R_sp
! rate of consumption of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_R_sc

! combined rate of production and consumption of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_SUM_R_s


! Numerical integration:
!---------------------------------------------------------------------//
! Previous time step's rate of change. Used for Adams-Bashforth
! time integration scheme.
 ! 1) particle mass
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dMdt_OLD
 ! 2) particle species mass
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dXdt_OLD
 ! 3) radius of unreacted core
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dRdt_OLD


! Interphase transfer variables.
!---------------------------------------------------------------------//
! Amount produced of gas phase species
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_gp
! Amount consumed of gas phase species
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_gc
! Net production (+) or consumption (-)  of gas phase species
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_SUM_R_g
! Amount of interphase mass transfer.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_PHASE
! Amount of gas phase enthalpy change due to phase change/chemical reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_HOR_g

! Phase change/reactive chemistry:
!---------------------------------------------------------------------//
! Actual number of Reactions
      INTEGER NO_OF_DES_RXNS

! Array linking all of the reaction data.
      TYPE(REACTION_BLOCK), DIMENSION(:), TARGET, ALLOCATABLE :: DES_Reaction


! Thermochemical data for discrese solids phases:
!---------------------------------------------------------------------//
      DOUBLE PRECISION DES_Thigh_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Tlow_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Tcom_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Ahigh_s(7, DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Alow_s(7, DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_HfrefoR_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_IC_PSrefoR(DIM_M, DIM_N_s)


      END MODULE DES_RXNS
