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

      LOGICAL FIRST_CALL

! DES - Species Equation
      LOGICAL DES_SPECIES_EQ(DIM_M)
      LOGICAL ANY_DES_SPECIES_EQ

! total number of discrete solids species for each phase
      INTEGER DES_NMAX(DIM_M)
! maximum number of species over all solids phases
      INTEGER MAX_DES_NMAX

      DOUBLE PRECISION OUTPUT_DATA_TIME

! species names
      CHARACTER(len=18) DES_SPECIES_NAME(DIM_N_ALL)

! molecular weight of discrete solids species
      DOUBLE PRECISION DES_MW_s (DIM_M, DIM_N_s)

! discrete solids species mass fractions (PARTICLES, 0:MAX_DES_NMAX))
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_X_s

! rate of production of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_R_sp
! rate of consumption of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_R_sc

! combined rate of production and consumption of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_SUM_R_s

! thermochemical data for discrese solids phases:
      DOUBLE PRECISION DES_Thigh_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Tlow_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Tcom_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Ahigh_s(7, DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_Alow_s(7, DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_HfrefoR_s(DIM_M, DIM_N_s)
      DOUBLE PRECISION DES_IC_PSrefoR(DIM_M, DIM_N_s)

! This indicates how a particle will respond to a chemical reaction.
!  1) [SHRINKING_CORE] : The outer radius of the particle reamins
!        constant while an unreacted core shrinks during the reaction.
!        In this model, a secondary radius (CORE_RAD) is tracked for
!        each particle to represent the radius of the core.
!  2) SHRINKING_PARTICLE : As components of the particle are consumed,
!        the size of the particle reduces to maintain a constant
!        density.
!  3) PROGRESSIVE_CONVERSION : The size of the particle remains constant
!        but the density of the particle is changed to reflect the 
!        consumption of the components of the particle.
      CHARACTER*64 REACTION_MODEL

! If using the shrinking core model, this term represents the radius
! of the unreacted core.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CORE_RAD
! Constnat density of the core for the shrinking core model. The
! Ro_Sol variable is adjusted to reflect a 'bluk density' of the 
! particle.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CORE_Rho


! Previous time step's rate of change. Used for Adams-Bashforth
! time integration scheme.
 ! 1) particle mass
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dMdt_OLD
 ! 2) particle species mass
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dXdt_OLD
 ! 3) radius of unreacted core
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dRdt_OLD




      END MODULE DES_RXNS
