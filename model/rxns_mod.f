      MODULE rxns
 
 
      Use param
      Use param1
 
 
!
!                      Rate of production of gas species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  R_gp 
!
!                      Rate of production of solids species
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  R_sp 
!
!                      Rate of consumption of gas species/X_g
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  RoX_gc 
!
!                      Rate of consumption of solids species/X_s
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  RoX_sc 
!
!                      Net production of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SUM_R_g 
!
!                      Net production of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  SUM_R_s 
!
!                      Rate of mass transfer from phase M to Phase L
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  R_phase 
!                        
!
!                      Molecular weight of all species: MW_g and MW_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MW_all 
!
!                      total number of species
      INTEGER          N_all
!
!                      Species name to index
      INTEGER, DIMENSION(:, :), ALLOCATABLE ::           SPECIES_ID2N
!
!                      Species index to name
      INTEGER, DIMENSION(:), ALLOCATABLE ::           SPECIES_N2IDg
      INTEGER, DIMENSION(:, :), ALLOCATABLE ::           SPECIES_N2IDs
!
!     Fraction of the mass source from RXn# that results in mass transfer
!     to "To phase #" from "From phase #."  Determined in check_data_09 and
!     used in rrates0.
!          R_temp(Rxn#, To phase #, From phase #)
!     e.g. R_temp(1,0,1) -  mass generation of gas phase from solids-1,
!          R_temp(1,0,2) -  mass generation of gas phase from solids-2,
!          R_temp(1,1,0) -  mass generation of solid-1 from gas = -R_temp(1,0,1)
!          R_temp(1,1,2) -  mass generation of solid-1 from solids-2.
!     Note, for example, that if gas is generated from solids-1 then
!     R_temp(1,0,1) > 0. The R-phase matrix is skew-symmetric and diagonal
!     elements are not needed. Only one of the two skew-symmetric elements
!     -- e.g., R_temp(1,0,1) or R_temp(1,1,0) -- needs to be specified.
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN, 0:DIM_M, 0:DIM_M) :: R_temp
!
!                      Species names
      CHARACTER*10     SPECIES_NAME(DIM_N_ALL)
!
!                      Reaction names
      CHARACTER*10, DIMENSION(DIMENSION_RXN) :: RXN_NAME
!
!                      whether rxn schemes and rate expressions specified
      LOGICAL, DIMENSION(DIMENSION_RXN) :: GOT_RXN, GOT_RATE 
!
!                      stoichiometry
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN, DIM_N_ALL) ::  STOICH
!
!                      stoichiometry x molecular weight
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN, DIM_N_ALL) ::  STOICHxMW
!
!                      Enthalpy change due to reaction.  The change is assigned
!                      to the phase identified by Rate_m4T.
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN) ::  Delta_H
!
!                      Which tempearure to use in rate expression?  Also,
!                      the phase to which the enthalpy change due to rxn is
!                      assigned.
!                      m = 0 => fluid phase, m = 1 => solids phase 1, . . .
      INTEGER, DIMENSION(DIMENSION_RXN) :: Rate_m4T
!
!                      Factors in rate expression:
!                      1. preexponential factor (A),
!                      2. Temperature exponent (n),
!                      3. Activation temperature (E/R)
!                      rate expression = A * T_m**n * exp(-(E/R) / T_m) * . . .
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN, 3)::  Rate_fac
!
!                      Exponents for concentartion dependence in the rate expression
      DOUBLE PRECISION, DIMENSION(DIMENSION_RXN, DIM_N_ALL) ::  Rate_exp
!
!                      Actual number of Reactions
      INTEGER          NO_OF_RXNS
!
 
 
!HPF$ align R_gp(:, *) with TT(:)
!HPF$ align R_sp(:, *, *) with TT(:)
!HPF$ align RoX_gc(:, *) with TT(:)
!HPF$ align RoX_sc(:, *, *) with TT(:)
!HPF$ align SUM_R_g(:) with TT(:)
!HPF$ align SUM_R_s(:, *) with TT(:)
!HPF$ align R_phase(:, *) with TT(:)

      END MODULE rxns
