C
C                      Rate of production of gas species
      DOUBLE PRECISION R_gp (DIMENSION_3, DIMENSION_N_g)
C
C                      Rate of production of solids species
      DOUBLE PRECISION R_sp (DIMENSION_3, DIMENSION_M, DIMENSION_N_s)
C
C                      Rate of consumption of gas species/X_g
      DOUBLE PRECISION RoX_gc (DIMENSION_3, DIMENSION_N_g)
C
C                      Rate of consumption of solids species/X_s
      DOUBLE PRECISION RoX_sc (DIMENSION_3, DIMENSION_M, DIMENSION_N_s)
C
C                      Net production of gas
      DOUBLE PRECISION SUM_R_g (DIMENSION_3)
C
C                      Net production of solids
      DOUBLE PRECISION SUM_R_s (DIMENSION_3, DIMENSION_M)
C
C                Mass transfer between phases.
C     Determine the mass g/(cm^3.s) transferred from one phase to the other.
C     e.g. R_temp(0,1) -  mass generation of gas phase from solids-1,
C          R_temp(0,2) -  mass generation of gas phase from solids-2,
C          R_temp(1,0) -  mass generation of solid-1 from gas = -R_temp(0,1)
C          R_temp(1,2) -  mass generation of solid-1 from solids-2.
C     Note, for example, that if gas is generated from solids-1 then
C     R_temp(0,1) > 0. The R-phase matrix is anti-symmetric. The diagonal 
C     elements are zero and are not needed.  Only one of the two anti-symmetric
C     elements -- e.g., R_temp(0,1) or R_temp(1,0) -- needs to be specified.
      DOUBLE PRECISION 
     &   R_temp (DIMENSION_RXN, 0:DIMENSION_M, 0:DIMENSION_M)
C
C                      Rate of mass transfer from phase M to Phase L
      DOUBLE PRECISION R_phase (DIMENSION_3, DIMENSION_LM+DIMENSION_M-1)
C
C                      Reaction names
      CHARACTER*10     RXN_NAME (DIMENSION_RXN)
C
C                      Species names
      CHARACTER*10     SPECIES_NAME (DIMENSION_N_all)
C
C                      Molecular weight of all species: MW_g and MW_s
      DOUBLE PRECISION MW_all (DIMENSION_N_all)
C
C                      total number of species
      INTEGER          N_all
C
C                      Species name to index
      INTEGER          SPECIES_ID2N(DIMENSION_N_all, 2)
C
C                      Species index to name
      INTEGER          SPECIES_N2IDg(DIMENSION_N_g)
      INTEGER          SPECIES_N2IDs(DIMENSION_M, DIMENSION_N_s)
C
C                      whether rxn schemes and rate expressions specified
      LOGICAL          GOT_RXN(DIMENSION_RXN), GOT_RATE (DIMENSION_RXN)
C
C                      stoichiometry
      DOUBLE PRECISION STOICH(DIMENSION_RXN, DIMENSION_N_all)
C
C                      stoichiometry x molecular weight
      DOUBLE PRECISION STOICHxMW(DIMENSION_RXN, DIMENSION_N_all)
C
C                      Enthalpy change due to reaction.  The change is assigned
C                      to the phase identified by Rate_m4T.
      DOUBLE PRECISION Delta_H(DIMENSION_RXN)
C
C                      Which tempearure to use in rate expression?  Also,
C                      the phase to which the enthalpy change due to rxn is
C                      assigned.
C                      m = 0 => fluid phase, m = 1 => solids phase 1, . . .
      Integer          Rate_m4T(DIMENSION_RXN)
C
C                      Factors in rate expression:
C                      1. preexponential factor (A),
C                      2. Temperature exponent (n),
C                      3. Activation temperature (E/R)
C                      rate expression = A * T_m**n * exp(-(E/R) / T_m) * . . .
      DOUBLE PRECISION Rate_fac(DIMENSION_RXN, 3)
C
C                      Exponents for concentartion dependence in the rate expression
      DOUBLE PRECISION Rate_exp(DIMENSION_RXN, DIMENSION_N_all)
C
C                      Actual number of Reactions
      INTEGER          NO_OF_RXNS
C
      COMMON /rxns_c/RXN_NAME, SPECIES_NAME
C
      COMMON /rxns_l/GOT_RXN, GOT_RATE
C
      COMMON /rxns_i/NO_OF_RXNS, SPECIES_ID2N, SPECIES_N2IDg,
     &               SPECIES_N2IDs,
     &               Rate_m4T, N_all
C
      COMMON /rxns_DP/
     &  R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g, SUM_R_s,
     &  R_phase, R_temp,
     &  STOICH, Rate_exp, Rate_fac, Delta_H, MW_all, STOICHxMW
