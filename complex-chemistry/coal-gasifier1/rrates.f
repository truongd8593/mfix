!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES(IER)                                            C
!  Purpose: Calculate gasification and combustion reaction rates       C
!                                                                      C
!  Author:  S. Venkatesan, M. Syamlal                 Date: 22-JUN-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Simplify MGAS chemistry; improve coal combustion rate      C
!           expression; add gas phase combustion reactions.            C
!  Author: M. Syamlal                                 Date: 15-JUL-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Separate char and coal reactions                           C
!  Author: M. Syamlal                                 Date: 24-AUG-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: Add sorbent reactions                                      C
!  Author: M. Syamlal                                 Date: 27-MAR-95  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: MFIX 2.0 changes, one solids phase                         C
!  Author: M. Syamlal                                 Date: 12-AUG-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!     Desai, P.R. and C.Y. Wen, "Computer Modeling of the MERC Fixed   C
!        Bed Gasifier," MERC/CR-78/3, March, 1978.                     C
!                                                                      C
!     Peters, N., "Premixed burning in diffusion flames -- The flame   C
!        zone model of Libby and Economos," Int. J. Heat Mass Transfer,C
!        Vol. 22, pp. 691-703, (1979).                                 C
!                                                                      C
!     Sergent, G.D. and I.W. Smith, "Combustion Rate of Bituminous     C
!        Coal Char in the Temperature Range 800 to 1700 K, Fuel,       C
!        Vol. 52, p.52 (1973).                                         C
!                                                                      C
!     Syamlal, M., and L.A. Bissett, "METC gasifier advanced           C
!       simulation (MGAS) model," DOE/METC-92/4108 (DE92001111),       C
!       January 1992                                                   C C
!     Syamlal, M., Derivation of shrinking core model, 1993            C
!                                                                      C
!     Wen, C.Y., H. Chen, and M. Onozaki, "User's Manual for Computer  C
!       Simulation and Design of the Moving Bed Coal Gasifier,"        C
!       DOE/MC/16474-1390 (DE83009533), January, 1982                  C
!                                                                      C
!     Westbrook, C.K., and F.L. Dryer, "Simplified mechanisms for the  C
!        oxidation of hydrocarbon fuels in flames," Combustion Science C
!        and technology, Vol. 27, pp. 31-43, (1981).                   C
!                                                                      C
!     Campbell, J.H., "The kinetics of decomposition of colorado oil   C
!        shale: II. Carbonate minerals," UCRL-52089 Part 2, National   C
!        Technical Information Service, Springfield, VA (1978).        C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables: TGSX, PO2, PCO, PCO2, PCH4, PH2, PH2O, DIFF,       C
!                   K_f, R_D, RXNA, K_a, TSX, TGX, RXNAF, RXNAB, EQ2,  C
!                   RXNB, RXNBF, RXNBB, EQ5, RXNC, EQ3, RXNCF, RXNCB,  C
!                   EQ6, RXND, RXNDF, RXNDB, A3, A4, A5, A6,           C
!                   RXNE, RXNEF, RXNEB, A7, A8, A9, RXNFF, RXNFB,      C
!                   RXNGF, RXNGB, RXNHF, RXNHB, RXNIF, RXNIB, K_r      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES(IER) 

      Use param
      Use param1
      Use fldvar
      Use geometry
      Use run
      Use indices
      Use physprop
      Use constant
      Use funits
      Use usr
      USE compar        !//d
      USE sendrecv      !// 400
 
      Use energy
      Use rxns
      
      IMPLICIT NONE
      
      INCLUDE  'usrnlst.inc'
!
!     Gas constant for O2 in cm^3.atm/g.K
      DOUBLE PRECISION R_O2
      PARAMETER (R_O2 = 82.06/32.)
!
!     Limit of maximum temperature (K)
      DOUBLE PRECISION MAX_TEMP
      PARAMETER (MAX_TEMP = 2500.0)
!
!     Limit of maximum sorbent temperature (K)
      DOUBLE PRECISION MAX_TSORB
      PARAMETER (MAX_TSORB = 1173.0)
!
!  Function subroutines
!
      LOGICAL COMPARE
!
!  Local Variables
!
!                      Local phase and cell ndex
      INTEGER          L, LM, IJK, M, N, IER
!
!                      Temperatures, Pressures, Proximate Analysis,
!                      Rate constants, activation energies,
!                      reaction rates
      DOUBLE PRECISION TGS1X, PO2, PCO, DIFF,&
                       PCO2, PCH4, PH2, PH2O, PATM, PATM_MW,&
                       K_f, R_D1, RXNA, K_a, RXNA1F, FAC,&
                       RXNA1B, CAR1, EP_s1,&
                       TGX, TS1X, EQ2, EQ3, EQ5, EQ6, A3, A4, A5,&
                       A6, A7, A8, A9,&
                       RXNB, RXNB1F, RXNB1B,&
                       RXNC, RXNC1F, RXNC1B,&
                       RXND, RXND1F, RXND1B,&
                       RXNE, RXNEF, RXNEB,&
                       RXNF0F, RXNF0B, RXNF1F, RXNF1B,&
                       RXNF2F, RXNF2B, RXNF3F, RXNF3B,&
                       RXNGF, RXNGB, RXNHF, RXNHB, RXNIF, RXNIB,&
                       RXNK1F, RXNK1B, RXNL1F,&
                       RXNL1B, EQCaO1,  K_r,&
                       X_coal1,TSORB1
      DOUBLE PRECISION R_tmp(0:MMAX, 0:MMAX)
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      R_tmp = UNDEFINED
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do private(ijk, R_tmp, L, LM, M, N, TGS1X, PO2, PCO, DIFF,&
!$omp&                       PCO2, PCH4, PH2, PH2O, PATM, PATM_MW,&
!$omp&                       K_f, R_D1, RXNA, K_a, RXNA1F, FAC,&
!$omp&                       RXNA1B, CAR1, EP_s1,&
!$omp&                       TGX, TS1X, EQ2, EQ3, EQ5, EQ6, A3, A4, A5,&
!$omp&                       A6, A7, A8, A9,&
!$omp&                       RXNB, RXNB1F, RXNB1B,&
!$omp&                       RXNC, RXNC1F, RXNC1B,&
!$omp&                       RXND, RXND1F, RXND1B,&
!$omp&                       RXNE, RXNEF, RXNEB,&
!$omp&                       RXNF0F, RXNF0B, RXNF1F, RXNF1B,&
!$omp&                       RXNF2F, RXNF2B, RXNF3F, RXNF3B,&
!$omp&                       RXNGF, RXNGB, RXNHF, RXNHB, RXNIF, RXNIB,&
!$omp&                       RXNK1F, RXNK1B, RXNL1F,&
!$omp&                       RXNL1B, EQCaO1,  K_r,&
!$omp&                       X_coal1,TSORB1)

      DO IJK = IJKSTART3, IJKEND3 
      
         R_gp(IJK, :) = ZERO
         RoX_gc(IJK, :) = ZERO
         R_sp(IJK, :, :) = ZERO
         RoX_sc(IJK, :, :) = ZERO
         SUM_R_G(IJK) = ZERO 
         HOR_G(IJK) = ZERO
         SUM_R_S(IJK, :) = ZERO 
         HOR_S(IJK, :) = ZERO 
	 R_PHASE(IJK, :) = ZERO
         D_p(1) = C(4)
         IF (FLUID_AT(IJK)) THEN 
!
!  User input is required in sections 1 through 4.
!
!1111111111111111111111111111111111111111111111111111111111111111111111111111111
!
!
! 1. Write the rates of various reactions:
!    Write the reaction rates for each of the reactions as RXNxF and RXNxB (both
!    quantities >= 0), where x identifies the reaction, F stands for forward
!    rate, and B stands for the backward rate.  The rates can be in
!    g-mole/(cm^3.s) or g/(cm^3.s).  For the sake of clarity, give the reaction
!    scheme and the units in a comment statement above the rate expression.
!    The volume (cm^3) is that of the computational cell.  Therefore, for
!    example, the rate term of a gas phase reaction will have a multiplicative
!    factor of epsilon. Note that X_g and X_s are mass fractions
!
!
!     The reaction rate expressions are taken mainly from
!     Wen et al. (1982), with the modifications proposed by Syamlal and
!     Bissett (1992).
!     MGAS chemistry was further modified:  all light hydrocarbons are lumped
!     to CH4, minor species H2S and NH3 are ignored, and combustion reactions
!     for H2, CH4, CO, and Tar are considered.
!
!     PAA  ash fraction -- Proximate analysis
!     PAFC Fixed carbon fraction -- Proximate analysis
!     EP_A  Ash layer void fraction
!
!              No   1   2    3    4   5    6   7    8
!     GAS Species  O2, CO, CO2, CH4, H2, H2O, N2, Tar
!
!                 No  1   2       3    4      5           6    7    8
!     Solids species FC, VM, Moist., Ash, CaCO3, CaMg(CO3)2, CaO, MgO

!  remove comments tagged !car to activate carbonate decomposition
!
 
          TGX   = MIN(MAX_TEMP, T_g(IJK))
          TS1X  = MIN(MAX_TEMP, T_s(IJK, 1))
          TSORB1= MIN(MAX_TSORB, T_s(IJK, 1))
          TGS1X = HALF * (TGX + TS1X)
!
!         COMPUTE PARTIAL PRESSURE of various gases IN ATM.  P_g in
!         dynes/cm^2
!
          PATM    = P_g(IJK) / 1013000.
          PATM_MW = PATM * MW_MIX_g(IJK)
          PO2     = PATM_MW * X_g(IJK, 1) / MW_g(1)
          PCO     = PATM_MW * X_g(IJK, 2) / MW_g(2)
          PCO2    = PATM_MW * X_g(IJK, 3) / MW_g(3)
          PCH4    = PATM_MW * X_g(IJK, 4) / MW_g(4)
          PH2     = PATM_MW * X_g(IJK, 5) / MW_g(5)
          PH2O    = PATM_MW * X_g(IJK, 6) / MW_g(6)
!
          EP_s1  = EP_s(IJK,1)
          X_coal1 = X_s(IJK, 1, 1) + X_s(IJK, 1, 2) + X_s(IJK, 1, 3)&
                 + X_s(IJK, 1, 4)
!
!         concentration of carbon, gmole/cc
!
          CAR1   = ROP_s(IJK,1) * X_s(IJK, 1, 1) / MW_s(1,1)
!
!         Combustion reactions
!
!
          RXNA1F = ZERO
          RXNF0F = ZERO
          RXNF1F = ZERO
          RXNF2F = ZERO
          RXNF3F = ZERO
          RXNA1B = ZERO
          RXNF0B = ZERO
          RXNF1B = ZERO
          RXNF2B = ZERO
          RXNF3B = ZERO
          IF ( PO2 .GT. ZERO ) THEN
!
!         a)   2C + O2 --> 2CO		g-mole/(cm^3.s)
!
!         Wen at al. (1982), Syamlal and Bissett (1992), Syamlal (1993)
!         Intrinsic rate from Desai and Wen (1978), originally from
!         Sergeant and Smith (1973).
!
            IF ( .NOT.COMPARE(EP_g(IJK), ONE) ) THEN
              IF(PAFC .NE. ZERO) THEN
                IF(X_s(IJK,1,4) .GT. ZERO) THEN
                  R_D1 = (X_s(IJK,1,1) * PAA&
                     / (X_s(IJK,1,4) * PAFC))**(1./3.)
                  R_D1 = MIN(ONE, R_D1)
                ELSE
                  R_D1 = ZERO
                ENDIF
              ELSE
                R_D1 = ZERO
              ENDIF
!
              DIFF = 4.26 * ((TGX/1800.)**1.75) / PATM
!

              IF(R_D1 .EQ. ZERO .OR. EP_s1 .EQ. ZERO) THEN
                RXNA1F = ZERO
              ELSE
                K_f = DIFF * N_sh(IJK, 1) / (D_p(1) * R_O2 * TGX)
                K_r = 8710. * EXP( -27000./1.987/TS1X) * R_D1*R_D1
!
                IF(R_D1 .GE. ONE) THEN
                  RXNA = ONE / (ONE / K_f + ONE / K_r)
                ELSE
                  K_a = 2. * DIFF * f_EP_A * R_D1&
                        / ( D_p(1) * ( ONE - R_D1 ) * R_O2 * TS1X )
                  RXNA = ONE / ( ONE /K_f + ONE/K_a + ONE/ K_r)
                ENDIF
                FAC = X_s(IJK,1,1)/(X_s(IJK,1,1) + 1.0E-6)
                RXNA1F = RXNA * PO2 * FAC * 6.0 * EP_s1&
                         / ( D_p(1) * 32.0 )
              ENDIF
            ENDIF
!
!           f0)  HYDROGEN COMBUSTION: 2H2 + O2 --> 2H2O  (mol/cm^3.s)
!           Peters (1979)
            IF(PH2 .GT. ZERO) THEN
              RXNF0F = 1.08E6 * EXP(-30000.0/(1.987*TGX)) * EP_g(IJK)&
                       * (RO_g(IJK)*X_g(IJK,1)/MW_g(1))&
                       * (RO_g(IJK)*X_g(IJK,5)/MW_g(5))
            ENDIF
!
!           f1)  METHANE COMBUSTION: CH4 + 2O2 --> CO2 + 2H2O
!           (mol/cm^3.s)   Westbrook and Dryer (1981)
            IF(PCH4 .GT. ZERO ) THEN
              RXNF1F = 6.7E12 * EXP(-48400.0/(1.987*TGX)) * EP_g(IJK)&
                * (RO_g(IJK)*X_g(IJK,1)/MW_g(1)) ** 1.3&
                * (RO_g(IJK)*X_g(IJK,4)/MW_g(4)) ** 0.2
            ENDIF
!
!           f2)  CO COMBUSTION: CO + 1/2O2 --> CO2 (mol/cm^3.s)
!           Westbrook and Dryer (1981)
            IF(PCO .GT. ZERO) THEN
              RXNF2F = 3.98E14 * EXP(-40000.0/(1.987*TGX)) * EP_g(IJK)&
                * (RO_g(IJK)*X_g(IJK,1)/MW_g(1)) ** 0.25&
                * (RO_g(IJK)*X_g(IJK,2)/MW_g(2))&
                * (RO_g(IJK)*X_g(IJK,6)/MW_g(6)) ** 0.5
            ENDIF
!
!           f3)  TAR COMBUSTION: Tar + f3_1 O2 --> f3_3 CO2 + f3_6 H2O
!           (mol/cm^3.s)  Westbrook and Dryer (1981), assuming that tar
!           is similar to C10H22.
            IF(X_g(IJK,8) .GT. ZERO) THEN
              RXNF3F = ZERO !3.8E11 * EXP(-30000.0/(1.987*TGX)) * EP_g(IJK)&
               ! * (RO_g(IJK)*X_g(IJK,1)/MW_g(1)) ** 1.5&
               ! * (RO_g(IJK)*X_g(IJK,8)/MW_g(8)) ** 0.25
            ENDIF
          ENDIF
!
!         Coal and Char phase reactions
!
          IF ( EP_s1 .GT. ZERO ) THEN
!
!           b) C + H2O --> CO + H2	g-mole/(cm^3.s)
!		Wen et al. (1982)
!
            EQ2 = EXP ( 17.2931 - 16326.1 / TGS1X )
            RXNB = AK2*EXP(-AE2/(1.987*TGS1X))*CAR1
            RXNB1F = RXNB *  PH2O
            RXNB1B = RXNB * PH2 * PCO / EQ2
!
!           c) C + CO2 --> 2CO		g-mole/(cm^3.s)
!		Wen et al. (1982)
!
            EQ5 = EXP ( 20.9238 - 20281.8 / TGS1X )
            RXNC = AK5*EXP(-AE5/(1.987*TGS1X))*CAR1
            RXNC1F = RXNC * PCO2
            RXNC1B = RXNC * PCO*PCO/EQ5
!
!           d) (1/2)C + H2 --> (1/2)CH4	g-mole/(cm^3.s)
!		Wen et al. (1982)
!
            EQ6 = EXP ( - 13.4289 + 10998.5 / TGS1X )
            RXND = EXP(  - 7.0869 - 8077.5 / TGS1X ) * CAR1
            RXND1F = RXND * PH2
            RXND1B = RXND * SQRT ( MAX(PCH4,ZERO) / EQ6 )
!
!           INITIAL STAGE KINETICS
!           g)    COAL MOISTURE  --> H2O  (ARBITRARY RATE EXPRESSION)
!        	g/(cm^3.s)
!
            RXNGF = AKM*EXP(-AEM/(1.987*TS1X))*ROP_s(IJK,1)*X_s(IJK,1,3)
            RXNGB = ZERO
!
!           h) VOLATILE MATTER  -->    TAR   +   GASES	    g/(cm^3.s)
!              (1.0)                (alphad)    (1-alphad )
!              (1-alphad = Sum of betad's )
!
            IF( ((X_s(IJK,1,2) / X_coal1) - VMSTAR(IJK)) .GT. ZERO)THEN
              RXNHF = AKD*EXP(-AED/(1.987*TS1X))*ROP_s(IJK,1)&
                      * (X_s(IJK,1,2) / X_coal1)
              RXNHB = AKD*EXP(-AED/(1.987*TS1X))*ROP_s(IJK,1)&
      	              * VMSTAR(IJK)
            ELSE
              RXNHF = ZERO
              RXNHB = ZERO
            ENDIF
!
!           k)   CaMg(CO3)_2 ---> CaCO3 + MgO + CO2	g-mole/(cm^3.s)
!
            RXNK1F = ZERO !2.E08 * EXP(-51000./(1.987*TSORB1)) * ROP_s(IJK,1)&
!car      	             * X_s(IJK, 1, 6)/MW_s(1,6)
            RXNK1B = ZERO
!car
!car            l)   CaCO3 ----> CaO + CO2		g-mole/(cm^3.s)	
!car
             RXNL1F = ZERO !1.3E10 * EXP(-55000./(1.987*TSORB1))&
! car     	              * ROP_s(IJK,1) * X_s(IJK, 1, 5)/MW_s(1,5)
!car             EQCaO1 = 1.03E08 * EXP(-21830./TSORB1)
             RXNL1B = ZERO !RXNL1F * PCO2 / EQCaO1&
!car                    * ( X_s(IJK, 1, 7)/(X_s(IJK, 1, 7)+1.e-4) )
           ELSE
             RXNB1F = ZERO
             RXNB1B = ZERO
             RXNC1F = ZERO
             RXNC1B = ZERO
             RXND1F = ZERO
             RXND1B = ZERO
             RXNGF  = ZERO
             RXNGB  = ZERO
             RXNHF  = ZERO
             RXNHB  = ZERO
             RXNK1F = ZERO
             RXNK1B = ZERO
             RXNL1F = ZERO
             RXNL1B = ZERO
           ENDIF
!
!          e)  WATER-GAS SHIFT REACTION: CO + H2O --> CO2 + H2
!                    	g-mole/(cm^3.s)
!          It is assumed that the shift reaction is catalyzed by char
!          and that the heat of reaction is to be added to char.
!          Wen et al. (1982)
           A3 = WG3 * 2.877E+05 * EXP(-27760.0/(1.987*TGS1X))
           EQ3 = EXP(-3.63061+3955.71/TGS1X)
           A4 = (PATM**(0.5-PATM/250.))/PATM/PATM
           A5 = EP_s1*PAA*RO_s(1)*EXP(-8.91+5553.0/TGS1X)
           RXNE = A3*A4*A5*EP_g(IJK)
           RXNEF = RXNE * PCO*PH2O
           RXNEB = RXNE * PCO2*PH2/EQ3
!
!         i) TAR  -->   CHAR    +  GASES		g/(cm^3.s)
!           (1.0)     (alphac)    (1-alphac = Sum of betac's)
!         It is assumed that tar cracking is catalyzed by char and that
!         the resulting carbon is deposited on char.
!
          IF(EP_s1 .GT. ZERO) THEN
            RXNIF = AKC*EXP(-AEC/(1.987*TGX))*ROP_g(IJK)*X_g(IJK, 8)
          ELSE
            RXNIF = ZERO
          ENDIF
          RXNIB = ZERO
!
!222222222222222222222222222222222222222222222222222222222222222222222222222222
!
! 2. Write the formation and consumption rates of various species:
!    Obtain the rates of formation and consumption of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.  All
!    the formation rates for gas species n are added to get R_gp (IJK, n).
!    All the consumption rates are added and then divided by X_g(IJK, n) to
!    get RoX_gc(IJK, n).  If X_g(IJK, n) is zero and species n is likely
!    to be consumed in a reaction then it is recommended that RoX_gc (IJK, n)
!    be initialized to the derivative of the consumption rate with respect to
!    X_g at X_g=0.
!    If the slope is not known analytically a small value such as 1.0e-9 may
!    instead be used.  A similar procedure is used for all the species in the
!    solids phases also.
!
!  GAS SPECIES
!  the gaseous species expressions have the forms:
!                                      R_gp(cell#, species#)
!                                      RoX_gc(cell#, species#)
!
!  (1) O2
      R_gp(IJK, 1) = ZERO
      IF(X_g(IJK, 1) .GT. ZERO) THEN
        RoX_gc(IJK, 1) = (RXNA1F + RXNF0F + 2. * RXNF1F&
                          + HALF * RXNF2F + f3_1 * RXNF3F) * MW_g(1)&
                          / X_g(IJK, 1)
      ELSE
        RoX_gc(IJK, 1) = 1.0e-9
      ENDIF
!
!  (2) CO
      R_gp(IJK, 2) = (2. * (RXNA1F ) + RXNB1F&
                      + 2. * (RXNC1F ) + RXNEB) * MW_g(2)&
                     + (RXNHF - RXNHB) * BETAD(2)&
                     + RXNIF * BETAC(2)
      IF(X_g(IJK, 2) .GT. ZERO) THEN
        RoX_gc(IJK, 2) = (RXNB1B + 2. * (RXNC1B )&
                          + RXNEF + RXNF2F) * MW_g(2) / X_g(IJK, 2)
      ELSE
        RoX_gc(IJK, 2) = 1.0e-9
      ENDIF
!
!  (3) CO2
      R_gp(IJK, 3) = (RXNC1B  + RXNEF + RXNF1F + RXNF2F&
                      + F3_3 * RXNF3F &
		      + RXNK1F&
                  + RXNL1F&
		     ) * MW_g(3)&
                     + (RXNHF - RXNHB) * BETAD(3)&
                     + RXNIF * BETAC(3)
      IF(X_g(IJK, 3) .GT. ZERO) THEN
        RoX_gc(IJK, 3) = (RXNC1F + RXNEB + RXNL1B )&
                         * MW_g(3) / X_g(IJK, 3)
      ELSE
        RoX_gc(IJK, 3) = 1.0e-9
      ENDIF
!
!  (4) CH4
      R_gp(IJK, 4) = HALF * (RXND1F ) * MW_g(4)&
                     + (RXNHF - RXNHB) * BETAD(4) + RXNIF * BETAC(4)
      IF(X_g(IJK, 4) .GT. ZERO) THEN
        RoX_gc(IJK, 4) = ((RXND1B ) * HALF + RXNF1F)* MW_g(4)&
                         / X_g(IJK, 4)
      ELSE
        RoX_gc(IJK, 4) = 1.0e-9
      ENDIF
!
!  (5) H2
      R_gp(IJK, 5) = (RXNB1F  + RXND1B  + RXNEF)&
                     * MW_g(5) + (RXNHF - RXNHB) * BETAD(5)&
                     + RXNIF * BETAC(5)
      IF(X_g(IJK, 5) .GT. ZERO) THEN
        RoX_gc(IJK, 5) = (RXNB1B  + RXND1F  + RXNEB&
                          + 2. * RXNF0F) * MW_g(5) / X_g(IJK, 5)
      ELSE
        RoX_gc(IJK, 5) = 1.0E-9
      ENDIF
!
!  (6) H2O
      R_gp(IJK, 6) = (RXNB1B + RXNEB + 2. * RXNF0F&
                      + 2. * RXNF1F + F3_6 * RXNF3F) * MW_g(6)&
                     + RXNGF + (RXNHF - RXNHB) * BETAD(6)&
                     + RXNIF * BETAC(6)
      IF(X_g(IJK, 6) .GT. ZERO) THEN
        RoX_gc(IJK, 6) = (RXNB1F  + RXNEF) * MW_g(6)&
                         / X_g(IJK, 6)
      ELSE
        RoX_gc(IJK, 6) = 1.0e-9
      ENDIF
!
!   (7) N2
!     inert here
      R_gp(IJK, 7) = ZERO
      RoX_gc(IJK, 7) = ZERO
!
!   (8) TAR
      R_gp(IJK, 8) = (RXNHF - RXNHB) * ALPHAD
      IF(X_g(IJK, 8) .GT. SMALL_NUMBER) THEN
        RoX_gc(IJK, 8) = (RXNF3F * MW_g(8) + RXNIF) / X_g(IJK, 8)
      ELSE
        RoX_gc(IJK, 8) = 1.0e-9
      ENDIF
!
!  SOLIDS SPECIES
!  COAL PSEUDO-SPECIES
!  the solid species' expressions have the form:
!                                   R_sp(cell#,solid phase#,solid species#)
!                                   RoX_sc(cell#,solid phase#,solid species#)
!
! (1)  CARBON
      R_sp(IJK, 1, 1) = (RXNB1B + RXNC1B + HALF * RXND1B) *&
                        MW_s(1, 1) +  RXNIF * ALPHAC
      IF(X_s(IJK, 1, 1) .GT. ZERO) THEN
        RoX_sc(IJK, 1, 1) = (2. * RXNA1F + RXNB1F + RXNC1F&
                             + HALF * RXND1F) * MW_s(1, 1)&
                             / X_s(IJK, 1, 1)
      ELSE
        RoX_sc(IJK, 1, 1) = 1.0e-7
      ENDIF
!
! (2)  VOLATILE MATTER
!
      R_sp(IJK, 1, 2) = RXNHB
      IF(X_s(IJK, 1, 2) .GT. ZERO) THEN
        RoX_sc(IJK, 1, 2) = RXNHF/X_s(IJK, 1, 2)
      ELSE
        RoX_sc(IJK, 1, 2) = 1.0e-7
      ENDIF
!
! (3)  MOISTURE
      R_sp(IJK, 1, 3) = RXNGB
      IF(X_s(IJK, 1, 3) .GT. ZERO) THEN
        RoX_sc(IJK, 1, 3) = RXNGF/X_s(IJK, 1, 3)
      ELSE
        RoX_sc(IJK, 1, 3) = 1.0e-7
      ENDIF
!
! (4) ASH
      R_sp(IJK, 1, 4) = ZERO
      RoX_sc(IJK,1,4) = ZERO
!
! (5) CaCO3
      R_sp(IJK, 1, 5) = (RXNK1F + RXNL1B )* MW_s(1,5)
     IF(X_s(IJK, 1, 5) .GT. ZERO) THEN
        RoX_sc(IJK,1,5) = RXNL1F * MW_s(1,5) / X_s(IJK, 1, 5)
      ELSE
        RoX_sc(IJK,1,5) = ZERO
      ENDIF
!
! (6) CaMg(CO3)2
 
      R_sp(IJK, 1, 6) = ZERO
      IF(X_s(IJK, 1, 6) .GT. ZERO) THEN
        RoX_sc(IJK,1,6) = RXNK1F * MW_s(1,6) / X_s(IJK, 1, 6)
      ELSE
        RoX_sc(IJK,1,6) = ZERO
      ENDIF
!
! (7) CaO
 
     R_sp(IJK, 1, 7) = RXNL1F * MW_s(1,7)
      IF(X_s(IJK, 1, 7) .GT. ZERO) THEN
        RoX_sc(IJK,1,7) = RXNL1B * MW_s(1,7) / X_s(IJK, 1, 7)
      ELSE
        RoX_sc(IJK,1,7) = ZERO
      ENDIF
!
! (8) MgO
      R_sp(IJK, 1, 8) = RXNK1F  * MW_s(1,8)
      RoX_sc(IJK,1,8) = ZERO
!
!
!3333333333333333333333333333333333333333333333333333333333333333333333333333333
!
! 3.  Determine the g/(cm^3.s) transferred from one phase to the other.
!          R_tmp(To phase #, From phase #)
!     e.g. R_tmp(0,1) -  mass generation of gas phase from solids-1,
!          R_tmp(0,2) -  mass generation of gas phase from solids-2,
!          R_tmp(1,0) -  mass generation of solid-1 from gas = -R_tmp(0,1)
!          R_tmp(1,2) -  mass generation of solid-1 from solids-2.
!     Note, for example, that if gas is generated from solids-1 then
!     R_tmp(0,1) > 0.
!     The R-phase matrix is skew-symmetric and diagonal elements are not needed.
!     Only one of the two skew-symmetric elements -- e.g., R_tmp(0,1) or
!     R_tmp(1,0) -- needs to be specified.
!
!
      R_tmp(0,1) =    RXNA1F * (2. * MW_s(1,1) )&
                    + (RXNB1F - RXNB1B) * MW_s(1,1)&
                    + (RXNC1F - RXNC1B) * MW_s(1,1)&
                    + (RXND1F - RXND1B) * (HALF * MW_s(1,1))&
                    +  RXNGF + (RXNHF - RXNHB)&
                    - RXNIF * ALPHAC  &
                   + (RXNK1F + RXNL1F - RXNL1B) * MW_g(3)
!
!
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
!
! 4.  Determine the heat of reactions in cal/(cm^3.s) at the
!     temperature T_g or T_s.  Note that for exothermic reactions
!     HOR_g or HOR_s should be negative. The assignment of heat of reaction
!     is user defined as it depends upon the microphysics near the interface,
!     which is averaged out in the multiphase flow equations.  For example,
!     heat of Reaction for the C + O2 reaction is split into parts;
!     CO formation is assigned to the solid phase and CO2 formation from CO to
!     the gas phase.
!
      HOR_g(IJK)    =  (-115596.0) * (RXNF0F - RXNF0B)&
                     + (-191759.0) * (RXNF1F - RXNF1B)&
                     +  (-67636.0) * (RXNF2F - RXNF2B)&
                     +  (HEATF3)   * (RXNF3F - RXNF3B)&
                     +       HEATC * (RXNIF - RXNIB)&
                     + (31000.0) * (RXNK1F )&
                    + (41000.0) * (RXNL1F  - RXNL1B )
      HOR_s(IJK, 1) =   (-52832.0) * (RXNA1F - RXNA1B)&
                     +   (31382.0) * (RXNB1F - RXNB1B)&
                     +   (41220.0) * (RXNC1F - RXNC1B)&
                     +   (-8944.5) * (RXND1F - RXND1B)&
                     +   (-9838.0) * (RXNEF - RXNEB)&
                     +     (540.5) *  RXNGF&
                     +       HEATD * (RXNHF - RXNHB)
!
!==============================================================================
!
!     No user input is required below this line
!-----------------------------------------------------------------------------
!   Determine g/(cm^3.s) of mass generation for each of the phases by adding
!   the reaction rates of all the individual species.

            SUM_R_G(IJK) = ZERO 
            IF (SPECIES_EQ(0)) THEN 
               IF (NMAX(0) > 0) THEN 
                  SUM_R_G(IJK) = SUM_R_G(IJK) + SUM(R_GP(IJK,:NMAX(0))-ROX_GC(&
                     IJK,:NMAX(0))*X_G(IJK,:NMAX(0))) 
               ENDIF 
            ENDIF 
!
            DO M = 1, MMAX 
               SUM_R_S(IJK,M) = ZERO 
               IF (SPECIES_EQ(M)) THEN 
                  IF (NMAX(M) > 0) THEN 
                     SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + SUM(R_SP(IJK,M,:NMAX(M))&
                        -ROX_SC(IJK,M,:NMAX(M))*X_S(IJK,M,:NMAX(M))) 
                  ENDIF 
               ENDIF 
            END DO 
	    
!
!
!     Store R_tmp values in an array.  Only store the upper triangle without
!     the diagonal of R_tmp array.
!
            DO L = 0, MMAX 
               DO M = L + 1, MMAX 
                  LM = L + 1 + (M - 1)*M/2 
                  IF (R_TMP(L,M) /= UNDEFINED) THEN 
                     R_PHASE(IJK,LM) = R_TMP(L,M) 
                  ELSE IF (R_TMP(M,L) /= UNDEFINED) THEN 
                     R_PHASE(IJK,LM) = -R_TMP(M,L) 
                  ELSE 
                     CALL START_LOG 
                     WRITE (UNIT_LOG, 1000) L, M 
                     CALL END_LOG 
                     CALL MFIX_EXIT
 !                    call mfix_exit(myPE)  
                  ENDIF 
               END DO 
            END DO 
	   
         ENDIF 
      END DO 
      
1000  FORMAT(/1X,70('*')//' From: RRATES',&
      /' Message: Mass transfer between phases ', I2, ' and ', I2,&
      ' (R_temp) not specified', /1X, 70('*')/)
      RETURN
      END SUBROUTINE RRATES
