!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES(IER)                                            C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!  Literature/Document References:                                     C
!                                                                      C
!     Desai, P.R. and C.Y. Wen, "Computer Modeling of the MER! Fixed   C
!        Bed Gasifier," MERC/CR-78/3, March, 1978.                     C
!                                                                      C
!     Syamlal, M., Rogers, W., and T.J. O'Brien, "MFIX Documentation   C
!        Theory Guide, DOE/METC-94/1004, NTIS/DE94000087, National     C
!        Technical Information Service, Springfield, VA, 1993.         C
!                                                                      C
!     Wen, C.Y., H. Chen, and M. Onozaki, "User's Manual for Computer  C
!       Simulation and Design of the Moving Bed Coal Gasifier,"        C
!       DOE/MC/16474-1390 (DE83009533), January, 1982                  C
!                                                                      C
!     Westbrook, C.K., and F.L. Dryer, "Simplified mechanisms for the  C
!        oxidation of hydrocarbon fuels in flames," Combustion Science C
!        and technology, Vol. 27, pp. 31-43, (1981).                   C
!                                                                      C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits 
      USE toleranc
      USE compar        !//d
      USE sendrecv      !// 400
      USE usr
      IMPLICIT NONE
      
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!
!     Gas constant for O2
      DOUBLE PRECISION R_O2
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Local phase and species indices
      INTEGER          L, LM, M, N

!                      cell index
      INTEGER          IJK
      
      DOUBLE PRECISION R_tmp(0:MMAX, 0:MMAX), rxn
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
!
!-----------------------------------------------
!
!  Function subroutines
!
      LOGICAL COMPARE
!
!  Local Variables
!
!                      Temperatures, Pressures, Proximate Analysis,
!                      Rate constants, activation energies,
!                      reaction rates
      DOUBLE PRECISION TGS1X, PO2, PCO, DIFF, &
                       PCO2, PATM, PATM_MW, &
                       K_f, R_D1, RXNA, K_a, RXNA1F, &
                       RXNA1B, CAR1, EP_s1,&
                       TGX, TS1X, EQ5, RXNB, RXNB1F, &
                       RXNB1B, RXNCF, RXNCB, K_r
!
      INCLUDE 'usrnlst.inc'
      INCLUDE 'species_indices.inc' 
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'

      R_O2 = GAS_CONST*9.86923E-7/32.  !cm^3.atm/g.K
      IER = 0
      R_tmp = UNDEFINED
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do private(ijk, R_tmp, L, LM, M, N)
      DO IJK = IJKSTART3, IJKEND3 
      
         IF (FLUID_AT(IJK)) THEN 
!
!
!  User input is required in sections 1 through 4.
!
!1111111111111111111111111111111111111111111111111111111111111111111111111111111
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
!     Hot Char  - solids-1
!     Cold Char - solids-2
!
!     PAA  -- Ash fraction from Proximate analysis
!     PAFC -- Fixed carbon fraction from Proximate analysis
!     f_EP_A -- Function of the ash-layer void fraction
!
!              No   1   2    3   4
!     GAS Species  O2, CO, CO2, N2
!
      TGX   = MIN(TMAX, T_g(IJK))
      TS1X  = MIN(TMAX, T_s(IJK,1))
      TGS1X = HALF * (TGX + TS1X)
 
!
!     COMPUTE PARTIAL PRESSURE of various gases IN ATM.  P_g in dynes/cm^2
!
      PATM    = P_g(IJK) / 1013000.
      PATM_MW = PATM * MW_MIX_g(IJK)
      PO2     = PATM_MW * X_g(IJK, O2) / MW_g(O2)
      PCO     = PATM_MW * X_g(IJK, CO) / MW_g(CO)
      PCO2    = PATM_MW * X_g(IJK, CO2) / MW_g(CO2)
!
      EP_s1  = EP_s(IJK,1)
!
!     concentration of carbon, gmole/cc
!
      CAR1   = ROP_s(IJK,1) * X_s(IJK, 1, FC) / MW_s(1,FC)
!
!  a)   COMBUSTION: 2C + O2 --> 2CO; g-mole/(cm^3.s)
!       Wen at al. (1982), Syamlal et al. (1993), Desai and Wen (1978)
!
      IF ( PO2 .GT. ZERO .AND. .NOT.COMPARE(EP_g(IJK), ONE) ) THEN
        IF(PAFC .NE. ZERO) THEN
          IF(X_s(IJK,1,FC) .GT. ZERO) THEN
            R_D1 = (X_s(IJK,1,FC) * PAA / (X_s(IJK,1,Ash) * PAFC))**(1./3.)
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
          K_f = DIFF * N_sh(IJK, 1) / (D_p0(1) * R_O2 * TGX)
          K_r = 8710. * EXP( -27000./1.987/TS1X) * R_D1*R_D1

          IF(R_D1 .GE. ONE) THEN
            RXNA = ONE / (ONE / K_f + ONE / K_r)
          ELSE
            K_a = 2. * DIFF * f_EP_A * R_D1 &
                / ( D_p0(1) * ( ONE - R_D1 ) * R_O2 * TS1X )
            RXNA = ONE / ( ONE /K_f + ONE/K_a + ONE/ K_r)
          ENDIF
          RXNA1F = RXNA * PO2 * 6.0 * EP_s1 / ( D_p0(1) * 32.0 )
        ENDIF
      ELSE
        RXNA1F = ZERO
      ENDIF
      RXNA1B = ZERO
!                                                            
! b)  CHAR-CO2 REACTION: C + CO2 --> 2CO;   g-mole/(cm^3.s)
!     Wen et al. (1982) rate was increased by a factor of 1000.
      IF ( EP_s1 .GT. ZERO ) THEN
        EQ5 = EXP ( 20.9238 - 20281.8 / TGS1X )                    
        RXNB = 1000.0 * 930.0*EXP(-45000./(1.987*TGS1X))*CAR1
        RXNB1F = RXNB * PCO2
        RXNB1B = RXNB * PCO*PCO/EQ5
      ELSE
        RXNB1F = ZERO
        RXNB1B = ZERO
      ENDIF
!
! c)  CO COMBUSTION: CO + 1/2O2 --> CO2; g-mol/(cm^3.s)
!      Westbrook and Dryer (1981)
      IF(PCO .GT. ZERO .AND. PO2 .GT. ZERO) THEN
        RXNCF =  3.98E14 * EXP(-40000.0/(1.987*TGX)) * EP_g(IJK) &
                * (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2)) ** 0.25 &
                * (RO_g(IJK)*X_g(IJK,CO)/MW_g(CO)) &
                * (RO_g(IJK)*0.1/18.) ** 0.5
!                * (RO_g(IJK)*X_g(IJK,H2O)/MW_g(H2O)) ** 0.5
      ELSE
        RXNCF = ZERO
      ENDIF
      RXNCB = ZERO
       
!  The following assignments are for writing out reaction rate data into .SPA
!  file and for visualizing the reaction rates with animate_mfix.  Also set in 
!  mfix.dat nRR to the number of reaction rates stored (7 in this case) and 
!  SPX_Dt, the intervals for writing .SPA file. The if statements are used
!  so that the user may change nRR from 0 to 7 in the datafile.
!  These arrays are not used for any calculations or for restart.
       if(nRR >= 1)ReactionRates(ijk, 1) = RXNA1F
       if(nRR >= 2)ReactionRates(ijk, 3) = RXNB1F - RXNB1B
       if(nRR >= 3)ReactionRates(ijk, 5) = RXNCF

!
!
!2222222222222222222222222222222222222222222222222222222222222222222222222222222
!
! 2. Write the formation and consumption rates of various species:
!    Obtain the rates of formation and consumption of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.
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
!

!  (1) O2
      R_gp(IJK, O2) = ZERO
      IF(X_g(IJK, O2) .GT. ZERO) THEN
        RoX_gc(IJK, O2) = (RXNA1F + HALF * RXNCF ) * MW_g(O2) &
                          / X_g(IJK, O2)
      ELSE
        RoX_gc(IJK, O2) = 1.0e-9
      ENDIF
!
!  (2) CO
      R_gp(IJK, CO) = (2. * (RXNA1F ) + 2. * (RXNB1F ) ) &
                      * MW_g(CO)
      IF(X_g(IJK, CO) .GT. ZERO) THEN
        RoX_gc(IJK, CO) = ( 2. * (RXNB1B ) + RXNCF) * MW_g(CO) &
                          / X_g(IJK, CO)
      ELSE
        RoX_gc(IJK, CO) = 1.0e-9
      ENDIF
!
!  (3) CO2
      R_gp(IJK, CO2) = (RXNB1B + RXNCF ) * MW_g(CO2)
      IF(X_g(IJK, CO2) .GT. ZERO) THEN
        RoX_gc(IJK, CO2) = (RXNB1F ) * MW_g(CO2) / X_g(IJK, CO2)
      ELSE
        RoX_gc(IJK, CO2) = 1.0e-9
      ENDIF
!
!   (7) N2
!     inert here
      R_gp(IJK, N2) = ZERO
      RoX_gc(IJK, N2) = ZERO

!
!  SOLIDS SPECIES
!
! (1)  CARBON
      R_sp(IJK, 1, FC) = (RXNB1B) * MW_s(1, FC)
      IF(X_s(IJK, 1, FC) .GT. ZERO) THEN
        RoX_sc(IJK, 1, FC) = (2. * RXNA1F + RXNB1F) * MW_s(1, FC) &
                             / X_s(IJK, 1, FC)
      ELSE
        RoX_sc(IJK, 1, FC) = 1.0e-7
      ENDIF
!
!
! (2) ASH
      R_sp(IJK, 1, Ash) = ZERO
      RoX_sc(IJK,1,Ash) = ZERO
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
      R_tmp(0,1) = RXNA1F *(2.*MW_s(1,FC))+ (RXNB1F-RXNB1B)*MW_s(1,FC)

!
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
!
! 4.  Determine the heat of reactions in cal/(cm^3.s) at the
!     temperature T_g or T_s1.  Note that for exothermic reactions
!     HOR_g (or HOR_s) will be negative. The assignment of heat of reaction
!     is user defined as it depends upon the microphysics near the interface,
!     which is averaged out in the multiphase flow equations.  For example,
!     heat of Reaction for the C + O2 reaction is split into parts;
!     CO formation is assigned to the solid phase and CO2 formation from CO to
!     the gas phase.
!     *** This section is no longer needed as the heats of reactions are  
!         calculated below.  If you need to override the automatic calculation, 
!         comment out the calculations below.   
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
	    ELSE
	      DO M = 1, MMAX
	        IF(R_tmp(0,M) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) + R_tmp(0,M)
		ELSEIF(R_tmp(M,0) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) - R_tmp(M,0)
		ENDIF
	      ENDDO 
            ENDIF 
!
            DO M = 1, MMAX 
               SUM_R_S(IJK,M) = ZERO 
               IF (SPECIES_EQ(M)) THEN 
                  IF (NMAX(M) > 0) THEN 
                     SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + SUM(R_SP(IJK,M,:NMAX(M))&
                        -ROX_SC(IJK,M,:NMAX(M))*X_S(IJK,M,:NMAX(M))) 
                  ENDIF 
	       ELSE
 	         DO L = 0, MMAX
	           IF(R_tmp(M,L) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) + R_tmp(M,L)
		   ELSEIF(R_tmp(L,M) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) - R_tmp(L,M)
		   ENDIF
	         ENDDO 
               ENDIF 
            END DO 
	    
!
!
!
            IF(ENERGY_EQ) THEN ! calculate heat of reactions only if energy eq. are solved
            HOR_G(IJK) = zero
            DO N = 1, NMAX(0)
	      rxn =  (R_gp(IJK, N) - RoX_gc(IJK, N) * X_g(IJK, N)) &
	              - SUM_R_g(IJK) * X_g(IJK, N)
              HOR_G(IJK) = HOR_G(IJK) + rxn * (HfrefoR_g(N)  + &
	                     (calc_ICpoR(T_G(IJK), Thigh_g(N), Tlow_g(N), &
			     Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N)) -  &
			      IC_PGrefoR(N)) )* GAS_CONST_cal / MW_g(N)
	                   
            END DO 
            IF (UNITS == 'SI') HOR_G(IJK) = 4183.925D0*HOR_G(IJK)    !in J/kg K
	    
            DO M = 1, MMAX 
              HOR_s(IJK, M) = zero
              DO N = 1, NMAX(M)
	        rxn =  (R_sp(IJK, M, N) - RoX_sc(IJK, M, N) * X_s(IJK, M, N)) &
	              - SUM_R_s(IJK, M) * X_s(IJK, M, N)
                HOR_s(IJK, M) = HOR_s(IJK, M) + rxn * (HfrefoR_s(M, N)  + &
	                     (calc_ICpoR(T_s(IJK, M), Thigh_s(M,N), Tlow_s(M,N), &
			     Tcom_s(M,N), Ahigh_s(1,M,N), Alow_s(1,M,N)) -  &
			      IC_PsrefoR(M,N)) )* GAS_CONST_cal / MW_s(M,N)
	                   
              END DO 
              IF (UNITS == 'SI') HOR_s(IJK, M) = 4183.925D0*HOR_s(IJK, M)    !in J/kg K
            END DO 
            ENDIF ! for energy_eq
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
                     IF(.not.DMP_LOG)call open_pe_log(ier)
                     WRITE (UNIT_LOG, 1000) L, M 
                     CALL END_LOG 
                     call mfix_exit(myPE)  
                  ENDIF 
               END DO 
            END DO 
	   
         ENDIF 
      END DO 
      
 1000 FORMAT(/1X,70('*')//' From: RRATES',/&
         ' Message: Mass transfer between phases ',I2,' and ',I2,&
         ' (R_tmp) not specified',/1X,70('*')/) 
      RETURN  
      END SUBROUTINE RRATES 
