
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                  ***********************************                 !
!                  *  PCCL gasification rates based  *                 !
!                  *  on 10% coal initial moisture.  *                 !
!                  ***********************************                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)


! Gas phase global variables:
!`````````````````````````````````````````````````````````````````````//
! Volume fraction.
      use fldvar, only: EP_g
! Species molecular weights
      use physprop, only: MW_g
! Mixture molecular weight
      use physprop, only: MW_MIX_g
! Pressure.
      use fldvar, only: P_g
! Density.
      use fldvar, only: RO_g
! Temperature.
      use fldvar, only: T_g
! Species mass fractions.
      use fldvar, only: X_g


! Solids phase global variables:
!`````````````````````````````````````````````````````````````````````//
! Particle diameter.
      use fldvar, only: D_p
! Species molecular weight.
      use physprop, only: MW_s
! Bulk density.
      use fldvar, only: ROP_s
! Material density.
      use fldvar, only: RO_S
! Temperature.
      use fldvar, only: T_s
! Species mass fractions.
      use fldvar, only: X_s
! Number of solids phases
      use param, only: DIMENSION_M


! Reaction global variables:
!`````````````````````````````````````````````````````````````````````//
! Number of reactions.
      use rxns, only: NO_OF_RXNS
! Size of Reaction rates arrays.
      use rxns, only: nRR
! Array for storing reaction rates for output
      use rxns, only: ReactionRates
! Reaction names for debugging.
      use parse, only: RXN_NAME

! Global parameters:
!`````````````````````````````````````````````````````````````````````//
! Gas constant (cal/mol.K)
      use constant, only: RGAS => GAS_CONST_cal
! Double percision aliases:
      use param1, only: ZERO
      use param1, only: SMALL_NUMBER
      use param1, only: HALF
      use param1, only: ONE

! Global user variables:
!`````````````````````````````````````````````````````````````````````//
! Sherwood number :: calculated in usr1.f
      use usr, only: N_sh

! User input: Proximate analysis
      use usr, only: PAC  ! Char
      use usr, only: PAV  ! Volatiles
      use usr, only: PAM  ! Moisture
      use usr, only: PAA  ! Ash

! Dummy variable required by usrnlst.inc (UNUSED)
      use usr, only: DUMMY_DP

      implicit none


! Passed arguments:
!`````````````````````````````````````````````````````````````````````//
! Fluid cell index
      INTEGER, INTENT(IN) :: IJK
! Array for storing reaction rates.
      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES


! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Index used for looping over all reactions.
      INTEGER L

! Specific gas consant for Oxygen (cm^3.atm/g.K)
      DOUBLE PRECISION, parameter :: RO2 = 2.564322d0

! Bounded Phase temperatues (K)
      DOUBLE PRECISION, parameter :: MAX_TEMP = 2.5d3
      DOUBLE PRECISION :: xTg   ! Gas

! Diffusion coefficient for O2 in N2. (cm^2/sec)
      DOUBLE PRECISION :: Diff
! Gas pressure:
      DOUBLE PRECISION :: Pg_atm     ! (atm)
      DOUBLE PRECISION :: PG_atmXMW  ! times mixture moleculre weight

! Gas phase molar concentrations: (g-mole/cm^3)
      DOUBLE PRECISION :: c_O2   ! Oxygen
      DOUBLE PRECISION :: c_CO   ! Carbon monoxide
      DOUBLE PRECISION :: c_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: c_CH4  ! Methane
      DOUBLE PRECISION :: c_H2   ! Hydrogen
      DOUBLE PRECISION :: c_H2O  ! Water Vapor
      DOUBLE PRECISION :: c_Tar  ! Tar

! Gas phase partial pressures: (atm)
      DOUBLE PRECISION :: p_O2   ! Oxygen
      DOUBLE PRECISION :: p_CO   ! Carbon monoxide
      DOUBLE PRECISION :: p_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: p_CH4  ! Methane
      DOUBLE PRECISION :: p_H2   ! Hydrogen
      DOUBLE PRECISION :: p_H2O  ! Water Vapor
      DOUBLE PRECISION :: p_Tar  ! Tar

! Gas phase mole fractions: (g-mol/g-mol Mix)
      DOUBLE PRECISION :: y_O2   ! Oxygen
      DOUBLE PRECISION :: y_CO   ! Carbon monoxide
      DOUBLE PRECISION :: y_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: y_CH4  ! Methane
      DOUBLE PRECISION :: y_H2   ! Hydrogen
      DOUBLE PRECISION :: y_H2O  ! Water Vapor
      DOUBLE PRECISION :: y_Tar  ! Tar

! Gas phase reaction limiters: (0,1]
      DOUBLE PRECISION :: r_O2   ! Oxygen
      DOUBLE PRECISION :: r_CO   ! Carbon monoxide
      DOUBLE PRECISION :: r_CO2  ! Carbon dioxide
      DOUBLE PRECISION :: r_CH4  ! Methane
      DOUBLE PRECISION :: r_H2   ! Hydrogen
      DOUBLE PRECISION :: r_H2O  ! Water Vapor
      DOUBLE PRECISION :: r_Tar  ! Tar

! These are local aliases for variables that need converted to CGS.
      DOUBLE PRECISION :: Pg                 ! Gas pressure
      DOUBLE PRECISION :: ROg                ! Gas Density

! Minimum amount of species required to facilitate a reaction.
      DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-6
! Minimum amount of inital species required to stage reactions.
      DOUBLE PRECISION, parameter :: PA_Limiter = 1.0d-4

! Reaction softener parameters.
      DOUBLE PRECISION :: rLimiter(NO_OF_RXNS)
      DOUBLE PRECISION, parameter :: r_pnt = 1.0d-4
      DOUBLE PRECISION, parameter :: r_exp = 3.0d0

      include 'species.inc'
      include 'usrnlst.inc'

! SI units to CGS:          Conversion     !   SI       |    CGS
      Pg   = P_g(IJK)      *  1.0d+1       !   Pa       |   Barye
      ROg  = RO_g(IJK)     *  1.0d-3       !   kg/m^3   |   g/cm^3


! Gas phase quantities:
!---------------------------------------------------------------------//
! Initialize gas phase molar concentrations and partial pressures.
      c_O2  = ZERO;  p_O2  = ZERO;  y_O2  = ZERO;  r_O2  = ONE;
      c_CO  = ZERO;  p_CO  = ZERO;  y_CO  = ZERO;  r_CO  = ONE;
      c_CO2 = ZERO;  p_CO2 = ZERO;  y_CO2 = ZERO;  r_CO2 = ONE;
      c_CH4 = ZERO;  p_CH4 = ZERO;  y_CH4 = ZERO;  r_CH4 = ONE;
      c_H2  = ZERO;  p_H2  = ZERO;  y_H2  = ZERO;  r_H2  = ONE;
      c_H2O = ZERO;  p_H2O = ZERO;  y_H2O = ZERO;  r_H2O = ONE;
      c_Tar = ZERO;  p_Tar = ZERO;  y_Tar = ZERO;  r_Tar = ONE;


! Calculate the bounded gas phase temperature (K)
      xTg  = min(MAX_TEMP, T_g(IJK))
! Gas pressure (atm)
      Pg_atm    = Pg / 101.325d4
! Gas pressure multipled by the molecular weight (atm)
      Pg_atmXmw = Pg_atm * MW_MIX_g(IJK)
! Calculate the diffusion coefficient for O2 in N2. Field, 1967.
      DIFF = 4.26d0 * ((xTg/1.8d3)**1.75d0) / Pg_atm

! Tar
      IF(X_g(IJK,Tar) .GT. c_Limiter) THEN
         c_Tar = ROg * X_g(IJK,Tar) / MW_g(Tar)       ! (g-mole/cm^3)
!        p_Tar = Pg_atmXmw * X_g(IJK,Tar) / MW_g(Tar) ! Not used
         r_Tar = (X_g(IJK,Tar)/(X_g(IJK,Tar) + r_pnt))**r_exp
      ENDIF
! Oxygen, O2
      IF(X_g(IJK,O2) .GT. c_Limiter) THEN
         c_O2 = ROg * X_g(IJK,O2) / MW_g(O2)          ! (g-mole/cm^3)
         p_O2 = Pg_atmXmw * X_g(IJK,O2) / MW_g(O2)    ! (atm)
         r_O2 = (X_g(IJK,O2)/(X_g(IJK,O2) + r_pnt))**r_exp
      ENDIF
! Hydrogen, H2
      IF(X_g(IJK,H2) .GT. c_Limiter) THEN
         c_H2 = ROg * X_g(IJK,H2) / MW_g(H2)          ! (g-mole/cm^3)
         p_H2 = Pg_atmXmw * X_g(IJK,H2) / MW_g(H2)    ! (atm)
         y_H2 = X_g(IJK,H2) * MW_MIX_g(IJK)/MW_g(H2)  ! (mol-H2/mol-Mix)
         r_H2 = (X_g(IJK,H2)/(X_g(IJK,H2) + r_pnt))**r_exp
      ENDIF
! Carbon monoxide, CO
      IF(X_g(IJK,CO) .GT. c_Limiter) THEN
         c_CO = ROg * X_g(IJK,CO) / MW_g(CO)           ! (g-mole/cm^3)
         p_CO = Pg_atmXmw * X_g(IJK,CO) / MW_g(CO)     ! (atm)
         y_CO = X_g(IJK,CO) * MW_MIX_g(IJK)/MW_g(CO)   ! (mol-CO/mol-Mix)
         r_CO = (X_g(IJK,CO)/(X_g(IJK,CO) + r_pnt))**r_exp
      ENDIF

! Water Vapor, H2O
      IF(X_g(IJK,H2O) .GT. c_Limiter) THEN
         c_H2O = ROg * X_g(IJK,H2O) / MW_g(H2O)         ! (g-mole/cm^3)
         p_H2O = Pg_atmXmw * X_g(IJK,H2O) / MW_g(H2O)   ! (atm)
         y_H2O = X_g(IJK,H2O) * MW_MIX_g(IJK)/MW_g(H2O) ! (mol-H2O/mol-Mix)
         r_H2O = (X_g(IJK,H2O)/(X_g(IJK,H2O) + r_pnt))**r_exp
      ENDIF

! Carbon dioxide, CO2
      IF(X_g(IJK,CO2) .GT. c_Limiter) THEN
!        c_CO2 = ROg * X_g(IJK,CO2) / MW_g(CO2)         ! Not used
         p_CO2 = Pg_atmXmw * X_g(IJK,CO2) / MW_g(CO2)   ! (atm)
         y_CO2 = X_g(IJK,CO2) * MW_MIX_g(IJK)/MW_g(CO2) ! (mol-CO2/mol-Mix)
         r_CO2 = (X_g(IJK,CO2)/(X_g(IJK,CO2) + r_pnt))**r_exp
      ENDIF
! Methane, CH4
      IF(X_g(IJK,CH4) .GT. c_Limiter) THEN
         c_CH4 = ROg * X_g(IJK,CH4) / MW_g(CH4)        ! (g-mole/cm^3)
!        p_CH4 = Pg_atmXmw * X_g(IJK,CH4) / MW_g(CH4)  ! Not used
         r_CH4 = (X_g(IJK,CH4)/(X_g(IJK,CH4) + r_pnt))**r_exp
      ENDIF


! Set up reaction limiters.
!---------------------------------------------------------------------//
      rLimiter = ONE

      rLimiter(CO_Combustion)     = min(r_O2,  r_CO)
      rLimiter(CH4_Combustion)    = min(r_O2,  r_CH4)
      rLimiter(H2_Combustion)     = min(r_O2,  r_H2)

!**********************************************************************!
!                                                                      !
!                 Homogeneous gas phase Reaction Rates                 !
!                                                                      !
!**********************************************************************!


! Tar-Cracking: Tar --> 1.65201*CO  + 0.452635*CO2 + ...
!---------------------------------------------------------------------//
! Ref: PCCL
      RATES(TarCracking) =  exp(-6.54d2/xTg) * EP_g(IJK) * c_Tar


! Tar-Combustion: Tar + 29.8O2 --> 20.8CO2 + 24.4H2O + 0.1SO3 + 0.8N2
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(TarCombustion) =  3.8d11 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
         (c_O2**1.5d0) * (c_Tar**0.25d0)


! CO Combustion: 2CO + O2 --> 2CO2
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(CO_Combustion) = 3.98d14 * exp(-2.013d4/xTg) * EP_g(IJK) * &
         (c_O2**0.25d0) * c_CO * (c_H2O**0.5d0)


! CH4 Combustion: CH4 + 2O2 --> CO2 + 2H2O
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)
      RATES(CH4_Combustion) = 6.7d12 * exp(-2.436d4/xTg) * EP_g(IJK) * &
         (c_O2**1.3d0) * (c_CH4**0.2d0)


! H2 Combustion: 2H2 + O2 --> 2H2O
!---------------------------------------------------------------------//
! Ref: Peters (1979)
      RATES(H2_Combustion) = 1.08d16 * exp(-1.51d4/xTg) * EP_g(IJK) *  &
         c_O2 * c_H2

!**********************************************************************!
!                         Unit Conversion                              !
!**********************************************************************!
! The reaction rates were calculated in CGS units (g-mol/sec.cm^3).
! Convert them to SI (kg-mol/sec.m^3).
      DO L=1, NO_OF_RXNS
         RATES(L) = RATES(L) * rLimiter(L) * 1.0d3
      ENDDO


!**********************************************************************!
!                         Save the Rates                               !
!**********************************************************************!
      DO L=1, min(nRR,NO_OF_RXNS)
         ReactionRates(IJK,L) = RATES(L)
      ENDDO

      RETURN
      END SUBROUTINE USR_RATES
