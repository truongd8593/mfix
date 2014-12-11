!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments: Write reaction rates in units of moles/sec (cgs and SI).  !
!                                                                      !
!  WARNING: Only discrete phase reactions should be specifed here.     !
!  Homogeneous gas phase reactions in DEM simulations must be given    !
!  in the continuum reaction hook (usr_rates.f).                       !
!                                                                      !
!  The call to usr_rates_des is made from inside a particle loop which !
!  is nested inside an IJK loop. Fluid grid calculations independent   !
!  of particle properties can be carried out in des/usr4_des.f to      !
!  reduce redundant calculations.                                      !
!                                                                      !
!  Example: Evaporation                                                !
!                                                                      !
!  mfix.dat input:                                                     !
!``````````````````````````````````````````````````````````````````````!
!    @(DES_RXNS)                                                       !
!      Evap { chem_eq = "Liquid --> Vapor" }                           !
!    @(DES_END)                                                        !
!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!  des/usr_rates_des.f input:                                          !
!``````````````````````````````````````````````````````````````````````!
!    Sa = Pi*(2.0d0*DES_RADIUS(NP))**2    ! Particle surface area      !
!    H2O_xfr = ...  ! An expression for mass transfer coefficient      !
!    Cmg_H2O = ...  ! Molar concentration grad of water vapor          !
!    DES_RATES(EVAP) = Sa * H2O_xfr * Cmg_H2O                          !
!``````````````````````````````````````````````````````````````````````!
!  * Species alias and reaction names given in the data file can be    !
!    used in reference to the reaction index in DES_RATES and a        !
!    species index in gas/solids phase variables.                      !
!                                                                      !
!  * Additional information is provided in section 4.11 of the code    !
!    Readme.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES_DES(NP, pM, IJK, DES_RATES)

      USE compar
      USE constant
      USE des_thermo
      USE des_rxns
      USE discretelement
      USE energy
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE rxns
      USE run
      use functions

! User input: Proximate analysis
      use usr, only: PAC  ! Char
      use usr, only: PAV  ! Volatiles
      use usr, only: PAM  ! Moisture
      use usr, only: PAA  ! Ash

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP  ! Global index of particle
      INTEGER, INTENT(IN) :: pM  ! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP

! Calculated reaction rates. (reacted moles per sec)
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)


! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Particle temperature (K), Mass (kg)
      DOUBLE PRECISION :: Tp, Mp


! Logicals used to 'order' heterogeneous reactions.
! Drying --> Pyrolysis --> Gasification & Combustion
      LOGICAL :: DRIED      ! Initial moisture is gone

! Drying specific variables:
!---------------------------------------------------------------------//
! The amount of energy required to vaporize all the particle's moisture.
      DOUBLE PRECISION :: reqd_Q
! The amount of energy availble to dry the particle.
      DOUBLE PRECISION :: aval_Q
! The amount of energy used to dry the particle.
      DOUBLE PRECISION :: used_Q
! Heat of Vaporization of water at 100 C
      DOUBLE PRECISION, PARAMETER :: Hv_H2O = 2259.95d3 ! (J/kg-H2O)

! Minimum amount of inital species required to stage reactions.
      DOUBLE PRECISION, parameter :: PA_Limiter = 1.0d-4
! Minimum amount of species required to facilitate a reaction.
      DOUBLE PRECISION, parameter :: c_Limiter = 1.0d-6


      INCLUDE '../species.inc'


! Initialization:
!`````````````````````````````````````````````````````````````````````//
! Set the particle temperature and mass.
      Tp = DES_T_s_NEW(NP)
      Mp = PMASS(NP)

! Logical set true if the majority of the intial moisture was driven
! from the solids. Used to suppress pyrolysis. The secondary logical 
! check is included for initially dry coal.
      DRIED = ((PAA * DES_X_s(NP,Moisture)) .LE.                       &
         (PAM * DES_X_s(NP,Ash) * 1.0d-3))  .OR.                       &
         (PAM .LT. PA_Limiter)

! Reaction rates:
!`````````````````````````````````````````````````````````````````````//

! Enact the Boiling Law.
      IF(.NOT.DRIED .AND. Tp >= 373.15d0) THEN
! Calculate the total amount of energy required to vaporize all of the
! moisture contained in the particle. (J/sec)
         reqd_Q = ((Mp/DTSOLID)*DES_X_s(NP,Moisture))*Hv_H2O
         aval_Q = max(0.0d0, Q_Source(NP))
! Take whatever energy is available to vaporize the moisture. (J/sec)
         used_Q = merge(reqd_Q, aval_Q, aval_Q > reqd_Q)
! Updated the energy source term to account for any loss.
         Q_Source(NP) = Q_Source(NP) - used_Q
! Calculate the drying rate based on the usable energy. (kmole/sec)
         DES_RATES(Drying) = used_Q/(MW_s(pM,Moisture)*Hv_H2O)
      ENDIF


! Pyrolysis:  Volatiles --> 
!---------------------------------------------------------------------//
! NOTE: Pyrolysis is suppressed if the majority of moisture has not
! been driven off.
         IF(DRIED) THEN
            IF(Tp > 373.0d0 .AND. DES_X_s(NP, Volatiles) > c_Limiter)  &
               DES_RATES(Pyrolysis) = 3.5930d4 * exp(-4.571d3/Tp) *    &
               (Mp*DES_X_s(NP,Volatiles)) / MW_s(pM,Volatiles)
         ENDIF

      RETURN
      END SUBROUTINE USR_RATES_DES
