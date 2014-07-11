!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_CHEMICAL_RXNS                                    !
!  Author: J.Musser                                   Date: 21-MAR-14  !
!                                                                      !
!  Purpose: Check chemical reactions specifications                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CHEMICAL_RXNS

! User defined reaction names from reaction blocks @(RXNS)
      use parse, only: RXN_NAME, DES_RXN_NAME
! Number of continuum reactions and data object
      use rxns, only: NO_OF_RXNS, REACTION
! Number of discrete reactions and data object
      use des_rxns, only: NO_OF_DES_RXNS, DES_REACTION
! User specifed species aliases used to specify reactions.
      use rxns, only: SPECIES_ALIAS_g, SPECIES_ALIAS_s

! Number of continuum solids
      use physprop, only: SMAX
! Number of discrete solids
      use discretelement, only: DES_MMAX
! Number of spcies comprising each phase
      use physprop, only: NMAX

! 
      use param1, only: UNDEFINED_I

      use rxn_com, only: checkSpeciesInc
      use rxn_com, only: checkDuplicateAliases


      use error_manager

      IMPLICIT NONE

! Error flag
      INTEGER :: IER

! Local representation of the number of solids phases.
      INTEGER :: lMMAX

! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_RXNS == UNDEFINED_I) NO_OF_RXNS = 0
      IF(NO_OF_DES_RXNS == UNDEFINED_I) NO_OF_DES_RXNS = 0

! If there are no chemical reactions, then skip this routine.
      IF(NO_OF_RXNS + NO_OF_DES_RXNS == 0) RETURN

      CALL INIT_ERR_MSG('CHECK_CHEMICAL_RXNS')

! Initialize the number of solids phases.
      lMMAX = SMAX + DES_MMAX

! Verify that the species aliases are unique.
      CALL checkDuplicateAliases(NMAX(0), SPECIES_ALIAS_g(:), &
         lMMAX, NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:))

! Verify that species aliases in the datafile match those in the 
! species.inc file.
      CALL checkSpeciesInc(NMAX(0), SPECIES_ALIAS_g(:), lMMAX,         &
         NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), NO_OF_RXNS, RXN_NAME(:), &
         NO_OF_DES_RXNS, DES_RXN_NAME)


      IF(NO_OF_RXNS > 0) THEN

         CALL CHECK_CHEMICAL_RXNS_COMMON(NO_OF_RXNS, REACTION)

      ENDIF



      CALL FINL_ERR_MSG

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_CHEMICAL_RXNS_COMMON                             !
!  Author: J.Musser                                   Date: 21-MAR-14  !
!                                                                      !
!  Purpose: Check chemical reactions specifications                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CHEMICAL_RXNS_COMMON(NUM_RXNS, RXN_BLOCK)

!      USE compar
      USE discretelement
      USE fldvar
!      USE funits
!      USE geometry
!      USE indices
!      USE mfix_pic
!      USE param
!      USE param1
      USE parse
      USE physprop
      USE run, only: ENERGY_EQ, SPECIES_EQ
      USE rxns

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NUM_RXNS
      TYPE(REACTION_BLOCK), TARGET, ALLOCATABLE :: RXN_BLOCK(:)


! loop/variable indices
      INTEGER L, lM, lN, M, N

      TYPE(REACTION_BLOCK), POINTER :: This

      DOUBLE PRECISION netMassTransfer(0:DIM_M)

      LOGICAL WARNED_USR(0:DIM_M)

! Indicates that a discrete phase model is active.
      LOGICAL foundDPM



      CALL INIT_ERR_MSG('CHECK_CHEMICAL_RXNS_COMMON')


! Allocate reaction blocks.
      allocate(RXN_BLOCK(NUM_RXNS))



! Initialize flag indicating if the user was already warned.
      WARNED_USR(:) = .FALSE.



! Loop over reaction data pulled from data file.
      DO L=1, NUM_RXNS

         This => RXN_BLOCK(L)

! Store the reaction name.
         This%Name = trim(RXN_NAME(L))

! This check should not be necessary. Pre-processing by make_mfix and
! reading the data file (PARSE_RXN) should have already caught any
! issues.
         IF(len_trim(This%Name) == 0) THEN
            WRITE(ERR_MSG, 1100) L
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1100 FORMAT('Error 1100: No reaction name identified for reaction ',  &
         I3,'.',/'This should have been caught during the parsing ',   &
         'routines.')

! Store the chemical equation.
         This%ChemEq = trim(RXN_CHEM_EQ(L))

! Verify that a chemical equation was given in the data file.
         IF(len_trim(This%ChemEq) == 0) THEN
            WRITE(*,1101) trim(This%Name)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: No chemical equation identified for ',       &
         'reaction ',A,'.',/'Please correct the mfix.dat file.')

! Take the data read from the data file and populate the reaction block.
         CALL setReaction(This, NMAX(0), SPECIES_ALIAS_g(:),lMMAX,     &
            NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), usrDH(L), usrfDH(L,:))

! If the energy equations are not being solved and a user provided
! heat of reaction is given, flag error and exit.
         IF(.NOT.ENERGY_EQ .AND. .NOT.This%Calc_DH) THEN
            WRITE(ERR_MSG,1200) trim(This%Name)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('Error 1200: Inconsistent user input. Energy equations ', &
         'are NOT being',/'solved and a user defined heat of reaction',&
         ' was detected',' for chemical',/' reaction ',A,'.',/'Please',&
         ' correct the mfix.dat file.')

! Skip empty reactions.
         IF(This%nSpecies == 0 .AND. This%nPhases == 0) THEN
            CYCLE

! Something went wrong while parsing the reaction. This is a sanity
! check and should never be true.
         ELSEIF((This%nPhases == 0 .AND. This%nSpecies /= 0) .OR.      &
            (This%nPhases /= 0 .AND. This%nSpecies == 0)) THEN

            WRITE(ERR_MSG,1201) trim(This%Name), This%nPhases,         &
               This%nSpecies
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1201 FORMAT('Error 1201: Illogical data returned from setReaction ',  &
         'for chemcial',/'reaction',1X,A,'.',//' Number of phases ',   &
         'identified: ',I3,/'Number of species identified: ',I3,//,    &
         'Please check the mfix.dat file.')

! Verify that the necessary information for each species in the reaction
! was defined.
         CALL checkThermoReqs(This, SPECIES_g, SPECIES_s, rDatabase,   &
            MW_G, MW_S, C_PG0, C_PS0)


! Verify Mass Balance (Mass of Reactants = Mass of Products)
!---------------------------------------------------------------------//
         IER = 0
         CALL checkMassBalance('CHECK_CHEMICAL_RXNS', This, &
            netMassTransfer(:), IER)
         IF(IER /= 0) THEN
            IF(DMP_LOG) CALL WRITE_RXN_SUMMARY(This, &
               SPECIES_ALIAS_g(:), SPECIES_ALIAS_s(:,:))
            CALL MFIX_EXIT(myPE)
         ENDIF

! Determine interphase exchanges
!---------------------------------------------------------------------//
         CALL calcInterphaseTxfr('CHECK_CHEMICAL_RXNS', This,   &
            netMassTransfer(:), ENERGY_EQ, SPECIES_EQ(:), &
            SPECIES_ALIAS_g(:), lMMAX, SPECIES_ALIAS_s(:,:))
      ENDDO

! Write a summary of the chemical reactions
!---------------------------------------------------------------------//
      DO L=1, NUM_RXNS
         This => RXN_BLOCK(L)
         CALL WRITE_RXN_SUMMARY(This, SPECIES_ALIAS_g(:), &
            SPECIES_ALIAS_s(:,:))
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_CHEMICAL_RXNS_COMMON




      END SUBROUTINE CHECK_CHEMICAL_RXNS

