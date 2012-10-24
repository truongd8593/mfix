!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_09                                          C
!  Purpose: Check chemical reactions specifications                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-MAR-98  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_09

      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE rxns
      USE indices
      USE funits 
      USE compar

      USE parse

      IMPLICIT NONE

! Error flag
      LOGICAL ERROR

! loop/variable indices
      INTEGER L, lM, lN, M, N, LL, MM

      INTEGER lS, lE

      DOUBLE PRECISION rSUM, pSUM
      DOUBLE PRECISION MWxStoich

      TYPE(REACTION_BLOCK), POINTER :: This

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

      DOUBLE PRECISION netMassTransfer(0:MMAX)

      LOGICAL lGas, lSolids

      INTEGER toPhase, toPhaseCount, mCount
      INTEGER fromPhase, fromPhaseCount
      INTEGER catPhase, cN

      INTEGER sprCount, sprPhase

      DOUBLE PRECISION lMW

      LOGICAL WARNED_USR(0:MMAX)

      DOUBLE PRECISION, PARAMETER :: massBalanceTol = 1.0d-3

! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_RXNS == 0) RETURN



      IF(CALL_ISAT) THEN
         IF(myPE == PE_IO) THEN
            write(*,1013)
            write(*,1000)
            write(UNIT_LOG,1013)
            write(UNIT_LOG,1000)
         ENDIF
         CALL MFiX_EXIT(myPE)


      ELSEIF(CALL_DI) THEN

         IF(myPE == PE_IO) THEN
            write(*,1014)
            write(*,1000)
            write(UNIT_LOG,1014)
            write(UNIT_LOG,1000)
         ENDIF
         CALL MFiX_EXIT(myPE)
      ENDIF

! Initialize flag indicating if the user was already warned.
      WARNED_USR(:) = .FALSE.

! Verify that the species aliases are unique.
      CALL checkDulpicateAliases()

! Verify that species aliases in the datafile match those in the 
! species.inc file.
      CALL checkSpeciesInc()

! Loop over reaction data pulled from data file.
      DO L=1, NO_OF_RXNS

         This => Reaction(L)

! Store the reaction name.
         This%Name = trim(RXN_NAME(L))
! This check should not be necessary. Pre-processing by make_mfix and
! reading the data file (PARSE_RXN) should have already caught any
! issues.
         IF(len_trim(This%Name) == 0) THEN
            IF(myPE == PE_IO) THEN
               write(*,1001) L
               write(*,1000)
               write(UNIT_LOG,1001) L
               write(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF
! Store the chemical equation.
         This%ChemEq = trim(RXN_CHEM_EQ(L))
! Verify that a chemical equation was given in the data file.
         IF(len_trim(This%ChemEq) == 0) THEN
            IF(myPE == PE_IO) THEN
               write(*,1002) trim(This%Name)
               write(*,1000)
               write(UNIT_LOG,1002) trim(This%Name)
               write(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF

! Take the data read from the data file and populate the reaction block.
         CALL setReaction(This, NMAX(0), SPECIES_ALIAS_g(:), MMAX,     &
            NMAX(1:MMAX), SPECIES_ALIAS_s(:,:), usrDH(L), usrfDH(L,:))

! If the energy equations are not being solved and a user provided
! heat of reaction is given, flag error and exit.
         IF(.NOT.ENERGY_EQ .AND. .NOT.This%Calc_DH) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,1012) trim(This%Name)
               WRITE(*,1000)
               WRITE(UNIT_LOG,1012) trim(This%Name)
               WRITE(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF

! Skip empty reactions.
         IF(This%nSpecies == 0 .AND. This%nPhases == 0) THEN
            CYCLE

! Something went wrong while parsing the reaction. This is a sanity
! check and should never be true.
         ELSEIF((This%nPhases == 0 .AND. This%nSpecies /= 0) .OR. &
            (This%nPhases /= 0 .AND. This%nSpecies == 0)) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,1007) trim(This%Name), This%nPhases, &
                  This%nSpecies
               WRITE(UNIT_LOG,1007) trim(This%Name), This%nPhases, &
                  This%nSpecies
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

! Initialize variables
         rSUM = ZERO
         pSUM = ZERO
         netMassTransfer(:) = ZERO
         sprCount = 0

! Verify that the molecular weights and stoichiometry are consistent and
! determine interphase mass exchanges.
         DO lN = 1, This%nSpecies
            M = This%Species(lN)%pMap
            N = This%Species(lN)%sMap

! Get the molecular weight.
            IF(M==0) THEN
! If the thermochemical database has not be read for this species and
! the molecular weight or specific heat coefficients are undefined,
! then read the database.
               IF(.NOT.rDatabase(M,N) .AND. (MW_g(N) == UNDEFINED .OR. &
                  This%Calc_DH)) THEN
                  IF(.NOT.WARNED_USR(0) .AND. myPE == PE_IO) THEN
                     WRITE(*,1003)
                     WRITE(UNIT_LOG,1003)
                     WARNED_USR(0) = .TRUE.
                  ENDIF
                  WRITE(*,1103) N
                  WRITE(UNIT_LOG,1103) N
                  CALL READ_DATABASE('TFM', 0, N, SPECIES_g(N), MW_g(N))
               ENDIF
               lMW = MW_g(N)
            ELSE
               IF(.NOT.rDatabase(M,N) .AND. (MW_s(M,N) == UNDEFINED    &
                  .OR. This%Calc_DH)) THEN
                  IF(.NOT.WARNED_USR(M) .AND. myPE == PE_IO) THEN
                     WRITE(*,1004) M
                     WRITE(UNIT_LOG,1004) M
                     WARNED_USR(M) = .TRUE.
                  ENDIF
                  WRITE(*,1104) M,N
                  WRITE(UNIT_LOG,1103) M, N
                  CALL READ_DATABASE('TFM', M, N, SPECIES_s(M,N), &
                     MW_s(M,N))
               ENDIF
               lMW = MW_s(M,N)
            ENDIF

! Multiply the molecular weight and stoichiometric coefficient.
            MWxStoich = lMW*This%Species(lN)%Coeff
! Store the molecular weight and MWxStoich
            This%Species(lN)%MWxStoich = MWxStoich
            This%Species(lN)%MW = lMW
! Calculate the net mass transfer for phase M.
!  0 : no interphase mass transfder
! >0 : gains mass from anther phase
! <0 : transfers mass to anther phase
            netMassTransfer(M) = netMassTransfer(M) + MWxStoich
! Calculate mass of reactants and products. Used to ensure mass balance.
            IF(MWxStoich < ZERO) THEN
               rSUM = rSUM - MWxStoich
               IF(M /= 0) THEN
                  sprCount = sprCount + 1
                  sprPhase = M
               ENDIF
            ELSE
               pSUM = pSUM + MWxStoich
            ENDIF
         ENDDO
! Verify that the mass of products equlas reactants: (Mass Balance)
         IF (.NOT.COMPARE(rSUM,pSUM)) THEN 
            IF(myPE == PE_IO) THEN
               WRITE(*,1005) trim(This%Name), rSUM, pSUM
               WRITE(UNIT_LOG,1005) trim(This%Name), rSUM, pSUM
               CALL WRITE_RXN_SUMMARY()
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
! Verify that there is at most one solids phase fule (reactant).
         IF(sprCount > 1) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,1006) trim(This%Name)
               WRITE(UNIT_LOG,1006) trim(This%Name)
            ENDIF
            CALL WRITE_RXN_SUMMARY()
            CALL MFIX_EXIT(myPE)
         ENDIF


! Determine interphase exchanges
!---------------------------------------------------------------------//
! Initialize interphase exchange terms.
         IF(Allocated(This%rPhase)) This%rPhase(:) = ZERO

! If there is only one phase referenced by the reaction, there there 
! should be no interphase mass transfer.
         IF(This%nPhases == 1) THEN
! Interphase mass transfer is set to zero. Small inconsistancies with
! molecular weights can resunt in a non-zero value for homogeneous
! reactions. Presumably, the mass balance check caught any major errors.
            This%rPhase(:) = ZERO
! Void interphase transfer flags.
            DO lN = 1, This%nSpecies
               M = This%Species(lN)%pMap
               This%Species(lN)%mXfr = M
            ENDDO
            This%Classification = "Homogeneous"
! This is a multiphase reaction. 
         ELSE
! Initialize.
            toPhaseCount = 0 
            fromPhaseCount = 0
            DO M = 0, MMAX
! Determine the number of phases with a net mass gain. Record the index
! of the last phase with a net mass gain.
               IF (netMassTransfer(M) > massBalanceTol) THEN 
                  toPhaseCount = toPhaseCount + 1 
                  toPhase = M
! Determine the number of phases with a net mass loss. Record the index
! index of the last phase with a net mass loss.
               ELSEIF(netMassTransfer(M) < -massBalanceTol) THEN 
                  fromPhaseCount = fromPhaseCount + 1 
                  fromPhase = M
               ENDIF 
            ENDDO

! Only one phase has a net mass gain.
            IF(toPhaseCount == 1) THEN 
! Interphase mass transfer flag.
               This%Classification = "Heterogeneous"
               DO M = 0, MMAX 
                  IF(M /= toPhase) THEN
                     IF (toPhase < M) THEN
                        LM = 1 + toPhase + ((M-1)*M)/2
                        This%rPhase(LM) = -netMassTransfer(M)
                     ELSE
                        LM = 1 + M + ((toPhase-1)*toPhase)/2
                        This%rPhase(LM) = netMassTransfer(M)
                     ENDIF

! Verify that if one phase's species equations are solved, that the 
! other phase's species equations are solved.
                     IF(abs(This%rPhase(LM)) > SMALL_NUMBER) THEN
                        IF((SPECIES_EQ(toPhase)             &
                              .AND. .NOT.SPECIES_EQ(M))     &
                           .OR. (.NOT.SPECIES_EQ(toPhase)   &
                              .AND. SPECIES_EQ(M))) THEN
                           IF(myPE == PE_IO) THEN
                              WRITE(*,1008)
                              WRITE(UNIT_LOG,1008)
                              IF(SPECIES_EQ(M)) THEN
                                 WRITE(*,1108) M, 'Solving'
                                 WRITE(*,1108) toPhase, 'Not Solving'
                                 WRITE(UNIT_LOG,1108) M, 'Solving'
                                 WRITE(UNIT_LOG,1108) toPhase, &
                                    'Not Solving'
                              ELSE
                                 WRITE(*,1108) toPhase, 'Solving'
                                 WRITE(*,1108) M, 'Not Solving'
                                 WRITE(UNIT_LOG,1108) M, 'Solving'
                                 WRITE(UNIT_LOG,1108) toPhase, &
                                    'Not Solving'
                              ENDIF
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY()
                           CALL MFiX_EXIT(myPE)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO 
! Set flags for enthalpy transfer associated with mass transfer.
               IF(ENERGY_EQ .AND. This%Calc_DH) THEN
                  DO lN = 1, This%nSpecies
                     M = This%Species(lN)%pMap
! The gas phase is referenced by the reaction. 
                     IF(M == 0) THEN
! The gas phase is the destination phase.
                        IF(toPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the 
! gas phase.
                           mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                           DO MM = 1, MMAX
                              LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phase MM.
                              IF(This%rPhase(LM) > 0) THEN
! Indicate that phase MM receives mass from phase M.
                                 This%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                                 This%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                                 mCount = mCount + 1
                              ENDIF
                           ENDDO
                           IF(mCount /= 1) THEN
                              IF(myPE == PE_IO) THEN
                                 WRITE(*,1009) trim(This%ChemEq)
                                 WRITE(*,1000)
                                 WRITE(UNIT_LOG,1009) trim(This%ChemEq)
                                 WRITE(UNIT_LOG,1000)
                              ENDIF
                              CALL WRITE_RXN_SUMMARY()
                              CALL MFiX_EXIT(myPE)
                           ENDIF
! A solids phase is the destination phase.
                        ELSE
! Since only one phase was detected with a net mass gain and the gas
! phase was detected as a source phase, then all the gas is assigned
! to the destination phase.
                           This%Species(lN)%mXfr = toPhase
! This variable is not used for gas/solids reactions.
                           This%Species(lN)%xXfr = ZERO
                        ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                     ELSE
! Check to see if phase M transfer mass to another solids phase.
                        DO LL = 1, MMAX-1
                           DO MM = LL + 1, MMAX
                              IF(M /= LL .AND. M /= MM) CYCLE
                              LM = LL + 1 + (MM-1)*MM/2
                              IF(This%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                              IF( M == LL .AND. &
                                 This%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                                 This%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                                 This%Species(lN)%xXfr = abs( &
                                    netMassTransfer(MM) / &
                                    netMassTransfer(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                              ELSEIF( M == MM .AND. &
                                 This%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                                 This%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                                 This%Species(lN)%xXfr = abs( &
                                    netMassTransfer(LL) / &
                                    netMassTransfer(MM))
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF ! Gas or Solids phase.
                  ENDDO ! Species Loop
               ENDIF ! Energy Equation
! If there is only one phase with a net mass loss, setup the array for
! interphase mass transfer.
            ELSEIF(fromPhaseCount == 1) THEN
               This%Classification = "Heterogeneous"
               DO M = 0, MMAX 
                  IF (M /= fromPhase) THEN
                     IF(M < fromPhase) THEN
                        LM = 1 + M + ((fromPhase-1)*fromPhase)/2
                        This%rPhase(LM) =  netMassTransfer(M)
                     ELSE
                        LM = 1 + fromPhase + ((M-1)*M)/2
                        This%rPhase(LM) = -netMassTransfer(M)
                     ENDIF

! Verify that if one phase's species equations are solved, that the 
! other phase's species equations are solved.
                     IF(abs(This%rPhase(LM)) > SMALL_NUMBER) THEN
                        IF((SPECIES_EQ(fromPhase)            &
                              .AND. .NOT.SPECIES_EQ(M))      &
                           .OR. (.NOT.SPECIES_EQ(fromPhase)  &
                              .AND. SPECIES_EQ(M))) THEN
                           IF(myPE == PE_IO) THEN
                              WRITE(*,1008)
                              WRITE(UNIT_LOG,1008)
                              IF(SPECIES_EQ(M)) THEN
                                 WRITE(*,1108) M, 'Solving'
                                 WRITE(*,1108) fromPhase, 'Not Solving'
                                 WRITE(UNIT_LOG,1108) M, 'Solving'
                                 WRITE(UNIT_LOG,1108) fromPhase, &
                                    'Not Solving'
                              ELSE
                                 WRITE(*,1108) toPhase, 'Solving'
                                 WRITE(*,1108) M, 'Not Solving'
                                 WRITE(UNIT_LOG,1108) fromPhase, &
                                    'Solving'
                                 WRITE(UNIT_LOG,1108) M, 'Not Solving'
                              ENDIF
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY()
                           CALL MFiX_EXIT(myPE)
                        ENDIF
                     ENDIF
                  ENDIF
               END DO 
! Set flags for enthalpy transfer associated with mass transfer.
               IF(ENERGY_EQ .AND. This%Calc_DH) THEN
                  DO lN = 1, This%nSpecies
                     M = This%Species(lN)%pMap
! Gas/solids reaction: Enthalpy transfer from mass transfer is only
! calculated from gas phase species.
                     IF(M == 0) THEN
! Gas phase is the source phase.
                        IF(fromPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the 
! gas phase.
                           mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                           DO MM = 1, MMAX
                              LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phases MM.
                              IF(This%rPhase(LM) < 0) THEN
! Indicate that phase MM receives mass from phase M.
                                 This%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                                 This%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                                 mCount = mCount + 1
                              ENDIF
                           ENDDO
                           IF(mCount /=1 ) THEN
                              IF(myPE == PE_IO) THEN
                                 WRITE(*,1009) trim(This%ChemEq)
                                 WRITE(*,1000)
                                 WRITE(UNIT_LOG,1009) trim(This%ChemEq)
                                 WRITE(UNIT_LOG,1000)
                              ENDIF
                              CALL WRITE_RXN_SUMMARY()
                              CALL MFiX_EXIT(myPE)
                           ENDIF
                        ELSE
! There can be only one solids phase fuel. Store the phase of the
! solids phase reactant.
                           This%Species(lN)%mXfr = fromPhase
! Mass fraction of transfered material.
! This variable is not currently used for gas/solids reactions.
                           This%Species(lN)%xXfr = ZERO
                        ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                     ELSE
! Check to see if phase M transfer mass to another solids phase.
                        DO LL = 1, MMAX-1
                           DO MM = LL + 1, MMAX
                              IF(M /= LL .AND. M /= MM) CYCLE
                              LM = LL + 1 + (MM-1)*MM/2
                              IF(This%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                              IF( M == LL .AND. &
                                 This%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                                 This%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                                 This%Species(lN)%xXfr = abs( &
                                    netMassTransfer(MM) /     &
                                    netMassTransfer(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                              ELSEIF( M == MM .AND. &
                                 This%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                                 This%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                                 This%Species(lN)%xXfr = abs( &
                                    netMassTransfer(LL) /      &
                                    netMassTransfer(MM))
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDIF ! Gas or Solids phase.
                  ENDDO ! Species Loop
               ENDIF ! Energy Equation

! If there are no phases with a net mass gain/loss, check to see if
! the reaction is turned off.
            ELSEIF(toPhaseCount == 0 .AND. fromPhaseCount == 0) THEN
! If the reaction is active, and there is no interphase mass transfer,
! classify the reaction as catalytic.
               IF(This%nPhases > 0) &
                  This%Classification = "Catalytic"
               This%rPhase(:)  = ZERO
! Set flags for enthalpy transfer associated with mass transfer.
               IF(ENERGY_EQ .AND. This%Calc_DH) THEN

! Identify the catalyst phase.
                  catPhase = -1
                  DO lN= 1, This%nSpecies
                     IF(COMPARE(This%Species(lN)%Coeff,ZERO)) THEN
                        IF(catPhase /= -1) THEN
                           IF(catPhase /= This%Species(lN)%pMap) THEN
                              IF(myPE == PE_IO) THEN
                                 WRITE(*,1009) trim(This%Name)
                                 WRITE(*,1000)
                                 WRITE(UNIT_LOG,1009) trim(This%Name)
                                 WRITE(UNIT_LOG,1000)
                              ENDIF
                              CALL WRITE_RXN_SUMMARY()
                              CALL MFiX_EXIT(myPE)
                           ENDIF
                        ELSE
                           catPhase = This%Species(lN)%pMap
                        ENDIF
                     ENDIF
                  ENDDO
! Verify that a catalyst phase was found.
                  IF(catPhase == -1) THEN
                     IF(myPE == PE_IO) THEN
                        WRITE(*,1010) 'catalyst', trim(This%Name)
                        WRITE(*,1000)
                        WRITE(UNIT_LOG,1010) 'catalyst', trim(This%Name)
                        WRITE(UNIT_LOG,1000)
                     ENDIF
                     CALL WRITE_RXN_SUMMARY()
                     CALL MFiX_EXIT(myPE)
                  ENDIF

! Identify the reactant phase.
                  toPhase = -1
                  DO lN = 1, This%nSpecies
                     IF(.NOT.COMPARE(This%Species(lN)%Coeff,ZERO)) THEN
                        IF(toPhase /= -1) THEN
                           IF(toPhase /= This%Species(lN)%pMap) THEN
                              IF(myPE == PE_IO) THEN
                                 WRITE(*,1009) trim(This%Name)
                                 WRITE(*,1000)
                                 WRITE(UNIT_LOG,1009) trim(This%Name)
                                 WRITE(UNIT_LOG,1000)
                              ENDIF
                              CALL WRITE_RXN_SUMMARY()
                              CALL MFiX_EXIT(myPE)
                           ENDIF
                        ELSE
                           toPhase = This%Species(lN)%pMap
                        ENDIF
                     ENDIF
                  ENDDO
! Verify that a reacting phase was found.
                  IF(toPhase == -1) THEN
                     IF(myPE == PE_IO) THEN
                        WRITE(*,1010) 'reacting', trim(This%Name)
                        WRITE(*,1000)
                        WRITE(UNIT_LOG,1010) 'reacting', trim(This%Name)
                        WRITE(UNIT_LOG,1000)
                     ENDIF
                     CALL WRITE_RXN_SUMMARY()
                     CALL MFiX_EXIT(myPE)
                  ENDIF

! Something when wrong.
                  IF(catPhase == toPhase) THEN
                     IF(myPE == PE_IO) THEN
                        WRITE(*,1011) trim(This%Name)
                        WRITE(*,1000)
                        WRITE(UNIT_LOG,1011) trim(This%Name)
                        WRITE(UNIT_LOG,1000)

                     ENDIF
                     CALL WRITE_RXN_SUMMARY()
                     CALL MFiX_EXIT(myPE)
! Gas/solid cataltyic reaction:
                  ELSEIF(toPhase == 0) THEN
                     DO lN = 1, This%nSpecies
                        IF(This%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                           This%Species(lN)%mXfr = catPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                           This%Species(lN)%xXfr = ZERO
                        ENDIF
                     ENDDO
                  ELSEIF(catPhase == 0) THEN
                     DO lN = 1, This%nSpecies
                        IF(This%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                           This%Species(lN)%mXfr = toPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                           This%Species(lN)%xXfr = ZERO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF ! Energy Equation
            ELSE
! Two or more phases have a net mass loss and two or more phases have
! a net mass gain. Therefore, the interphase mass transfer cannot be
! concluded.
               IF(myPE == PE_IO) THEN
                  CALL WRITE_RXN_SUMMARY()
                  WRITE(*,1009) trim(This%ChemEq)
                  WRITE(*,1000)
                  WRITE(UNIT_LOG,1009) trim(This%ChemEq)
                  WRITE(UNIT_LOG,1000)
               ENDIF
               CALL MFiX_EXIT(myPE)
            ENDIF
         ENDIF
      ENDDO

! Loop over reaction data pulled from data file.
      DO L=1, NO_OF_RXNS
         This => Reaction(L)
         CALL WRITE_RXN_SUMMARY()
      ENDDO


      RETURN

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1001 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1001:',      &
         ' No reaction name identified for reaction number ',I4)

 1002 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1002:',      &
         ' No chemical equation identified for chemcial ',/            &
         ' reaction ',A,'.')

 1003 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/ ' Message 1003:'     &
         ' Molecular weight and/or thermochemical data is undefined',/ &
         ' for a gas phase reactant or product. Thus, the',            &
         ' thermochemical database',/' will be used to gather the',    &
         ' necessary data.',/1X,70('*')/)

 1103 FORMAT(/'  Searching thermochemical databases for gas phase',    &
         ' species ',I2)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/ ' Message 1004:'     &
         ' Molecular weight and/or thermochemical data is undefined',/ &
         ' for a solids phase ',I2,' reactant/product. Thus, the',     &
         ' thermochemical',/' database will be used to gather the',    &
         ' necessary data.',/1X,70('*')/)

 1104 FORMAT(/'  Searching thermochemical databases for solids phase ',&
         I2,', species ',I2)


 1005 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/                      &
         ' Error 1005: Stoichiometry is not consistent with molecular',&
         ' weights',/' for reaction ',A,'.',/' Mass of reactants: ',   &
         F12.4,/' Mass of products:  ',F12.4,/1X,70('*')/)

 1006 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/                      &
         ' Error 1006: More than one solids phase fules was detected.',&
         ' Unable to',/' determine solids/solids heat of reaction',    &
         ' unambiguously for',/' reaction ',A,'.',/1X,70('*')/)

 1007 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1007:',      &
         ' Illogical data returned from setReaction for chemcial',/    &
         ' reaction',1X,A,'.',//' Number of phases identified: ',I3,/  &
         ' Number of species identified: ',I3,//' Fatal error!',       &
         ' Calling MFiX_Exit.',/1X,70('*')/)

 1008 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1008:',      &
         ' A chemical reaction or phase change was detected between',/ &
         ' a phases solving species equations and another phase not',  &
         ' solving',/' species equations.',/)

 1108 FORMAT(' Phase ',I2,': ',A,' species equations.')


 1009 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1009:',      &
         ' Reaction complexity exceeds implementation capabilities.',/ &
         ' Unable to determine unambiguously interphase heat or mass', &
         ' transfer.',//' Reaction: ',A,//' Consider splitting the',   &
         ' chemcial reaction equation into two or more',/' separate',  &
         ' equations. The same reaction rate calculated in usr_rates',/&
         ' can be used for the multiple reactions to ensure mass')


 1010 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1010:',      &
         ' Unable to determine ',A,' phase for catalytic reaction'/    &
         1X,A,'.')

 1011 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1011:',      &
         ' Unable to distinguish catalyst phase from reacting phase',/ &
         ' for catalytic reaction ',A,'.')

 1012 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/' Error 1012:',       &
         ' Inconsistent user input. Energy equations are NOT being',/  &
         ' solved and a user defined heat of reaction was detected',   &
         ' for chemical',/' reaction ',A,'.')

 1013 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/' Error 1013:',       &
         ' Inconsistent user input. ISAT for chemical reactions',/     &
         ' require a different input method thaN non-ISAT reactions.')


 1014 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/' Error 1014:',       &
         ' Inconsistent user input. Direct Integration for chemical',/ &
         ' reactions require a different input method thaN non-ISAT',  &
         ' reactions.')

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkDulpicateAlaises()                              !
!                                                                      !
!  Purpose: Loop through species in all phases and ensure that no two  !
!  entries are the same. ***Warning*** Species aliases that were not   !
!  specified are skipped. Non-specified aliases are treated later in   !
!  parse_mod.f/mapAliases.                                             !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkDulpicateAliases()

      IMPLICIT NONE

      INTEGER lM1, lN1
      INTEGER lM2, lN2

      CHARACTER(len=32) lSA1, lSA2

! Set variables for error messages.
      lM1 = 0
      lM2 = 0

! Compare gas phase aliases.
      DO lN1 = 1, NMAX(0)
         lSA1 = SPECIES_ALIAS_g(lN1)
         IF(len_trim(lSA1) == 0) CYCLE
         DO lN2=lN1+1,NMAX(0)
            lSA2 = SPECIES_ALIAS_g(lN2)
            IF(len_trim(lSA2) == 0) CYCLE
            IF(compareAliases(lSA1, lSA2)) GoTo 100
         ENDDO
! Compare gas and solids phase aliases.
         DO lM2 = 1, MMAX
            DO lN2 = 1, NMAX(lM2)
               lSA2 = SPECIES_ALIAS_s(lM2,lN2)
               IF(len_trim(lSA2) == 0) CYCLE
               IF(compareAliases(lSA1, lSA2)) GoTo 100
            ENDDO
         ENDDO
      ENDDO
! Compare aliaes between solids phases
      DO lM1 = 1, MMAX
         DO lN1 = 1, NMAX(lM1)
            lSA1 = SPECIES_ALIAS_s(lM1,lN1)
            IF(len_trim(lSA1) == 0) CYCLE
! Self phase comparison.
            lM2 = lM1
            DO lN2=lN1+1,NMAX(lM2)
               lSA2 = SPECIES_ALIAS_s(lM2,lN2)
               IF(len_trim(lSA2) == 0) CYCLE
               IF(compareAliases(lSA1, lSA2)) GoTo 100
            ENDDO
! Compare with other phases.
            DO lM2 = lM1+1, MMAX
               DO lN2 = 1, NMAX(lM2)
                  lSA2 = SPECIES_ALIAS_s(lM2,lN2)
                  IF(len_trim(lSA2) == 0) CYCLE
                  IF(compareAliases(lSA1, lSA2)) GoTo 100
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN

 100  IF(myPE == PE_IO) THEN
         WRITE(*,1001)
         WRITE(*,1101) lM1, lN1, lSA1
         WRITE(*,1101) lM2, lN2, lSA2
         WRITE(*,1201)

         WRITE(UNIT_LOG,1001)
         WRITE(UNIT_LOG,1101) lM1, lN1, lSA1
         WRITE(UNIT_LOG,1101) lM2, lN2, lSA2
         WRITE(UNIT_LOG,1201)

         CALL MFIX_EXIT(myPE)
      ENDIF



 1001 FORMAT(/1X,70('*')/' From: CHECK_DATA_09 -->',                   &
         ' checkDulpicateAliases',/' Error 1001: Non-unique species',  &
         ' aliases. Species aliases must be unique',/' so that',       &
         ' chemcial equation entries can be linked to a specific',     &
         ' phase.',//' Please refer to the Readme file for specifying',&
         ' chemcial reactions.'/)

 1101 FORMAT(' Phase: ',I2,', Species: ',I3,' - Alias: ',A)

 1201 FORMAT(1X,70('*')/)


      END SUBROUTINE checkDulpicateAliases

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkSpeciesInc()                                    !
!                                                                      !
!  Purpose: Loop through the species.inc file and verify that the      !
!  match those provided in the datafile.                               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkSpeciesInc()

! Input/Output status.
      INTEGER IOS
! File unit.
      INTEGER, PARAMETER :: FUNIT = 167
! Full path to Burcat and Ruscic database
      CHARACTER(len=256) FILENAME
      CHARACTER(len=128) INPUT
! Loop counter for checking all three locations for data
      INTEGER LC

      INTEGER POS

      INTEGER lIndex
      CHARACTER(len=64) lName

! Length of noncomment string
      INTEGER LINE_LEN
! Integer function which returns COMMENT_INDEX
      INTEGER, EXTERNAL :: SEEK_COMMENT 
! Blank line function
      LOGICAL, EXTERNAL :: BLANK_LINE

! Full path to model directory.
      INCLUDE 'mfix_directory_path.inc'

! Initialize loop counter.
      LC = 0
      DO
         LC = LC + 1
! First check the model directory for the species.inc file.
         IF(LC == 1) THEN
            FILENAME = trim(MFIX_PATH)//'/species.inc'
            OPEN(UNIT=FUNIT,FILE=trim(FILENAME),STATUS='OLD',IOSTAT=IOS)
! Cycle on opening error
            IF(IOS /= 0) CYCLE
            IF(myPE == PE_IO) THEN
               WRITE(*,1000) '/model/species.inc'         ! (screen)
               WRITE(UNIT_LOG,1000) '/model/species.inc'  ! (log file)
            ENDIF
! Check the local run directory.
	        ELSEIF(LC == 2) THEN
            OPEN(UNIT=FUNIT,FILE='species.inc',STATUS='OLD',IOSTAT= IOS)
! Cycle on opening error
	           IF(IOS /= 0) CYCLE
            IF(myPE == PE_IO) THEN
               WRITE(*,1000) 'species.inc'         ! (screen)
               WRITE(UNIT_LOG,1000) 'species.inc'  ! (log file)
            ENDIF
        	ELSE
! No species.inc file was located.
            IF(myPE == PE_IO) THEN
               WRITE(*,1004)
               WRITE(UNIT_LOG,1004)
            ENDIF
            EXIT
            CLOSE(FUNIT)
        	ENDIF

         REWIND(FUNIT)
         READ_LP: DO
            READ(FUNIT,"(A)",IOSTAT=IOS) INPUT
            IF(IOS > 0) THEN
               WRITE(*,1001) trim(adjustl(FILENAME))
               CALL MFiX_EXIT(myPE)
            ELSEIF(IOS<0)THEN
! All entries have been processed.
               CLOSE(FUNIT)
               EXIT
            ENDIF
! Clean up the input.
            LINE_LEN = SEEK_COMMENT(INPUT,LEN(INPUT)) - 1 
            CALL REMOVE_COMMENT (INPUT, LINE_LEN + 1, LEN(INPUT)) 
            CALL MAKE_UPPER_CASE (INPUT, LINE_LEN) 
            CALL REPLACE_TAB (INPUT, LINE_LEN) 
! Skip empty entires.
            IF(LINE_LEN <= 0) CYCLE
            IF(BLANK_LINE(INPUT)) CYCLE

            POS = INDEX(INPUT,"INTEGER, PARAMETER ::")
            IF(POS /= 0) THEN
               INPUT = INPUT((POS + 21):)
            ELSE
               CYCLE
            ENDIF

            POS = INDEX(INPUT,"=")
            IF(POS /= 0) THEN
! Store the chemical equation.
               WRITE(lName,"(A)",IOSTAT=IOS) &
                  trim(adjustl(INPUT(:(POS-1))))
               IF(IOS /= 0) THEN
                  IF(myPE == PE_IO) THEN
                     WRITE(*,1002) 'name', trim(INPUT)
                     WRITE(UNIT_LOG,1002) 'name', trim(INPUT)
                  ENDIF
                  CALL MFiX_EXIT(myPE)
               ENDIF
               READ(INPUT((POS+1):),*,IOSTAT=IOS) lIndex
               IF(IOS /= 0) THEN
                  IF(myPE == PE_IO) THEN
                     WRITE(*,1002) 'index', trim(INPUT)
                     WRITE(UNIT_LOG,1002) 'index', trim(INPUT)
                  ENDIF
                  CALL MFiX_EXIT(myPE)
               ENDIF

! Match against what was provided in the datafile:
! Gas phase species aliases.
               DO lN=1, NMAX(0)
                  IF(compareAliases(SPECIES_ALIAS_g(lN), &
                     lName, lN, lIndex)) CYCLE READ_LP
               ENDDO
! Solids phase species aliases.
               DO lM = 1, MMAX
                  DO lN=1, NMAX(lM)
                     IF(compareAliases(SPECIES_ALIAS_s(lM, lN), &
                        lName, lN, lIndex)) CYCLE READ_LP
                  ENDDO
               ENDDO
! Reaction Names
               DO lN=1, NO_OF_RXNS
                  IF(compareAliases(RXN_NAME(lN), &
                     lName, lN, lIndex)) CYCLE READ_LP
               ENDDO
! No match was made.
               IF(myPE == PE_IO) THEN
                  WRITE(*,1003) trim(lName)
                  WRITE(UNIT_LOG,1002) trim(lName)
               ENDIF
               CALL MFiX_EXIT(myPE)

            ENDIF

         ENDDO READ_LP

         CLOSE(FUNIT)
         EXIT

      ENDDO

 1000 FORMAT(/2X,'Verifying reaction aliases in ',A)

 1001 FORMAT(/1X,70('*'),' From: CHECK_DATA_09 --> checkSpeciesInc',/  &
         ' Error 1001: There was a problem reading file: ',A,/         &
         1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: CHECK_DATA_09 --> checkSpeciesInc',/  &
         ' Error 1002: Unable to obtain alias ',A,' from species.inc', &
         ' file.',//' INPUT: ',A,//1X,70('*')/)

 1003 FORMAT(/1X,70('*')/' From: CHECK_DATA_09 --> checkSpeciesInc',/  &
         ' Error 1003: A match could not be made for an entry in the', &
         ' species.inc',/' file. (',A,')',//' If the reaction names',  &
         ' or species aliases were changed in the data file,',/        &
         ' recomplie the code using make_mfix to correct the issue.',/ &
         1X,70('*')/)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DATA_09 --> checkSpeciesInc',/  &
         ' Warning 1004: Unable to locate original species.inc file.', &
         ' No',/' verification of mfix.dat species aliases or',        &
         ' reaction names can be',/' preformed.',/1X,70('*')/)

      END SUBROUTINE checkSpeciesInc


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: compareAlaises()                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION compareAliases(lS1, lS2, N1, N2)

      CHARACTER(len=*), INTENT(IN) :: lS1, lS2

      INTEGER, OPTIONAL, INTENT(IN) :: N1, N2

      CALL MAKE_UPPER_CASE (lS1, len(lS1)) 
      CALL MAKE_UPPER_CASE (lS2, len(lS2)) 

      compareAliases = .FALSE.
      IF(trim(lS1) == trim(lS2)) compareAliases = .TRUE.

      IF(.NOT.compareAliases) RETURN

      IF(PRESENT(N1) .AND. PRESENT(N2)) THEN
         IF(N1 == N2) THEN
            compareAliases = .TRUE.
         ELSE
            compareAliases = .FALSE.
         ENDIF
      ENDIF

      RETURN
      END FUNCTION compareAliases



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: WRITE_RXN_SUMMARY                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE WRITE_RXN_SUMMARY

      CHARACTER*72, OUTPUT
      CHARACTER*72, full, divided, empty

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

      empty = ''
      CALL WRITE0(empty)

      full = ''
      WRITE(full,2000)

      divided = ''
      WRITE(divided,2005) 

! Lead bar
      CALL WRITE0(full)
! Reaction Nmae
      OUTPUT = ''      
      WRITE(OUTPUT, 2001)trim(This%Name)
      OUTPUT(72:72) = '|'
      CALL WRITE0(OUTPUT)

! Row Divider
      CALL WRITE0(full)

      OUTPUT = ''
      WRITE(OUTPUT, 2002)trim(This%ChemEq(1:54))
      OUTPUT(72:72) = '|'
      CALL WRITE0(OUTPUT)

      CALL WRITE0(full)

      IF(This%nSpecies > 0) THEN

         OUTPUT = ''      
         WRITE(OUTPUT, 2007)trim(This%Classification)
         OUTPUT(72:72) = '|'
         CALL WRITE0(OUTPUT)
! Row Divider
         CALL WRITE0(full)


         WRITE(OUTPUT,2003); CALL WRITE0(OUTPUT)
         WRITE(OUTPUT,2004); CALL WRITE0(OUTPUT)
         CALL WRITE0(divided)
      ENDIF


      DO lN = 1, This%nSpecies

         M = This%Species(lN)%pMap
         N = This%Species(lN)%sMap

         WRITE(OUTPUT,2006)

         IF(M == 0) THEN
            lS = (9-int(len_trim(SPECIES_ALIAS_g(N))/2))
            lE = lS + len_trim(SPECIES_ALIAS_g(N))
            OUTPUT(lS:lE) = trim(SPECIES_ALIAS_g(N))
            WRITE(OUTPUT(32:35),"(A)") 'Gas'
         ELSE
            lS = (9-int(len_trim(SPECIES_ALIAS_s(M,N))/2))
            lE = lS + len_trim(SPECIES_ALIAS_s(M,N))
            OUTPUT(lS:lE) = trim(SPECIES_ALIAS_s(M,N))
            WRITE(OUTPUT(30:36),"(A,I2)") 'Solid',M
         ENDIF
         WRITE(OUTPUT(43:44),"(I2)") N
         WRITE(OUTPUT(51:60),"(F9.4)") This%Species(lN)%MW

         IF(COMPARE(This%Species(lN)%Coeff, ZERO)) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") ZERO
            WRITE(OUTPUT(63:71),"(A)") 'Catalyst'
         ELSEIF(This%Species(lN)%Coeff < ZERO) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") -This%Species(lN)%Coeff
            WRITE(OUTPUT(63:71),"(A)") 'Reactant'
         ELSE
            WRITE(OUTPUT(17:26),"(F9.4)")  This%Species(lN)%Coeff
            WRITE(OUTPUT(63:70),"(A)") 'Product'
         ENDIF
         CALL WRITE0(OUTPUT)
         CALL WRITE0(divided)

      ENDDO

      CALL WRITE0(empty)

      RETURN


 2000 FORMAT(2X,'|',68('-'),'|')

 2001 FORMAT(2X,'| Name: ',A)
 2002 FORMAT(2x,'| Chemical Eq: ',A)

 2003 FORMAT('  | Species  |   Stoich    |         | Species |',       &
              ' Molecular  |          |')

 2004 FORMAT('  |  Alias   |   Coeff.    |  Phase  |  Index  |',       &
              '   Weight   |   Type   |')


 2005 FORMAT(2X,'|',10('-'),'|',13('-'),'|',9('-'),'|',9('-'),'|',     &
             12('-'),'|',10('-'),'|')

 2006 FORMAT(2X,'|',10(' '),'|',13(' '),'|',9(' '),'|',9(' '),'|',     &
             12(' '),'|',10(' '),'|')


 2007 FORMAT(2X,'| Classification: ',A)

      END SUBROUTINE WRITE_RXN_SUMMARY


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: WRITE0                                                  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE WRITE0(LINE)
      CHARACTER(len=*), INTENT(IN) :: LINE

      IF(myPE == PE_IO) THEN
         WRITE(*,*) trim(LINE)
         WRITE(UNIT_LOG,*) trim(LINE)
      ENDIF

      END SUBROUTINE WRITE0

      END SUBROUTINE CHECK_DATA_09
