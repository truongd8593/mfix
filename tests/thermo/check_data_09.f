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

      USE compar
      USE discretelement
      USE fldvar
      USE funits 
      USE geometry
      USE indices
      USE mfix_pic
      USE param 
      USE param1 
      USE parse
      USE physprop
      USE run
      USE rxns


      IMPLICIT NONE

! Error flag
      INTEGER IER

! loop/variable indices
      INTEGER L, lM, lN, M, N

      TYPE(REACTION_BLOCK), POINTER :: This

      DOUBLE PRECISION netMassTransfer(0:DIM_M)

      LOGICAL WARNED_USR(0:DIM_M)

! Local representation of the number of solids phases.
! TFM or Hybrid lMMAX = MMAX, otherwise lMMAX = 0
      INTEGER lMMAX

! Indicates that a discrete phase model is active.
      LOGICAL foundDPM

! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_RXNS == 0) RETURN

! Initialize the number of solids phases. This forces all heterogeneous
! reactions with a DPM model to use a different reaction block.
      lMMAX = MMAX
      IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_HYBRID) lMMAX = 0
      IF(DISCRETE_ELEMENT .AND. MPPIC) lMMAX = 0

! Flag that a DPM model is in use.
      foundDPM = .FALSE.
      IF(DISCRETE_ELEMENT .OR. MPPIC) foundDPM = .TRUE.

! Initialize flag indicating if the user was already warned.
      WARNED_USR(:) = .FALSE.

! Verify that the species aliases are unique.
      CALL checkDulpicateAliases('CHECK_DATA_09', NMAX(0), &
         SPECIES_ALIAS_g(:),lMMAX, NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:))

! Verify that species aliases in the datafile match those in the 
! species.inc file.
      CALL checkSpeciesInc('CHECK_DATA_09',NMAX(0),SPECIES_ALIAS_g(:), &
         lMMAX, NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), NO_OF_RXNS,       &
         RXN_NAME(:), 'species.inc', foundDPM)

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
         CALL setReaction(This, NMAX(0), SPECIES_ALIAS_g(:),lMMAX,     &
            NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), usrDH(L), usrfDH(L,:))

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
                  rDatabase(M,N) = .TRUE.
               ENDIF
               This%Species(lN)%MW = MW_g(N)
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
                  rDatabase(M,N) = .TRUE.
               ENDIF
               This%Species(lN)%MW = MW_s(M,N)
            ENDIF
         ENDDO

! Verify Mass Balance (Mass of Reactants = Mass of Products)
!---------------------------------------------------------------------//
         IER = 0
         CALL checkMassBalance('CHECK_DATA_09', This, &
            netMassTransfer(:), IER)
         IF(IER /= 0) THEN
            IF(DMP_LOG) CALL WRITE_RXN_SUMMARY(This, &
               SPECIES_ALIAS_g(:), SPECIES_ALIAS_s(:,:))
            CALL MFIX_EXIT(myPE)
         ENDIF

! Determine interphase exchanges
!---------------------------------------------------------------------//
         CALL calcInterphaseTxfr('CHECK_DATA_09', This,   &
            netMassTransfer(:), ENERGY_EQ, SPECIES_EQ(:), &
            SPECIES_ALIAS_g(:), lMMAX, SPECIES_ALIAS_s(:,:))
      ENDDO

! Write a summary of the chemical reactions
!---------------------------------------------------------------------//
      OPEN(678,FILE='POST_Thermo.dat', status='NEW')
      DO L=1, NO_OF_RXNS
         This => Reaction(L)
         CALL WRITE_RXN_SUMMARY(This, SPECIES_ALIAS_g(:), &
            SPECIES_ALIAS_s(:,:))
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

 1007 FORMAT(//1X,70('*')/' From: CHECK_DATA_09',/' Error 1007:',      &
         ' Illogical data returned from setReaction for chemcial',/    &
         ' reaction',1X,A,'.',//' Number of phases identified: ',I3,/  &
         ' Number of species identified: ',I3,//' Fatal error!',       &
         ' Calling MFiX_Exit.',/1X,70('*')/)


 1012 FORMAT(/1X,70('*')/' From: CHECK_DATA_09',/' Error 1012:',       &
         ' Inconsistent user input. Energy equations are NOT being',/  &
         ' solved and a user defined heat of reaction was detected',   &
         ' for chemical',/' reaction ',A,'.')

      END SUBROUTINE CHECK_DATA_09
