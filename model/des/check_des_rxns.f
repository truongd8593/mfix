!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_RXNS                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_RXNS

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      Use compar
      USE des_rxns
      Use des_thermo
      Use discretelement
      Use funits  
      Use run
      Use param
      Use param1
      Use parse
      Use physprop
      USE rxns
      USE stiff_chem, only : STIFF_CHEMISTRY

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Dummy loop indices for solids phase and species
      INTEGER L, M, N, lM, lN, llN

! Number of processors used. (DES reactive chemistry is currently limited
! to serial runs!)
      INTEGER CHECK_MPI
! Error Flag
      INTEGER IER

      TYPE(REACTION_BLOCK), POINTER :: This

      DOUBLE PRECISION netMassTransfer(0:DIM_M)

! Logical indicating that an error message has been sent to the user.
      LOGICAL WARNED_USR(0:DIM_M)

      LOGICAL lSpeciesEq(0:DIM_M)


! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_DES_RXNS == 0) RETURN

! Verify the reaction model.

      IF(trim(REACTION_MODEL) == 'VARIABLE_DENSITY') THEN
         RM_VARIABLE_DENSITY = .TRUE.
         RM_SHRINKING_PARTICLE = .FALSE.
      ELSEIF(trim(REACTION_MODEL) == 'SHRINKING_PARTICLE') THEN
         RM_VARIABLE_DENSITY = .FALSE.
         RM_SHRINKING_PARTICLE = .TRUE.
      ELSE
         IF(DMP_LOG) THEN
            WRITE(*,1009) trim(REACTION_MODEL)
            WRITE(UNIT_LOG,1009) trim(REACTION_MODEL)
         ENDIF
         CALL MFiX_EXIT(myPE)
      ENDIF

      ALLOCATE( DES_Reaction(NO_OF_DES_RXNS) )

! Check the number of processors. DES reactive chemistry is currently 
! limited to serial runs.
      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE (*, 9001)
         WRITE (UNIT_LOG, 9001)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Stiff chemistry solver is a TFM reaction model not for DES.
      IF(STIFF_CHEMISTRY) THEN
         IF(DMP_LOG) THEN
            write(*,9003)
            write(*,1000)
            write(UNIT_LOG,9003)
            write(UNIT_LOG,1000)
         ENDIF
         CALL MFiX_EXIT(myPE)
      ENDIF

! Initialize flag indicating if the user was already warned.
      WARNED_USR(:) = .FALSE.

! Verify that the species aliases are unique.
      CALL checkDulpicateAliases('CHECK_DES_RXNS', NMAX(0),   &
         SPECIES_ALIAS_g(:),DES_MMAX, DES_NMAX_s(1:DES_MMAX), &
         DES_SPECIES_ALIAS_s(:,:))

! Verify that species aliases in the datafile match those in the 
! species.inc file.
      CALL checkSpeciesInc('CHECK_DES_RXNS',NMAX(0),SPECIES_ALIAS_g(:),&
         DES_MMAX, DES_NMAX_s(1:DES_MMAX), DES_SPECIES_ALIAS_s(:,:),   &
         NO_OF_DES_RXNS, DES_RXN_NAME(:), 'species.inc', .FALSE.)

! Loop over reaction data pulled from data file.
      DO L=1, NO_OF_DES_RXNS

         This => DES_Reaction(L)

! Store the reaction name.
         This%Name = trim(DES_RXN_NAME(L))
! This check should not be necessary. Pre-processing by make_mfix and
! reading the data file (PARSE_RXN) should have already caught any
! issues.
         IF(len_trim(This%Name) == 0) THEN
            IF(DMP_LOG) THEN
               write(*,1001) L
               write(*,1000)
               write(UNIT_LOG,1001) L
               write(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF
! Store the chemical equation.
         This%ChemEq = trim(DES_RXN_CHEM_EQ(L))
! Verify that a chemical equation was given in the data file.
         IF(len_trim(This%ChemEq) == 0) THEN
            IF(DMP_LOG) THEN
               write(*,1002) trim(This%Name)
               write(*,1000)
               write(UNIT_LOG,1002) trim(This%Name)
               write(UNIT_LOG,1000)
            ENDIF
         CALL MFiX_EXIT(myPE)
         ENDIF

! Take the data read from the data file and populate the reaction block.
         CALL setReaction(This, NMAX(0), SPECIES_ALIAS_g(:), DES_MMAX, &
            DES_NMAX_s(1:DES_MMAX), DES_SPECIES_ALIAS_s(:,:),          &
            DES_usrDH(L), DES_usrfDH(L,:))
! Skip empty reactions.
         IF(This%nSpecies == 0 .AND. This%nPhases == 0) THEN
            CYCLE
! Something went wrong while parsing the reaction. This is a sanity
! check and should never be true.
         ELSEIF((This%nPhases == 0 .AND. This%nSpecies /= 0) .OR. &
            (This%nPhases /= 0 .AND. This%nSpecies == 0)) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1007) trim(This%Name), This%nPhases, &
                  This%nSpecies
               WRITE(UNIT_LOG,1007) trim(This%Name), This%nPhases, &
                  This%nSpecies
            ENDIF
            CALL MFIX_EXIT(myPE)
! Homogeneous gas phase reactions must be given in a TFM reaction block.
         ELSEIF(This%nPhases == 1) THEN
! This is a quick check to make sure all the phases match. If they this
! check fails, then something when wrong while parsing the reaction.
            DO lN = 1, This%nSpecies - 1
               M = This%Species(lN)%pMap
               DO llN = lN + 1,This%nSpecies
                  IF(M /= This%Species(llN)%pMap) THEN
                     IF(DMP_LOG) THEN
                        WRITE(*,1010) trim(This%Name)
                        WRITE(UNIT_LOG,1010) trim(This%Name)
                     ENDIF
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDDO
            ENDDO
            IF(M == 0 .AND. DMP_LOG) THEN
               WRITE(UNIT_LOG,1011) trim(This%Name)
               WRITE(UNIT_LOG,1000)
               WRITE(UNIT_LOG,1011) trim(This%Name)
               WRITE(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF

! If the energy equations are not being solved and a user provided
! heat of reaction is given, flag error and exit.
         IF(.NOT.This%Calc_DH) THEN
            IF(.NOT.ENERGY_EQ .OR. .NOT.DES_ENERGY_EQ) THEN
               IF(DMP_LOG) THEN
                  WRITE(*,1008) trim(This%Name)
                  WRITE(*,1000)
                  WRITE(UNIT_LOG,1008) trim(This%Name)
                  WRITE(UNIT_LOG,1000)
               ENDIF
               CALL MFiX_EXIT(myPE)
            ENDIF
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
                  IF(.NOT.WARNED_USR(0) .AND. DMP_LOG) THEN
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
               IF(.NOT.DES_rDatabase(M,N) .AND. &
                  (DES_MW_s(M,N) == UNDEFINED .OR. This%Calc_DH)) THEN
                  IF(.NOT.WARNED_USR(M) .AND. DMP_LOG) THEN
                     WRITE(*,1004) M
                     WRITE(UNIT_LOG,1004) M
                     WARNED_USR(M) = .TRUE.
                  ENDIF
                  WRITE(*,1104) M,N
                  WRITE(UNIT_LOG,1103) M, N
                  CALL READ_DATABASE('DEM', M, N, DES_SPECIES_s(M,N),  &
                     DES_MW_s(M,N))
                  DES_rDatabase(M,N) = .TRUE.
               ENDIF
               This%Species(lN)%MW = DES_MW_s(M,N)
            ENDIF
         ENDDO

! Verify Mass Balance (Mass of Reactants = Mass of Products)
!---------------------------------------------------------------------//
         IER = 0
         CALL checkMassBalance('CHECK_DES_RXNS', This, &
            netMassTransfer(:), IER)
         IF(IER /= 0) THEN
            IF(DMP_LOG) CALL WRITE_RXN_SUMMARY(This, &
               SPECIES_ALIAS_g(:), DES_SPECIES_ALIAS_s(:,:))
            CALL MFIX_EXIT(myPE)
         ENDIF

! Determine interphase exchanges
!---------------------------------------------------------------------//
! Construct a temp array with logicals for species equations.
         lSpeciesEq(0) = SPECIES_EQ(0)
         lSpeciesEq(1:DIM_M) = DES_SPECIES_EQ(1:DIM_M)

         CALL calcInterphaseTxfr('CHECK_DES_RXNS', This,  &
            netMassTransfer(:), ENERGY_EQ, lSpeciesEq(:), &
            SPECIES_ALIAS_g(:), DES_MMAX, DES_SPECIES_ALIAS_s(:,:))
      ENDDO

! Write a summary of the chemical reactions
!---------------------------------------------------------------------//
      DO L=1, NO_OF_DES_RXNS
         This => DES_Reaction(L)
         CALL WRITE_RXN_SUMMARY(This, SPECIES_ALIAS_g(:), &
            DES_SPECIES_ALIAS_s(:,:))
      ENDDO

      RETURN

! Error Messages
!---------------------------------------------------------------------//

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1001 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1001:',     &
         ' No reaction name identified for reaction number ',I4)

 1002 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1002:',     &
         ' No chemical equation identified for chemcial ',/            &
         ' reaction ',A,'.')

 1003 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/ ' Message 1003:'    &
         ' Molecular weight and/or thermochemical data is undefined',/ &
         ' for a gas phase reactant or product. Thus, the',            &
         ' thermochemical database',/' will be used to gather the',    &
         ' necessary data.',/1X,70('*')/)

 1103 FORMAT(/'  Searching thermochemical databases for gas phase',    &
         ' species ',I2)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/ ' Message 1004:'    &
         ' Molecular weight and/or thermochemical data is undefined',/ &
         ' for a discrete solids phase ',I2,' reactant/product. Thus,',&
         ' the',/' thermochemical database will be used to gather the',&
         ' necessary data.',/1X,70('*')/)

 1104 FORMAT(/'  Searching thermochemical databases for discrete',     &
         ' solids phase ',I2,', species ',I2)


 1005 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1005:',      &
         ' Stoichiometry is not consistent with molecular weights',/   &
         ' for reaction ',A,'.',/' Mass of reactants: ',F12.4,/        &
         ' Mass of products:  ',F12.4,/1X,70('*')/)

 1006 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1006:',      &
         ' More than one solids phase fules was detected. Unable to'/, &
         ' determine solids/solids heat of reaction unambiguously for',&
         /' reaction ',A,'.',/1X,70('*')/)


 1007 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1007:',     &
         ' Illogical data returned from setReaction for chemcial',/    &
         ' reaction',1X,A,'.',//' Number of phases identified: ',I3,/  &
         ' Number of species identified: ',I3,//' Fatal error!',       &
         ' Calling MFiX_Exit.',/1X,70('*')/)


 1008 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1008:',      &
         ' Inconsistent user input. Energy equations are NOT being',/  &
         ' solved and a user defined heat of reaction was detected',   &
         ' for chemical',/' reaction ',A,'.')

 1009 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1009:',     &
         ' Unknown reaction model: ',A,//' Available options:',/       &
         3X,'VARIABLE_DENSITY - Constant particle diameter (default)',/&
         3X,'SHRINKING_PARTICLE - Constant particle density.',/        &
         1X,70('*')//)

 1010 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1010:',     &
         ' Illogical data returned from setReaction for chemcial',/    &
         ' reaction',1X,A,'.',//' A single phase was reported but two',&
         ' or more species are associated',/' with different phases!'//&
         ' Fatal error! Calling MFiX_Exit.',/1X,70('*')/)

 1011 FORMAT(//1X,70('*')/' From: CHECK_DES_RXNS',/' Error 1011:',     &
         ' A homogeneous gas phase reaction was detected in the DEM',/ &
         ' reaction block. Homogeneous gas phase reactions must be',   &
         ' specified in',/' the continuum reaction block [@(RXNS)...@',&
         '(END)].',/' Reaction: ',A,//' Only gas-particle reactions',  &
         ' should appear in the DEM reaction block as',/               &
         ' calculations are preformed on a per-particle basis.')

 9001 FORMAT(/1X,70('*')/ ' From: CHECK_DES_RXNS',/' Error 9001:',     &
         ' DES reactive chemistry modules are limited to serail runs.',&
         /' Set nodesi, nodesj, and nodesk to 1 in the mfix.dat file.',&
         /1X,70('*')/)

 9003 FORMAT(/1X,70('*')/' From: CHECK_DES_RXNS',/' Error 9003:',      &
         ' Inconsistent user input. STIFF_CHEMISTRY for chemical',/    &
         ' reactions is not defined for DES.')

      END SUBROUTINE CHECK_DES_RXNS
