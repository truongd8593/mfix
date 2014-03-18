      MODULE RXN_COM
 
      Use param
      Use param1
      USE compar
      Use funits
 
! The following data types are used to group chemical reaction data.
!-----------------------------------------------------------------------

! Species belong to PHASE_ associated with a particular reaction.
      TYPE SPECIES_
! A link between the reacting species' arbitrary index and the 
! associated phase index in MFiX.
         INTEGER pMap
! A link between the reacting species' arbitrary index and the 
! associated species index in MFiX.
         INTEGER sMap
! Stoichiometric coefficient of the species from chemical equation.
         DOUBLE PRECISION Coeff
! Molecular weight
         DOUBLE PRECISION MW
! Fractional mass transfer
         DOUBLE PRECISION xXfr
! Index indicating enthalpy transfer associated with mass transfer.
         INTEGER mXfr
! Molecular weight of speices multipling the stoichiometric coefficient
         DOUBLE PRECISION MWxStoich
      END TYPE SPECIES_

! Grouping of reaction information.
      TYPE REACTION_BLOCK
! Name of reaction construct from data file.
         CHARACTER*32 Name
! User defined chemical equation from data file.
         CHARACTER*512 ChemEq
! Reaction classification: Homogeneous, Heterogeneous, Catalytic.
         CHARACTER*16 Classification
! Indicates if the automated heat of reaction is to be calculated (T) or
! if the user has supplied a heat of reaction (F).
         LOGICAL Calc_DH
! Number of phases associated with the reaction.
         INTEGER nPhases
! Number of species associated with the reaction.
         INTEGER nSpecies
! User-specified heat of reaction split among phases by fracDH
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: HoR
! Interphase mass transfer.
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rPhase
! Reactant/Product information
         TYPE(SPECIES_), DIMENSION(:), ALLOCATABLE :: Species

      END TYPE REACTION_BLOCK

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
      SUBROUTINE checkDulpicateAliases(CALLER, lNg, lSAg, &
         lMMx, lNs, lSAs)

      IMPLICIT NONE
! Calling routine.
      CHARACTER(len=*), INTENT(IN) :: CALLER
! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs

! Loop indicles.
      INTEGER lM1, lN1  ! Phase, Species
      INTEGER lM2, lN2  ! Phase, Species

      CHARACTER(len=32) lSA1, lSA2

! Set variables for error messages.
      lM1 = 0
      lM2 = 0

! Compare gas phase aliases.
      DO lN1 = 1, lNg
         lSA1 = lSAg(lN1)
         IF(len_trim(lSA1) == 0) CYCLE
         DO lN2=lN1+1,lNg
            lSA2 = lSAg(lN2)
            IF(len_trim(lSA2) == 0) CYCLE
            IF(compareAliases(lSA1, lSA2)) GoTo 100
         ENDDO
! Compare gas and solids phase aliases.
         DO lM2 = 1, lMMx
            DO lN2 = 1, lNs(lM2)
               lSA2 = lSAs(lM2,lN2)
               IF(len_trim(lSA2) == 0) CYCLE
               IF(compareAliases(lSA1, lSA2)) GoTo 100
            ENDDO
         ENDDO
      ENDDO
! Compare aliases between solids phases
      DO lM1 = 1, lMMx
         DO lN1 = 1, lNs(lM1)
            lSA1 = lSAs(lM1,lN1)
            IF(len_trim(lSA1) == 0) CYCLE
! Self phase comparison.
            lM2 = lM1
            DO lN2=lN1+1,lNs(lM2)
               lSA2 = lSAs(lM2,lN2)
               IF(len_trim(lSA2) == 0) CYCLE
               IF(compareAliases(lSA1, lSA2)) GoTo 100
            ENDDO
! Compare with other phases.
            DO lM2 = lM1+1, lMMx
               DO lN2 = 1, lNs(lM2)
                  lSA2 = lSAs(lM2,lN2)
                  IF(len_trim(lSA2) == 0) CYCLE
                  IF(compareAliases(lSA1, lSA2)) GoTo 100
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN

 100  IF(DMP_LOG) THEN
         WRITE(*,1001) trim(CALLER)
         WRITE(*,1101) lM1, lN1, lSA1
         WRITE(*,1101) lM2, lN2, lSA2
         WRITE(*,1201)

         WRITE(UNIT_LOG,1001) trim(CALLER)
         WRITE(UNIT_LOG,1101) lM1, lN1, lSA1
         WRITE(UNIT_LOG,1101) lM2, lN2, lSA2
         WRITE(UNIT_LOG,1201)

         CALL MFIX_EXIT(myPE)
      ENDIF


 1001 FORMAT(/1X,70('*')/' From: ',A,' --> RXN_COM -->',               &
         ' checkDulpicateAliases',/' Error 1001: Non-unique species',  &
         ' aliases. Species aliases must be unique',/' so that',       &
         ' chemical equation entries can be linked to a specific',     &
         ' phase.',//' Please refer to the Readme file for specifying',&
         ' chemical reactions.'/)

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
      SUBROUTINE checkSpeciesInc(CALLER, lNg, lSAg, lMMx, lNs, lSAs, &
         lNRxn, lRNames, lFile, lfDPM)

      IMPLICIT NONE

! Calling routine.
      CHARACTER(len=*), INTENT(IN) :: CALLER
! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! Number of reactions
      INTEGER, INTENT(IN) :: lNRxn
! Reaction Names (aliases)
      CHARACTER(len=32), DIMENSION(DIMENSION_RXN), INTENT(IN) :: lRNames

! Generated species include file.
      CHARACTER(len=*), INTENT(IN) :: lFile
! Flag when a discrete phase model is active.
      LOGICAL, INTENT(IN) :: lfDPM

! Input/Output status.
      INTEGER IOS
! File unit.
      INTEGER, PARAMETER :: FUNIT = 167
! Full path to Burcat and Ruscic database
      CHARACTER(len=256) FILENAME
      CHARACTER(len=128) INPUT
! Loop counters
      INTEGER LC  ! file locations
      INTEGER lM, lN  ! Phase, Species/Reactions

! String position
      INTEGER POS

      INTEGER lIndex
! Loop inde

      CHARACTER(len=64) lName

! Flag to give DPM message
      LOGICAL :: lDPM_Flag
! Length of noncomment string
      INTEGER LINE_LEN
! Integer function which returns COMMENT_INDEX
      INTEGER, EXTERNAL :: SEEK_COMMENT 
! Blank line function
      LOGICAL, EXTERNAL :: BLANK_LINE

! Full path to model directory.
      INCLUDE 'mfix_directory_path.inc'

! Initialize
      lDPM_Flag = .FALSE.
      LC = 0  ! Loop counter.

      DO
         LC = LC + 1
! Check the local run directory.
         IF(LC == 2) THEN
            FILENAME = trim(MFIX_PATH)//'/'//trim(lFile)
            OPEN(UNIT=FUNIT,FILE=trim(FILENAME),STATUS='OLD',IOSTAT=IOS)
! Cycle on opening error
            IF(IOS /= 0) CYCLE
            IF(DMP_LOG) THEN
               WRITE(*,1000) '/model/'//trim(lFile)         ! (screen)
               WRITE(UNIT_LOG,1000) '/model/'//trim(lFile)  ! (log file)
            ENDIF

!Check the model directory for the species.inc file.
	        ELSEIF(LC == 1) THEN
            FILENAME = trim(lFile)
            OPEN(UNIT=FUNIT,FILE=trim(FILENAME),STATUS='OLD',IOSTAT=IOS)
! Cycle on opening error
	           IF(IOS /= 0) CYCLE
            IF(DMP_LOG) THEN
               WRITE(*,1000) trim(lFile)         ! (screen)
               WRITE(UNIT_LOG,1000) trim(lFile)  ! (log file)
            ENDIF
        	ELSE
! No species.inc file was located.
            IF(DMP_LOG) THEN
               WRITE(*,1004) trim(CALLER)
               WRITE(UNIT_LOG,1004) trim(CALLER)
            ENDIF
            EXIT
            CLOSE(FUNIT)
        	ENDIF

         REWIND(FUNIT)
         READ_LP: DO
            READ(FUNIT,"(A)",IOSTAT=IOS) INPUT
            IF(IOS > 0) THEN
               WRITE(*,1001) trim(CALLER), trim(adjustl(FILENAME))
               CALL MFiX_EXIT(myPE)
            ELSEIF(IOS<0)THEN
! All entries have been processed.
               CLOSE(FUNIT)
! Give DPM Message if needed.
               IF(lDPM_Flag .AND. DMP_LOG) WRITE(*,1005) trim(CALLER)
               EXIT
            ENDIF
! Clean up the input.
            LINE_LEN = SEEK_COMMENT(INPUT,LEN(INPUT)) - 1 
            CALL REMOVE_COMMENT(INPUT, LINE_LEN + 1, LEN(INPUT)) 
            CALL MAKE_UPPER_CASE(INPUT, LINE_LEN) 
            CALL REPLACE_TAB(INPUT, LINE_LEN) 
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
                  IF(DMP_LOG) THEN
                     WRITE(*,1002) trim(CALLER), 'name', trim(INPUT)
                     WRITE(UNIT_LOG,1002) trim(CALLER), 'name', &
                        trim(INPUT)
                  ENDIF
                  CALL MFiX_EXIT(myPE)
               ENDIF
               READ(INPUT((POS+1):),*,IOSTAT=IOS) lIndex
               IF(IOS /= 0) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1002) trim(CALLER), 'index', trim(INPUT)
                     WRITE(UNIT_LOG,1002) trim(CALLER), 'index', &
                        trim(INPUT)
                  ENDIF
                  CALL MFiX_EXIT(myPE)
               ENDIF

! Match against what was provided in the datafile:
! Gas phase species aliases.
               DO lN=1, lNg
                  IF(compareAliases(lSAg(lN), lName, lN, lIndex)) &
                     CYCLE READ_LP
               ENDDO
! Solids phase species aliases.
               DO lM = 1, lMMx
                  DO lN=1, lNs(lM)
                     IF(compareAliases(lSAs(lM,lN), lName, lN, lIndex))&
                        CYCLE READ_LP
                  ENDDO
               ENDDO
! Reaction Names
               DO lN=1, lNRxn
                  IF(compareAliases(lRNames(lN), lName, lN, lIndex)) &
                     CYCLE READ_LP
               ENDDO
! No match was made.
               IF(lfDPM .AND. trim(CALLER)=='CHECK_DATA_09') THEN
! If this routine was invoked by CHECK_DATA_09, ignore the error as
! the species could be associated with a discrete solids phase.
                  lDPM_Flag = .TRUE.
! This is a TFM only run. Since a match could not be made, flag this
! as an error and exit.
               ELSE
                  IF(DMP_LOG) THEN
                     WRITE(*,1003)trim(CALLER), trim(lName)
                     IF(lfDPM) WRITE(*,1103)
                     WRITE(*,1203)
                     WRITE(UNIT_LOG,1003) trim(CALLER), trim(lName)
                     IF(lfDPM) WRITE(UNIT_LOG,1103)
                     WRITE(UNIT_LOG,1203)
                  ENDIF
                  CALL MFiX_EXIT(myPE)
               ENDIF

            ENDIF

         ENDDO READ_LP

         CLOSE(FUNIT)
         EXIT

      ENDDO

 1000 FORMAT(/2X,'Verifying reaction aliases in ',A)

 1001 FORMAT(/1X,70('*'),/' From: ',A,' --> RXN_COM --> checkSpeciesInc'&
         ,/' Error 1001: There was a problem reading file: ',A,/       &
         1X,70('*')/)

 1002 FORMAT(/1X,70('*'),/' From: ',A,' --> RXN_COM --> checkSpeciesInc'&
         ,/' Error 1002: Unable to obtain alias ',A,' from species.inc'&
         ,' file.',//' INPUT: ',A,//1X,70('*')/)

 1003 FORMAT(/1X,70('*'),/' From: ',A,' --> RXN_COM --> checkSpeciesInc'&
         ,/' Error 1003: A match could not be made for an entry in the'&
         ,' species.inc',/' file. (',A,')',//' If the reaction names', &
         ' or species aliases were changed in the data file,',/        &
         ' recomplie the code using make_mfix to correct the issue.')

 1103 FORMAT(/' Note: A discrete phase model is being used',           &
         ' (DEM/MPPIC/Hybrid). Verify',/' that all gas/DPM solids',    &
         ' heterogeneous reactions are specified within',/' a DPM ',   &
         ' reaction block [@(DPM_RXNS)...@(END)].')

 1203 FORMAT(1X,70('*'))

 1004 FORMAT(/1X,70('*'),/' From: ',A,' --> RXN_COM --> checkSpeciesInc'&
         ,/' Warning 1004: Unable to locate original species.inc file.'&
         ,' No',/' verification of mfix.dat species aliases or',       &
         ' reaction names can be',/' preformed.',/1X,70('*')/)

 1005 FORMAT(/1X,70('*'),/' From: ',A,' --> RXN_COM --> checkSpeciesInc'&
         ,/' Message 1005: One or more species in the species.inc',    &
         ' file were not',/' matched to any gas or continuous solids', &
         ' phase species. Error detection',/' is being deferred to',   &
         ' check_des_rxns.f.',/1X,70('*')/)

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

      IMPLICIT NONE

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
      SUBROUTINE WRITE_RXN_SUMMARY(RxN, lSAg, lSAs)

      IMPLICIT NONE

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(IN) :: RxN

! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs


      CHARACTER*72, OUTPUT
      CHARACTER*72, full, divided, empty

      CHARACTER*32 lSP

      INTEGER lN, M, N
      INTEGER lS, lE

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

      empty = ''
      CALL WRITE_RS0(empty)

      full = ''
      WRITE(full,2000)

      divided = ''
      WRITE(divided,2005) 

! Lead bar
      CALL WRITE_RS0(full)
! Reaction Nmae
      OUTPUT = ''      
      WRITE(OUTPUT, 2001)trim(RxN%Name)
      OUTPUT(72:72) = '|'
      CALL WRITE_RS0(OUTPUT)

! Row Divider
      CALL WRITE_RS0(full)

      OUTPUT = ''
      WRITE(OUTPUT, 2002)trim(RxN%ChemEq(1:54))
      OUTPUT(72:72) = '|'
      CALL WRITE_RS0(OUTPUT)

      CALL WRITE_RS0(full)

      IF(RxN%nSpecies > 0) THEN

         OUTPUT = ''      
         WRITE(OUTPUT, 2007)trim(RxN%Classification)
         OUTPUT(72:72) = '|'
         CALL WRITE_RS0(OUTPUT)
! Row Divider
         CALL WRITE_RS0(full)

         WRITE(OUTPUT,2003); CALL WRITE_RS0(OUTPUT)
         WRITE(OUTPUT,2004); CALL WRITE_RS0(OUTPUT)
         CALL WRITE_RS0(divided)
      ENDIF


      DO lN = 1, RxN%nSpecies

         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap

         WRITE(OUTPUT,2006)

         IF(M == 0) THEN
            IF(len_trim(lSAg(N)) > 8) THEN
               lSP = lSAg(N)
               OUTPUT(5:13) = lSP(1:8)
            ELSE
              lS = (9-int(len_trim(lSAg(N))/2))
              lE = lS + len_trim(lSAg(N))
               OUTPUT(lS:lE) = trim(lSAg(N))
            ENDIF
            WRITE(OUTPUT(32:35),"(A)") 'Gas'
         ELSE
            IF(len_trim(lSAs(M,N)) > 8) THEN
               lSP = lSAs(M,N)
               OUTPUT(5:13) = lSP(1:8)
            ELSE
               lS = (9-int(len_trim(lSAs(M,N))/2))
               lE = lS + len_trim(lSAs(M,N))
               OUTPUT(lS:lE) = trim(lSAs(M,N))
            ENDIF
            WRITE(OUTPUT(30:36),"(A,I2)") 'Solid',M
         ENDIF
         WRITE(OUTPUT(43:44),"(I2)") N
         WRITE(OUTPUT(51:60),"(F9.4)") RxN%Species(lN)%MW

         IF(COMPARE(RxN%Species(lN)%Coeff, ZERO)) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") ZERO
            WRITE(OUTPUT(63:71),"(A)") 'Catalyst'
         ELSEIF(RxN%Species(lN)%Coeff < ZERO) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") -RxN%Species(lN)%Coeff
            WRITE(OUTPUT(63:71),"(A)") 'Reactant'
         ELSE
            WRITE(OUTPUT(17:26),"(F9.4)")  RxN%Species(lN)%Coeff
            WRITE(OUTPUT(63:70),"(A)") 'Product'
         ENDIF
         CALL WRITE_RS0(OUTPUT)
         CALL WRITE_RS0(divided)

      ENDDO

      CALL WRITE_RS0(empty)

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

      contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: WRITE_RS0                                               !
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
      SUBROUTINE WRITE_RS0(LINE)

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: LINE

      IF(DMP_LOG) THEN
         WRITE(*,*) trim(LINE)
         WRITE(UNIT_LOG,*) trim(LINE)
      ENDIF

      END SUBROUTINE WRITE_RS0


      END SUBROUTINE WRITE_RXN_SUMMARY



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: checkMassBalance                                        !
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
      SUBROUTINE checkMassBalance(CALLER, RxN, lnMT, IER)

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: CALLER

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN

      DOUBLE PRECISION, INTENT(OUT) :: lnMT(0:DIM_M)
      INTEGER, INTENT(OUT) :: IER

      INTEGER M, N, lN ! Phase, Species
      DOUBLE PRECISION rSUM, pSUM
      DOUBLE PRECISION MWxStoich

      INTEGER sprCount, sprIndex

      DOUBLE PRECISION, PARAMETER :: massBalanceTol = 1.0d-3

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

! Initialize variables
      IER = 0
      rSUM = ZERO
      pSUM = ZERO
      lnMT(:) = ZERO
      sprCount = 0

! Verify that the molecular weights and stoichiometry are consistent and
! determine interphase mass exchanges.
      DO lN = 1, RxN%nSpecies
         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap

! Multiply the molecular weight and stoichiometric coefficient.
         MWxStoich = RxN%Species(lN)%MW * RxN%Species(lN)%Coeff
         RxN%Species(lN)%MWxStoich = MWxStoich
! Calculate the net mass transfer for phase M.
!  0 : no interphase mass transfder
! >0 : gains mass from anther phase
! <0 : transfers mass to anther phase
         lnMT(M) = lnMT(M) + MWxStoich
! Calculate mass of reactants and products. Used to ensure mass balance.
         IF(MWxStoich < ZERO) THEN
            rSUM = rSUM - MWxStoich
            IF(M /= 0) THEN
               sprCount = sprCount + 1
               IF(sprCount == 1) THEN
                  sprIndex = M
! Verify that there is at most one solids phase fule (reactant).
               ELSEIF( M /= sprIndex) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1002) trim(CALLER), trim(RxN%Name)
                     WRITE(UNIT_LOG,1002) trim(CALLER), trim(RxN%Name)
                     IER = 1
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            pSUM = pSUM + MWxStoich
         ENDIF
      ENDDO
! Verify that the mass of products equlas reactants: (Mass Balance)
      IF (.NOT.COMPARE(rSUM,pSUM)) THEN 
         IF(DMP_LOG) THEN
            WRITE(*,1001) trim(CALLER), trim(RxN%Name), rSUM, pSUM
            WRITE(UNIT_LOG,1001) trim(CALLER), trim(RxN%Name), rSUM,pSUM
            IER = 1
         ENDIF
      ENDIF

      RETURN

! Error Messages
!---------------------------------------------------------------------//

 1001 FORMAT(/1X,70('*')/' From: ',A,' --> RXN_COM -->',               &
         ' checkMassBalance',/' Error 1001: Stoichiometry is not',     &
         ' consistent with molecular weights',/' for reaction ',A,'.',/&
         ' Mass of reactants: ',F12.4,/' Mass of products:  ',F12.4,/  &
         1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: ',A,' --> RXN_COM -->',               &
         ' checkMassBalance',/' Error 1002: More than one solids',     &
         ' phase fules was detected. Unable to',/' determine solids/', &
         'solids heat of reaction unambiguously for',/' reaction ',A,  &
         '.',/1X,70('*')/)

      END SUBROUTINE checkMassBalance



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: calcInterphaseTxfr                                      !
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
      SUBROUTINE calcInterphaseTxfr(CALLER, RxN, lnMT, lEEq, lSEq, &
         lSAg, lMMx, lSAs)

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: CALLER

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN

      DOUBLE PRECISION, INTENT(IN) :: lnMT(0:DIM_M)
! Energy equation flag
      LOGICAL, INTENT(IN) :: lEEq
! Gas/Solids Species Eq Flag
      LOGICAL, INTENT(IN) :: lSEq(0:DIM_M)
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs

      INTEGER toPhase, toPhaseCount, mCount
      INTEGER fromPhase, fromPhaseCount
      INTEGER catPhase

      INTEGER sprPhase


      INTEGER M, MM, LL
      INTEGER lM, lN

      DOUBLE PRECISION, PARAMETER :: massBalanceTol = 1.0d-3

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

! Initialize interphase exchange terms.
      IF(Allocated(RxN%rPhase)) RxN%rPhase(:) = ZERO

! If there is only one phase referenced by the reaction, there there 
! should be no interphase mass transfer.
      IF(RxN%nPhases == 1) THEN
! Interphase mass transfer is set to zero. Small inconsistancies with
! molecular weights can resunt in a non-zero value for homogeneous
! reactions. Presumably, the mass balance check caught any major errors.
         RxN%rPhase(:) = ZERO
! Void interphase transfer flags.
         DO lN = 1, RxN%nSpecies
            M = RxN%Species(lN)%pMap
            RxN%Species(lN)%mXfr = M
         ENDDO
         RxN%Classification = "Homogeneous"
! This is a multiphase reaction. 
      ELSE
! Initialize.
         toPhaseCount = 0 
         fromPhaseCount = 0
         DO M = 0, lMMx
! Determine the number of phases with a net mass gain. Record the index
! of the last phase with a net mass gain.
            IF (lnMT(M) > massBalanceTol) THEN 
               toPhaseCount = toPhaseCount + 1 
               toPhase = M
! Determine the number of phases with a net mass loss. Record the index
! index of the last phase with a net mass loss.
            ELSEIF(lnMT(M) < -massBalanceTol) THEN 
               fromPhaseCount = fromPhaseCount + 1 
               fromPhase = M
            ENDIF 
         ENDDO

! Only one phase has a net mass gain.
         IF(toPhaseCount == 1) THEN 
! Interphase mass transfer flag.
            RxN%Classification = "Heterogeneous"
            DO M = 0, lMMx 
               IF(M /= toPhase) THEN
                  IF (toPhase < M) THEN
                     LM = 1 + toPhase + ((M-1)*M)/2
                     RxN%rPhase(LM) = -lnMT(M)
                  ELSE
                     LM = 1 + M + ((toPhase-1)*toPhase)/2
                     RxN%rPhase(LM) = lnMT(M)
                  ENDIF

! Verify that if one phase's species equations are solved, that the 
! other phase's species equations are solved.

                  IF(abs(RxN%rPhase(LM)) > SMALL_NUMBER) THEN
                     IF((lSEq(toPhase) .AND. .NOT.lSEq(M)) .OR. &
                        (.NOT.lSEq(toPhase) .AND. lSEq(M))) THEN
                        IF(DMP_LOG) THEN
                           WRITE(*,1001) trim(CALLER)
                           WRITE(UNIT_LOG,1001) trim(CALLER)
                           IF(lSEq(M)) THEN
                              WRITE(*,1101) M, 'Solving'
                              WRITE(*,1101) toPhase, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) toPhase, &
                                 'Not Solving'
                           ELSE
                              WRITE(*,1101) toPhase, 'Solving'
                              WRITE(*,1101) M, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) toPhase, &
                                 'Not Solving'
                           ENDIF
                           WRITE(*,1000)
                           WRITE(UNIT_LOG,1000)
                        ENDIF
                        CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
                        CALL MFiX_EXIT(myPE)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO 

! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN
               DO lN = 1, RxN%nSpecies
                  M = RxN%Species(lN)%pMap
! The gas phase is referenced by the reaction. 
                  IF(M == 0) THEN
! The gas phase is the destination phase.
                     IF(toPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the 
! gas phase.
                        mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                        DO MM = 1, lMMx
                           LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phase MM.
                           IF(RxN%rPhase(LM) > 0) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                              RxN%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                              mCount = mCount + 1
                           ENDIF
                        ENDDO
                        IF(mCount /= 1) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), &
                              lSAs(:,:))
                           CALL MFiX_EXIT(myPE)
                        ENDIF

! A solids phase is the destination phase.
                     ELSE
! Since only one phase was detected with a net mass gain and the gas
! phase was detected as a source phase, then all the gas is assigned
! to the destination phase.
                        RxN%Species(lN)%mXfr = toPhase
! This variable is not used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                  ELSE
! Check to see if phase M transfer mass to another solids phase.
                     DO LL = 1, lMMx-1
                        DO MM = LL + 1, lMMx
                           IF(M /= LL .AND. M /= MM) CYCLE
                           LM = LL + 1 + (MM-1)*MM/2
                           IF(RxN%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                           IF( M == LL .AND. &
                              RxN%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                              RxN%Species(lN)%xXfr =  &
                                 abs(lnMT(MM) / lnMT(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                           ELSEIF( M == MM .AND. &
                              RxN%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                              RxN%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(LL) / lnMT(MM))
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF ! Gas or Solids phase.
               ENDDO ! Species Loop
            ENDIF ! Energy Equation
! If there is only one phase with a net mass loss, setup the array for
! interphase mass transfer.
         ELSEIF(fromPhaseCount == 1) THEN
            RxN%Classification = "Heterogeneous"
            DO M = 0, lMMx 
               IF (M /= fromPhase) THEN
                  IF(M < fromPhase) THEN
                     LM = 1 + M + ((fromPhase-1)*fromPhase)/2
                     RxN%rPhase(LM) =  lnMT(M)
                  ELSE
                     LM = 1 + fromPhase + ((M-1)*M)/2
                     RxN%rPhase(LM) = -lnMT(M)
                  ENDIF

! Verify that if one phase's species equations are solved, that the 
! other phase's species equations are solved.
                  IF(abs(RxN%rPhase(LM)) > SMALL_NUMBER) THEN
                     IF((lSEq(fromPhase) .AND. .NOT.lSEq(M)) .OR.   &
                        (.NOT.lSEq(fromPhase) .AND. lSEq(M))) THEN
                        IF(DMP_LOG) THEN
                           WRITE(*,1001) trim(CALLER)
                           WRITE(UNIT_LOG,1001) trim(CALLER)
                           IF(lSEq(M)) THEN
                              WRITE(*,1101) M, 'Solving'
                              WRITE(*,1101) fromPhase, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) fromPhase, &
                                 'Not Solving'
                           ELSE
                              WRITE(*,1101) toPhase, 'Solving'
                              WRITE(*,1101) M, 'Not Solving'
                              WRITE(UNIT_LOG,1101) fromPhase, 'Solving'
                              WRITE(UNIT_LOG,1101) M, 'Not Solving'
                           ENDIF
                           WRITE(*,1000)
                           WRITE(UNIT_LOG,1000)
                        ENDIF
                        CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
                        CALL MFiX_EXIT(myPE)
                     ENDIF
                  ENDIF
               ENDIF
            END DO 

! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN
               DO lN = 1, RxN%nSpecies
                  M = RxN%Species(lN)%pMap
! Gas/solids reaction: Enthalpy transfer from mass transfer is only
! calculated from gas phase species.
                  IF(M == 0) THEN
! Gas phase is the source phase.
                     IF(fromPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the 
! gas phase.
                        mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                        DO MM = 1, lMMx
                           LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phases MM.
                           IF(RxN%rPhase(LM) < 0) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                              RxN%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                              mCount = mCount + 1
                           ENDIF
                        ENDDO
                        IF(mCount /=1 ) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),  &
                              lSAs(:,:))
                           CALL MFiX_EXIT(myPE)
                        ENDIF
                     ELSE
! There can be only one solids phase fuel. Store the phase of the
! solids phase reactant.
                        RxN%Species(lN)%mXfr = fromPhase
! Mass fraction of transfered material.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                  ELSE
! Check to see if phase M transfer mass to another solids phase.
                     DO LL = 1, lMMx-1
                        DO MM = LL + 1, lMMx
                           IF(M /= LL .AND. M /= MM) CYCLE
                           LM = LL + 1 + (MM-1)*MM/2
                           IF(RxN%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                           IF( M == LL .AND. &
                              RxN%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(MM) / lnMT(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                           ELSEIF( M == MM .AND. &
                              RxN%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                              RxN%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(LL) / lnMT(MM))
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
            IF(RxN%nPhases > 0) RxN%Classification = "Catalytic"
            RxN%rPhase(:)  = ZERO
! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN

! Identify the catalyst phase.
               catPhase = -1
               DO lN= 1, RxN%nSpecies
                  IF(COMPARE(RxN%Species(lN)%Coeff,ZERO)) THEN
                     IF(catPhase /= -1) THEN
                        IF(catPhase /= RxN%Species(lN)%pMap) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),    &
                              lSAs(:,:))
                           CALL MFiX_EXIT(myPE)
                        ENDIF
                     ELSE
                        catPhase = RxN%Species(lN)%pMap
                     ENDIF
                  ENDIF
               ENDDO
! Verify that a catalyst phase was found.
               IF(catPhase == -1) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1003) trim(CALLER), 'catalyst', &
                        trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1003) trim(CALLER), &
                        'catalyst', trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
                  CALL MFiX_EXIT(myPE)
               ENDIF

! Identify the reactant phase.
               toPhase = -1
               DO lN = 1, RxN%nSpecies
                  IF(.NOT.COMPARE(RxN%Species(lN)%Coeff,ZERO)) THEN
                     IF(toPhase /= -1) THEN
                        IF(toPhase /= RxN%Species(lN)%pMap) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), &
                              lSAs(:,:))
                           CALL MFiX_EXIT(myPE)
                        ENDIF
                     ELSE
                        toPhase = RxN%Species(lN)%pMap
                     ENDIF
                  ENDIF
               ENDDO
! Verify that a reacting phase was found.
               IF(toPhase == -1) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1003) trim(CALLER), 'reacting', &
                        trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1003) trim(CALLER), 'reacting', &
                        trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
                  CALL MFiX_EXIT(myPE)
               ENDIF

! Something when wrong.
               IF(catPhase == toPhase) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1004) trim(CALLER), trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1004) trim(CALLER),trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
                  CALL MFiX_EXIT(myPE)
!Gas/solid cataltyic reaction:
               ELSEIF(toPhase == 0) THEN
                  DO lN = 1, RxN%nSpecies
                     IF(RxN%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                        RxN%Species(lN)%mXfr = catPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
                  ENDDO
               ELSEIF(catPhase == 0) THEN
                  DO lN = 1, RxN%nSpecies
                     IF(RxN%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                        RxN%Species(lN)%mXfr = toPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF ! Energy Equation
         ELSE
! Two or more phases have a net mass loss and two or more phases have
! a net mass gain. Therefore, the interphase mass transfer cannot be
! concluded.
            IF(DMP_LOG) THEN
               CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:))
               WRITE(*,1002) trim(CALLER), trim(RxN%ChemEq)
               WRITE(*,1000)
               WRITE(UNIT_LOG,1002) trim(CALLER), trim(RxN%ChemEq)
               WRITE(UNIT_LOG,1000)
            ENDIF
            CALL MFiX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN

! Error Messages
!---------------------------------------------------------------------//

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1001 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1001: A chemical reaction or', &
         ' phase change was detected between',/' a phases solving',    &
         ' species equations and another phase not solving',/          &
         ' species equations.',/)

 1101 FORMAT(' Phase ',I2,': ',A,' species equations.')

 1002 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1002: Reaction complexity',    &
         ' exceeds implementation capabilities.',/' Unable to',        &
         ' determine unambiguously interphase heat or mass transfer.', &
         //' Reaction: ',A,//' Consider splitting the chemical',       &
         ' reaction equation into two or more',/' separate equations.',&
         ' The same reaction rate calculated in usr_rates',/' can be', &
         ' used for the multiple reactions to ensure mass')

 1003 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1003: Unable to determine ',A, &
         ' phase for catalytic reaction'/1X,A,'.')

 1004 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1004: Unable to distinguish',  &
         ' catalyst phase from reacting phase',/' for catalytic',      &
         ' reaction ',A,'.')

      END SUBROUTINE calcInterphaseTxfr


      END MODULE RXN_COM
