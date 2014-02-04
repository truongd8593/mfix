!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_GAS_PHASE                                        !
!  Purpose: Check the gas phase input section                          !
!                                                                      !
!  Author: P.Nicoletti                                Date: 02-DEC-91  !
!          J.Musser                                   Date: 01-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GAS_PHASE


! Global Variables:
!---------------------------------------------------------------------//
      USE compar
      USE param 
      USE param1 
      USE physprop
      USE funits 
      USE run
      USE indices
      USE rxns

! Global Parameters:
!---------------------------------------------------------------------//

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables::
!---------------------------------------------------------------------//
      INTEGER :: N 

! Flag that the energy equations are solved and constant gas phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPG

! Flag that the average molecular weight (MW_AVG) and constant gas
! phase density are undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL MWg_ROg

! Flag that the gas phase species equations are solved and the 
! molecular weight for a species is not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWg

! Flag that the user was already warned why a call to the thermo-
! chemical database is being made.
      LOGICAL WARNED_USR

! Flag indicating that the thermochemical database header was output 
! to the screen. (Miminize messages)
      LOGICAL thermoHeader


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GAS_PHASE")


! CHECK MU_g0
      IF (MU_G0 <= ZERO) THEN 
         WRITE(ERR_MSG,1002) 'MU_G0', MU_G0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! CHECK K_g0
      IF (K_G0 < ZERO) THEN 
         WRITE(ERR_MSG,1002) 'K_G0', K_G0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! CHECK C_pg0
      IF (C_PG0 < ZERO) THEN 
         WRITE(ERR_MSG,1002) 'C_PG0', C_PG0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! CHECK DIF_g0
      IF (DIF_G0 < ZERO) THEN
         WRITE(ERR_MSG,1002) 'DIF_g0', DIF_g0
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 




 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',E14.6,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')



! Legacy checks for species equations. If values for the new input
! method are found, MFIX exits.
      IF(USE_RRATES) THEN
         IF(SPECIES_EQ(0)) THEN

            IF(NMAX(0) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000)'NMAX(0)','specified'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

            ELSEIF(NMAX_G /= UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000)'NMAX_g', 'undefined'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

            ELSEIF(NMAX(0) > DIMENSION_N_G) THEN
               WRITE(ERR_MSG,1003)'NMAX(0)', NMAX(0)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            NMAX(0) = 1
         ENDIF
      ENDIF

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be 'A,/'when ',    &
         'USE_RRATES is .TRUE.'/,'Please correct the mfix.dat file')



      IF(NMAX_g /= UNDEFINED_I)THEN
! Verify that the number of species is within range.
         IF(NMAX_g > DIMENSION_N_G) THEN
            WRITE(ERR_MSG,1053)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE
! Copy the new keyword entry into the runtime variable.
            NMAX(0) = NMAX_g
         ENDIF
! Neither of the variables for given in the data file.
      ELSE
         IF(SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG,1052)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE
            NMAX(0) = 1
            NMAX_g = 1
         ENDIF
      ENDIF


! CHECK MW_AVG
      IF (SPECIES_EQ(0)) THEN
! MW_AVG is defined and the gas phase species equations are solved, then
! the user specified average molecular weight is ignored. The gas phase
! mixture molecular weight (MW_MIX_g) is used instead.
         IF (MW_AVG /= UNDEFINED) THEN 
            WRITE (ERR_MSG, 1410)
            CALL FLUSH_ERR_MSG
            MW_AVG = UNDEFINED 
         ENDIF 

 1410 FORMAT('Message: MW_AVG not needed when solving species equations')

      ELSE 
! When the species equations are not solved and the gas phase is
! compressible, verify that the user provided average molecular weight
! has a physical value. (This does not include the case where MW_AVG
! is UNDEFINED.)
         IF (RO_G0 == UNDEFINED) THEN 
            IF (MW_AVG <= ZERO) THEN 
               WRITE(ERR_MSG, 1002) 'MW_AVG', MW_AVG
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ELSE 
! Gas density for incompressible flows must be positive.
            IF (RO_G0 < ZERO) THEN
               WRITE(ERR_MSG, 1002) 'RO_G0', RO_G0
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
! Incompressible simulations do not need MW_AVG. Notify the user that 
! the provided data is ignored.
            IF (MW_AVG /= UNDEFINED)THEN
               WRITE(ERR_MSG, 1400)
               CALL FLUSH_ERR_MSG
            ENDIF

 1400 FORMAT('Message: MW_AVG not needed when RO_g0 is specified.')

         ENDIF 
      ENDIF


! Check MW_g
      IF(.NOT.USE_RRATES) THEN
! Initialize flag indicating the database was read for a species.
         rDatabase(0,:) = .FALSE.
! Flag indicating if the user was already warned.
         WARNED_USR = .FALSE.
! Flag indicating the search header was already written.
         thermoHeader = .FALSE.
! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
         EEQ_CPG = .FALSE.
         IF(ENERGY_EQ .AND. C_PG0 == UNDEFINED) EEQ_CPG = .TRUE.
! Flag that the average molecular weight and constant gas density are
! undefined.
         MWg_ROg = .FALSE.
         IF(MW_AVG == UNDEFINED) THEN
            DO N=1,NMAX(0)
               IF(MW_g(N) == UNDEFINED) MWg_ROg = .TRUE.
            ENDDO
         ENDIF
         IF(MWg_ROg .AND. RO_G0==UNDEFINED) THEN
            MWg_ROg = .TRUE.
         ELSE
            MWg_ROg = .FALSE.
         ENDIF

! Loop over the gas phase species
         DO N = 1, NMAX(0)
            SEQ_MWg = .FALSE.
            IF(SPECIES_EQ(0) .AND. MW_G(N)==UNDEFINED) SEQ_MWg = .TRUE.

! If the gas phase species equations are solved, or the gas phase is
! compressible and the user has not supplied the average molecular 
! weight, or the energy equations are solved and individual species
! molecular wegiths are not given, then try and get them from the 
! thermochemical database.
! A final thermochemical check is preformed in check_data_09. If neither
! of the above conditions result in species data being read from the
! database AND a particular species is referenced by a chemical equation
! then a call to read_database is forced.
            IF(EEQ_CPG  .OR. MWg_ROg .OR. SEQ_MWg) THEN
! Notify the user of the reason the thermochemical database is used.
               IF(.NOT.WARNED_USR) THEN
                  IF(EEQ_CPG) THEN
                     WRITE(ERR_MSG,1054)
                     CALL FLUSH_ERR_MSG
                  ENDIF
                  IF(SEQ_MWg) THEN
                     WRITE(ERR_MSG,1055)
                     CALL FLUSH_ERR_MSG
                  ENDIF
                  IF(MWg_ROg) THEN
                     WRITE(ERR_MSG,1056)
                     CALL FLUSH_ERR_MSG
                  ENDIF
                  WARNED_USR = .TRUE.
               ENDIF
! Flag that the species name is not provided.
               IF(SPECIES_g(N) == UNDEFINED_C) THEN
                  WRITE(ERR_MSG,1057)N
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
! Read the database.
               IF(.NOT.thermoHeader) THEN
                  WRITE(ERR_MSG,1058)
                  CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
                  thermoHeader = .TRUE.
               ENDIF
               WRITE(ERR_MSG,1059) N, trim(SPECIES_g(N))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               CALL READ_DATABASE('TFM', 0, N, SPECIES_g(N),MW_g(N))
! Flag variable to stating that the database was read.
               rDatabase(0,N) = .TRUE.
            ENDIF 
         ENDDO
! Flag the legacy variable.
         DATABASE_READ = .TRUE.

! Verify that no molecular weight data was given for species that do
! not exist.
         DO N = NMAX(0) + 1, DIMENSION_N_G 
            IF (MW_G(N) /= UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'MW_g defined for N > NMAX(0)', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1501) N 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
! Legacy checks for MW_g
      ELSE
         IF (SPECIES_EQ(0) .OR. &
            MW_AVG==UNDEFINED .AND. RO_G0==UNDEFINED) THEN 

            DO N = 1, NMAX(0)
               IF(MW_G(N) == UNDEFINED) THEN
                  CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                     'Value of MW_g not specified', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1500) N 
! No need to abort since MW will be read from database
!               CALL ERROR_ROUTINE (' ', ' ', 1, 3)                
               ENDIF 
            ENDDO 
            DO N = NMAX(0) + 1, DIMENSION_N_G 
               IF (MW_G(N) /= UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                     'MW_g defined for N > NMAX(0)', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1501) N 
                  CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
               ENDIF 
            ENDDO 
         ENDIF 
      ENDIF


! Finalize the error manager
      CALL FINL_ERR_MSG

      RETURN  


 1050 FORMAT('Message: NMAX(0) is specied for the gas phase. This is', &
         ' a legacy',/'variable and NMAX_g should be used. Copying',   &
         ' NMAX(0) to NMAX_g.')

 1051 FORMAT('Message: NMAX_g and NMAX(0) are both given for the gas', &
         ' phase and do',/' not match. NMAX(0) is a legacy variable',  &
         ' and is not required.',/' Please correct the data file.')

 1052 FORMAT(' Message: The number of gas speices (NMAX_g) is not',    &
         ' specified. Please',/' correct the data file.')

 1053 FORMAT('Message: The number of gas species (NMAX_g) is too ',    &
         'large. Please',/' correct the data file.')

 1054 FORMAT('Message: The energy equations are being solved (ENERGY_',&
         'EQ) and the',/' constant gas specific heat is undefined',    &
         ' (C_PG0). Thus, the thermo-',/' chemical database will be',  &
         ' used to gather specific heat data on the',/' individual',   &
         ' gas phase species.')

 1055 FORMAT('Message: Gas phase species equations are being solved, ',&
         'and one or more',/'species molecular weights are undefined.',&
         ' Thus, the thermochemical',/'database will be used to',      &
         ' gather molecular weight data on the gas',/'phase species.')

 1056 FORMAT(//1X,70('*')/' From: CHECK_DATA_05',/                     &
         ' Message: MW_AVG and RO_G0 are undefined. Thus, the',        &
         ' thermochemical',/' database will be used to gather',        &
         ' molecular weight data on the gas',/' phase species.',/      &
         1X,70('*')/)

 1057 FORMAT('Message: Gas phase species ',I2,' name (SPECIES_g) is',  &
         ' undefined.',/'Please correct the data file.')

 1058 FORMAT(/'  Searching thermochemical databases for gas phase',    &
         ' species data')

 1059 FORMAT(/2x,'>',I3,': Species: ',A)



 1500 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' not specified') 
 1501 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' specified') 

      END SUBROUTINE CHECK_GAS_PHASE
