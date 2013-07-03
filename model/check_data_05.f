!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_05                                          C
!  Purpose: check the gas phase input section                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MU_g0, MW_AVG, RO_g0                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_05 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE compar
      USE param 
      USE param1 
      USE physprop
      USE funits 
      USE run
      USE indices
      USE rxns

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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

!-----------------------------------------------

      IF(.NOT.USE_RRATES) THEN
         IF(NMAX_g == UNDEFINED_I .AND. NMAX(0) /= UNDEFINED_I) THEN
! The number of gas phase species was given in NMAX. Warn the user
! to use the correct variable. Copy the old variable entry to the new
! variable name for later pre-processing and continue.
            IF(DMP_LOG) THEN
               WRITE(*,1050)myPE
               WRITE(UNIT_LOG,1050)myPE
            ENDIF
! Verify that the number of species is within range.
            IF(NMAX(0) > DIMENSION_N_G .AND. myPE == PE_IO) THEN
               WRITE(*,1053) myPE
               WRITE(UNIT_LOG,1053) myPE
               CALL MFIX_EXIT(myPE)
            ELSE
! Copy the legacy entry into the new variable.
               NMAX_g = NMAX(0)
            ENDIF
! If for whatever reason, the number of species are given in both
! variables, make sure they match.
         ELSEIF(NMAX(0) /= UNDEFINED_I .AND. NMAX_g /= UNDEFINED_I) THEN
            IF(NMAX(0) /= NMAX_g) THEN
               IF(myPE .EQ. PE_IO) THEN
                  WRITE(*,1051)myPE
                  WRITE(UNIT_LOG,1051)myPE
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
! Notify the user of the correct variable name.
            WRITE(*,1050)myPE
            WRITE(UNIT_LOG,1050)myPE
         ELSEIF(NMAX_g /= UNDEFINED_I .AND. NMAX(0) == UNDEFINED_I) THEN
! Verify that the number of species is within range.
            IF(NMAX_g > DIMENSION_N_G .AND. myPE == PE_IO) THEN
               WRITE(*,1053) myPE
               WRITE(UNIT_LOG,1053) myPE
               CALL MFIX_EXIT(myPE)
            ELSE
! Copy the new keyword entry into the runtime variable.
               NMAX(0) = NMAX_g
            ENDIF
! Neither of the variables for given in the data file.
         ELSE
            IF(SPECIES_EQ(0)) THEN
               IF(myPE .EQ. PE_IO) THEN
                  WRITE(*,1052) myPE
                  WRITE(UNIT_LOG,1052) myPE
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ELSE
               NMAX(0) = 1
               NMAX_g = 1
            ENDIF
         ENDIF
      ELSE
! Legacy check.
         IF (NMAX(0) == UNDEFINED_I) THEN 
            IF (SPECIES_EQ(0)) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'Number of gas species (NMAX(0)) not specified',0,2) 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ELSE 
               NMAX(0) = 1 
            ENDIF 
         ELSEIF (NMAX(0) > DIMENSION_N_G) THEN 
            CALL ERROR_ROUTINE ('CHECK_DATA_05',&
               'NMAX(0) is too large',0,2)
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) NMAX(0), DIMENSION_N_G 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF
      ENDIF


! CHECK MW_AVG
      IF (SPECIES_EQ(0)) THEN
! MW_AVG is defined and the gas phase species equations are solved, then
! the user specified average molecular weight is ignored. The gas phase
! mixture molecular weight (MW_MIX_g) is used instead.
         IF (MW_AVG /= UNDEFINED) THEN 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1410) 
            MW_AVG = UNDEFINED 
         ENDIF 
      ELSE 
! When the species equations are not solved and the gas phase is
! compressible, verify that the user provided average molecular weight
! has a physical value. (This does not include the case where MW_AVG
! is UNDEFINED.)
         IF (RO_G0 == UNDEFINED) THEN 
            IF (MW_AVG <= ZERO) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05',&
                  'MW_AVG is unphysical',0,2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) MW_AVG 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ELSE 
! Gas density for incompressible flows must be positive.
            IF (RO_G0 < ZERO) THEN 
               CALL ERROR_ROUTINE ('CHECK_DATA_05', &
                  'Value of RO_g0 is unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1300) RO_G0 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF
! Incompressible simulations do not need MW_AVG. Notify the user that 
! the provided data is ignored.
            IF (MW_AVG /= UNDEFINED .AND. DMP_LOG)WRITE (UNIT_LOG, 1400) 
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
                  IF(EEQ_CPG .AND. myPE .EQ. PE_IO) THEN
                     WRITE(*,1054)
                     WRITE(UNIT_LOG,1054)
                  ENDIF
                  IF(SEQ_MWg .AND. myPE .EQ. PE_IO) THEN
                     WRITE(*,1055)
                     WRITE(UNIT_LOG,1055)
                  ENDIF
                  IF(MWg_ROg .AND. myPE .EQ. PE_IO) THEN
                     WRITE(*,1056)
                     WRITE(UNIT_LOG,1056)
                  ENDIF
                  WARNED_USR = .TRUE.
               ENDIF
! Flag that the species name is not provided.
               IF(SPECIES_g(N) == UNDEFINED_C) THEN
                  IF(myPE .EQ. PE_IO) THEN
                     WRITE(*,1057)N
                     WRITE(UNIT_LOG,1057)N
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF
! Read the database.
               IF(myPE .EQ. PE_IO) THEN
                  IF(.NOT.thermoHeader) THEN
                     WRITE(*,1058)
                     WRITE(UNIT_LOG,1058)
                     thermoHeader = .TRUE.
                  ENDIF
                  WRITE(*,1059) N, trim(SPECIES_g(N))
                  WRITE(UNIT_LOG,1059) N, trim(SPECIES_g(N))
               ENDIF
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

! CHECK MU_g0
      IF (MU_G0 <= ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05', &
            'MU_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) MU_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK K_g0
      IF (K_G0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05', &
            'K_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1110) K_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK C_pg0
      IF (C_PG0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05',&
            'C_pg0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1120) C_PG0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK DIF_g0
      IF (DIF_G0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_05',&
            'DIF_g0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1130) DIF_G0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 


      RETURN  

 1000 FORMAT(1X,/,1X,'NMAX(0) in mfix.dat = ',I3,/,1X,&
         'DIMENSION_N_g in param.inc = ',I3) 
 1050 FORMAT(//1X,70('*')/' (PE ',I6,'): From: CHECK_DATA_05',/        &
         ' Message: NMAX(0) is specied for the gas phase. This is a',  &
         ' legacy',/' variable and NMAX_g should be used. Copying',    &
         ' NMAX(0) to NMAX_g.',/1X,70('*')/)

 1051 FORMAT(//1X,70('*')/' (PE ',I6,'): From: CHECK_DATA_05',/        &
         ' Message: NMAX_g and NMAX(0) are both given for the gas',    &
         ' phase and do',/' not match. NMAX(0) is a legacy variable',  &
         ' and is not required.',/' Please correct the data file.',    &
         /1X,70('*')/)

 1052 FORMAT(//1X,70('*')/' (PE ',I6,'): From: CHECK_DATA_05',/        &
         ' Message: The number of gas speices (NMAX_g) is not',        &
         ' specified. Please',/' correct the data file.',/1X,70('*')/)

 1053 FORMAT(//1X,70('*')/' (PE ',I6,'): From: CHECK_DATA_05',/        &
         ' Message: The number of gas species (NMAX_g) is too large!', &
         ' Please',/' correct the data file.',/1X,70('*')/)

 1054 FORMAT(//1X,70('*')/' From: CHECK_DATA_05',/                     &
         ' Message: The energy equations are being solved (ENERGY_EQ)',&
         ' and the',/' constant gas specific heat is undefined',       &
         ' (C_PG0). Thus, the thermo-',/' chemical database will be',  &
         ' used to gather specific heat data on the',/' individual',   &
         ' gas phase species.',/1X,70('*')/)

 1055 FORMAT(//1X,70('*')/' From: CHECK_DATA_05',/                     &
         ' Message: Gas phase species equations are being solved, and',&
         ' one or more',/' species molecular weights are undefined.',  &
         ' Thus, the thermochemical',/' database will be used to',     &
         ' gather molecular weight data on the gas',/' phase species.',&
         /1X,70('*')/)

 1056 FORMAT(//1X,70('*')/' From: CHECK_DATA_05',/                     &
         ' Message: MW_AVG and RO_G0 are undefined. Thus, the',        &
         ' thermochemical',/' database will be used to gather',        &
         ' molecular weight data on the gas',/' phase species.',/      &
         1X,70('*')/)

 1057 FORMAT(/1X,70('*')/' From: CHECK_DATA_05',/                      &
         ' Message: Gas phase species ',I2,' name (SPECIES_g) is',     &
         ' undefined.',/' Please',' correct the data file.',           &
         /1X,70('*')/)

 1058 FORMAT(/'  Searching thermochemical databases for gas phase',    &
         ' species data')

 1059 FORMAT(2x,'>',I3,': Species: ',A)

 1100 FORMAT(1X,/,1X,'MU_g0   in mfix.dat = ',G12.5) 
 1110 FORMAT(1X,/,1X,'K_g0   in mfix.dat = ',G12.5) 
 1120 FORMAT(1X,/,1X,'C_pg0   in mfix.dat = ',G12.5) 
 1130 FORMAT(1X,/,1X,'DIF_g0   in mfix.dat = ',G12.5)          
 1200 FORMAT(1X,/,1X,'MW_AVG in mfix.dat = ',G12.5) 
 1300 FORMAT(1X,/,1X,'RO_g0   in mfix.dat = ',G12.5) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_DATA_05',/&
         ' Message: Since RO_g0 is specified MW_AVG will be ignored',/, 1X&
         ,70('*')/) 
 1410 FORMAT(/1X,70('*')//' From: CHECK_DATA_05',/&
         ' Message: Since gas phase rxns are specified MW_AVG',&
         ' will be ignored',/1X,70('*')/) 
 1500 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' not specified') 
 1501 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' specified') 

      END SUBROUTINE CHECK_DATA_05 
