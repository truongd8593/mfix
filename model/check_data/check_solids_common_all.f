!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_ALL                                 !
!  Purpose: Check the solid phase input that is common to all solids   !
!  phase models.                                                       !
!                                                                      !
!    ****** DO NOT PUT MODEL SPECIFIC CHECKS IN THIS ROUTINE ******    !
!                                                                      !
!  Use the companion routines for checks specific to a particular      !
!  solids phase model:                                                 !
!                                                                      !
!    > CHECK_CONTINUUM_SOLDIS :: TRM solids phase model)               !
!    > CHECK_DES_SOLDIS       :: DEM solids phase model)               !
!    > CHECK_MPPIC_SOLDIS     :: MPPIC solids phase model)             !
!                                                                      !
!                                                                      !
!  Author: P.Nicoletti                               Date: 02-DEC-91   !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_ALL

      USE param 
      USE param1 
      USE run
      USE indices
      USE physprop
      USE constant
      USE discretelement
      USE funits 
      USE mfix_pic
      USE compar
      USE rxns

      use error_manager

      IMPLICIT NONE

      INTEGER :: LC, N

      INTEGER :: MMAX_L

!-----------------------------------------------


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_ALL")

! Set the number of solids phases to be checked.
      MMAX_L = SMAX + DES_MMAX

! Check D_p0
      DO LC = 1, MMAX_L
         IF(D_P0(LC) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) iVar('D_p0',LC)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(D_P0(LC) <= ZERO)THEN
            WRITE(ERR_MSG, 1001) iVar('D_p0',LC), iVal(D_P0(LC))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      DO LC = MMAX_L+1, DIM_M
         IF(D_P0(LC) /= UNDEFINED)THEN
            WRITE(ERR_MSG,1002) 'D_p0'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDDO 

! Check K_s0
      IF (K_S0 < ZERO) THEN 
         WRITE(ERR_MSG, 1001) 'K_s0', iVal(K_s0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check C_ps0
      IF (C_PS0 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'C_ps0', iVal(C_ps0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Check the input specifications for solids species.
      IF(USE_RRATES)THEN
         CALL CHECK_SOLIDS_SPECIES_LEGACY(MMAX_L)
      ELSE
         CALL CHECK_SOLIDS_SPECIES(MMAX_L)
      ENDIF


! Check the solids density input parameters.
      CALL CHECK_SOLIDS_DENSITY(MMAX_L)

! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: Too many',A,' values ',       &
         'specified',/'Please correct the mfix.dat file.')


      END SUBROUTINE CHECK_SOLIDS_COMMON_ALL


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES                                     !
! Purpose: Clear solids phase species data from continuum variables    !
! and notify the user.                                                 !
!                                                                      !
! Author: J. Musser                                  Date: 02-NOV-12   !
! Reviewer:                                          Date:             !
!                                                                      !
! Literature/Document References: MFiX Readme                          !
!                                                                      !
! Variables modified: NMAX(1:DIM_M), NMAX_s(:), SPECIES_S(:,:),        !
!                     SPECIES_ALIAS_S(:,:)                             !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES(MMAX_LL)

      USE compar
      USE funits 
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE RUN

      use error_manager

      IMPLICIT NONE


! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Flag that the user was already warned why a call to the thermo-
! chemical database is being made.
      LOGICAL WARNED_USR

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPS

! Flag that the solids phase species equations are solved and the 
! molecular weight for a species are not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWs

      LOGICAL thermoHeader

      INTEGER :: LC, N


      CALL INIT_ERR_MSG("CHECK_SOLIDS_SPECIES")


! Reconcile the new species input method with the legacy input method.
      DO LC=1, MMAX_LL
         IF(SPECIES_EQ(LC)) THEN

            IF(NMAX_s(LC) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,1000) iVar('NMAX_S',LC)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX_s(LC) > DIM_N_S) THEN
               WRITE(ERR_MSG,1001) trim(iVar('NMAX_s',LC)),            &
                  trim(iVal(NMAX_s(LC)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSE
               NMAX(LC) = NMAX_s(LC)
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            NMAX(LC) = merge(1, NMAX_S(LC), NMAX_S(LC) == UNDEFINED_I)
         ENDIF
      ENDDO

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be 'A,/'when ',    &
         'USE_RRATES is .TRUE.'/,'Please correct the mfix.dat file')


! Flag indicating if the user was already warned.
      WARNED_USR = .FALSE.
! Flag indicating the search header was already written.
      thermoHeader = .FALSE.
! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
      EEQ_CPS = .FALSE.
      IF(ENERGY_EQ .AND. C_PS0 == UNDEFINED) EEQ_CPS = .TRUE.

! Check MW_s if solids species are present    
      DO LC = 1, MMAX_LL
! Initialize flag indicating the database was read for a species.
         rDatabase(LC,:) = .FALSE.
         DO N = 1, NMAX(LC)
! Legacy code may use 'hard coded' functions for specific heat and
! therefore should avoid the follow check.
! Flag that the solids phase species equations are solved and the 
! molecular weight for a species are not given in the data file.
            SEQ_MWs = (SPECIES_EQ(LC) .AND. MW_S(LC,N) == UNDEFINED)

! The thermodynamic data base provides the specific heat, molecular
! weights, and heat of formation.
! Check the thermodynamic database if:
!   1) the energy equation is solved and constant solids phase specific
!      heat isn't given. (EEQ_CPS = .TRUE.)
!   2) the species energy equation is solved and the molecular weights
!      for the solids phase species are not given. (SEQ_MWs = .TRUE.)
! A final thermochemical check is preformed in check_data_09. If neither
! of the above conditions result in species data being read from the
! database AND a particular species is referenced by a chemical equation
! then a call to read_database is forced.
            IF(EEQ_CPS  .OR. SEQ_MWs) THEN
! Notify the user of the reason the thermochemical database is used.
               IF(.NOT.WARNED_USR) THEN

                  IF(EEQ_CPS) THEN
                     WRITE(ERR_MSG, 2001)
                     CALL FLUSH_ERR_MSG
                  ENDIF
 2001 FORMAT('Message 2001: The energy equations are being solved ',   &
         '(ENERGY_EQ) and the',/'specified constant solids specific ', &
         'heat is undefined (C_PS0). Thus,',/'the thermochemical data',&
         'base will be used to gather specific heat data',/'on the ',  &
         'individual soids phase species.')

                  IF(SEQ_MWs) THEN
                     WRITE(ERR_MSG,2002)LC
                     CALL FLUSH_ERR_MSG
                  ENDIF
 2002 FORMAT('Message 2002: Solids phase',I2,' species equations are ',&
         'being solved, and one',/'or more species molecular weights ',&
         'are undefined. Thus, the thermo-',/'chemical database will ',&
         'be used to gather molecular weight data on the',/'solids ',  &
         'phase speicies.')


! Set a flag to prevent the user from getting the same message over
! and over again.
                  WARNED_USR = .TRUE.
               ENDIF
! Flag that the species name is not provided.
               IF(SPECIES_s(LC,N) == UNDEFINED_C) THEN
                  WRITE(ERR_MSG,2003) LC, N
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
 2003 FORMAT('Message 2003: Solids phase ',I2,' species ',I2,' name', &
         ' (SPECIES_s) is undefined.',/'Please correct the data file.')

               CALL READ_DATABASE('TFM', LC, N, SPECIES_s(LC,N),      &
                  MW_S(LC,N))
! Flag variable to stating that the database was read.
               rDatabase(LC,N) = .TRUE.
! Flag the legacy variable to prevent re-reading the database.
               DATABASE_READ = .TRUE.
            ENDIF
         ENDDO ! Loop over species


! Verify that no additional species information was given.
         DO N = NMAX(LC) + 1, DIM_N_S 
            IF (MW_S(LC,N) /= UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_SOLIDS_COMMON', &
                  'MW_s defined for N > NMAX(m)', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDDO ! Loop over solids phases

      CALL FINL_ERR_MSG

 1410 FORMAT(1X,/,1X,'Solids phase = ',I2,'   Species = ',I3) 


      RETURN
      END SUBROUTINE CHECK_SOLIDS_SPECIES


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES_LEGACY                              !
! Purpose: Clear solids phase species data from continuum variables    !
! and notify the user.                                                 !
!                                                                      !
! Author: J. Musser                                  Date: 03-FEB-14   !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY(MMAX_LL)

      USE compar
      USE funits 
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE RUN

      use error_manager


      IMPLICIT NONE

! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

      INTEGER :: LC, N


      CALL INIT_ERR_MSG("CHECK_SOLIDS_SPECIES_LEGACY")

! Reconcile the new species input method with the legacy input method.
      DO LC=1, MMAX_LL
         IF(SPECIES_EQ(LC)) THEN

! Legacy checks for species equations.
            IF(NMAX_s(LC) /= UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) iVar('NMAX_s',LC), 'undefined'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX(LC) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) iVar('NMAX',LC),'specified'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX(LC) > DIM_N_S) THEN
               WRITE(ERR_MSG,1001) iVar('NMAX',LC), iVal(NMAX(LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            IF(NMAX(LC) == UNDEFINED) NMAX(LC) = 1
         ENDIF
      ENDDO

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be 'A,/'when ',    &
         'USE_RRATES is .TRUE.'/,'Please correct the mfix.dat file')


! Check MW_s if solids species are present    
      DO LC = 1, MMAX_LL
! Initialize flag indicating the database was read for a species.
         DO N = 1, NMAX(LC)
            IF (MW_S(LC,N) == UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_SOLIDS_COMMON', &
                  'Species molecular weight undefined', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
            ENDIF
         ENDDO ! Loop over species
         DO N = NMAX(LC) + 1, DIM_N_S 
            IF (MW_S(LC,N) /= UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('CHECK_SOLIDS_COMMON', &
                  'MW_s defined for N > NMAX(m)', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDDO ! Loop over solids phases


      CALL FINL_ERR_MSG


 1410 FORMAT(1X,/,1X,'Solids phase = ',I2,'   Species = ',I3) 

      RETURN
      END SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY


 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON                                     !
!  Purpose: check the solid phase input section                        !
!                                                                      !
!  Author: P.Nicoletti                               Date: 02-DEC-91   !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DENSITY(MMAX_LL)

      USE param 
      USE param1 
      USE run
      USE indices
      USE physprop
      USE constant
      USE discretelement
      USE funits 
      USE mfix_pic
      USE compar
      USE rxns

      use error_manager

      IMPLICIT NONE

! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL


      INTEGER :: LC, N


! Flag that a no fatal error was detected.
      LOGICAL :: PASSED
! Flag that a fatal error was detected.
      LOGICAL :: FAILED, WARNED_USR


! Solids phase density calculation
      DOUBLE PRECISION, EXTERNAL :: EOSS0


      CALL INIT_ERR_MSG("CHECK_SOLIDS_DENSITY")


! Initialize the flag for variable solids density.
      SOLVE_ROs = .FALSE.


      DO LC = 1, MMAX_LL


! Set the flag for variable solids density if any of the input parameters
! are specified.
         IF(INERT_SPECIES(LC) /= UNDEFINED_I) SOLVE_ROs(LC) = .TRUE.
         DO N=1, DIM_N_s
            IF(RO_Xs0(LC,N) /= UNDEFINED) SOLVE_ROs(LC) = .TRUE.
            IF(X_S0(LC,N) /= UNDEFINED) SOLVE_ROs(LC) = .TRUE.
         ENDDO


! Verify that one -and only one- solids density model is in use.
         IF(RO_S0(LC) == UNDEFINED .AND. .NOT.SOLVE_ROs(LC)) THEN
            WRITE(ERR_MSG, 1101) LC
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1101: No solids density information for phase ',  &
         I2,'.',/'Please correct the mfix.dat file.')

! Check if the constant solids phase density is physical.
         IF(RO_S0(LC) /= UNDEFINED .AND. SOLVE_ROs(LC)) THEN
            WRITE(ERR_MSG, 1100) LC
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1101 FORMAT('Error 1101: Conflicting solids density input specified ',&
         'for solids',/'phase ',I2,'. Constant solids density ',       &
         'specified (RO_s0) along with one',/'or more of the variable',&
         ' solids density parameters:',/'RO_Xs0, X_s0, INERT_SPECIES.',&
         /'Please correct the mfix.dat file.')

         ENDIF



! Check that there is sufficient information to calculate a variable
! solids density. The first round of checks only determine that a user
! has specified an input value. Checks on the values follows.
! Initialize the error flag.

! Check physical restrictions on variable solids density input.
         IF(SOLVE_ROs(LC)) THEN

! Verify that the inert species index is defined and in range.
            IF(INERT_SPECIES(LC) == UNDEFINED_I) THEN

               WRITE(ERR_MSG,1000) iVar('INERT_SPECIES',LC)
               CALL FLUSH_ERR_MSG
               FAILED = .TRUE.


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')


            ELSEIF(INERT_SPECIES(LC) < 1 .OR.                          &
               INERT_SPECIES(LC) > NMAX(LC)) THEN
               IF(DMP_LOG) THEN
                  IF(.NOT.WARNED_USR) THEN
                     WRITE(*,1305)
                     WRITE(UNIT_LOG,1305)
                     WARNED_USR = .TRUE.
                  ENDIF
                  WRITE(*,1307) LC, NMAX(LC)
                  WRITE(UNIT_LOG,1307) LC, NMAX(LC)
               ENDIF
               FAILED = .TRUE.
            ENDIF
            DO N=1, NMAX(LC)
               IF(RO_Xs0(LC,N) == UNDEFINED) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1305)
                        WRITE(UNIT_LOG,1305)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1308)'RO_Xs0',LC, N
                     WRITE(UNIT_LOG,1308)'RO_Xs0',LC, N
                  ENDIF
                  FAILED = .TRUE.
               ELSEIF(RO_Xs0(LC,N) < ZERO) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1305)
                        WRITE(UNIT_LOG,1305)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1309)'RO_Xs0', LC, N, RO_Xs0(LC,N)
                     WRITE(UNIT_LOG,1309)'RO_Xs0', LC, N, RO_Xs0(LC,N)
                  ENDIF
                  FAILED = .TRUE.
               ENDIF

               IF(X_s0(LC,N) == UNDEFINED) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1305)
                        WRITE(UNIT_LOG,1305)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1308)'X_s0',LC, N
                     WRITE(UNIT_LOG,1308)'X_s0',LC, N
                  ENDIF
                  FAILED = .TRUE.
               ELSEIF(X_s0(LC,N) < ZERO .OR. X_s0(LC,N) >= ONE) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1305)
                        WRITE(UNIT_LOG,1305)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1309)'X_s0', LC, N, X_s0(LC,N)
                     WRITE(UNIT_LOG,1309)'X_s0', LC, N, X_s0(LC,N)
                  ENDIF
                  FAILED = .TRUE.
               ENDIF
            ENDDO
            IF(FAILED) THEN
               IF(DMP_LOG) THEN
                  WRITE(*,9999)
                  WRITE(UNIT_LOG,9999)
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF



! Check physical restrictions on constant solids density input.
         ELSE
            IF(RO_S0(LC) <= ZERO) THEN
               WRITE(ERR_MSG,1101) iVar('RO_s0',LC), iVal(RO_s0(LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF


         ENDIF


! Verify that the inert species mass fraction is non-zero.
            IF(X_s0(LC,INERT_SPECIES(LC)) == ZERO) THEN
               IF(DMP_LOG) THEN
                  WRITE(*,1310) LC, INERT_SPECIES(LC)
                  WRITE(UNIT_LOG,1310) LC, INERT_SPECIES(LC)
                  WRITE(*,9999)
                  WRITE(UNIT_LOG,9999)
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF

! All of the information for variable solids density has been verified
! as of this point. Calculate and store the baseline density.
            BASE_ROs(LC) = EOSS0(LC)

         ENDIF ! SOLVE_RO_s
      ENDDO



! Ensure no values are given for undefined solids phases.
      FAILED = .FALSE.
      WARNED_USR = .FALSE.
      DO LC = MMAX_LL + 1, DIM_M
         IF (RO_S0(LC) /= UNDEFINED) THEN
            IF(DMP_LOG) THEN
               IF(.NOT.WARNED_USR) THEN
                  WRITE(*,1312) MMAX_LL
                  WRITE(UNIT_LOG,1312) MMAX_LL
                  WARNED_USR = .TRUE.
               ENDIF
               WRITE(*,1313) LC, RO_S0(LC)
               WRITE(UNIT_LOG,1313) LC, RO_S0(LC)
            ENDIF
            FAILED = .TRUE.
         ENDIF
         IF (INERT_SPECIES(LC) /= UNDEFINED_I) THEN
            IF(DMP_LOG) THEN
               IF(.NOT.WARNED_USR) THEN
                  WRITE(*,1312) MMAX_LL
                  WRITE(UNIT_LOG,1312)
                  WARNED_USR = .TRUE.
               ENDIF
               WRITE(*,1314) LC, INERT_SPECIES(LC)
               WRITE(UNIT_LOG,1314) LC, INERT_SPECIES(LC)
            ENDIF
            FAILED = .TRUE.
         ENDIF
         DO N=1, DIM_N_s
            IF(RO_Xs0(LC,N) /= UNDEFINED .AND. RO_Xs0(LC,N) /=ZERO) THEN
               IF(DMP_LOG) THEN
                  IF(.NOT.WARNED_USR) THEN
                     WRITE(*,1312) mmax_ll         ! Screen Error
                     WRITE(UNIT_LOG,1312) MMAX_LL  ! Log Error
                     WARNED_USR = .TRUE.
                  ENDIF
                  WRITE(*,1315) 'RO_XS0',LC,N, RO_XS0(LC,N)
                  WRITE(UNIT_LOG,1315) 'RO_XS0',LC,N, RO_XS0(LC,N)
               ENDIF
               FAILED = .TRUE.
            ENDIF
            IF(X_S0(LC,N) /= UNDEFINED .AND. X_S0(LC,N) /= ZERO) THEN
               IF(DMP_LOG) THEN
                  IF(.NOT.WARNED_USR) THEN
                     WRITE(*,1312) MMAX_LL         ! Screen Error
                     WRITE(UNIT_LOG,1312) MMAX_LL  ! Log Error
                     WARNED_USR = .TRUE.
                  ENDIF
                  WRITE(*,1315) 'X_S0',LC,N, X_S0(LC,N)
                  WRITE(UNIT_LOG,1315) 'X_S0',LC,N, X_S0(LC,N)
               ENDIF
               FAILED = .TRUE.
            ENDIF
         ENDDO
         IF(FAILED)THEN
            IF(DMP_LOG) THEN
               WRITE(*,9999)
               WRITE(UNIT_LOG,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDDO 





! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN  

 1300 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1300:'      &
         ' No solids density information for phase ',I2,'.')


 1301 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1301:'      &
         ' Unphysical solids density (RO_s0) for phase ',I2,'.')

 1302 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1310:'      &
         ' Conflicting solids phase density values defined:',I2,'.',2/,&
         '  > Constant density provided: RO_s0(',I2,') = ',g11.5)

 1303 FORMAT('  > Variable density parameter provided: ',A,'(',I2,     &
         ') = ',I3)

 1304 FORMAT('  > Variable density parameter provided: ',A,'(',I2,     &
         ',',I3,') = ',g11.5)


 1305 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1305:'      &
         ' One or more invalid variable solid density parameters:',/)
  
 1306 FORMAT(3x,'> INERT_SPECIES(',I2,') is UNDEFINED.')

 1307 FORMAT(3x,'> INERT_SPECIES(',I2,') is out of range :: [1,',I3,']')

 1308 FORMAT(3x,'> ',A,'(',I2,',',I3,') is UNDEFINED.')

 1309 FORMAT(3x,'> ',A,'(',I2,',',I3,') is unphysical :: ',g11.5)

 1310 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1310:'      &
         ' Invalid baseline inert speices mass fraction.',/            &
         ' The inert spcies mass fraction must be greater than zero.',/&
         ' Phase ',I2,' Inert Species: ',I3,'  X_s0 = 0.0')

 1312 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1305:'      &
         ' One or more parameters defined with solid phase index',/    &
         ' exceeding MMAX = ',I2,/)

 1313 FORMAT('  > RO_s0(',I2,') = ',g11.5)

 1314 FORMAT('  > INERT_SPECIES(',I2,') = ',I3)

 1315 FORMAT('  > ',A,'(',I2,',',I3,') = ',g11.5)


 1400 FORMAT(1X,/,1X,'RO_S0(',I2,') = ',G12.5,/,1X,'MMAX in mfix = ',I2,/) 



 9999 FORMAT(/' Please refer to the Readme file on the required input',&
         ' and make',/' the necessary corrections to the data file.',  &
         /1X,70('*')//)

      END SUBROUTINE CHECK_SOLIDS_DENSITY
