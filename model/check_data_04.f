!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_04                                           C
!  Purpose: check the solid phase input section                        C
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
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_04 

!-----------------------------------------------
! Modules
!-----------------------------------------------
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

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: LC, N
      CHARACTER*85 LONG_STRING      

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

! Flag indicating that the thermochemical database header was output 
! to the screen. (Miminize messages)
      LOGICAL thermoHeader

! Flag that a no fatal error was detected.
      LOGICAL :: PASSED
! Flag that a fatal error was detected.
      LOGICAL :: FAILED

! Solids phase density calculation
      DOUBLE PRECISION, EXTERNAL :: EOSS0
!-----------------------------------------------


      IF ((DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_HYBRID) & 
      .OR. (DISCRETE_ELEMENT .AND. MPPIC) ) THEN
! override possible user settings on the following continuum flags

! Set close_packed to true to prevent possible issues stemming from the
! pressure correction equation.  Specifically, if closed_packed is false
! then a mixture pressure correction equation is invoked and this is not
! correctly setup for DEM.  To do so would require ensuring that
! 1) the solids phase continuum quantities used in these equations are
!    correctly set based on their DEM counterparts and 
! 2) the pressure correction coefficients for such solids phases are 
!    also calculated (currently these calculations are turned off 
!    when using DEM)
         CLOSE_PACKED(:) = .TRUE. 
         MOMENTUM_X_EQ(1:DIM_M) = .FALSE.
         MOMENTUM_Y_EQ(1:DIM_M) = .FALSE.
         MOMENTUM_Z_EQ(1:DIM_M) = .FALSE. 

! Check continuum solids phase species data. Clear anything that was
! specified and warn the user.
         CALL CHECK_04_TFM_DEM()

         RETURN
      ENDIF

! Check MMAX
      IF (MMAX<0 .OR. MMAX>DIM_M) THEN 
         CALL ERROR_ROUTINE ('check_data_04', &
            'MMAX not specified or unphysical', 0, 2) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000) MMAX, DIM_M 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! Check D_p0
! Only need to check for real phases (for GHD theory)      
      DO LC = 1, SMAX 
         IF (D_P0(LC)<ZERO .OR. D_P0(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'D_p0 not specified or unphysical', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1100) LC, D_P0(LC) 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 
      DO LC = SMAX + 1, DIM_M 
         IF (D_P0(LC) /= UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_data_04', &
               'too many D_p0 values specified', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1200) LC, D_P0(LC), MMAX
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 

      IF (MODEL_B) THEN 
         DO LC = 1, MMAX 
            IF (.NOT.CLOSE_PACKED(LC)) THEN 
               CALL ERROR_ROUTINE ('check_data_04', ' ', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1420) LC
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDIF 

      IF (MMAX > 0) THEN 
         IF (C_E == UNDEFINED) CALL ERROR_ROUTINE ('check_data_04', &
            'Coefficient of restitution (C_e) not specified', 1, 1) 
         IF (C_F==UNDEFINED .AND. MMAX>=2 .AND. KT_TYPE .EQ. UNDEFINED_C) &
            CALL ERROR_ROUTINE ('check_data_04',&
               'Coefficient of friction (C_f) not specified',1,1) 

         IF ((FRICTION .OR. SCHAEFFER) .AND. (PHI == UNDEFINED)) &
            CALL ERROR_ROUTINE ('check_data_04', &
               'Angle of internal friction (Phi) not specified',1,1)
            LONG_STRING = 'Angle of wall-particle friction (Phi_w) &
               &not specified'
         IF ((FRICTION .OR. JENKINS) .AND. (PHI_W == UNDEFINED)) &
             CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)
      ENDIF 

! Check EP_star
      IF (MMAX > 0) THEN 
         IF (EP_STAR<ZERO .OR. EP_STAR>ONE) &
            CALL ERROR_ROUTINE ('check_data_04',&
            'Value of EP_star is unphysical', 1, 1) 
      ENDIF 

! Yu_Standish and Fedors_Landel correlations are used with more than one solids phase
      LONG_STRING = 'MMAX must be >= 2 for Yu_Standish or &
         &Fedors_Landel correlations'
      IF (SMAX < 2 .AND. (YU_STANDISH .OR. FEDORS_LANDEL)) &
         CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)

! Fedors_Landel correlations is limited to a binary mixture of powders
      IF (SMAX > 2 .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04', &
            'Fedors_Landel requires MMAX = 2 ', 1, 1)

! Must choose between Yu_Standish and Fedors_Landel correlations, can't use both.
      LONG_STRING = 'Cannot use both Yu_Standish and Fedors_Landel &
         &correlations'
      IF (YU_STANDISH .AND. FEDORS_LANDEL) &
         CALL ERROR_ROUTINE ('check_data_04',TRIM(LONG_STRING),1,1)
      IF (SIGM_BLEND) THEN
         TANH_BLEND = .FALSE. ! Setting tanh blending to be false
      ENDIF

! Define restitution coefficient matrix      
      DO LC = 1, SMAX 
         DO N = 1, SMAX
            IF(r_p(LC,N) == UNDEFINED) r_p(LC,N) = C_e
            r_p(N,LC) = r_p(LC,N) ! just need to define r_p(1,2) and r_p(2,1) will be set.
         ENDDO
      ENDDO

! Check NMAX_s: The number of solids phase species should be defined
! in NMAX_s. However, for legacy access, (USE_RRATES = .T.) the number
! of species will appear in NMAX(m). 
! --> NMAX(M) must be defined if species_eq(M) true. 
! --> Only need to check species for real phases (for GHD theory) 
      DO LC = 1, MMAX
! If the legacy connection for reaction rates is not being used, check
! to see if the number of species was provided under NMAX(M).
         IF(.NOT.USE_RRATES) THEN
            IF(NMAX_s(LC) == UNDEFINED_I) THEN
               IF(NMAX(LC) /= UNDEFINED_I) THEN
! The number of solids phase species was given in NMAX. Warn the user
! to use the correct variable. Copy the old variable entry to the new
! variable name for later pre-processing and continue.
                  IF(DMP_LOG) THEN
                     WRITE(*,1056)LC
                     WRITE(UNIT_LOG,1056)LC
                  ENDIF
                  NMAX_s(LC) = NMAX(LC)
               ENDIF
! If for whatever reason, the number of species are given in both
! variables, make sure they match.
            ELSEIF(NMAX(LC) /= UNDEFINED_I .AND. &
               NMAX_s(LC) /= UNDEFINED_I) THEN
! The two varaibles do not match. Flag the error and exit.
               IF(NMAX(LC) /= NMAX_s(LC)) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1055)LC
                     WRITE(UNIT_LOG,1055)LC
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF
            ELSE
! Copy the new keyword entry into the runtime variable.
               NMAX(LC) = NMAX_s(LC)
            ENDIF
         ENDIF
! If the number of solids phase species wasn't given in either keyword,
! Give the error and exit.
         IF (NMAX(LC) == UNDEFINED_I) THEN 
            IF (SPECIES_EQ(LC)) THEN 
               CALL ERROR_ROUTINE ('check_data_04', &
                  'Number of species not specified', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1045) LC 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ELSE 
               NMAX(LC) = 1 
            ENDIF 
         ENDIF 
         IF (NMAX(LC) > DIM_N_S) THEN 
            CALL ERROR_ROUTINE ('check_data_04', 'NMAX is too large', 0, 2) 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1050) LC, NMAX, DIM_N_S 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      ENDDO 

! Flag indicating if the user was already warned.
      WARNED_USR = .FALSE.
! Flag indicating the search header was already written.
      thermoHeader = .FALSE.
! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
      EEQ_CPS = .FALSE.
      IF(ENERGY_EQ .AND. C_PS0 == UNDEFINED) EEQ_CPS = .TRUE.

! Check MW_s if solids species are present    
      DO LC = 1, SMAX
! Initialize flag indicating the database was read for a species.
         rDatabase(LC,:) = .FALSE.
         DO N = 1, NMAX(LC)
! Legacy code may use 'hard coded' functions for specific heat and
! therefore should avoid the follow check.
            IF(.NOT.USE_RRATES) THEN
! Flag that the solids phase species equations are solved and the 
! molecular weight for a species are not given in the data file.
               SEQ_MWs = .FALSE.
               IF(SPECIES_EQ(LC) .AND. MW_S(LC,N) == UNDEFINED)        &
                  SEQ_MWs = .TRUE.

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
                     IF(EEQ_CPS .AND. DMP_LOG) THEN
                        WRITE(*,1058)
                        WRITE(UNIT_LOG,1058)
                     ENDIF
                     IF(SEQ_MWs .AND. DMP_LOG) THEN
                        WRITE(*,1059)LC
                        WRITE(UNIT_LOG,1059)LC
                     ENDIF
! Set a flag to prevent the user from getting the same message over
! and over again.
                     WARNED_USR = .TRUE.
                  ENDIF
! Flag that the species name is not provided.
                  IF(SPECIES_s(LC,N) == UNDEFINED_C) THEN
                     IF(DMP_LOG) THEN
                        WRITE(*,1060) LC, N
                        WRITE(UNIT_LOG,1060) LC, N
                     ENDIF
                     CALL MFIX_EXIT(myPE)
                  ENDIF
! Read the database.
                  IF(DMP_LOG) THEN
                     IF(.NOT.thermoHeader) THEN
                        WRITE(*,1061) LC
                        WRITE(UNIT_LOG,1061) LC
                        thermoHeader = .TRUE.
                     ENDIF
                     WRITE(*,1062) N, trim(SPECIES_s(LC,N))
                     WRITE(UNIT_LOG,1062) N, trim(SPECIES_s(LC,N))
                  ENDIF
                  CALL READ_DATABASE('TFM', LC, N, SPECIES_s(LC,N),   &
                     MW_S(LC,N))
! Flag variable to stating that the database was read.
                  rDatabase(LC,N) = .TRUE.
! Flag the legacy variable to prevent re-reading the database.
                  DATABASE_READ = .TRUE.
               ENDIF
            ELSE
! This is a legacy check.
               IF (MW_S(LC,N) == UNDEFINED) THEN 
                  CALL ERROR_ROUTINE ('check_data_04', &
                     'Species molecular weight undefined', 0, 2) 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
               ENDIF
            ENDIF ! USE_RRATES
         ENDDO ! Loop over species
         DO N = NMAX(LC) + 1, DIM_N_S 
            IF (MW_S(LC,N) /= UNDEFINED) THEN 
               CALL ERROR_ROUTINE ('check_data_04', &
                  'MW_s defined for N > NMAX(m)', 0, 2) 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1410) LC, N 
               CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
            ENDIF 
         ENDDO 
      ENDDO ! Loop over solids phases

! Check RO_s/0
      DO LC = 1, SMAX
! Check if the constant solids phase density is physical.
         IF(RO_S0(LC) /= UNDEFINED) THEN
            SOLVE_ROs(LC) = .FALSE.
! Check if the constant solids phase density is physical.
            IF(RO_S0(LC) <= ZERO) THEN
               IF(DMP_LOG) THEN
                  WRITE(*,1301) LC;  WRITE(*,9999)
                  WRITE(UNIT_LOG,1301) LC;  WRITE(UNIT_LOG,9999)
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF
! Ensure that no variable solids density is given.
            FAILED = .FALSE.
            WARNED_USR = .FALSE.
            IF (INERT_SPECIES(LC) /= UNDEFINED_I) THEN
               IF(DMP_LOG) THEN
                  IF(.NOT.WARNED_USR) THEN
                     WRITE(*,1302) LC, LC, RO_S0(LC)
                     WRITE(UNIT_LOG,1302) LC, LC, RO_S0(LC)
                     WARNED_USR = .TRUE.
                  ENDIF
                  WRITE(*,1303) 'INERT_SPECIES',LC,INERT_SPECIES(LC)
                  WRITE(UNIT_LOG,1303) 'INERT_SPECIES', LC,            &
                     INERT_SPECIES(LC)
               ENDIF
               FAILED = .TRUE.
            ENDIF
            DO N=1, DIM_N_s
               IF(RO_Xs0(LC,N) /= UNDEFINED .AND.                      &
                  RO_Xs0(LC,N) /=ZERO) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1302) LC, LC, RO_S0(LC)
                        WRITE(UNIT_LOG,1303) LC, LC, RO_S0(LC)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1304) 'RO_XS0',LC,N, RO_XS0(LC,N)
                     WRITE(UNIT_LOG,1304) 'RO_XS0',LC,N, RO_XS0(LC,N)
                  ENDIF
                  FAILED = .TRUE.
               ENDIF
               IF(X_S0(LC,N) /= UNDEFINED .AND. X_S0(LC,N) /= ZERO) THEN
                  IF(DMP_LOG) THEN
                     IF(.NOT.WARNED_USR) THEN
                        WRITE(*,1302) LC, LC, RO_S0(LC)
                        WRITE(UNIT_LOG,1302) LC, LC, RO_S0(LC)
                        WARNED_USR = .TRUE.
                     ENDIF
                     WRITE(*,1304) 'X_S0',LC,N, X_S0(LC,N)
                     WRITE(UNIT_LOG,1304) 'X_S0',LC, N, X_S0(LC,N)
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

! Check if any variable density parameters are provided.
         ELSE
            IF(INERT_SPECIES(LC) /= UNDEFINED_I) SOLVE_ROs(LC) = .TRUE.
            DO N=1, DIM_N_s
               IF(RO_Xs0(LC,N) /= UNDEFINED) SOLVE_ROs(LC) = .TRUE.
               IF(X_S0(LC,N) /= UNDEFINED) SOLVE_ROs(LC) = .TRUE.
            ENDDO
         ENDIF

! No solids density was provided.
         IF(RO_S0(LC) == UNDEFINED .AND. .NOT.SOLVE_ROs(LC)) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1301) LC; WRITE(*,9999)
               WRITE(UNIT_LOG,1301) LC; WRITE(UNIT_LOG,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

! Check that there is sufficient information to calculate a variable
! solids density. The first round of checks only determine that a user
! has specified an input value. Checks on the values follows.
! Initialize the error flag.
         IF(SOLVE_ROs(LC)) THEN
            FAILED = .FALSE.
            WARNED_USR = .FALSE.
! Verify that the inert species index is defined and in range.
            IF(INERT_SPECIES(LC) == UNDEFINED_I) THEN
               IF(DMP_LOG) THEN
                  IF(.NOT.WARNED_USR) THEN
                     WRITE(*,1305)
                     WRITE(UNIT_LOG,1305)
                     WARNED_USR = .TRUE.
                  ENDIF
                  WRITE(*,1306) LC, NMAX(LC)
                  WRITE(UNIT_LOG,1306) LC, NMAX(LC)
               ENDIF
               FAILED = .TRUE.
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
      DO LC = SMAX + 1, DIM_M
         IF (RO_S0(LC) /= UNDEFINED) THEN
            IF(DMP_LOG) THEN
               IF(.NOT.WARNED_USR) THEN
                  WRITE(*,1312) SMAX
                  WRITE(UNIT_LOG,1312) SMAX
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
                  WRITE(*,1312) SMAX
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
                     WRITE(*,1312) SMAX         ! Screen Error
                     WRITE(UNIT_LOG,1312) SMAX  ! Log Error
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
                     WRITE(*,1312) SMAX         ! Screen Error
                     WRITE(UNIT_LOG,1312) SMAX  ! Log Error
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






! CHECK MU_s0
      IF (MU_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'MU_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1500) MU_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK K_s0
      IF (K_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'K_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1510) K_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      
! CHECK C_ps0
      IF (C_PS0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04', &
            'C_ps0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1520) C_PS0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

! CHECK DIF_s0
      IF (DIF_S0 < ZERO) THEN 
         CALL ERROR_ROUTINE ('CHECK_DATA_04',&
            'DIF_s0 value is unphysical', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1530) DIF_S0 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 


      RETURN  

 1000 FORMAT(1X,/,1X,'MMAX        in  mfix.dat = ',I6,/,1X,&
         'DIM_M in  param.inc  = ',I6,/) 
 1045 FORMAT(1X,/,1X,'NMAX is not specified for solids phase',I2) 
 1050 FORMAT(1X,/,1X,'NMAX(',I2,')   in  mfix.dat = ',I6,/,1X,&
         'DIM_N_s in  param.inc  = ',I6,/) 

 1055 FORMAT(/1X,70('*')/' From: CHECK_DATA_04',/                      &
         ' Message: NMAX_s and NMAX are both given for solids phase ', &
         I2,' in the',/' data file and do not match. NMAX is a legacy',&
         ' variable and is not',/' required. Please correct the data', &
         ' file.',/1X,70('*')/)

 1056 FORMAT(/1X,70('*')/' From: CHECK_DATA_04',/                      &
         ' Message: NMAX is specified for solids phase ',I2,'. This',  &
         ' is a legacy',/' variable, and NMAX_s should be used.',      &
         ' Copying NMAX to NMAX_s.',/1X,70('*')/)

 1058 FORMAT(/1X,70('*')/' From: CHECK_DATA_04',/                      &
         ' Message: The energy equations are being solved (ENERGY_EQ)',&
         ', and the',/' specified constant solids specific heat is',   &
         ' undefined (C_PS0). Thus,',/' the thermochemical database',  &
         ' will be used to gather specific heat data',/' on the',      &
         ' individual soids phase species.',/1X,70('*')/)

 1059 FORMAT(/1X,70('*')/' From: CHECK_DATA_04',/                      &
         ' Message: Solids phase ',I2,' species equations are being',  &
         ' solved, and one',/' or more species molecular weights are', &
         ' undefined. Thus, the thermo-',/' chemical database will be',&
         ' used to gather molecular weight data on the',/' solids',    &
         ' phase speicies.',/1X,70('*')/)

 1060 FORMAT(/1X,70('*')/'  From: CHECK_DATA_04',/                     &
         ' Message: Solids phase ',I2,' species ',I2,' name',          &
         ' (SPECIES_s) is undefined.',/' Please correct the data file.'&
         ,/1X,70('*')/)

 1061 FORMAT(/'  Searching thermochemical databases for solids phase ',&
         I2,' species data')

 1062 FORMAT(/2x,'>',I3,': Species: ',A)


 1100 FORMAT(1X,/,1X,'D_p0(',I2,') in mfix.dat = ',G12.5) 
 1200 FORMAT(1X,/,1X,'D_p0(',I2,') = ',G12.5,/,1X,'MMAX in mfix = ',I2,/) 

 1300 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1300:'      &
         ' No solids density information for phase ',I2,'.')

 1301 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1301:'      &
         ' Unphysical solids density (RO_s0) for phase ',I2,'.')

 1302 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1310:'      &
         ' Conflicting solids phase density values defined:',I2,'.',2/,&
         '  > Constant density provided: RO_s0(',I2,') = ',g11.5)

 1303 FORMAT('  > Variable density parameter provided: ',A,'(',I2,     &
         ') = ',I3)

 1304 FORMAT('  > Variable density parameter provided: ',A,'(',I2,     &
         ',',I3,') = ',g11.5)


 1305 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1305:'      &
         ' One or more invalid variable solid density parameters:',/)
  
 1306 FORMAT(3x,'> INERT_SPECIES(',I2,') is UNDEFINED.')

 1307 FORMAT(3x,'> INERT_SPECIES(',I2,') is out of range :: [1,',I3,']')

 1308 FORMAT(3x,'> ',A,'(',I2,',',I3,') is UNDEFINED.')

 1309 FORMAT(3x,'> ',A,'(',I2,',',I3,') is unphysical :: ',g11.5)

 1310 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1310:'      &
         ' Invalid baseline inert speices mass fraction.',/            &
         ' The inert spcies mass fraction must be greater than zero.',/&
         ' Phase ',I2,' Inert Species: ',I3,'  X_s0 = 0.0')

 1312 FORMAT(//1X,70('*')/' From: CHECK_DATA_04',/,' Error 1305:'      &
         ' One or more parameters defined with solid phase index',/    &
         ' exceeding MMAX = ',I2,/)

 1313 FORMAT('  > RO_s0(',I2,') = ',g11.5)

 1314 FORMAT('  > INERT_SPECIES(',I2,') = ',I3)

 1315 FORMAT('  > ',A,'(',I2,',',I3,') = ',g11.5)


 1400 FORMAT(1X,/,1X,'RO_S0(',I2,') = ',G12.5,/,1X,'MMAX in mfix = ',I2,/) 


 1410 FORMAT(1X,/,1X,'Solids phase = ',I2,'   Species = ',I3) 
 1420 FORMAT(1X,/,1X,'Solids phase = ',I2,' is not Close_Packed.',/,&
         ' With Model B all solids phases should have that property')

 1500 FORMAT(1X,/,1X,'MU_s0   in mfix.dat = ',G12.5)
 1510 FORMAT(1X,/,1X,'K_s0   in mfix.dat = ',G12.5)
 1520 FORMAT(1X,/,1X,'C_ps0   in mfix.dat = ',G12.5)
 1530 FORMAT(1X,/,1X,'DIF_s0   in mfix.dat = ',G12.5)


 9999 FORMAT(/' Please refer to the Readme file on the required input',&
         ' and make',/' the necessary corrections to the data file.',  &
         /1X,70('*')//)


      END SUBROUTINE CHECK_DATA_04 


!----------------------------------------------------------------------!
! Subroutine: CHECK_04_TFM_DEM                                         !
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
      SUBROUTINE CHECK_04_TFM_DEM()

      USE compar
      USE funits 
      USE param 
      USE param1 
      USE physprop
      USE rxns

      IMPLICIT NONE

! Flag indicating if the user was warned.
      LOGICAL WARNED_USR0, WARNED_USR1
! Loop indices
      INTEGER lM, lN ! Phase, Species

! Flag indicating if the user was warned.
      WARNED_USR0 = .FALSE.
      DO lM = 1, DIM_M
         WARNED_USR1 = .FALSE.
! Check if continuum phase variables are specified.
         IF(NMAX_s(lM) /= UNDEFINED_I .OR. &
            NMAX(lM) /= UNDEFINED_I) THEN
            IF(.NOT.WARNED_USR0 .AND. DMP_LOG) THEN
               WRITE(*,1001)
               WRITE(*,1101) lM
               WRITE(UNIT_LOG,1001)
               WRITE(UNIT_LOG,1101) lM
            ELSEIF(.NOT.WARNED_USR1 .AND. DMP_LOG) THEN
               WRITE(*,1101) lM
               WRITE(UNIT_LOG,1101) lM
            ENDIF
! Clear the entry.
            NMAX(lM) = 0
            NMAX_s(lM) = 0
! Set message flags.
            WARNED_USR0 = .TRUE.
            WARNED_USR1 = .TRUE.
         ELSE
! In a DEM simulation, MMAX can be used to specify the number of
! discrete solids phases. These values are set to zero for to aid 
! routines that loop over solids phase species. 
            NMAX(lM) = 0
            NMAX_s(lM) = 0
         ENDIF

         DO lN = 1, DIM_N_S
! Verify that a species name was not provided.
            IF(SPECIES_s(lM,lN) /= UNDEFINED_C) THEN
               IF(.NOT.WARNED_USR0 .AND. DMP_LOG) THEN
                  WRITE(*,1001)
                  WRITE(*,1101) lM
                  WRITE(UNIT_LOG,1001)
                  WRITE(UNIT_LOG,1101) lM
               ELSEIF(.NOT.WARNED_USR1 .AND. DMP_LOG) THEN
                  WRITE(*,1101) lM
                  WRITE(UNIT_LOG,1101) lM
               ENDIF
! Clear the entry.
               SPECIES_s(lM,lN) = UNDEFINED_C
! Set message flags.
               WARNED_USR0 = .TRUE.
               WARNED_USR1 = .TRUE.
            ENDIF
! Verify that a species alias was not provided.
            IF(SPECIES_ALIAS_s(lM,lN) /= UNDEFINED_C) THEN
               IF(.NOT.WARNED_USR0 .AND. DMP_LOG) THEN
                  WRITE(*,1001)
                  WRITE(*,1101) lM
                  WRITE(UNIT_LOG,1001)
                  WRITE(UNIT_LOG,1101) lM
               ELSEIF(.NOT.WARNED_USR1 .AND. DMP_LOG) THEN
                  WRITE(*,1101) lM
                  WRITE(UNIT_LOG,1101) lM
               ENDIF
! Clear the entry.
               SPECIES_ALIAS_s(lM,lN) = UNDEFINED_C
! Set message flags.
               WARNED_USR0 = .TRUE.
               WARNED_USR1 = .TRUE.
           ENDIF
         ENDDO
      ENDDO

      IF(WARNED_USR0 .AND. DMP_LOG) THEN
         WRITE(*,1201)
         WRITE(UNIT_LOG,1201)
      ENDIF

 1001 FORMAT(/1X,70('*')/' From: CHECK_DATA_04 --> CHECK_04_TFM_DEM',/ &
         ' Warning 1001: Inconsistent input in the data file. A',      &
         ' discrete element',/' model (DEM/MPPIC) is used and species',&
         ' information for the following',/' TFM solids phases are',   &
         ' specified.'/)

 1101 FORMAT(' TFM Solids Phase: ',I2)

 1201 FORMAT(/' The species information for the above solids phases',  &
         ' is being cleared.',/' Please refer to the Readme file for', &
         ' information on specifying discrete',/' phase species.',/    &
         1X,70('*')/)

      END SUBROUTINE CHECK_04_TFM_DEM
