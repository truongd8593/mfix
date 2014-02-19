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
!    > CHECK_SOLIDS_CONTINUUM :: TRM solids phase model)               !
!    > CHECK_SOLIDS_DES       :: DEM solids phase model)               !
!    > CHECK_SOLIDS_MPPIC     :: MPPIC solids phase model)             !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_ALL


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Use legacy rrates implementation.
      use rxns, only: USE_RRATES
! Number of continuum solids phases.
      use physprop, only: SMAX
! Number of disctete (DEM/MPPIC) solids.
      use discretelement, only: DES_MMAX
! User specified: Constant solids specific heat.
      use physprop, only: C_PS0
! User specified: Constant solids thermal conductivity.
      use physprop, only: K_S0
! User specified: Initial solids diameter.
      use physprop, only: D_P0


! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M
! Parameter constants
      use param1, only: UNDEFINED, ZERO


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: M
! Total number of all solids
      INTEGER :: MMAX_L

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_ALL")

! Set the number of solids phases to be checked.
      MMAX_L = SMAX + DES_MMAX

! Check D_p0
      DO M = 1, MMAX_L
         IF(D_P0(M) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('D_p0',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(D_P0(M) <= ZERO)THEN
            WRITE(ERR_MSG, 1001) trim(iVar('D_p0',M)), iVal(D_P0(M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      DO M = MMAX_L+1, DIM_M
         IF(D_P0(M) /= UNDEFINED)THEN
            WRITE(ERR_MSG,1002) trim(iVar('D_p0',M))
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

! Check solids drag model selection.
      CALL CHECK_SOLIDS_DRAG

! Check the solids density input parameters.
      CALL CHECK_SOLIDS_DENSITY(MMAX_L)

! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.'&
         ,/'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_ALL



!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES                                     !
! Purpose: Check solids species input.                                 !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DRAG

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve energy equations
      use run, only: ENERGY_EQ
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! User specifed drag type, as string
      use run, only: DRAG_TYPE
! User drag type as ENUM
      use run, only: DRAG_TYPE_ENUM
! Possible DRAG_TYPE_ENUM values:
      use run, only: SYAM_OBRIEN
      use run, only: GIDASPOW
      use run, only: GIDASPOW_PCF
      use run, only: GIDASPOW_BLEND
      use run, only: GIDASPOW_BLEND_PCF
      use run, only: WEN_YU
      use run, only: WEN_YU_PCF
      use run, only: KOCH_HILL
      use run, only: KOCH_HILL_PCF
      use run, only: BVK
      use run, only: HYS

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phase species.
      USE param, only: DIM_N_s
! Parameter constants.
      use param1, only: UNDEFINED, UNDEFINED_I, UNDEFINED_C

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DRAG")


      SELECT CASE(trim(adjustl(DRAG_TYPE)))

      CASE ('SYAM_OBRIEN')
         DRAG_TYPE_ENUM = SYAM_OBRIEN

      CASE ('GIDASPOW')
         DRAG_TYPE_ENUM = GIDASPOW

      CASE ('GIDASPOW_PCF')
         DRAG_TYPE_ENUM = GIDASPOW_PCF

      CASE ('GIDASPOW_BLEND')
         DRAG_TYPE_ENUM = GIDASPOW_BLEND

      CASE ('GIDASPOW_BLEND_PCF')
         DRAG_TYPE_ENUM = GIDASPOW_BLEND_PCF

      CASE ('WEN_YU')
         DRAG_TYPE_ENUM = WEN_YU

      CASE ('WEN_YU_PCF')
         DRAG_TYPE_ENUM = WEN_YU_PCF

      CASE ('KOCH_HILL')
         DRAG_TYPE_ENUM = KOCH_HILL

      CASE ('KOCH_HILL_PCF')
         DRAG_TYPE_ENUM = KOCH_HILL_PCF

      CASE ('BVK')
         DRAG_TYPE_ENUM = BVK

      CASE ('HYS')
         DRAG_TYPE_ENUM = HYS

      CASE DEFAULT
         
      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_DRAG


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES                                     !
! Purpose: Check solids species input.                                 !
!                                                                      !
! Author: J. Musser                                  Date: 07-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve energy equations
      use run, only: ENERGY_EQ
! Flag: Solve species equations
      use run, only: SPECIES_EQ

! Flag: Database for phase X was read for species Y
      use rxns, only: rDatabase
! Solids phase species database names.
      use rxns, only: SPECIES_s
! Solids phase molecular weights.
      use physprop, only: MW_s
! Number of solids phase species.
      use physprop, only: NMAX, NMAX_s
! User specified: Constant solids specific heat
      use physprop, only: C_PS0

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phase species.
      USE param, only: DIM_N_s
! Parameter constants.
      use param1, only: UNDEFINED, UNDEFINED_I, UNDEFINED_C

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Subroutine Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//
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

! Loop counters.
      INTEGER :: M, N


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_SPECIES")


! Reconcile the new species input method with the legacy input method.
      DO M=1, MMAX_LL
         IF(SPECIES_EQ(M)) THEN
            IF(NMAX_s(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,1000) iVar('NMAX_s',M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX_s(M) > DIM_N_S) THEN
               WRITE(ERR_MSG,1001) trim(iVar('NMAX_s',M)),            &
                  trim(iVal(NMAX_s(M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSE
               NMAX(M) = NMAX_s(M)
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            NMAX(M) = merge(1, NMAX_s(M), NMAX_s(M) == UNDEFINED_I)
         ENDIF

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
         EEQ_CPS = (ENERGY_EQ .AND. C_PS0 == UNDEFINED)
         IF(EEQ_CPS)THEN
            WRITE(ERR_MSG,2000)
            CALL FLUSH_ERR_MSG
         ENDIF

 2000 FORMAT('Message: 2000 The energy equations are being solved ',   &
         '(ENERGY_EQ) and',/'the constant solids specific heat is ',   &
         'undefined (C_PS0). Thus, the',/'thermochemical database ',   &
         'will be used to gather specific heat data on',/'the ',       &
         'individual gas phase species.')

         SEQ_MWs = .FALSE.
         DO N=1,NMAX(M)
            IF(MW_s(M,N) == UNDEFINED) THEN
               IF(SPECIES_EQ(M)) SEQ_MWs = .TRUE.
            ENDIF
         ENDDO

         IF(SEQ_MWs) THEN
            WRITE(ERR_MSG, 2001) M
            CALL FLUSH_ERR_MSG
         ENDIF

 2001 FORMAT('Message 2001: One or more species molecular weights are',&
         ' undefined and',/'the solids phase species equations are ',  &
         'solved (SOLVE_EQ(',I2,')). The',/'thermochemical database ', &
         'will be used in an attempt to gather missing',/'molecular ', &
         'weight data.')

! Initialize flag indicating the database was read for a species.
         rDatabase(M,:) = .FALSE.

         IF(EEQ_CPS .OR. SEQ_MWs) THEN

            WRITE(ERR_MSG, 3000) M
            CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 3000 FORMAT('Message 3000: Searching thermochemical databases for ',  &
         'solids phase',I3,/'species data.',/'  ')

            DO N = 1, NMAX(M)

! Notify the user of the reason the thermochemical database is used.
               IF(EEQ_CPS .OR. MW_s(M,N) == UNDEFINED) THEN

! Flag that the species name is not provided.
                  IF(SPECIES_s(M,N) == UNDEFINED_C) THEN
                     WRITE(ERR_MSG,1000) trim(iVar('SPECIES_s',M,N))
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ENDIF
! Update the log files.
                  WRITE(ERR_MSG, 3001) N, trim(SPECIES_s(M,N))
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
                  3001 FORMAT(/2x,'>',I3,': Species: ',A)

                  CALL READ_DATABASE('TFM', M, N, SPECIES_s(M,N),      &
                     MW_S(M,N))
! Flag variable to stating that the database was read.
                  rDatabase(M,N) = .TRUE.
               ENDIF
            ENDDO ! Loop over species
            CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
         ENDIF

! Verify that no additional species information was given.
         DO N = NMAX(M) + 1, DIM_N_S 
            IF(MW_S(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1002) trim(iVar('MW_s',M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDDO 
      ENDDO ! Loop over solids phases

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.'&
         ,/'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_SPECIES


!----------------------------------------------------------------------!
! Subroutine: CHECK_SOLIDS_SPECIES_LEGACY                              !
! Purpose: These are legacy checks for using rrates.f to specify       !
! chemcial reactions.                                                  !
!                                                                      !
! Author: J. Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! Solids phase molecular weights.
      use physprop, only: MW_s
! Number of solids phase species.
      use physprop, only: NMAX, NMAX_s
! Flag: Database was read. (legacy)
      use physprop, only: DATABASE_READ

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of gas phase species.
      USE param, only: DIM_N_s
! Constants.
      USE param1, only: UNDEFINED_I, UNDEFINED, ZERO

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: M, N


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_SPECIES_LEGACY")

! Reconcile the new species input method with the legacy input method.
      DO M=1, MMAX_LL
         IF(SPECIES_EQ(M)) THEN

! Legacy checks for species equations.
            IF(NMAX_s(M) /= UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) trim(iVar('NMAX_s',M)), 'undefined'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,2000) trim(iVar('NMAX',M)), 'specified'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(NMAX(M) > DIM_N_S) THEN
               WRITE(ERR_MSG,1002) trim(iVar('NMAX',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
         ELSE
            IF(NMAX(M) == UNDEFINED) NMAX(M) = 1
         ENDIF
      ENDDO

! Check MW_s if solids species are present    
      DO M = 1, MMAX_LL
! Initialize flag indicating the database was read for a species.
         DO N = 1, NMAX(M)
            IF(MW_S(M,N) == UNDEFINED) THEN
               WRITE(ERR_MSG, 2000) trim(iVar('MW_s',M,N)), 'specified'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO ! Loop over species
         DO N = NMAX(M) + 1, DIM_N_S 
            IF(MW_S(M,N) /= UNDEFINED) THEN 
               WRITE(ERR_MSG,1002) trim(iVar('MW_s',M,N))
            ENDIF 
         ENDDO 
      ENDDO ! Loop over solids phases

! Set the legacy database flag. (Also in check_gas_phase.f)
      DATABASE_READ = .FALSE.

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.'&
         ,/'Please correct the mfix.dat file.')

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be 'A,/'when ',    &
         'USE_RRATES is .TRUE.'/,'Please correct the mfix.dat file')

      END SUBROUTINE CHECK_SOLIDS_SPECIES_LEGACY


 
!----------------------------------------------------------------------!
!  Subroutine: CHECK_SOLIDS_DENSITY                                    !
!  Purpose: check the solid phase density input                        !
!                                                                      !
!  Author: J.Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_SOLIDS_DENSITY(MMAX_LL)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve variable solids density.
      use run, only: SOLVE_ROs
! User specified: constant solids density
      use physprop, only: RO_s0
! Calculated baseline variable solids density.
      use physprop, only: BASE_ROs
! Baseline species densities
      use physprop, only: RO_Xs0
! Baseline species mass fractions.
      use physprop, only: X_s0
! Index of inert solids species
      use physprop, only: INERT_SPECIES
! Number of solids phase species.
      use physprop, only: NMAX

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE
      use param1, only: UNDEFINED, UNDEFINED_I
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of solids species.
      use param, only: DIM_N_s

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Arguments:
!---------------------------------------------------------------------//
! Total number of solids phases
      INTEGER, intent(in) :: MMAX_LL

! External functions:
!---------------------------------------------------------------------//
! Solids phase density caMulation
      DOUBLE PRECISION, EXTERNAL :: EOSS0

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: M, N


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DENSITY")

! Initialize the flag for variable solids density.
      SOLVE_ROs = .FALSE.

! Check each solids phase.
      DO M = 1, MMAX_LL

! Set the flag for variable solids density if any of the input parameters
! are specified.
         IF(INERT_SPECIES(M) /= UNDEFINED_I) SOLVE_ROs(M) = .TRUE.
         DO N=1, DIM_N_s
            IF(RO_Xs0(M,N) /= UNDEFINED) SOLVE_ROs(M) = .TRUE.
            IF(X_S0(M,N) /= UNDEFINED) SOLVE_ROs(M) = .TRUE.
         ENDDO


! Verify that one -and only one- solids density model is in use.
         IF(RO_S0(M) == UNDEFINED .AND. .NOT.SOLVE_ROs(M)) THEN
            WRITE(ERR_MSG, 1101) M
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1101: No solids density information for phase ',  &
         I2,'.',/'Please correct the mfix.dat file.')

! Check if the constant solids phase density is physical.
         ELSEIF(RO_S0(M) /= UNDEFINED .AND. SOLVE_ROs(M)) THEN
            WRITE(ERR_MSG, 1100) M
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1101 FORMAT('Error 1101: Conflicting solids density input specified ',&
         'for solids',/'phase ',I2,'. Constant solids density ',       &
         'specified (RO_s0) along with one',/'or more of the variable',&
         ' solids density parameters:',/'RO_Xs0, X_s0, INERT_SPECIES.',&
         /'Please correct the mfix.dat file.')

         ENDIF

! Check physical restrictions on variable solids density input.
         IF(SOLVE_ROs(M)) THEN
! Check INERT_SPECIES
            IF(INERT_SPECIES(M) == UNDEFINED_I) THEN
               WRITE(ERR_MSG,1000) trim(iVar('INERT_SPECIES',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(INERT_SPECIES(M) < 1 .OR.                          &
               INERT_SPECIES(M) > NMAX(M)) THEN
               WRITE(ERR_MSG, 1001) trim(iVar('INERT_SPECIES',M)),    &
                  trim(iVal(INERT_SPECIES(M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            DO N=1, NMAX(M)
! Check RO_Xs0
               IF(RO_Xs0(M,N) == UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('RO_Xs0',M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ELSEIF(RO_Xs0(M,N) < ZERO) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('RO_Xs0',M,N)),        &
                     trim(iVal(RO_xs0(M,N)))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
! Check X_s0
               IF(X_s0(M,N) == UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('X_s0',M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ELSEIF(X_s0(M,N) < ZERO .OR. X_s0(M,N) >= ONE) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('X_s0',M,N)),        &
                     trim(iVal(X_s0(M,N))) 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO

! Check for input overflow.
            DO N=NMAX(M)+1, DIM_N_s
               IF(RO_Xs0(M,N) /= UNDEFINED) THEN
                  WRITE(ERR_MSG,1002) trim(iVar('RO_Xs0',M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               IF(X_s0(M,N) /= UNDEFINED) THEN
                  WRITE(ERR_MSG,1002) trim(iVar('X_s0',M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO

! Check X_s0(Inert)
            IF(X_s0(M,INERT_SPECIES(M)) == ZERO) THEN
               WRITE(ERR_MSG,1102) M, INERT_SPECIES(M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1102 FORMAT('Error 1102: Invalid baseline inert speices mass',/       &
         'fraction. The inert spcies mass fraction must be greater ',  &
         'than zero.',/' Phase ',I2,' Inert Species: ',I3,' X_s0 = 0.0')

! All of the information for variable solids density has been verified
! as of this point. Calculate and store the baseline density.
            BASE_ROs(M) = EOSS0(M)

! Check physical restrictions on constant solids density input.
         ELSEIF(RO_S0(M) <= ZERO) THEN
            WRITE(ERR_MSG,1101) iVar('RO_s0',M), iVal(RO_s0(M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

      ENDDO

! Check for input overflow.
      DO M = MMAX_LL+1, DIM_M
         IF(RO_S0(M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1002) trim(iVar('RO_s0',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(INERT_SPECIES(M) /= UNDEFINED_I) THEN
            WRITE(ERR_MSG,1002) trim(iVar('INERT_SPECIES',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         DO N=1, DIM_N_s
            IF(RO_Xs0(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1002) trim(iVar('RO_Xs0',M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(X_s0(M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1002) trim(iVar('X_s0',M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDDO

! Finalize the error messges
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal input: ',A,' specified out of range.'&
         ,/'Please correct the mfix.dat file.')

 1300 FORMAT(//1X,70('*')/' From: CHECK_SOLIDS_COMMON',/,' Error 1300:'      &
         ' No solids density information for phase ',I2,'.')

      END SUBROUTINE CHECK_SOLIDS_DENSITY
