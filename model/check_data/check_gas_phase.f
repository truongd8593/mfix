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



! Check the input specifications for gas species.
      IF(USE_RRATES)THEN
         CALL CHECK_GAS_SPECIES_LEGACY
      ELSE
         CALL CHECK_GAS_SPECIES
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



 1057 FORMAT('Message: Gas phase species ',I2,' name (SPECIES_g) is',  &
         ' undefined.',/'Please correct the data file.')



 1500 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' not specified') 
 1501 FORMAT(1X,/,1X,'MW_g for gas species ',I3,' specified') 

      END SUBROUTINE CHECK_GAS_PHASE



!----------------------------------------------------------------------!
! Subroutine: CHECK_GAS_SPECIES                                     !
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
      SUBROUTINE CHECK_GAS_SPECIES

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
      INTEGER  :: MMAX_LL


! Flag that the user was already warned why a call to the thermo-
! chemical database is being made.
      LOGICAL WARNED_USR

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
! If true, a call to the thermochemical database is made.
      LOGICAL EEQ_CPG

! Flag that the solids phase species equations are solved and the 
! molecular weight for a species are not given in the data file.
! If true, a call to the thermochemical database is made.
      LOGICAL SEQ_MWg

      LOGICAL MWg_ROg

!      LOGICAL thermoHeader

      INTEGER :: N


      CALL INIT_ERR_MSG("CHECK_GAS_SPECIES")


! Reconcile the new species input method with the legacy input method.
      IF(SPECIES_EQ(0)) THEN

         IF(NMAX_g == UNDEFINED_I) THEN
            WRITE(ERR_MSG,1000) 'NMAX_g'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(NMAX_g > DIM_N_G) THEN
            WRITE(ERR_MSG,1001) 'NMAX_g', iVal(NMAX_g)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE
            NMAX(0) = NMAX_g
         ENDIF
! Set the number of species to one if the species equations are not solved and
! the number of species is not specified.
      ELSE
         NMAX(0) = merge(1, NMAX_g, NMAX_g == UNDEFINED_I)
      ENDIF

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

! Flag that the energy equations are solved and specified solids phase
! specific heat is undefined.
      EEQ_CPG = (ENERGY_EQ .AND. C_PG0 == UNDEFINED)
      IF(EEQ_CPG) THEN
         WRITE(ERR_MSG,2000)
         CALL FLUSH_ERR_MSG
      ENDIF

 2000 FORMAT('Message: 2000 The energy equations are being solved ',   &
         '(ENERGY_EQ) and',/'the constant gas specific heat is ',      &
         'undefined (C_PG0). Thus, the thermo-',/'chemical database ', &
         'will be used to gather specific heat data on the',/          &
         'individual gas phase species.')

      MWg_ROg = .FALSE.
      SEQ_MWg = .FALSE.
      IF(MW_AVG == UNDEFINED) THEN
         DO N=1,NMAX(0)
            IF(MW_g(N) == UNDEFINED) THEN
               IF(RO_G0 == UNDEFINED) MWg_ROg = .TRUE.
               IF(SPECIES_EQ(0)) SEQ_MWg = .TRUE.
            ENDIF
         ENDDO
      ENDIF

      IF(MWg_ROg) THEN
         WRITE(ERR_MSG, 2001)
         CALL FLUSH_ERR_MSG
      ENDIF

 2001 FORMAT('Message 2001: MW_AVG and RO_G0 are undefined and one or',&
         ' more species',/'molecular weights are undefined. The therm',&
         'ochemical database will be',/'used in an attempt to gather ',&
         'missing molecular weight data.')

      IF(MWg_ROg) THEN
         WRITE(ERR_MSG, 2002)
         CALL FLUSH_ERR_MSG
      ENDIF

 2002 FORMAT('Message 2002: One or more species molecular weights are',&
         ' undefined and',/'the gas phase species equations are being',&
         ' solved (SOLVE_EQ(0)). The',/'thermochemical database will ',&
         'be used in an attempt to gather missing',/'molecular weight',&
         ' data.')

      IF(EEQ_CPG  .OR. SEQ_MWg .OR. MWg_ROg) THEN
! Initialize flag indicating the database was read for a species.
         rDatabase(0,:) = .FALSE.

         WRITE(ERR_MSG, 3000)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 3000 FORMAT('Message 3000: Searching thermochemical databases for ',  &
         'gas phase',/'species data.',/'  ')


         DO N = 1, NMAX(0)
! Notify the user of the reason the thermochemical database is used.
! Flag that the species name is not provided.
            IF(SPECIES_g(N) == UNDEFINED_C) THEN
               WRITE(ERR_MSG,1000) trim(iVar('SPECIES_g',N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
! Update the log files.
            WRITE(ERR_MSG, 3001) N, trim(SPECIES_g(N))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            3001 FORMAT(2x,'>',I3,': Species: ',A)
! Read the database.
            CALL READ_DATABASE('TFM', 0, N, SPECIES_g(N), MW_g(N))
! Flag variable to stating that the database was read.
            rDatabase(0,N) = .TRUE.
! Flag the legacy variable to prevent re-reading the database.
            DATABASE_READ = .TRUE.
         ENDDO ! Loop over species

      ENDIF

! Verify that no additional species information was given.
      DO N = NMAX(0) + 1, DIM_N_G
         IF(MW_G(N) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('MW_g',N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_GAS_SPECIES


!----------------------------------------------------------------------!
! Subroutine: CHECK_GAS_SPECIES_LEGACY                                 !
! Purpose: These are legacy checks for using rrates.f to specify       !
! chemcial reactions.                                                  !
!                                                                      !
! Author: J. Musser                                  Date: 03-FEB-14   !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_GAS_SPECIES_LEGACY


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Solve species equations
      use run, only: SPECIES_EQ
! Gas phase molecular weights.
      use physprop, only: MW_g
! Number of gas phase species.
      use physprop, only: NMAX, NMAX_g


! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of gas phase species.
      USE param, only: DIM_N_g
! Constants.
      USE param1, only: UNDEFINED_I, UNDEFINED, ZERO


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter.
      INTEGER :: N


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GAS_SPECIES_LEGACY")


! Reconcile the new species input method with the legacy input method.
      IF(SPECIES_EQ(0)) THEN
! Legacy checks for species equations.
         IF(NMAX_g /= UNDEFINED_I) THEN
            WRITE(ERR_MSG,2000) 'NMAX_g', 'undefined'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(NMAX(0) == UNDEFINED_I) THEN
            WRITE(ERR_MSG,2000) trim(iVar('NMAX',0)), 'specified'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(NMAX(0) > DIM_N_G) THEN
            WRITE(ERR_MSG,1001) trim(iVar('NMAX',0)), iVal(NMAX(0))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! Set the number of species to one if the species equations are not
! solved and the number of species is not specified.
      ELSE
         IF(NMAX(0) == UNDEFINED_I) NMAX(0) = 1
      ENDIF

! Check MW_g if solids species are present    
      DO N = 1, NMAX(0)
         IF(MW_G(N) == UNDEFINED) THEN 
            WRITE(ERR_MSG,2000)trim(iVar('MW_g',N)), 'specified'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(MW_G(N) <= ZERO) THEN 
            WRITE(ERR_MSG,1001)trim(iVar('MW_g',N)), iVal(MW_G(N)) 
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO ! Loop over species
      DO N = NMAX(0) + 1, DIM_N_G
         IF(MW_G(N) /= UNDEFINED) THEN 
            WRITE(ERR_MSG,1001)trim(iVar('MW_g',N)), iVal(MW_G(N)) 
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO 

      CALL FINL_ERR_MSG

      RETURN

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 2000 FORMAT('Error 2000: Invalid input. ',A,' must be 'A,/'when ',    &
         'USE_RRATES is .TRUE.'/,'Please correct the mfix.dat file')

      END SUBROUTINE CHECK_GAS_SPECIES_LEGACY

