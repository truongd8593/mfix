!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_THERMO                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_THERMO

      Use compar
      Use des_thermo
      Use discretelement
      Use funits  
      Use interpolation
      Use physprop
      Use run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER I, M ! dummy loop index

! Number of processors used. (DES heat transfer is currently limited
! to serial runs!)
      INTEGER CHECK_MPI


! The size of the thermodynamic neighborhoods needed for particle-fluid-
! particle conduction and radiation.
      DOUBLE PRECISION SZ_PFP, SZ_RAD
!-----------------------------------------------      
! Functions
!-----------------------------------------------      



!------------------------------------------------------------------------->>>>
      OUTPUT_DATA_TIME = 0.0d0
!-------------------------------------------------------------------------<<<<


      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_THERMO ---------->'

      IF(DES_ENERGY_EQ)THEN

! Loop over all solids phases
         DO M=1,MMAX
! Verify that the thermal conductivity values are physical.
            IF(DES_C_ps0(M) .LT. ZERO)THEN
               IF(DMP_LOG) THEN
                  WRITE(UNIT_LOG,1001)'DES_C_ps0','unphysical',M
                  WRITE(*,1001)'DES_C_ps0','unphysical',M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF
! Verify that the thermal conductivity values are physical and defined.
            IF(DES_K_s0(M) .LT. ZERO)THEN
               IF(DMP_LOG) THEN
                  WRITE(UNIT_LOG,1001)'DES_K_s0','unphysical',M
                  WRITE(*,1001)'DES_K_s0','unphysical',M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ELSEIF(DES_K_s0(M) == UNDEFINED)THEN
               IF(DMP_LOG) THEN
                  WRITE(UNIT_LOG,1001)'DES_K_s0','undefined',M
                  WRITE(*,1001)'DES_K_s0','undefined',M
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO

! Convection Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         IF(DES_CONV_EQ) THEN
! Ensure that the continuum energy equations will be solved.
            IF(.NOT.ENERGY_EQ)THEN
               IF(DMP_LOG) WRITE(*,1000)
               CALL MFIX_EXIT(myPE)
            ENDIF

            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES convection model will be solved.'

            SELECT CASE(TRIM(DES_CONV_CORR))
! Verify the selected convective heat transfer coefficient model
               CASE ('RANZ_1952')
                  ! Ranz, W.E. and Marshall, W.R., "Frication and transfer 
                  ! coefficients for single particles and packed beds," 
                  ! Chemical Engineering Science, Vol. 48, No. 5, pp 247-253,
                  ! 1925
! The Ranz and Marshall correlation needs no additional input from the
! user.
                  IF(DMP_LOG) WRITE(*,'(9X,A,A)') &
                     'Heat transfer coefficient',&
                     ' model: Ranz and Marshall (1952)'

! If the heat transfer coefficient correlation provided by the user does
! not match one of the models outlined above, flag the error and exit.
               CASE DEFAULT
                  IF(DMP_LOG) THEN
                     WRITE(*,'(6X,A)') &
                        'INVALID DES CONVECTION MODEL (DES_CONV_CORR):'
                     WRITE(*,'(6X,A)')'The available models include:'
                     WRITE(*,'(9X,A)')'RANZ_1952 - Ranz and Marshall (1952)'
                  ENDIF
                  CALL MFIX_EXIT
            END SELECT


         ELSE
            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES convection model will NOT be solved.'
         ENDIF ! END convection model checks


! Conduction Equations:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         IF(DES_COND_EQ) THEN
            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES conduction model will be solved.'

! Particle-particle Conduction:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            IF(DES_COND_EQ_PP) THEN
               DO M=1,MMAX
! Verify that a thermal conductivity value is specified for each solids
! phase and that the value is physical.
                  IF(DES_K_S0(M) == UNDEFINED)THEN
                     IF(DMP_LOG) THEN
                        WRITE(UNIT_LOG,1001)'DES_K_s0','undefined',M
                        WRITE(*,1001)'DES_K_s0','undefined',M
                     ENDIF
                     CALL MFIX_EXIT(myPE)
                  ELSEIF(DES_K_S0(M) .LT. ZERO) THEN
                     IF(DMP_LOG) THEN
                        WRITE(UNIT_LOG,1001)'DES_K_s0','unphysical',M
                        WRITE(*,1001)'DES_K_s0','unphysical',M
                     ENDIF
                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDDO
               IF(DMP_LOG) WRITE(*,'(6X,A)')&
                  'Solving the particle-particle conduction model.'
! A secondary neighborhood is established surrounding a particle to 
! determine the neighboring particles for possible particle-particle
! heat transfer.
               FIND_THERMO_NBRHD = .TRUE.
            ELSE
               IF(DMP_LOG) WRITE(*,'(5X,A)')&
                  'NOT Solving the particle-particle conduction model.'
            ENDIF

! Particle-fluid-particle Conduction:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            IF(DES_COND_EQ_PFP) THEN
! Verify that either the simulation is coupled with a gas-phase
! simulation or that K_g0 has been provided in the dat file.
               IF(K_G0 == UNDEFINED .AND. &
                  .NOT.DES_CONTINUUM_COUPLED)THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1002)
                     WRITE(UNIT_LOG,1002)
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF

! Set the default value for the minimum distance separating particles'
! surfaces.
               IF(DES_MIN_COND_DIST == UNDEFINED)THEN
                  DES_MIN_COND_DIST = 1.0D-04 ! cm
                  IF (UNITS == 'SI') DES_MIN_COND_DIST = &
                     DES_MIN_COND_DIST/100.0  ! m
               ENDIF

               IF(DMP_LOG) WRITE(*,'(6X,A)')&
                  'Solving the particle-fluid-particle conduction model.'
! A secondary neighborhood is established surrounding a particle to 
! determine the neighboring particles for possible particle-fluid-
! particle heat transfer.
               FIND_THERMO_NBRHD = .TRUE.
            ELSE
               IF(DMP_LOG) WRITE(*,'(6X,A)')&
                  'NOT Solving the particle-fluid-particle conduction model.'
            ENDIF

         ELSE 
            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES conduction model will NOT be solved.'
! Set the logicals for the conduction sub-models to false.
! This should not be needed, but is done for precaution.
            DES_COND_EQ_PFP = .FALSE.
            DES_COND_EQ_PP  = .FALSE.

         ENDIF ! END conduction model checks

! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         IF(DES_RADI_EQ) THEN
            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES radiation model will be solved.'
! A secondary neighborhood is established surrounding a particle to 
! determine the neighboring particles for possible radiative 
! heat transfer.
            FIND_THERMO_NBRHD = .TRUE.
! Set the value of the Stefan-Boltzman Constant based on the untis
            IF(UNITS == 'SI')THEN
!              W/((m^2).K^4)
               SB_CONST = 5.6704d0*(10.0d0**(-8))
            ELSE
!              cal/((cm^2).sec.K^4)
               SB_CONST = 1.355282d0*(10.0d0**(-12))
            ENDIF

            DO M=1,MMAX
! Verify that a emmisivity value is specified for each solids pase
               IF(DES_Em(M) == UNDEFINED)THEN
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1001)'DES_Em','undefined',M
                     WRITE(*,1001)'DES_Em','undefined',M
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ELSEIF(DES_Em(M) .LT. ZERO .OR. &
                  DES_Em(M) .GT. ONE) THEN
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1001)'DES_Em','unphysical',M
                     WRITE(*,1001)'DES_Em','unphysical',M
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDDO
         ELSE
            IF(DMP_LOG) WRITE(*,'(6X,A)') &
               'DES radiation model will NOT be solved.'
         ENDIF ! END radiation model checks

! Calculate the size of the thermodynamic neighbor (if needed)
         IF( FIND_THERMO_NBRHD) THEN
            SZ_PFP = ZERO
! If solving conductive  heat transfer, the minimum neighborhood size
! is [1.3]*DES_RADIUS by default. A buffer of 0.1 is added to ensure
! that the neighbor is seen as soon as possible. The size of the fluid
! lens for particle-fluid-particle heat transfer can be adjusted by 
! specifying a different value in the mfix.dat file.
            IF(DES_COND_EQ_PP .OR. DES_COND_EQ_PFP) &
               SZ_PFP = MAX_RADIUS * ( 2.05d0 + FLPC )
            SZ_RAD = ZERO
! If radiation heat transfer is solved, the radius of the radiative
! heat transfer domain is [3.0]*DES_RADIUS by default. The size of the
! radiative heat transfer domain on the particle may be adjusted by
! specifying a new value in the mfix.dat file.
            IF(DES_RADI_EQ) &
               SZ_RAD = MAX_RADIUS * RDPC
            NBRHD_SZ = MAX( SZ_PFP, SZ_RAD )
            IF(DMP_LOG) WRITE(*,'(6X,A,F10.6)') &
               'Radius of thermodynamic neighborhood is ',NBRHD_SZ
         ENDIF ! FIND_THERMO_NBRHD

! Check the number of processors. DES reactive chemistry is currently 
! limited to serial runs.
         CHECK_MPI = NODESI * NODESJ * NODESK
         IF(CHECK_MPI.NE.1) THEN
            WRITE (UNIT_LOG, 1003)
            CALL MFIX_EXIT(myPE)
         ENDIF

      ELSE ! DES_ENERGY_EQ

         IF(DMP_LOG) WRITE(*,'(6X,A)') &
            'DES energy equations are not being solved.'

! Reinitialize the heat transfer logicals to false.
         DES_CONV_EQ = .FALSE.
         DES_COND_EQ = .FALSE.
         DES_RADI_EQ = .FALSE.

         DES_COND_EQ_PFP = .FALSE.
         DES_COND_EQ_PP  = .FALSE.

      ENDIF ! DES_ENERGY_EQ

      IF(DMP_LOG) WRITE(*,'(1X,A)')&
         '<---------- END CHECK_DES_THERMO ----------'

      RETURN

 1000 FORMAT(/1X,70('*')/, ' From: CHECK_DES_THERMO',/, ' Message: ',&
         'The DES convection model requires that the energy ', &
         'equations',/1X,'(ENERGY_EQ) for the continuum phase must ',&
         'also be solved.',/1X,70('*')/)

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_THERMO',/, ' Message: ',&
         A,' is ',A,' for solids phase ',I2,'. Check mfix.dat.', &
         /1X,70('*')/)

 1002 FORMAT(/1X,70('*')/,' From: CHECK_DES_THERMO',/,' Message:',     &
      ' The particle-fluid-particle model requires that either a'/,    &
      ' simulation be coupled with a gas-phase simultion or that',     &
      ' K_g0 be',/,' provided in the mfix.dat file. Check mfix.dat.',  &
      /1X,70('*')/)

 1003 FORMAT(/1X,70('*')/ ' From: CHECK_DES_THERMO',/' Message: ',&
         'DES heat transfer modules are currently limited to serial',/&
         ' runs. Set nodesi, modesj, and nodesk to 1 in the mfix.dat',&
         ' file.',/1X,70('*')/)

      END SUBROUTINE CHECK_DES_THERMO
