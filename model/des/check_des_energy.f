!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_THERMO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_ENERGY

!-----------------------------------------------
! Modules
!-----------------------------------------------
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

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_ENERGY ---------->'

      IF(.NOT.DES_ENERGY_EQ)THEN
! Reinitialize the heat transfer logicals to false.
         DES_CONV_EQ     = .FALSE.
         DES_COND_EQ     = .FALSE.
         DES_COND_EQ_PFP = .FALSE.
         DES_COND_EQ_PP  = .FALSE.
         DES_RADI_EQ     = .FALSE.

         IF(DMP_LOG) WRITE(*,'(6X,A)') &
            'DES energy equations are not being solved.'

         IF(DMP_LOG) WRITE(*,'(1X,A)')&
            '<---------- END CHECK_DES_THERMO ----------'
         RETURN
      ENDIF

! Check the number of processors. DES reactive chemistry is currently 
! limited to serial runs.
      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE (UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)
      ENDIF

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

! Verify the selected convective heat transfer coefficient model
         SELECT CASE(TRIM(DES_CONV_CORR))
! Ranz, W.E. and Marshall, W.R., "Frication and transfer coefficients
! for single particles and packed beds,"  Chemical Engineering Science,
! Vol. 48, No. 5, pp 247-253, 1952.
            CASE ('RANZ_1952')
! The Ranz and Marshall correlation needs no additional input.
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
                  WRITE(*,'(9X,A)') &
                     'RANZ_1952 - Ranz and Marshall (1952)'
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
            IF(DMP_LOG) WRITE(*,'(6X,A)')&
               'Solving the particle-particle conduction model.'
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
         ELSE
            IF(DMP_LOG) WRITE(*,'(6X,A,A)')&
               'NOT Solving the particle-fluid-particle ',&
               'conduction model.'
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
! Set the value of the Stefan-Boltzman Constant based on the untis
         IF(UNITS == 'SI')THEN
            SB_CONST = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
         ELSE
            SB_CONST = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
         ENDIF
      ELSE
         IF(DMP_LOG) WRITE(*,'(6X,A)') &
            'DES radiation model will NOT be solved.'
      ENDIF ! END radiation model checks


      IF(DMP_LOG) WRITE(*,'(1X,A)')&
         '<---------- END CHECK_DES_THERMO ----------'

      RETURN

 1000 FORMAT(/1X,70('*')/, ' From: CHECK_DES_THERMO',/, ' Message: ',&
         'The DES convection model requires that the energy ', &
         'equations',/1X,'(ENERGY_EQ) for the continuum phase must ',&
         'also be solved.',/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/,' From: CHECK_DES_THERMO',/,' Message:',     &
      ' The particle-fluid-particle model requires that either a'/,    &
      ' simulation be coupled with a gas-phase simultion or that',     &
      ' K_g0 be',/,' provided in the mfix.dat file. Check mfix.dat.',  &
      /1X,70('*')/)

 1003 FORMAT(/1X,70('*')/ ' From: CHECK_DES_THERMO',/' Message: ',&
         'DES heat transfer modules are currently limited to serial',/&
         ' runs. Set nodesi, modesj, and nodesk to 1 in the mfix.dat',&
         ' file.',/1X,70('*')/)

 1004 FORMAT(/1X,70('*')/ ' From: CHECK_DES_THERMO',/' Message: ',&
         'urrently, the DEM heat transfer models are not capable of',/&
         ' being restarted.',/1X,70('*')/)

      END SUBROUTINE CHECK_DES_ENERGY
