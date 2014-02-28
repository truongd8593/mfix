!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_DES_SOLIDS                                     !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES


! Global Variables:
!---------------------------------------------------------------------//
! NONE

! Global Parameters:
!---------------------------------------------------------------------//
! NONE

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! NONE

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES")


!      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
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
!         CLOSE_PACKED(:) = .TRUE. 
!         MOMENTUM_X_EQ(1:DIM_M) = .FALSE.
!         MOMENTUM_Y_EQ(1:DIM_M) = .FALSE.
!         MOMENTUM_Z_EQ(1:DIM_M) = .FALSE. 

! Check continuum solids phase species data. Clear anything that was
! specified and warn the user.
!         CALL CHECK_04_TFM_DEM()
!         RETURN
!      ENDIF




! Check settings on cluster identification
!      IF(DES_CALC_CLUSTER) THEN
!         IF(CLUSTER_LENGTH_CUTOFF .EQ. UNDEFINED) THEN
!            IF(DMP_LOG) WRITE(UNIT_LOG,1101)
!            IF(DMP_LOG) WRITE(*,1101)
!            CALL MFIX_EXIT(myPE)
!         ENDIF
!         IF(FACTOR_RLM < &
!            1.d0+CLUSTER_LENGTH_CUTOFF/(2.d0*MAX_RADIUS)) THEN
!            IF(DMP_LOG) WRITE(UNIT_LOG,1102)
!            IF(DMP_LOG) WRITE(*,1102)
!            CALL MFIX_EXIT(myPE)
!         ENDIF
!      ENDIF


      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_SOLIDS_DES




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_DES_ENERGY                                  !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES_ENERGY

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use compar, only: NODESI, NODESJ, NODESK


      use run, only: ENERGY_EQ
      use run, only: UNITS
      use des_thermo, only: DES_CONV_CORR_ENUM
      use des_thermo, only: RANZ_1952


      Use discretelement
      Use funits  
      Use interpolation
      Use physprop
      Use des_thermo


      use error_manager

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


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES_ENERGY")


      IF(.NOT.ENERGY_EQ)THEN
! Reinitialize the heat transfer logicals to false.
         DES_CONV_EQ     = .FALSE.
         DES_COND_EQ     = .FALSE.
         DES_COND_EQ_PFP = .FALSE.
         DES_COND_EQ_PP  = .FALSE.
         DES_RADI_EQ     = .FALSE.

         CALL FINL_ERR_MSG

         RETURN
      ENDIF



! Check the number of processors. DES reactive chemistry is currently 
! limited to serial runs.
      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE(ERR_MSG, 2000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2000 FORMAT('Error 2000: Currently, simulations with discrete solids',&
        ' heat transfer',/'modules are limited to serial runs. Please',&
        ' correct the mfix.dat file.')



! Gas/Solids convection:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify the selected convective heat transfer coefficient model
      SELECT CASE(TRIM(DES_CONV_CORR))
! Ranz, W.E. and Marshall, W.R., "Frication and transfer coefficients
! for single particles and packed beds,"  Chemical Engineering Science,
! Vol. 48, No. 5, pp 247-253, 1952.
      CASE ('RANZ_1952')
         DES_CONV_CORR_ENUM = RANZ_1952
! If the heat transfer coefficient correlation provided by the user does
! not match one of the models outlined above, flag the error and exit.
      CASE DEFAULT
         WRITE(ERR_MSG,1001)'DES_CONV_CORR', trim(DES_CONV_CORR)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')




! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      IF(DES_RADI_EQ) THEN

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






! Conduction Equations:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      IF(DES_COND_EQ) THEN

! Set the default value for the minimum distance separating particles'
! surfaces.
         IF(DES_MIN_COND_DIST == UNDEFINED)THEN
            DES_MIN_COND_DIST = 1.0D-04 ! cm
            IF (UNITS == 'SI') DES_MIN_COND_DIST = &
               DES_MIN_COND_DIST/100.0  ! m
         ENDIF
      ENDIF





      CALL FINL_ERR_MSG


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


 1004 FORMAT(/1X,70('*')/ ' From: CHECK_DES_THERMO',/' Message: ',&
         'urrently, the DEM heat transfer models are not capable of',/&
         ' being restarted.',/1X,70('*')/)

      END SUBROUTINE CHECK_SOLIDS_DES_ENERGY

