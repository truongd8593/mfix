!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke gas/solids coupled simulation.
      USE discretelement, only: DES_CONTINUUM_COUPLED
! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only: GENER_PART_CONFIG
! Runtime Flag: Interpolate mean field quantities.
      USE discretelement, only: DES_INTERP_ON
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Invoke TFM/DEM hybrid model.
      USE discretelement, only: DES_CONTINUUM_HYBRID
! Runtime Flag: Invoke DEM cluster detection.
      USE discretelement, only: DES_CALC_CLUSTER
! Runtime Flag: Utilize cutcell geometry.
      USE cutcell, only: CARTESIAN_GRID
! Runtime Flag: Interpolate DEM field quanties.
      USE discretelement, only: DES_INTERP_MEAN_FIELDS
! Runtime Flag: Invoke gas/solids coupled simulation.
      USE discretelement, only: DES_CONTINUUM_COUPLED

! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! TFM solids phase diameters and densities. (DEM default)
      USE physprop, only: D_p0, RO_s0

! User specified integration method.
      USE discretelement, only: DES_INTG_METHOD
      USE discretelement, only: INTG_ADAMS_BASHFORTH
      USE discretelement, only: INTG_EULER
! User specified neighbor search method.
      USE discretelement, only: DES_NEIGHBOR_SEARCH
! User specified data out format (VTP, TecPlot)
      USE discretelement, only: DES_OUTPUT_TYPE
! Max/Min particle radi
      USE discretelement, only: MAX_RADIUS, MIN_RADIUS
! Max distance between two particles in a cluster.
      USE discretelement, only: CLUSTER_LENGTH_CUTOFF
! Runtime Flag: Periodic boundaries
      USE discretelement, only: DES_PERIODIC_WALLS
      USE discretelement, only: DES_PERIODIC_WALLS_X
      USE discretelement, only: DES_PERIODIC_WALLS_Y
      USE discretelement, only: DES_PERIODIC_WALLS_Z

! Subroutine access.
      USE desgrid, only: DESGRID_CHECK

      use physprop, only: MMAX

      USE run, only: MOMENTUM_X_EQ
      USE run, only: MOMENTUM_Y_EQ
      USE run, only: MOMENTUM_Z_EQ

      USE physprop, only: CLOSE_PACKED

      USE mpi_utility


! Global Parameters:
!---------------------------------------------------------------------//
!      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, lM  ! Solids phase Index


      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE")


      DES_D_p0 = UNDEFINED
      DES_RO_s = UNDEFINED

      MAX_RADIUS = -UNDEFINED
      MIN_RADIUS =  UNDEFINED


      DO M=MMAX+1, MMAX+DES_MMAX
! Copy of the input keyword values into discrete solids arrays. We may be
! able to remove the DES_ specific variables moving forward.
         DES_D_p0(M) = D_p0(M)
         DES_RO_s(M) = RO_s0(M)
! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks        
         MAX_RADIUS = MAX(MAX_RADIUS, 0.5d0*DES_D_P0(M))
         MIN_RADIUS = MIN(MIN_RADIUS, 0.5d0*DES_D_P0(M))
      ENDDO


! Set close_packed to true to prevent possible issues stemming from the
! pressure correction equation.  Specifically, if closed_packed is false
! then a mixture pressure correction equation is invoked and this is not
! correctly setup for DEM.  To do so would require ensuring that
! 1) the solids phase continuum quantities used in these equations are
!    correctly set based on their DEM counterparts and 
! 2) the pressure correction coefficients for such solids phases are 
!    also calculated (currently these calculations are turned off 
!    when using DEM)
      CLOSE_PACKED((MMAX+1):DIM_M) = .TRUE.


! Turn off the 'continuum' equations for discrete solids if the user
! specified them.  We could make use of these flags.
      MOMENTUM_X_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Y_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Z_EQ((MMAX+1):DIM_M) = .FALSE. 

! Derive periodicity from cyclic boundary flags.
      DES_PERIODIC_WALLS_X = CYCLIC_X .OR. CYCLIC_X_PD
      DES_PERIODIC_WALLS_Y = CYCLIC_Y .OR. CYCLIC_Y_PD
      DES_PERIODIC_WALLS_Z = CYCLIC_Z .OR. CYCLIC_Z_PD

      DES_PERIODIC_WALLS = (DES_PERIODIC_WALLS_X .OR.                  &
        DES_PERIODIC_WALLS_Y .OR. DES_PERIODIC_WALLS_Z)

! Overwrite user's input in case of DEM (no fluid)
      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
         DES_INTERP_ON = .FALSE.

      ELSE
! MPPIC and DES/CUTCELL simulations require that the mean feilds be
! interpolated without regard to user specifications.
         IF(MPPIC.OR.CARTESIAN_GRID) DES_INTERP_MEAN_FIELDS = .TRUE.

! DES_INTERP_MEAN_FIELDS is invoked if des_interp_on is true to remain
! consistent with previous implementations.
         IF(DES_INTERP_ON)  DES_INTERP_MEAN_FIELDS= .TRUE.
      ENDIF
            



! Check for valid neighbor search option.
      SELECT CASE(DES_NEIGHBOR_SEARCH)
      CASE (1) ! N-Square
      CASE (2)
         WRITE(ERR_MSG,2001) 2, 'QUADTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (3)
         WRITE(ERR_MSG,2001) 3, 'OCTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (4) ! Grid based
      CASE DEFAULT
         WRITE(ERR_MSG,2001) DES_NEIGHBOR_SEARCH,'UNKNOWN'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2001 FORMAT('Error 2001:Invalid DES_NEIGHBOR_SEARCH method: ',I2,1X,  &
         A,/'Please correct the mfix.dat file.')

      END SELECT


! Check the output file format 
      IF(DES_OUTPUT_TYPE == UNDEFINED_C) DES_OUTPUT_TYPE = 'PARAVIEW'
      SELECT CASE(trim(DES_OUTPUT_TYPE))
      CASE ('PARAVIEW')
      CASE ('TECPLOT')
      CASE DEFAULT
         WRITE(ERR_MSG,2010) trim(DES_OUTPUT_TYPE)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2010 FORMAT('Error 2010:Invalid DES_OUTPUT_TYPE: ',A,/'Please '       &
         'correct the mfix.dat file.')

      END SELECT


! Check for valid integration method
      SELECT CASE(trim(DES_INTG_METHOD))
      CASE ('EULER')
         INTG_EULER = .TRUE.
         INTG_ADAMS_BASHFORTH = .FALSE.
         !DES_INTG_METHOD_ENUM = 1
      CASE ('ADAMS_BASHFORTH')
         INTG_EULER = .FALSE.
         INTG_ADAMS_BASHFORTH = .TRUE.
         !DES_INTG_METHOD_ENUM = 2
      CASE DEFAULT
         WRITE(ERR_MSG,2020) trim(DES_INTG_METHOD)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2020 FORMAT('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please '      &
         'correct the mfix.dat file.')

      END SELECT



! Set flags for energy equations
      CALL CHECK_SOLIDS_COMMON_DISCRETE_ENERGY

! Check thermodynamic properties of discrete solids.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_THERMO


! Check geometry constrains.
!      CALL CHECK_DES_GEOMETRY

! Check settings on cohesion model
!      CALL CHECK_DES_COHESION

! Check settings for collision models
!      IF(.NOT.MPPIC) CALL CHECK_DES_COLLISION

! Check TFM/DEM Hybrid model settings.
!      IF (DES_CONTINUUM_HYBRID) CALL CHECK_DES_HYBRID
! Check settings for particle generation.
!      IF(GENER_PART_CONFIG) CALL CHECK_DES_PCONFIG
! Check quantities related to MPPIC      
!      IF(MPPIC) CALL CHECK_DES_MPPIC

! Check coupling settings.
!      CALL CHECK_DES_COUPLING
           
! the entire checking and setting up indices for desgridsearch
! moved to desgrid_mod to accomodate parallelization
! this is now conducted regardless of neighbor search option      
!      CALL DESGRID_CHECK


      CALL FINL_ERR_MSG


      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY                      !
!                                                                      !
!  Purpose: Check input parameters for solving discrete solids phase   !
!  energy equations.  Only DEM simulations (neither hybrid nor MPPIC)  !
!  can invoke particle-particle heat transfer. Therefore, checks for   !
!  those functions are reseved for later.                              !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY


! Global Variables:
!---------------------------------------------------------------------//
      use run, only: ENERGY_EQ
      use run, only: UNITS

      use discretelement, only: DES_MMAX

      use physprop, only: MMAX
      Use compar, only: NODESI, NODESJ, NODESK

      use des_thermo, only: DES_CONV_CORR
      use des_thermo, only: DES_CONV_CORR_ENUM
      use des_thermo, only: RANZ_1952

      use des_thermo, only: DES_COND_EQ
      use des_thermo, only: DES_COND_EQ_PFP
      use des_thermo, only: DES_COND_EQ_PP

      Use des_thermo, only: SB_CONST
      Use des_thermo, only: DES_Em

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M
! Number of processors used. (Heat transfer is limited to serial runs!)
      INTEGER :: CHECK_MPI


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_ENERGY")


      IF(.NOT.ENERGY_EQ)THEN
! Reinitialize the heat transfer logicals to false.
         DES_COND_EQ     = .FALSE.
         DES_COND_EQ_PFP = .FALSE.
         DES_COND_EQ_PP  = .FALSE.

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


! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify that a emmisivity value is specified for each solids pase
      DO M = MMAX+1, MMAX+DES_MMAX
         IF(DES_Em(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_Em',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
     ENDDO


! Set the value of the Stefan-Boltzman Constant based on the untis
      IF(UNITS == 'SI')THEN
         SB_CONST = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
      ELSE
         SB_CONST = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
      ENDIF


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_THERMO                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use compar
      Use des_thermo
      Use des_rxns
      Use discretelement
      Use funits  
      Use interpolation
      Use param1
      Use physprop
      Use run

      use physprop, only: MMAX

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      

! Loop index
      INTEGER M, N ! Phase, Species

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
!-----------------------------------------------      


! Set the flag identifying that at least one of the species equations
! are being solved.
      ANY_DES_SPECIES_EQ = &
         ANY(SPECIES_EQ((MMAX+1):(MMAX+DES_MMAX)))

! Initialize MAX_DES_NMAX. This is the maximum number of species found
! over all solids phases. This is used for allocating DES_X_s
      MAX_DES_NMAX = 1

! Loop over all solids phases
      DO M = MMAX+1, MMAX+DES_MMAX

! Copy of the input keyword values into discrete solids arrays. We may be
! able to remove the DES_ specific variables moving forward.
!         DES_C_ps0(M) = C_ps0
!         DES_K_s0(M)  = K_s0

         DES_NMAX_s(M) = NMAX(M)

         MAX_DES_NMAX = MAX(MAX_DES_NMAX, DES_NMAX_s(M))

!        DES_SPECIES_s(M,N),    
!        DES_MW_S(M,N)

      ENDDO ! DES_MMAX


      RETURN

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_THERMO',/,' Error 1001: ',&
         A,' is ',A,' for discrete solids phase ',I2,/' Please',       &
         ' correct the data file.',/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Message 1002:',  &
         ' The discrete phase energy equations are being solved',/     &
         ' (DES_ENERGY_EQ) and the specified constant solids specific',&
         ' heat is',/' undefined (DES_C_PS0). Thus, the',              &
         ' thermochemical database will be used',/' to gather',        &
         ' specific heat data on the individual soids phase species.',/&
         1X,70('*')/)

 1003 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Message 1003:',  &
         ' Discrete solids phase ',I2,' species equations are being',/ &
         ' solved, and one or more species molecular weights are',     &
         ' undefined. Thus,',/' the thermochemical database will be',  &
         ' used to gather molecular weight',/' data on the solids',    &
         ' phase speicies.',/1X,70('*')/)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DES_THERMO',/' Error 1004:',    &
         ' Discrete solids phase ',I2,' species ',I3,' name',          &
         ' (DES_SPECIES_s)',/' is undefined. Please correct the data', &
         ' file. ',/1X,70('*')/)

 1005 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Error 1005:',  &
         ' DES_NMAX_s is too large for discrete phase ',I2,'.',        &
         ' The maximum',/' number of species is ',I3,'.',/1X,70('*')/)

 1006 FORMAT(/1X,70('*')/, ' From: CHECK_DES_RXNS',/, ' Error 1006:',  &
         ' Number of species for discrete phase ',I2,' is not',        &
         ' specified.',/1X,70('*')/)


 1100 FORMAT(/'  Searching thermochemical databases for discrete',     &
         ' solids phase ',I2,', species ',I2)

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO
