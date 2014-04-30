!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
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


      use physprop, only: MMAX

      USE run, only: MOMENTUM_X_EQ
      USE run, only: MOMENTUM_Y_EQ
      USE run, only: MOMENTUM_Z_EQ

      use run, only: RUN_TYPE
      use discretelement, only: GENER_PART_CONFIG

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
            

! Overrite for restart cases.
      IF(TRIM(RUN_TYPE) .NE. 'NEW') GENER_PART_CONFIG = .FALSE.

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
      CALL CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY


! Check TFM/DEM Hybrid model settings.
!      IF (DES_CONTINUUM_HYBRID) CALL CHECK_DES_HYBRID

! Check settings for particle generation.
!      IF(GENER_PART_CONFIG) CALL CHECK_DES_PCONFIG
           


      CALL FINL_ERR_MSG


      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose: Check input parameters for solving discrete solids phase   !
!  energy equations.  Only DEM simulations (neither hybrid nor MPPIC)  !
!  can invoke particle-particle heat transfer. Therefore, checks for   !
!  those functions are reseved for later.                              !
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
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO

      use des_rxns, only: ANY_DES_SPECIES_EQ
      use des_rxns, only: MAX_DES_NMAX
      use run, only: SPECIES_EQ
      use discretelement, only: DES_MMAX
      use des_rxns, only: DES_NMAX_s
      use physprop, only: NMAX
      use physprop, only: MMAX

      use error_manager

      IMPLICIT NONE


! Loop index
      INTEGER M, N ! Phase, Species



!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_THERMO")


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


      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY                   !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check user input data                                      !
!                                                                      !
!  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

!-----------------------------------------------
! Modules 
!-----------------------------------------------      
      USE geometry, only: COORDINATES
      USE geometry, only: DO_I, DO_J, DO_K
      USE geometry, only: NO_I, NO_J, NO_K
      USE geometry, only: ZLENGTH

! Flag: Use DES E-L model
      use discretelement, only: DISCRETE_ELEMENT
      USE discretelement, only: DES_CONTINUUM_COUPLED
      USE discretelement, only: MAX_RADIUS

! flag to tell if using stl represenation in discrete models 
      USE discretelement, only: USE_STL_DES
! flag to tell if using CG 
      USE cutcell, only: cartesian_grid
! flag to tell if using stl represenation in CG 
      USE cutcell, only: use_stl 

      use param1, only: UNDEFINED_I

!flag to force conversion of regular bounding box to 
!triangular facets for particle-wall interactions in discrete models 
      use discretelement, only: des_convert_box_to_facets

! Flag: Use cohesion
      use discretelement, only: USE_COHESION

! Flag: Use MPPIC E-L model
      use mfix_pic, only: MPPIC

      use error_manager
 
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
      DOUBLE PRECISION :: MIN_DEPTH

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY")


! DEM/MPPIC is restriced to CARTESIAN coordinates.
      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error: 1100: DES and MPPIC models only support ',        &
         'CARTESIAN coordinates.')


! Check dimension. This is redundant with check_data_03.
      IF(NO_I .OR. NO_J) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1200 FORMAT('Error 1200: Illegal geometry for DEM/MPPIC. 2D ',        &
         'simulations are',/'restricted to the XY plane. Please ',     &
         'correct the mfix.dat file.')


      IF(DES_CONTINUUM_COUPLED)THEN
! Check that the depth of the simulation exceeds the largest particle
! to ensure correct calculation of volume fraction. This is important
! for coupled simulations.
         MIN_DEPTH = 2.0d0*MAX_RADIUS
         IF(ZLENGTH < MIN_DEPTH)THEN
            WRITE(ERR_MSG, 1300)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1300 FORMAT('Error 1300: The maximum particle diameter exceeds the ', &
         'simulation',/'depth (ZLENGTH). Please correct the mfix.dat ',&
         'file.')


      IF(CARTESIAN_GRID.and.discrete_element.and..not.use_stl) then 
         write(err_msg, '((A, A,/, A, /, A, /))')& 
         'cartesian grid and discrete modeds (DEM or PIC) only', & 
         'work with stl representation for walls', &
         'Quadrics and polygons are no longer supported for discrete models', &
         'Switch to STL representation for using PIC or DEM models with cartesian grid'
         CALL FLUSH_ERR_MSG(abort = .true.)
      endif

      IF(CARTESIAN_GRID.and.USE_STL.and..not.use_stl_des) then 
         write(err_msg, '(3(A,/))')'Detected cartesian grid using STL representation', & 
         'USE_STL_DES either not specified or set as false in the input file', &
         'Forcing the discrete model to use STL for particle-wall interactions as well'
         CALL FLUSH_ERR_MSG
         USE_STL_DES = .true.
      endif
      
      IF(MPPIC.AND..NOT.USE_STL_DES) THEN 
         write(err_msg, '(3(A,/))') & 
         'PIC model detected but USE_STL_DES left undefined or specified as false', &
         'Particle-wall interactions in PIC model resoved as triangle-parcel', & 
         'Forcing USE_STL_DES to true for PIC model. Bounding box will be converted to facets.'
         CALL FLUSH_ERR_MSG
         USE_STL_DES = .true.
      ENDIF


      IF(.not.cartesian_grid.and.use_stl_des.and..not.des_convert_box_to_facets) then 
         write(err_msg, '(3(A,/))')'USE_STL_DES detected for particle-wall interactions', & 
         'but DES_CONVERT_BOX_TO_FACETS not set to true to convert bounding box to facets', &
         'Set DES_CONVERT_BOX_TO_FACETS to true and re-run'         
         CALL FLUSH_ERR_MSG(abort = .true.)
      endif
      
      IF(use_stl_des.and.use_cohesion) then 
         write(err_msg, '(3(A,/))') & 
         'The cohesion force model has not been implemented in new', &
         'routines for STL facet based particle-wall interactions', & 
         'This will be restored shortly. Sorry :('
         CALL FLUSH_ERR_MSG(abort = .true.)
      endif
      

! Verify that there are no internal obstacles.
!      IF(.NOT.CARTESIAN_GRID) THEN
!         DO K = KSTART1, KEND1 
!         DO J = JSTART1, JEND1
!         DO I = ISTART1, IEND1 
!            IJK  = FUNIJK(I,J,K)
!            IF(.NOT.FLUID_AT(IJK)) THEN
!               WRITE(ERR_MSG,1400)
!               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
!            ENDIF
!         ENDDO
!         ENDDO
!         ENDDO
!      ENDIF

! 1400 FORMAT('Error 1400: DES simulations cannot have defined ',       &
!         'internal obstacles.',/'Please correct the data file')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY
