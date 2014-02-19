!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!  Author: P.Nicoletti                                Date: 27-NOV-91  !
!          J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_RUN_CONTROL 


! Global Variables:
!---------------------------------------------------------------------//
! New or restart
      USE run, only: RUN_TYPE
! Brief description of simulation.
      USE run, only: DESCRIPTION
! Simulation units: SI, CGS
      USE run, only: UNITS
! Simulation start/stop times.
      USE run, only: TIME, TSTOP
! Time step size, one over time step size.
      USE run, only: DT, ODT
! Flag: Use K-Epsilon turbulence model.
      USE run, only: K_EPSILON
! Turbulence lenghth scale.
      use constant, only: L_SCALE0
! Flag: Solve species eq.
      USE run, only: SPECIES_EQ, ANY_SPECIES_EQ
! Radial distribution function
      USE run, only: RDF_TYPE
! Number of solids phases.
      USE physprop, only: MMAX, SMAX
! Number of scalar equations to solve
      USE scalars, only: NSCALAR, phase4scalar
! Phase associated with scalar transport.
      USE scalars, only: phase4scalar
! Number of discrete solids phases.
      use discretelement, only: DES_MMAX

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED, UNDEFINED_C
      USE param1, only: ONE, ZERO

! Global Module proceedures:
!---------------------------------------------------------------------//
      USE error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Dummy integer
      INTEGER :: N


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_RUN_CONTROL")


! Clear out the run description if not specified.
      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' ' 

! Verify UNITS input.
      IF(UNITS == UNDEFINED_C) THEN
         WRITE(ERR_MSG,1000) 'UNITS'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF((UNITS /= 'CGS') .AND. (UNITS /= 'SI')) THEN
         WRITE(ERR_MSG,1001) 'UNITS', UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify that DT is valid.
      IF (DT < ZERO) THEN
         WRITE(ERR_MSG,1002) 'DT', DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Steady-state simulation.
      ELSEIF(DT == UNDEFINED .OR. DT == ZERO) THEN
         ODT = ZERO 
         TIME = ZERO 

! Transient simulation.
      ELSE
! Calculate one over the initial timestep.
         ODT = ONE/DT 
! Verify the remaining time settings.
         IF (TIME == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TIME'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TSTOP == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TSTOP'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TIME < ZERO) THEN
            WRITE(ERR_MSG,1002)'TIME', TIME
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TSTOP < ZERO) THEN
            WRITE(ERR_MSG,1002) 'TSTOP', TSTOP
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF

! Verify the run type.
      IF(.NOT.(RUN_TYPE=='NEW' .OR. RUN_TYPE=='RESTART_1'              &
         .OR. RUN_TYPE=='RESTART_2')) THEN
         WRITE(ERR_MSG,1001) 'RUN_TYPE', RUN_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Turbulence model:
      IF (K_Epsilon .AND. L_SCALE0 /= ZERO) THEN
         WRITE(ERR_MSG,2001) 
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 2001 FORMAT('Error 2001: Cannot set K_EPSILON = .T. and specify ',    &
         'L_SCALE0 /= ZERO')
      ENDIF

! Set variable ANY_SPECIES_EQ
      N = max(SMAX, DES_MMAX)
      ANY_SPECIES_EQ = any(SPECIES_EQ(:N))
     
! Check phase specification for Scalars
      DO N = 1, NScalar
         IF(Phase4Scalar(N) < 0 .OR. Phase4Scalar(N) > MMAX) THEN
            WRITE(ERR_MSG,1002) iVar('Phase4Scalar',N), phase4scalar(N)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Checks requred for the subgird drag models.
      CALL CHECK_SUBGRID_MODEL

! Clear the error manager
      CALL FINL_ERR_MSG
      

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_RUN_CONTROL 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SUBGRID_MODEL                                     !
!  Purpose: Check the subgrid drag model interactions.                 !
!                                                                      !
!  Author: J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SUBGRID_MODEL

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Specify friction model (Schaeffer model/Princeton model)
      USE run, only: FRICTION
! Flag: Solve granular energy eq
      USE run, only: GRANULAR_ENERGY
! Flag: SOlve K-Epsilong Eq.
      USE run, only: K_EPSILON
! Flag: Impose a mean shear on flow field.
      USE run, only: SHEAR
! Flag: Invoke Schaeffer and KT-Theory blending
      USE run, only: BLENDING_STRESS
! User specifed drag model
      USE run, only: DRAG_TYPE
! Ratio of filter size to computational cell size
      USE run, only: FILTER_SIZE_RATIO
! User specifed subgrid model: IGCI or MILIOLI
      USE run, only: SUBGRID_TYPE
! Flag: Include wall effect term
      USE run, only: SUBGRID_WALL
! Initial turbulcence length scale
      use constant, only: L_SCALE0
! Specularity coefficient for particle-wall collisions
      use constant, only: PHIP
! Flag: Use cartesian grid model
      USE cutcell, only : CARTESIAN_GRID
! Flag: Use discrete element solids model
      use discretelement, only: DISCRETE_ELEMENT
! Flag: Use MP-PIC solids model
      use mfix_pic, only: MPPIC

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: ZERO, UNDEFINED_C


! Global Module proceedures:
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! NONE

! If the modles are not being used, return.
      IF(SUBGRID_TYPE == UNDEFINED_C .AND. .NOT.SUBGRID_WALL) RETURN


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SUBGRID_MODEL")


      IF(SUBGRID_TYPE == UNDEFINED_C .AND. SUBGRID_WALL) THEN
         WRITE(ERR_MSG,2011)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2011 FORMAT('Error 2011: Invalid input. SUBGRID_WALL cannot be used ',&
          'without',/'specifying a SUBGRID_TYPE.',/'Please correct ',  &
          'the mfix.dat file.')

      IF(SUBGRID_TYPE /= 'IGCI' .AND. SUBGRID_TYPE /= 'MILIOLI') THEN
         WRITE(ERR_MSG,1001) 'SUBGRID_TYPE', SUBGRID_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(DRAG_TYPE /= 'WEN_YU')THEN
         WRITE(ERR_MSG, 2012)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2012 FORMAT('Error 2012: Invalid input. WEN_YU is the only DRAG_TYPE',&
          ' available',/'when using the SUBGRID model.',/'Please ',    &
          'correct the mfix.dat file.')

      IF(DISCRETE_ELEMENT .OR. MPPIC) THEN
         WRITE(ERR_MSG, 2013)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2013 FORMAT('Error 2013: Invalid input. The SUBRID model is not ',    &
          'available',/'with discrete solids phases.',/'Please ',      &
          'correct the mfix.dat file.')

! Impose the subgrid limitations.
      IF(FILTER_SIZE_RATIO <= ZERO) THEN
         WRITE(ERR_MSG, 1002)'FILTER_SIZE_RATIO', FILTER_SIZE_RATIO
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(GRANULAR_ENERGY) THEN
         WRITE(ERR_MSG, 2010) 'GRANULAR_ENERGY', 'FALSE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(K_EPSILON) THEN
         WRITE(ERR_MSG, 2010) 'K_EPSILON', 'FALSE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(BLENDING_STRESS) THEN
         WRITE(ERR_MSG, 2010) 'BLENDING_STRESS', 'FALSE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(FRICTION) THEN
         WRITE(ERR_MSG, 2010) 'FRICTION', 'FALSE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(SHEAR) THEN
         WRITE(ERR_MSG, 2010) 'SHEAR', 'FALSE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(PHIP /= ZERO) THEN
         WRITE(ERR_MSG, 2010) 'PHIP', 'ZERO'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ENDIF

      IF(SUBGRID_WALL .AND. .NOT.CARTESIAN_GRID) THEN
         WRITE(ERR_MSG, 2010) 'CARTESIAN_GRID', 'TRUE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2010 FORMAT('Error 2010: Invalid input. ',A,' must be 'A,/'when ',    &
         'using the SUBGRID model.'/,'Please correct the mfix.dat',    &
         ' file.')

      CALL FINL_ERR_MSG

      RETURN


 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SUBGRID_MODEL
