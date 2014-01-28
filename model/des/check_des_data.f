!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DES_DATA                                          C
!  Purpose: Check user input data                                      C
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 10-Nov-08  C
!  Reviewer:                                                           C
!  Comments: Most all user's input data are checked here               C
!  Revision: Some of the checks made in des_allocate_arrays are        C
!            moved here. In addition write statments are now made to   C
!            log file                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DES_DATA

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
! Runtime Flag: Printing screen messages. PE_IO Only.
      USE discretelement, only: PRINT_DES_SCREEN

! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Number of solids phases.
      USE physprop, only: MMAX
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
! DES Neighborhood size.
      USE discretelement, only: FACTOR_RLM
! Max/Min particle radi
      USE discretelement, only: MAX_RADIUS, MIN_RADIUS
! Max distance between two particles in a cluster.
      USE discretelement, only: CLUSTER_LENGTH_CUTOFF
! Subroutine access.
      USE desgrid, only: DESGRID_CHECK
! File unit for LOG messages.
      USE funits, only: UNIT_LOG

      USE mpi_utility

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M  ! Solids phase Index


! Set flag for screen messages.
      PRINT_DES_SCREEN = (myPE.EQ.pe_IO)
      IF(DMP_LOG) WRITE(UNIT_LOG,1001) PRINT_DES_SCREEN
      IF(DMP_LOG) WRITE(*,1001) PRINT_DES_SCREEN

      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
! MMAX, D_p0 and RO_s0 are to be strictly associated with the continuum
! model and are no longer to be used in the DEM. However, when not using
! DES_CONTINUUM_HYBRID the current DEM code still assumes that the user
! specifies MMAX, D_P0 and RO_S0.  These are then linked to their
! respective DES variables.  Note that valid values of MMAX, D_P0 and 
! RO_s0 are ensured by check_data_04.
         DES_MMAX = MMAX
         DO M = 1, DES_MMAX
            DES_D_p0(M) = D_p0(M)
            DES_RO_s(M) = RO_s0(M)
         ENDDO
      ENDIF

      IF (DES_MMAX == UNDEFINED_I) THEN
         IF (DMP_LOG) WRITE(UNIT_LOG, 1090)
         CALL MFIX_EXIT(myPE)
      ELSE
! The following checks on DES_mmax, density and diameter are similar to
! those in check_data_04.  If .not.des_continuum_hybrid then these
! checks are basically redundant, however, they will facilitate a move
! to make DEM operate separately from MFIX's check_data_04.  
! ---------------------------------------------------------------->>>    
! Check mmax: 
         IF (DES_MMAX<0 .OR. DES_MMAX>DIM_M) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG, 1070)
            CALL MFIX_EXIT(myPE)
         ENDIF
! Check for valid diameter values. 
! Valid diameters are needed: 1) when using gener_part_config the
! diameters are used to identify the number of particles in a given
! solids phase and 2) to identify which solids phase each particle
! belongs in (this sorting is conducted in particles_in_cell)
         DO M = 1,DES_MMAX
            IF (DES_D_P0(M)<ZERO .OR. DES_D_P0(M)==UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1043)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = DES_MMAX+1, DIM_M
            IF (DES_D_P0(M) /= UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1044)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! Check for valid density values
         DO M = 1,DES_MMAX
            IF (DES_RO_S(M)<ZERO .OR. DES_RO_S(M)==UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1071)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
         DO M = DES_MMAX+1, DIM_M
            IF (DES_RO_S(M) /= UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG, 1072)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO

      ENDIF   ! end if/else(des_mmax==undefined_i)
! End checks on DES_MMAX, density and diameter      
! ----------------------------------------------------------------<<<

! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks        
      MAX_RADIUS = ZERO
      MIN_RADIUS = LARGE_NUMBER
      DO M = 1,DES_MMAX
         MAX_RADIUS = MAX(MAX_RADIUS, 0.5d0*DES_D_P0(M))
         MIN_RADIUS = MIN(MIN_RADIUS, 0.5d0*DES_D_P0(M))
      ENDDO

! Overwrite user's input in case of DEM (no fluid)
      IF(.NOT.DES_CONTINUUM_COUPLED) DES_INTERP_ON = .FALSE.
      
! Check for valid neighbor search option      
      IF(DES_NEIGHBOR_SEARCH.EQ.2)  then 
         DES_NEIGHBOR_SEARCH = 4 
         IF(DMP_LOG) WRITE(unit_log, 1046) 2, 'QUADTREE'
      ELSEIF (DES_NEIGHBOR_SEARCH.EQ.3) THEN
         DES_NEIGHBOR_SEARCH = 4 
         IF(DMP_LOG) WRITE(unit_log, 1046) 3, 'OCTREE'
      ENDIF
      
      IF(DES_NEIGHBOR_SEARCH.EQ.1) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 1, 'N-SQUARE'
      ELSEIF(DES_NEIGHBOR_SEARCH.EQ.4) THEN
         IF(DMP_LOG) WRITE(unit_log, 1045) 4, 'GRID BASED'
      ELSE
         IF(DMP_LOG) WRITE(UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)         
      ENDIF

! Check the output file format 
      IF (DES_OUTPUT_TYPE /= UNDEFINED_C) THEN
         IF(TRIM(DES_OUTPUT_TYPE) /= 'TECPLOT') THEN
            IF(DMP_LOG) WRITE (*, 1038)
            CALL MFIX_EXIT(myPE)   
         ENDIF
      ENDIF

! Check for valid integration method
      IF (TRIM(DES_INTG_METHOD) == 'ADAMS_BASHFORTH') THEN
         INTG_ADAMS_BASHFORTH = .TRUE.
         INTG_EULER = .FALSE.
      ELSEIF(TRIM(DES_INTG_METHOD) == 'EULER') THEN
         INTG_ADAMS_BASHFORTH = .FALSE.
         INTG_EULER = .TRUE.
! stop if the specified integration method is unavailable
      ELSE
         IF(DMP_LOG) WRITE (UNIT_LOG, 1034)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Check settings on cluster identification
      IF(DES_CALC_CLUSTER) THEN
         IF(CLUSTER_LENGTH_CUTOFF .EQ. UNDEFINED) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1101)
            IF(DMP_LOG) WRITE(*,1101)
            CALL MFIX_EXIT(myPE)
         ENDIF
         IF(FACTOR_RLM < &
            1.d0+CLUSTER_LENGTH_CUTOFF/(2.d0*MAX_RADIUS)) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1102)
            IF(DMP_LOG) WRITE(*,1102)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

! Check geometry constrains.
      CALL CHECK_DES_GEOMETRY
! Set flags for energy equations
      CALL CHECK_DES_ENERGY
! Check thermodynamic properties of discrete solids.
      CALL CHECK_DES_THERMO
! Check settings on cohesion model
      CALL CHECK_DES_COHESION
! Check settings for collision models
      IF(.NOT.MPPIC) CALL CHECK_DES_COLLISION
! Check TFM/DEM Hybrid model settings.
      IF (DES_CONTINUUM_HYBRID) CALL CHECK_DES_HYBRID
! Check settings for particle generation.
      IF(GENER_PART_CONFIG) CALL CHECK_DES_PCONFIG
! Check quantities related to MPPIC      
      IF(MPPIC) CALL CHECK_DES_MPPIC

! Check coupling settings.
      CALL CHECK_DES_COUPLING
           
! the entire checking and setting up indices for desgridsearch
! moved to desgrid_mod to accomodate parallelization
! this is now conducted regardless of neighbor search option      
      CALL DESGRID_CHECK


      RETURN

 1001 FORMAT(/2X,'FULL_LOG = ',L3,' FOR DES ON IO PROC') 

 1003 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_NEIGHBOR_SEARCH in mfix.dat',/10X,&
         'Must be > 0 or < 5',/1X,70('*')/)




 1034 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_INTG_METHOD in mfix.dat.  Must be',&
         /10X,'EULER (default/undefined) or ADAMS_BASHFORTH',&
         /1X,70('*')/)

  
 1038 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'The only option for DES_OUTPUT_DATA is TECPLOT',&
          /1X,70('*')/)

         
 1043 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'D_P0 must be defined and >0 in mfix.dat for M=1,',&
         'MMAX.',/1X,70('*')/)
 1044 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Too many D_P0 are defined in mfix.dat for given ',&
         'MMAX.',/1X,70('*')/)

 1045 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_NEIGHBOR_SEARCH set to ', I2, ' ',A,/1X,70('*')/)

 1046    FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_NEIGHBOR_SEARCH specified as ', I2, ' ',A,/, &
         'forced to grid based neighbor search  (4)', /1X,70('*')/)

 1054 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
          'does not equal D_P0(1).',/1X,70('*'))

 1055 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'WARNING: zlength or dz(1) is used to calculate the ',&
          'number',/10X,'of particles in the 2D simulation when ',&
          'GENER_PART_CONFIG is T and DIMN = 2.',/1X,70('*'))

 1056     FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'DES_EPS_XSTART, DES_EPS_YSTART, and DES_EPS_ZSTART are obsolete flags.', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))

 1057     FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
          'Gener_part_config is true but using deprecated flags', /, &
          'VOL_FRAC is an obsolete flag', /, &
          'IC region is now based on usual IC_flags', /, & 
          'Delete these flags from mfix.dat and restart', /, &
          'Exitting',/1X,70('*'))

 
 1070 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'MMAX not specified or unphysical in mfix.dat.' ,&
         /1X,70('*')/)

 1071 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'RO_s must be defined and >0 in mfix.dat for M=1,',&
         'MMAX.',/1X,70('*')/)

 1072 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Too many RO_s are defined in mfix.dat for given ',&
         'MMAX.',/1X,70('*')/)

 1090 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_MMAX not defined in mfix.dat. It must be defined',/10X,&
         'when using DES_CONTINUUM_HYBRID',/1X,70('*')/)

 1101 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'ERROR: CLUSTER_LENGTH_CUTOFF must be defined when ',/10X,&
         'using function DES_CALC_CLUSTER.',/1X,70('*')/)          

 1102 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: CLUSTER_LENGTH_CUTOFF outside of neighbor ',/10X,&
         'search distance. Increase FACTOR_RLM.',/1X,70('*')/) 

         END SUBROUTINE CHECK_DES_DATA
