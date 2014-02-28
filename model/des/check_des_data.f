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

 1101 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'ERROR: CLUSTER_LENGTH_CUTOFF must be defined when ',/10X,&
         'using function DES_CALC_CLUSTER.',/1X,70('*')/)          

 1102 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: CLUSTER_LENGTH_CUTOFF outside of neighbor ',/10X,&
         'search distance. Increase FACTOR_RLM.',/1X,70('*')/) 

         END SUBROUTINE CHECK_DES_DATA
