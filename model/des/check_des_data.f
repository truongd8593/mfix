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
! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only: GENER_PART_CONFIG
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Invoke TFM/DEM hybrid model.
      USE discretelement, only: DES_CONTINUUM_HYBRID
! Runtime Flag: Printing screen messages. PE_IO Only.
      USE discretelement, only: PRINT_DES_SCREEN
! Subroutine access.
      USE desgrid, only: DESGRID_CHECK
! File unit for LOG messages.
      USE funits, only: UNIT_LOG

      USE mpi_utility

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! None.


! Set flag for screen messages.
      PRINT_DES_SCREEN = (myPE.EQ.pe_IO)
      IF(DMP_LOG) WRITE(UNIT_LOG,1001) PRINT_DES_SCREEN
      IF(DMP_LOG) WRITE(*,1001) PRINT_DES_SCREEN

! Check TFM/DEM Hybrid model settings.
      IF (DES_CONTINUUM_HYBRID) CALL CHECK_DES_HYBRID
! Check settings for particle generation.
      IF(GENER_PART_CONFIG) CALL CHECK_DES_PCONFIG

! the entire checking and setting up indices for desgridsearch
! moved to desgrid_mod to accomodate parallelization
! this is now conducted regardless of neighbor search option      
      CALL DESGRID_CHECK


      RETURN

 1001 FORMAT(/2X,'FULL_LOG = ',L3,' FOR DES ON IO PROC') 


      END SUBROUTINE CHECK_DES_DATA
