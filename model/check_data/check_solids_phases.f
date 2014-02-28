!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_PHASES                                     !
!  Purpose: Driver routine for calls to solids phase checks.           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_PHASES


! Global Variables:
!---------------------------------------------------------------------//
! Runtime flag specifying TFM solids
      use run, only: TFM_SOLIDS
! Runtime flag specifying DEM solids
      use run, only: DEM_SOLIDS
! Runtime flag specifying MPPIC solids
      use run, only: PIC_SOLIDS

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_PHASES")

! Checks common to all solids models.
      CALL CHECK_SOLIDS_COMMON_ALL

! Checks common to discrete solids phases (DEM, MPPIC).
      IF(DEM_SOLIDS .OR. PIC_SOLIDS) &
         CALL CHECK_SOLIDS_COMMON_DISCRETE

! Checks specific to the particular solids phase.
      IF(TFM_SOLIDS) CALL CHECK_SOLIDS_CONTINUUM
      IF(DEM_SOLIDS) CALL CHECK_SOLIDS_DES
      IF(PIC_SOLIDS) CALL CHECK_SOLIDS_MPPIC

      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_SOLIDS_PHASES
