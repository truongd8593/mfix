!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_MPPIC                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MPPIC


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
!      use geometry, only: IMAX 
!      use geometry, only: JMAX
!      use geometry, only: KMAX
! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K

! Global Parameters:
!---------------------------------------------------------------------//
!      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MPPIC")




      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_SOLIDS_MPPIC
