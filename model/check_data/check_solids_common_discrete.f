!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


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
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE")


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




      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE
