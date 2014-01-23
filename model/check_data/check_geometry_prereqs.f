!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY_PREREQS                                  !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY_PREREQS


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: IMAX 
      use geometry, only: JMAX
      use geometry, only: KMAX
! Runtime flag specifying 2D simulations
      use geometry, only: NO_K

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY_PREREQS")

! Verify that the domain decomposition was specifed.
      IF(IMAX == UNDEFINED_I .OR. JMAX == UNDEFINED_I .OR.             &
         (.NOT.NO_K .AND. KMAX == UNDEFINED_I) ) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF



      CALL FINL_ERR_MSG

      RETURN  


 1000 FORMAT('Error 1000: IMAX or JMAX or KMAX not specified in ',     &
          'mfix.dat')

      END SUBROUTINE CHECK_GEOMETRY_PREREQS
