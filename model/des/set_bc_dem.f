!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM

      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE mfix_pic
      use mpi_utility

      use bc

      use error_manager

      IMPLICIT NONE


      CALL INIT_ERR_MSG("SET_BC_DEM")

! The variable PARTICLES should already be set by this point if using
! gener_part_config option      
      IF(PARTICLES == UNDEFINED_I .AND. MAX_PIS /= UNDEFINED_I)THEN
         PARTICLES = 0

      ELSEIF(PARTICLES == UNDEFINED_I .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1200 FORMAT('Error 1200: Either PARTICLES or MAX_PIS must specified.',&
         'Please correct the mfix.dat file.')

      ELSEIF(PARTICLES == 0 .AND. MAX_PIS == UNDEFINED_I) THEN
         WRITE(ERR_MSG, 1201)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1201 FORMAT('Error 1201: MAX_PIS must be specified in the mfix.dat ', &
         'file if',/' there are no initial particles (PARTICLES = 0).')

      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
! Inlet/outlet for MPPIC are based off the regular mfix declarations, 
! and so DEM_BCMI could still be zero.
      IF(PARTICLES == 0 .AND. DEM_BCMI == 0) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG
      ENDIF

 1202 FORMAT('WARNING 1202: The system is initiated with no particles',&
         ' and no',/'solids inlet was detected.')


! Check MAX_PIS requirements
      IF(DEM_BCMI == 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         MAX_PIS = PARTICLES
      ELSEIF(DEM_BCMI /= 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(ERR_MSG, 1203)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1203 FORMAT('Error 1203: The maximum number of particles (MAX_PIS) ', &
         'must be',/'spedified if a DEM inlet is specifed. Please ',   &
         'correct the mfix.dat ',/'file.')

      IF(DEM_BCMI > 0) CALL SET_BC_DEM_MI
      IF(DEM_BCMO > 0) CALL SET_BC_DEM_MO

! Set the flag that one or more DEM MI/MO exists.
      DEM_MIO = (DEM_BCMI /= 0 .OR. DEM_BCMO /= 0)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_DEM
