!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_COHESION                                     !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check/set parameters for DES cohesion models.              !
!                                                                      !
!  Comments: Original code moved from CHECK_DES_COHESION               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_COHESION

! Global Variables:
!---------------------------------------------------------------------//

! Runtime Flag: Invoke a cohesion model for DES simulation
      use discretelement, only: USE_COHESION
! Largest discrete particle diameter.
      use discretelement, only: MAX_RADIUS
! Runtime Flag: Invoke Square well model.
      use discretelement, only: SQUARE_WELL
! Runtime Flag: Invoke Van der Waals model.
      use discretelement, only: VAN_DER_WAALS
! User specified parameter for increase neighbor search area.
      use discretelement, only: FACTOR_RLM
! File unit for LOG messages.
      use funits, only: UNIT_LOG

! User specified parameters for square-well model.
      use discretelement, only: MASTER_WELL_DEPTH
      use discretelement, only: MASTER_WALL_WELL_DEPTH
      use discretelement, only: RADIUS_RATIO
      use discretelement, only: WALL_RADIUS_RATIO

! User specified parameters for Van der Waals model.
      use discretelement, only: VDW_INNER_CUTOFF
      use discretelement, only: VDW_OUTER_CUTOFF
      use discretelement, only: HAMAKER_CONSTANT
      use discretelement, only: WALL_VDW_INNER_CUTOFF
      use discretelement, only: WALL_VDW_OUTER_CUTOFF
      use discretelement, only: WALL_HAMAKER_CONSTANT

      use mpi_utility

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Local integer error flag.
      INTEGER :: IER
! Neighborhood size for Van der Waals force.
      DOUBLE PRECISION :: VDW_NEIGHBORHOOD



! Override the following settings if cohesion not used.
      IF(.NOT.USE_COHESION) THEN
         SQUARE_WELL = .FALSE.
         VAN_DER_WAALS = .FALSE.
         WALL_VDW_OUTER_CUTOFF = ZERO
         RETURN
      ENDIF

! Verify that only one cohesion model is specified.
      IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
         IF(myPE == PE_IO) WRITE(*, 1000)
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         CALL MFIX_EXIT(myPE)            

! Verify that at a cohesion model is specified.
      ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
         IF(myPE == PE_IO) WRITE(*, 1001)
         IF(DMP_LOG) WRITE(UNIT_LOG,1001)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Squarewell model checks
      IF (SQUARE_WELL) THEN

! Warn that the Square Well model may not function properly.
         IF(myPE == PE_IO) WRITE(*, 1100)
         IF(DMP_LOG) WRITE(UNIT_LOG,1100)

         IER=0
         IF(MASTER_WELL_DEPTH .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'MASTER_WELL_DEPTH'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'MASTER_WELL_DEPTH'
            IER = 1
         ENDIF
         IF(MASTER_WALL_WELL_DEPTH .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'MASTER_WALL_WELL_DEPTH'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'MASTER_WALL_WELL_DEPTH'
            IER = 1
         ENDIF

         IF(RADIUS_RATIO .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'RADIUS_RATIO'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'RADIUS_RATIO'
            IER = 1

         ENDIF
         IF(WALL_RADIUS_RATIO .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1101) 'WELL_RADIUS_RATIO'
            IF(DMP_LOG) WRITE(UNIT_LOG,1101) 'WELL_RADIUS_RATIO'
            IER = 1
         ENDIF
         IF(IER /= 0) CALL MFIX_EXIT(myPE)
      ENDIF


! Van der Waals model checks.
      IF (VAN_DER_WAALS) THEN

         IER = 0
         IF (VDW_INNER_CUTOFF .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'VDW_INNER_CUTOFF'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'VDW_INNER_CUTOFF'
            IER = 1
         ENDIF

         IF(VDW_OUTER_CUTOFF .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'VDW_OUTER_CUTOFF'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'VDW_OUTER_CUTOFF'
            IER = 1
         ENDIF

         IF(HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'HAMAKER_CONSTANT'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'HAMAKER_CONSTANT'
            IER = 1
         ENDIF

         IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED)THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'WALL_VDW_INNER_CUTOFF'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'WALL_VDW_INNER_CUTOFF'
            IER = 1
         ENDIF

         IF (WALL_VDW_OUTER_CUTOFF .EQ. UNDEFINED)THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'WALL_VDW_OUTER_CUTOFF'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'WALL_VDW_OUTER_CUTOFF'
            IER = 1
         ENDIF

         IF(WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            IF(myPE == PE_IO) WRITE(*, 1201) 'WALL_HAMAKER_CONSTANT'
            IF(DMP_LOG) WRITE(UNIT_LOG,1201) 'WALL_HAMAKER_CONSTANT'
            IER = 1
         ENDIF

         VDW_NEIGHBORHOOD = 1.0d0 + (VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS))
         IF (FACTOR_RLM < VDW_NEIGHBORHOOD) THEN
            IF(myPE == PE_IO) WRITE(*, 1202)
            IF(DMP_LOG) WRITE(UNIT_LOG,1202)
            IER = 1
         ENDIF

         IF(IER /= 0) CALL MFIX_EXIT(myPE)
      ENDIF                


      RETURN


 1000 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Error 1000:',  &
         ' Cannot use SQUARE_WELL and VAN_DER_WAALS cohesion',/        &
         ' models simultaneously.',/1X,70('*')/)            

 1001 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Error 1001:',  &
         ' A cohesion model was not selected. Specify one of the',/    &
         ' available models in the mfix.dat file.',/1X,70('*')/)


!<-------------------- Square-well model messages. ------------------->!


 1100 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Warning 1100:',&
         ' The square-well cohesion model has not been verified',/     &
         ' with the current DEM implemenation.',/1X,70('*')/)

 1101 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Error 1101:',  &
         ' Missing input data for square-well cohesion model.',/       &
         ' Input parameter ',A,' is UNDEFINED.',/1X,70('*')/)


!<------------------- Van der Waals model messages. ----------------->!


 1201 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Error 1201:',  &
         ' Missing input data for Van der Waals cohesion model.',/     &
         ' Input parameter ',A,' is UNDEFINED.',/1X,70('*')/)

 1202 FORMAT(/1X,70('*')/' From: CHECK_DES_COHESION',/' Error 1202: ', &
         ' VDW_OUTER_CUTOFF outside of the neighbor search distance.',/&
         ' Increase FACTOR_RLM to increase the search distance.',/     &
         1X,70('*')/)

      END SUBROUTINE CHECK_DES_COHESION
