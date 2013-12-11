!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_COHESION                                     !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check/set parameters for DES cohesion models.              !
!                                                                      !
!  Comments: Original code moved from CHECK_DES_DATA                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_COHESION

      use discretelement
      use mpi_utility

      IMPLICIT NONE

      IF(USE_COHESION) THEN

         IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1081)
            CALL MFIX_EXIT(myPE)            
         ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG,1082)
            CALL MFIX_EXIT(myPE)
         ELSEIF (SQUARE_WELL) THEN     
! note that square_well cohesion code may not function correctly with
! all aspects of the current DEM code so proceed at your own risk. 
            IF(DMP_LOG) WRITE(UNIT_LOG, 1080)
            IF (MASTER_WELL_DEPTH .EQ. UNDEFINED .AND. &
                MASTER_WALL_WELL_DEPTH .EQ. UNDEFINED .AND. &
                RADIUS_RATIO .EQ. UNDEFINED .AND. &
                WALL_RADIUS_RATIO .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1083)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF (VAN_DER_WAALS) THEN
            IF (VDW_INNER_CUTOFF .EQ. UNDEFINED .AND. &
                VDW_OUTER_CUTOFF .EQ. UNDEFINED .AND. &
                HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1084)
               CALL MFIX_EXIT(myPE)         
            ENDIF
            IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED .AND. &
                WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1085)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF (FACTOR_RLM < &
                1.d0+VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS)) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1100)
               IF(DMP_LOG) WRITE(*,1100)
            ENDIF
         ENDIF                
      ELSE
! override any of the following settings if cohesion not used              
          WALL_VDW_OUTER_CUTOFF = ZERO
          SQUARE_WELL = .FALSE.
          VAN_DER_WAALS = .FALSE.
      ENDIF

      RETURN


 1080 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'USE_COHESION has not been verified with the current',/10X,&
         'DEM code so proceed at your own risk.',/1X,70('*')/)

 1081 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Cannot use both SQUARE_WELL and VAN_DER_WAALS',/10X,&
         'cohesion models'/1X,70('*')/)            

 1082 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'A cohesion model has not been selected'/1X,70('*')/)

 1083 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH',/10X,&
         'RADIUS_RATIO and/or WALL_RADIUS_RATIO are undefined ',&
         'but must be',/10X,'defined when SQUARE_WELL=T',/1X,70('*')/)

 1084 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'VDW_INNER_CUTOFF, VDW_OUTER_CUTOFF, HAMAKER_CONSTANT',/10X,&
         'are undefined but must be defined when ',&
         'VAN_DER_WAALS=T',/1X,70('*')/)

 1085 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WALL_VDW_INNER_CUTOFF, WALL_HAMAKER_CONSTANT are ',&
         'undefined',/10X, 'but must be defined when ',&
         'VAN_DER_WAALS=T',/1X,70('*')/)

 1100 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: VDW_OUTER_CUTOFF outside of neighbor ',/10X,&
         'search distance. Increase FACTOR_RLM.',/1X,70('*')/) 

      END SUBROUTINE CHECK_DES_COHESION
