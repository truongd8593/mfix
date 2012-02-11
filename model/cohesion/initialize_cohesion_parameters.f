!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INITIALIZE_COHESION_PARAMETERS                         C
!  Purpose: This module gives initial values to all                    C
!           cohesion parameters                                        C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE INITIALIZE_COHESION_PARAMETERS

!------------------------------------------------
! Modules      
!------------------------------------------------      
      USE constant
      USE discretelement
      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
      INTEGER :: I,J
!------------------------------------------------

      IF(COHESION_DEBUG.gt.0)THEN
        WRITE(*,*) '**START INITIALIZE COHESION PARAMETERS'
      ENDIF

    
      IF (SQUARE_WELL) THEN
         DO I=1, PARTICLES+2*DIMN
            IS_LINKED(I)=0
            LINKS(I,1)=0
            DO J=2,MAXNEIGHBORS
              LINKS(I,J)=-1
            ENDDO
            WELL_WIDTH(I)=RADIUS_RATIO*DES_RADIUS(I)
            WELL_DEPTH(I)=MASTER_WELL_DEPTH
            IF(I.GT.PARTICLES)THEN ! Walls
               WELL_WIDTH(I)=WALL_RADIUS_RATIO*MAX_RADIUS
               WELL_DEPTH(I)=MASTER_WALL_WELL_DEPTH
               PMASS(I)=99999999.d0
            ENDIF
         ENDDO
      ENDIF   ! square_well

      IF (VAN_DER_WAALS) THEN
! Surface energy set so that force stays constant at inner cut off
         SURFACE_ENERGY=HAMAKER_CONSTANT/&
            (24.d0*Pi*VDW_INNER_CUTOFF*VDW_INNER_CUTOFF)
         WALL_SURFACE_ENERGY=WALL_HAMAKER_CONSTANT/&
            (24.d0*Pi*WALL_VDW_INNER_CUTOFF*WALL_VDW_INNER_CUTOFF)
      ENDIF   ! end if van_der_waals


      MAX_PART_IN_GRID=3


      IF(COHESION_DEBUG.gt.0)THEN
        WRITE(*,*) '**END INITIALIZE COHESION PARAMETERS'
      ENDIF


      END SUBROUTINE INITIALIZE_COHESION_PARAMETERS
