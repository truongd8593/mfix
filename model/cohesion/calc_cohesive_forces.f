!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COHESIVE_FORCES                                   C
!  Purpose: Module to calculate cohesive forces using:                 C
!       a) square-well potential                                       C
!       b) Hamaker van der Waals model                                 C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_COHESIVE_FORCES

      USE discretelement


!-----Debugging
!      IF(CALLED.eq.COH_DEBUG_STEP)THEN
!        COHESION_DEBUG=.TRUE.
!      ELSE 
!         COHESION_DEBUG=.FALSE.
!      END IF

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *, '**START COHESIVE CALCS**'
      END IF

!-----Update search grids used to identify cohesive particle interactions
      CALL UPDATE_SEARCH_GRIDS

!-----COHESION MODELS

!-----Square-well potential
      IF(SQUARE_WELL)THEN
         CALL CALC_SQUARE_WELL      
      END IF

!-----van der Waals
      IF(VAN_DER_WAALS)THEN
        CALL CALC_VAN_DER_WAALS
      END IF


!-----Liquid bridging (Future work)



      IF(COHESION_DEBUG.gt.0)THEN
        PRINT *, '**END COHESIVE CALCS**'
      END IF
 
      END SUBROUTINE CALC_COHESIVE_FORCES
