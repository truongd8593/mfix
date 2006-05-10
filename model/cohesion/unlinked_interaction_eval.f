!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UNLINKED_INTERACTION_EVAL(I,J)                         C
!  Purpose: Module to determine type of interaction for unlinked       C
!           particles                                                  C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE UNLINKED_INTERACTION_EVAL(I,J)

!-----MODULES USED
      USE discretelement
      USE run
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,K
      DOUBLE PRECISION R(3)
      DOUBLE PRECISION VEL_DIF(3)
      DOUBLE PRECISION DIST
      DOUBLE PRECISION REL_VEL
      DOUBLE PRECISION CONST_A
      DOUBLE PRECISION CONST_B
      DOUBLE PRECISION MASS_NORM

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START UNLINKED EVALUATION'
      END IF

!-----Calculate relative positions and velocities
      DIST=0
      REL_VEL=0
      DO K = 1, DIMN
        R(K) = DES_POS_NEW(J,K)-DES_POS_NEW(I,K)
        IF(DES_PERIODIC_WALLS) THEN
          IF(DES_PERIODIC_WALLS_X) THEN
             IF(K.eq.1) THEN
                IF(R(K).gt.(EX2-WX1)) THEN
                    R(K)=R(K)-(EX2-WX1)
                END IF
                IF(R(K).lt.-(EX2-WX1)) THEN
                    R(K)=R(K)+(EX2-WX1)
                END IF
             END IF
          END IF
          IF(DES_PERIODIC_WALLS_Y) THEN
             IF(K.eq.2) THEN
                IF(R(K).gt.(TY2-BY1)) THEN
                    R(K)=R(K)-(TY2-BY1)
                END IF
                IF(R(K).lt.-(TY2-BY1)) THEN
                    R(K)=R(K)+(TY2-BY1)
                END IF
             END IF
          END IF
          IF(DES_PERIODIC_WALLS_Z) THEN
             IF(K.eq.3) THEN
                IF(R(K).gt.(NZ2-SZ1)) THEN
                    R(K)=R(K)-(NZ2-SZ1)
                END IF
                IF(R(K).lt.-(NZ2-SZ1)) THEN
                    R(K)=R(K)+(NZ2-SZ1)
                END IF
             END IF
          END IF
        END IF
        VEL_DIF(K) = DES_VEL_NEW(I,K)-DES_VEL_NEW(J,K)
        REL_VEL=REL_VEL+R(K)*VEL_DIF(K)
        DIST=DIST+(R(K))**2
      END DO
      DIST=sqrt(DIST)

      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *, 'part 1=',I,'pos',DES_POS_NEW(I,1),DES_POS_NEW(I,2)
         PRINT *, 'part 2=',J,'pos',DES_POS_NEW(J,1),DES_POS_NEW(J,2)
         PRINT *, 'REL VEL=', REL_VEL,'DIST=',DIST&
            ,'WELL_WIDTH=',WELL_WIDTH(I), WELL_WIDTH(J)
      END IF

!-----Are particles inside square well and approaching?
      IF(REL_VEL.gt.0)THEN
        IF(DIST.lt.(WELL_WIDTH(I)+WELL_WIDTH(J)))THEN
           CALL CALC_APP_COH_FORCE(I,J)
           CALL ADD_PART_TO_LINK_LIST(I,J)
        END IF  !  INSIDE WELL CHECK
      END IF  ! REL VEL CHECK

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END UNLINKED EVALUATION'
      END IF

      END SUBROUTINE UNLINKED_INTERACTION_EVAL
