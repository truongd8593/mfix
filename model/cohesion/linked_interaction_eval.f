!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LINKED_INTERACTION_EVAL(I,J)                           C
!  Purpose: Module to determine type of interaction for                C
!           linked particles                                           C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE LINKED_INTERACTION_EVAL(I,J)

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,K
      INTEGER CHECK_AGG
      DOUBLE PRECISION R(NDIM)
      DOUBLE PRECISION VEL_DIF(NDIM)
      DOUBLE PRECISION DIST
      DOUBLE PRECISION REL_VEL
      DOUBLE PRECISION CONST_A
      DOUBLE PRECISION CONST_B
      DOUBLE PRECISION MASS_NORM
      DOUBLE PRECISION COMBINED_WELL_DEPTH

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START LINKED EVALUATION'
      END IF



!-----Calculate relative positions and velocities
      DIST=0
      REL_VEL=0
      DO K = 1, DIMN
        R(K) = DES_POS_NEW(K,J)-DES_POS_NEW(K,I)
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
        VEL_DIF(K) = DES_VEL_NEW(K,I)-DES_VEL_NEW(K,J)
        REL_VEL=REL_VEL+R(K)*VEL_DIF(K)
        DIST=DIST+(R(K))**2
      END DO
      DIST=sqrt(DIST)
      COMBINED_WELL_DEPTH=(WELL_DEPTH(I)+WELL_DEPTH(J))/2.0

      IF(I.gt.PARTICLES)THEN
        COMBINED_WELL_DEPTH=WELL_DEPTH(I)
      END IF
      IF(J.gt.PARTICLES)THEN
        COMBINED_WELL_DEPTH=WELL_DEPTH(J)
      END IF


!----Are particles outside the square well and departing?
      IF(REL_VEL.lt.0)THEN
        IF(DIST.gt.(WELL_WIDTH(I)+WELL_WIDTH(J)))THEN
          CONST_A=REL_VEL/DIST
          MASS_NORM=PMASS(I)*PMASS(J)/(PMASS(I)+PMASS(J))
          CONST_B=CONST_A**2-2*COMBINED_WELL_DEPTH/MASS_NORM
          IF(CONST_B.gt.0)THEN
            CALL CALC_ESC_COH_FORCE(I,J)
            CALL REMOVE_PART_FROM_LINK_LIST(I,J)
          ELSE
            CALL CALC_CAP_COH_FORCE(I,J)
          END IF
        END IF
      END IF

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END LINKED EVALUATION'
      END IF

      END SUBROUTINE LINKED_INTERACTION_EVAL
