!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_ESC_COH_FORCE(I,J)                                C
!  Purpose: This module will calculate the attractive force            C
!           imposed on particles by an escaping-cohesive               C
!           interaction                                                C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE CALC_ESC_COH_FORCE(I,J)

!-----MODULES USED
      USE discretelement
      USE run
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,K
      DOUBLE PRECISION R(NDIM)
      DOUBLE PRECISION VEL_DIF(NDIM)
      DOUBLE PRECISION DIST
      DOUBLE PRECISION REL_VEL
      DOUBLE PRECISION CONST_A
      DOUBLE PRECISION CONST_B
      DOUBLE PRECISION MASS_NORM
      DOUBLE PRECISION NORM_REL_VEL
      DOUBLE PRECISION IMPULSE
      DOUBLE PRECISION COMBINED_WELL_DEPTH
      LOGICAL DEBUG_WRITE

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START ESCAPE CALC'
      END IF

      IF(J.gt.PARTICLES+DIMN*2)THEN
          STOP
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
        DIST=DIST+(R(K))**2
        VEL_DIF(K) = DES_VEL_NEW(K,I)-DES_VEL_NEW(K,J)
        REL_VEL=REL_VEL+R(K)*VEL_DIF(K)
      END DO
      DIST=sqrt(DIST)
      NORM_REL_VEL=REL_VEL/DIST
      MASS_NORM=PMASS(I)*PMASS(J)/(PMASS(I)+PMASS(J))
      COMBINED_WELL_DEPTH=(WELL_DEPTH(I)+WELL_DEPTH(J))/2.0

      IF(I.gt.PARTICLES)THEN
        COMBINED_WELL_DEPTH=WELL_DEPTH(I)
      END IF
      IF(J.gt.PARTICLES)THEN
        COMBINED_WELL_DEPTH=WELL_DEPTH(J)
      END IF


!-----Calculate momentum impulse transferred during cohesive interaction
      CONST_A=NORM_REL_VEL
      CONST_B=sqrt(CONST_A*CONST_A-2*COMBINED_WELL_DEPTH/MASS_NORM)
      IMPULSE=MASS_NORM*(CONST_A+CONST_B)


!-----Convert momentum impulse into force on particle
      DO K=1,DIMN
         FC(K,I)=FC(K,I)-IMPULSE*R(K)/(DIST*DTSOLID)
         FC(K,J)=FC(K,J)+IMPULSE*R(K)/(DIST*DTSOLID)
      END DO


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END ESCAPE CALC'
      END IF

      END SUBROUTINE CALC_ESC_COH_FORCE
