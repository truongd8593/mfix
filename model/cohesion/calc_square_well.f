!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_SQUARE_WELL                                       C
!  Purpose: Module to calculate cohesive forces between particles      C
!    using a square-well potential                                     C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C



      SUBROUTINE CALC_SQUARE_WELL
      
!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,N
      INTEGER LINK
      INTEGER CHECK_LINK
      INTEGER X_INDEX, Y_INDEX, Z_INDEX
      INTEGER Z_START, Z_END
      INTEGER XX, YY, ZZ, KK
      LOGICAL NOT_CHECKED(PARTICLES+6)
      INTEGER CHECKED_PARTICLES_LIST(100)
      LOGICAL NOT_WALL
      INTEGER MAX_CHECK
      LOGICAL ALREADY_EXISTS

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START SQUARE WELL CALCS'
      END IF

!!The checked particles list is necessary to guard against redundant
!!  calculations with the wall since a particle will neighbor three
!!  wall cells if it is near a wall.  This list also guards against
!!  a linked particle moving outside of the search area.

      DO J=1,PARTICLES+2*DIMN
         NOT_CHECKED(J)=.TRUE.
      END DO

      IF(DIMN.eq.3)THEN
         MAX_CHECK=26*MAX_PART_IN_GRID
      ELSE
         MAX_CHECK=8*MAX_PART_IN_GRID
      END IF

      DO J=1, MAX_CHECK
         CHECKED_PARTICLES_LIST(J)=-1
      END DO

!-----START LOOP OVER ALL PARTICLES
      DO I = 1, PARTICLES
        N=1

        IF(COHESION_DEBUG.gt.1)THEN
           PRINT *,'Particle=',I
	   PRINT *,'NEIGHBORS=',(NEIGHBOURS(I,J), J=1, MAXNEIGHBORS)
	   PRINT *,'LINKS=',(LINKS(I,J), J=1, MAXNEIGHBORS)
        END IF

!-------Check all neighboring grids for particles

        IF(DIMN.eq.2)THEN
           Z_START=1
           Z_END=1
        ELSE
           Z_START=PART_GRID(I,3)-1
           Z_END=PART_GRID(I,3)+1
        END IF

      
        DO X_INDEX=PART_GRID(I,1)-1,PART_GRID(I,1)+1
           DO Y_INDEX=PART_GRID(I,2)-1,PART_GRID(I,2)+1
             DO Z_INDEX=Z_START,Z_END

             XX=X_INDEX
             YY=Y_INDEX
             ZZ=Z_INDEX
             NOT_WALL=.TRUE.

             IF(XX.eq.0)THEN !West Wall
               IF(DES_PERIODIC_WALLS_X) THEN
                 XX=SEARCH_GRIDS(1)
               ELSE
                  IF(NOT_CHECKED(PARTICLES+1))THEN !!guard against multiple intrxns 
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+1)
                    NOT_CHECKED(PARTICLES+1)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+1
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(YY.eq.0)THEN !Bottom WAll
               IF(DES_PERIODIC_WALLS_Y) THEN
                 YY=SEARCH_GRIDS(2)
               ELSE
                  IF(NOT_CHECKED(PARTICLES+2))THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+2)
                    NOT_CHECKED(PARTICLES+2)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+2
                    N=N+1
                    NOT_WALL=.FALSE.
                  END IF
               END IF
             END IF

             IF(XX.eq.SEARCH_GRIDS(1)+1)THEN !East Wall
               IF(DES_PERIODIC_WALLS_X) THEN
                  XX=1
               ELSE
                  IF(NOT_CHECKED(PARTICLES+3))THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+3)
                    NOT_CHECKED(PARTICLES+3)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+3
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(YY.eq.SEARCH_GRIDS(2)+1)THEN !Top Wall
               IF(DES_PERIODIC_WALLS_Y) THEN
                 YY=0
               ELSE
                  IF(NOT_CHECKED(PARTICLES+4))THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+4)
                    NOT_CHECKED(PARTICLES+4)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+4
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(DIMN.EQ.3) THEN
             IF(ZZ.eq.0)THEN !Bottom WAll
               IF(DES_PERIODIC_WALLS_Z) THEN
                 YY=SEARCH_GRIDS(3)
               ELSE
                  IF(NOT_CHECKED(PARTICLES+5))THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+5)
                    NOT_CHECKED(PARTICLES+5)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+5
                    N=N+1
                    NOT_WALL=.FALSE.
                  END IF
               END IF
             END IF

             IF(ZZ.eq.SEARCH_GRIDS(3)+1)THEN !Top Wall
               IF(DES_PERIODIC_WALLS_Z) THEN
                 ZZ=0
               ELSE
                  IF(NOT_CHECKED(PARTICLES+6))THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,PARTICLES+6)
                    NOT_CHECKED(PARTICLES+6)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+6
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF
             END IF

             IF(COHESION_DEBUG.gt.2)THEN
               PRINT *,'****XX=',XX
               PRINT *,'****YY=',YY
             END IF

             IF((PART_IN_GRID(XX,YY,ZZ,1).gt.0).AND.NOT_WALL)THEN
               DO KK=2,PART_IN_GRID(XX,YY,ZZ,1)+1               

                 J=PART_IN_GRID(XX,YY,ZZ,KK)

                 IF(J.gt.PARTICLES+DIMN*2)THEN
                    PRINT *, 'STOP 3'
                    STOP
                 END IF 

                 IF(I.gt.J)THEN
                    IF(NOT_CHECKED(J))THEN
                      NOT_CHECKED(J)=.FALSE.
                      CHECKED_PARTICLES_LIST(N)=J
                      N=N+1

                      LINK = CHECK_LINK(I,J)
                      IF(LINK.eq.1)THEN
                         CALL LINKED_INTERACTION_EVAL(I,J)
                      ELSE
                         CALL UNLINKED_INTERACTION_EVAL(I,J)
                      END IF
                    END IF
                  END IF

                END DO ! Loop over particles in grid
             END IF ! Is there a particle in grid and it is not a wall

             END DO ! Loop in z-direction
           END DO ! Loop in y-direction
        END DO ! Loop in x-direction        

!-------Check any linked particles outside neighbor range
        IF(COHESION_DEBUG.gt.1)THEN
           PRINT *,'****START CHECK OUTSIDE RANGE'
        END IF

        IF(LINKS(I,1).gt.0)THEN
          N=2
          J=LINKS(I,N)

          IF(J.gt.PARTICLES+DIMN*2)THEN
            PRINT *, 'STOP 2'
            STOP
          END IF 

          DO WHILE(J.gt.0)
            IF((I.gt.J).OR.(J.gt.PARTICLES))THEN
               IF(NOT_CHECKED(J))THEN
                  IF(J.gt.PARTICLES)THEN
                    CALL CHECK_SW_WALL_INTERACTION(I,J)
                  ELSE
                    CALL LINKED_INTERACTION_EVAL(I,J)
                  END IF
               END IF
            END IF
            N=N+1
            J=LINKS(I,N)
          END DO
        END IF
      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *,'****END CHECK OUTSIDE RANGE'
      END IF

!-----Reset checked particle list
      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *,'****START RESET CHECKED LIST'
      END IF
      N=1
      J=CHECKED_PARTICLES_LIST(N)

      IF(J.gt.PARTICLES+DIMN*2)THEN
          PRINT *, 'STOP 1'
          STOP
      END IF 

      DO WHILE(J.gt.0)
        NOT_CHECKED(CHECKED_PARTICLES_LIST(N))=.TRUE.
        CHECKED_PARTICLES_LIST(N)=-1
        N=N+1
        J=CHECKED_PARTICLES_LIST(N)
      END DO
      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *,'****END RESET CHECKED LIST'
      END IF

      END DO !End Loop over all particles


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END SQUARE WELL CALCS'
      END IF

      END SUBROUTINE CALC_SQUARE_WELL


