!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_VAN_DER_WAALS                                     C
!  Purpose: This subroutine calculates the attractive force between    C
!      particles consistent with the Hamaker van der Waals             C
!      model (Seville, Processing of Particulate Solids, p105)         C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE CALC_VAN_DER_WAALS


!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I, J, N, K, II
      INTEGER X_INDEX, Y_INDEX, Z_INDEX
      INTEGER XX, YY, ZZ, KK
      INTEGER Z_START, Z_END
      DOUBLE PRECISION RADIUS, R(NDIM)
      DOUBLE PRECISION FORCE, DIST
      LOGICAL NOT_CHECKED(PARTICLES+6)
      INTEGER CHECKED_PARTICLES_LIST(50)
      LOGICAL NOT_WALL

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START VAN DER WAALS CALCS'
      END IF


!!The checked particles list is necessary to guard against redundant
!!  calculations with the wall since a particle will neighbor three
!!  wall cells if it is near a wall.

      DO J=1,PARTICLES+2*DIMN
         NOT_CHECKED(J)=.TRUE.
      END DO

      DO J=1,50
         CHECKED_PARTICLES_LIST(J)=-1
      END DO


!-----START LOOP OVER ALL PARTICLES
      DO I = 1, PARTICLES
        N=1


!-------Check all neighboring grids for particles

        IF(DIMN.eq.2)THEN
           Z_START=1
           Z_END=1
        ELSE
           Z_START=PART_GRID(3,I)-1
           Z_END=PART_GRID(3,I)+1
        END IF

        DO X_INDEX=PART_GRID(1,I)-1,PART_GRID(1,I)+1
           DO Y_INDEX=PART_GRID(2,I)-1,PART_GRID(2,I)+1
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
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+1)
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
                  IF(NOT_CHECKED(PARTICLES+2))THEN !!guard against multiple intrxns 
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+2)
                    NOT_CHECKED(PARTICLES+2)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+2
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(XX.eq.SEARCH_GRIDS(1)+1)THEN !East Wall
               IF(DES_PERIODIC_WALLS_X) THEN
                 XX=0
               ELSE
                  IF(NOT_CHECKED(PARTICLES+3))THEN !!guard against multiple intrxns 
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+3)
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
                  IF(NOT_CHECKED(PARTICLES+4))THEN !!guard against multiple intrxns 
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+4)
                    NOT_CHECKED(PARTICLES+4)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+4
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF
           
             IF(ZZ.eq.0)THEN !Bottom WAll
               IF(DES_PERIODIC_WALLS_Z) THEN
                 ZZ=SEARCH_GRIDS(3)
               ELSE
                  IF(NOT_CHECKED(PARTICLES+5))THEN !!guard against multiple intrxns 
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+5)
                    NOT_CHECKED(PARTICLES+5)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+5
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(ZZ.eq.SEARCH_GRIDS(3)+1)THEN !Top Wall
               IF(DES_PERIODIC_WALLS_Z) THEN
                 ZZ=0
               ELSE
                  IF(NOT_CHECKED(PARTICLES+6))THEN !!guard against multiple intrxns 
                    CALL CHECK_VDW_WALL_INTERACTION(I,PARTICLES+6)
                    NOT_CHECKED(PARTICLES+6)=.FALSE.
                    CHECKED_PARTICLES_LIST(N)=PARTICLES+6
                    N=N+1
                  END IF
               END IF
               NOT_WALL=.FALSE.
             END IF

             IF(COHESION_DEBUG.gt.2)THEN
               PRINT *,'****XX=',XX
               PRINT *,'****YY=',YY
             END IF


             IF((PART_IN_GRID(XX,YY,ZZ,1).gt.0).AND.NOT_WALL)THEN
               DO KK=2,PART_IN_GRID(XX,YY,ZZ,1)+1
               J=PART_IN_GRID(XX,YY,ZZ,KK)

               IF(I.gt.J)THEN
                  IF(NOT_CHECKED(J))THEN

      RADIUS=0.0
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
        RADIUS=RADIUS+R(K)*R(K)
      END DO

                     NOT_CHECKED(J)=.FALSE.
                     CHECKED_PARTICLES_LIST(N)=J
                     N=N+1
                     RADIUS=SQRT(RADIUS)
                     DIST=RADIUS-DES_RADIUS(I)-DES_RADIUS(J)
                     IF(DIST.lt.VDW_OUTER_CUTOFF)THEN
                        IF(DIST.gt.VDW_INNER_CUTOFF)THEN
                           FORCE=HAMAKER_CONSTANT*DES_RADIUS(I)*DES_RADIUS(J)/&
                                (6*DIST*DIST*(DES_RADIUS(I)+DES_RADIUS(J)))
                     ELSE

                           FORCE=4*3.14*SURFACE_ENERGY*DES_RADIUS(I)*DES_RADIUS(J)/&
                                (DES_RADIUS(I)+DES_RADIUS(J))
                     END IF                     
                        DO II=1,DIMN
                           FC(II,I)=FC(II,I)+R(II)/RADIUS*FORCE
                           FC(II,J)=FC(II,J)-R(II)/RADIUS*FORCE
                        END DO

                     END IF   !! is particle within cutoff ?
                  END IF   !!  has particle already been checked?
               END IF    !!guard against double counting

               END DO    !!Loop over particles in grid
             END IF    !!is there particle in grid ?

             END DO     !! z-loop     
           END DO     !! y-loop
         END DO    !! x-loop 

!-----Reset checked particle list
      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *,'****START RESET CHECKED LIST'
      END IF
      N=1
      J=CHECKED_PARTICLES_LIST(N)
      DO WHILE(J.gt.0)
        NOT_CHECKED(CHECKED_PARTICLES_LIST(N))=.TRUE.
        CHECKED_PARTICLES_LIST(N)=-1
        N=N+1
        J=CHECKED_PARTICLES_LIST(N)
      END DO
      IF(COHESION_DEBUG.gt.1)THEN
         PRINT *,'****END RESET CHECKED LIST'
      END IF
 
      END DO  !! particle-loop


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END VAN DER WAALS CALCS'
      END IF

      END SUBROUTINE CALC_VAN_DER_WAALS

