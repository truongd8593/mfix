!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES(PARTS)                                 C
!  Purpose: DES - allocating DES arrays                                C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MAKE_ARRAYS_DES(PARTS)
      
      USE discretelement
      IMPLICIT NONE
      
      INTEGER LN, K, PARTS, FAC
      DOUBLE PRECISION  DUMMY


      IF(ZONES.EQ.1) THEN
         OPEN(UNIT=10, FILE='pos_input.dat', STATUS='OLD')

         DO LN = 1, PARTICLES

            READ (10, *) DUMMY, (DES_POS_OLD(K,LN), K=1,NDIM), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,NDIM)

            
            IF(LN.LT.PARTICLES/5) THEN
               FAC = 60
            ELSE IF(LN.LT.2*PARTICLES/5) THEN
               FAC = 50
            ELSE IF(LN.LT.3*PARTICLES/5) THEN
               FAC = 60
            ELSE IF(LN.LT.4*PARTICLES/5) THEN
               FAC = 50
            ELSE IF(LN.LT.PARTICLES) THEN
               FAC = 60
            END IF

            PR(LN) = DES_RADIUS(LN) - DES_RADIUS(LN)/FAC
            RO_Sol(LN) = ROs
            DES_POS_OLD(1,LN) = DES_POS_OLD(1,LN) + 0.45
            DES_POS_OLD(3,LN) = RADIUS_EQ

            OMEGA_OLD(1,LN) = 0.0
            OMEGA_OLD(2,LN) = 0.0
            OMEGA_OLD(3,LN) = 0.0

         END DO

      ELSE IF (ZONES.GT.1) THEN

         OPEN(UNIT=10, FILE='zone1.dat', STATUS='OLD')
         DO LN = 1, ZN1 
            READ (10, *) DUMMY, (DES_POS_OLD(K,LN), K=1,NDIM), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,NDIM)
            DES_RADIUS(LN) = DES_RADIUS(LN) - DES_RADIUS(LN)/100
            DO K = 1, 3
               DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 1.0
               IF(K.EQ.2) THEN
                  DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 0.001
               END IF
            END DO
            RO_Sol(LN) = ROs
            OMEGA_OLD(1,LN) = 0.0
            OMEGA_OLD(2,LN) = 0.0
            OMEGA_OLD(3,LN) = 0.0
         END DO

         IF (ZONES.GE.2) THEN
            OPEN(UNIT=10, FILE='zone2.dat', STATUS='OLD')

            DO LN = ZN1 + 1, ZN1 + ZN2
               READ (10, *) DUMMY, (DES_POS_OLD(K,LN), K=1,NDIM), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,NDIM)
               DES_RADIUS(LN) = DES_RADIUS(LN) - DES_RADIUS(LN)/100
               DO K = 1, 3
                  DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 1
                  IF(K.EQ.2) THEN
                     DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 0.015
                  END IF
               END DO
               RO_Sol(LN) = ROs
               OMEGA_OLD(1,LN) = 0.0
               OMEGA_OLD(2,LN) = 0.0
               OMEGA_OLD(3,LN) = 0.0
            END DO
         END IF

         IF(ZONES.GE.3) THEN
            OPEN(UNIT=10, FILE='zone3.dat', STATUS='OLD')

            DO LN = ZN1 + ZN2 + 1, ZN1 + ZN2 + ZN3
               READ (10, *) DUMMY, (DES_POS_OLD(K,LN), K=1,NDIM), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,NDIM)
               DES_RADIUS(LN) = DES_RADIUS(LN) - DES_RADIUS(LN)/100
               DO K = 1, 3
                  DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 1
                  IF(K.EQ.2) THEN
                     DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 0.04
                  END IF
               END DO
               RO_Sol(LN) = ROs
               OMEGA_OLD(1,LN) = 0.0
               OMEGA_OLD(2,LN) = 0.0
               OMEGA_OLD(3,LN) = 0.0
            END DO
         END IF

         IF(ZONES.GE.4) THEN
            OPEN(UNIT=10, FILE='zone4.dat', STATUS='OLD')

            DO LN = ZN1 + ZN2 + ZN3 + 1, ZN1 + ZN2 + ZN3 + ZN4
               READ (10, *) DUMMY, (DES_POS_OLD(K,LN), K=1,NDIM), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,NDIM)
               DES_RADIUS(LN) = DES_RADIUS(LN) - DES_RADIUS(LN)/100
               DO K = 1, 3
                  DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 1.0
                  IF(K.EQ.2) THEN
                     DES_POS_OLD(K,LN) = DES_POS_OLD(K,LN) + 0.005
                  END IF
               END DO
               OMEGA_OLD(1,LN) = 0.0
               OMEGA_OLD(2,LN) = 0.0
               OMEGA_OLD(3,LN) = 0.0
            END DO
         END IF
      END IF

      DO LN = 1, PARTICLES
         DO K = 1,DIMN
            DES_POS_NEW(K,LN) = DES_POS_OLD(K,LN)
            DES_VEL_NEW(K,LN) = DES_VEL_OLD(K,LN)
            OMEGA_NEW(K,LN) = OMEGA_OLD(K,LN)
         END DO
      END DO

!     DO LN =1, PARTICLES
!     PRINT *,(DES_POS_OLD(K,LN), K=1,DIMN), DES_RADIUS(LN), (DES_VEL_OLD(K,LN),K=1,DIMN)
!     END DO

      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES 


