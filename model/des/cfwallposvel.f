!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLPOSVEL(L, STIME, I)
!  Purpose:  DES -Calculate the position and velocity of wall particle C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLPOSVEL(L, STIME, I)

      Use discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE

      INTEGER L, I 
      DOUBLE PRECISION STIME, A, OMEGA_W, F, COSOMEGAT, SINOMEGAT
      DOUBLE PRECISION DES_R
      
!     
!---------------------------------------------------------------------
!     Assigning wall position and velocity
!---------------------------------------------------------------------
!     

      A = 0
      OMEGA_W = 0.0
      IF(DES_F.NE.0.0) THEN
         OMEGA_W = 2*22*DES_F/7
         IF(UNITS == "CGS") THEN
            A = DES_GAMMA*980/(OMEGA_W*OMEGA_W)
         ELSE
            A = DES_GAMMA*9.81/(OMEGA_W*OMEGA_W)
         END IF
      END IF

      SINOMEGAT = SIN(OMEGA_W*STIME)
      COSOMEGAT = COS(OMEGA_W*STIME)

      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = 0.0
      END IF

      DES_WALL_VEL(1,I) = 0.0
      DES_WALL_VEL(2,I) = 0.0
      DES_WALL_VEL(3,I) = 0.0

      IF(I.EQ.1) THEN
         DES_WALL_POS(1,I) = WX1 - DES_R	
         DES_WALL_POS(2,I) = DES_POS_NEW(2,L)
         DES_WALL_POS(3,I) = DES_POS_NEW(3,L)
         WALL_NORMAL(1,1) = -1.0
         WALL_NORMAL(2,1) = 0.0
         WALL_NORMAL(3,1) = 0.0

      ELSE IF(I.EQ.2) THEN
         DES_WALL_POS(1,I) = EX2 + DES_R
         DES_WALL_POS(2,I) = DES_POS_NEW(2,L)
         DES_WALL_POS(3,I) = DES_POS_NEW(3,L)
         WALL_NORMAL(1,2) = 1.0
         WALL_NORMAL(2,2) = 0.0
         WALL_NORMAL(3,2) = 0.0

      ELSE IF(I.EQ.3) THEN
         DES_WALL_POS(1,I) = DES_POS_NEW(1,L)
         DES_WALL_POS(2,I) = BY1 - DES_R + (A*SINOMEGAT)
         DES_WALL_POS(3,I) = DES_POS_NEW(3,L)
         DES_WALL_VEL(2,I) =  A*OMEGA_W*COSOMEGAT
         WALL_NORMAL(1,3) = 0.0
         WALL_NORMAL(2,3) = -1.0
         WALL_NORMAL(3,3) = 0.0

      ELSE IF(I.EQ.4) THEN
         DES_WALL_POS(1,I) = DES_POS_NEW(1,L)
         DES_WALL_POS(2,I) = TY2 + DES_R
         DES_WALL_POS(3,I) = DES_POS_NEW(3,L)
         WALL_NORMAL(1,4) = 0.0
         WALL_NORMAL(2,4) = 1.0
         WALL_NORMAL(3,4) = 0.0

      ELSE IF(I.EQ.5) THEN
         DES_WALL_POS(1,I) = DES_POS_NEW(1,L)
         DES_WALL_POS(2,I) = DES_POS_NEW(2,L)
         DES_WALL_POS(3,I) = SZ1 - DES_R
         WALL_NORMAL(1,5) = 0.0
         WALL_NORMAL(2,5) = 0.0
         WALL_NORMAL(3,5) = -1.0

      ELSE IF(I.EQ.6) THEN
         DES_WALL_POS(1,I) = DES_POS_NEW(1,L)
         DES_WALL_POS(2,I) = DES_POS_NEW(2,L)
         DES_WALL_POS(3,I) = NZ2 + DES_R
         WALL_NORMAL(1,6) = 0.0
         WALL_NORMAL(2,6) = 0.0
         WALL_NORMAL(3,6) = 1.0
      END IF
      
!     IF(CALLED.GT.1095) THEN
!     IF((L.EQ.1).AND.(I.EQ.3)) THEN
!     PRINT *, DES_F, DES_GAMMA, OMEGA_W, A, STIME, SINOMEGAT, COSOMEGAT
!     PRINT *, L, I, (DES_POS_NEW(K,L),K=1,DIMN), (DES_WALL_POS(K,I),K=1,DIMN)
!     PRINT *, (DES_VEL_NEW(K,L),K=1,DIMN), (DES_WALL_VEL(K,I), K=1,DIMN)
!     STOP
!     END IF
!     END IF

      RETURN
      END SUBROUTINE CFWALLPOSVEL



