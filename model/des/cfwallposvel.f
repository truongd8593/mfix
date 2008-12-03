!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLPOSVEL(L, I)
!>  Purpose:  DES -Calculate the position and velocity of wall particle 
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLPOSVEL(L, I)

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
      DOUBLE PRECISION A, OMEGA_W, F, COSOMEGAT, SINOMEGAT
      DOUBLE PRECISION DES_R, OOMEGAW2
      
!     
!---------------------------------------------------------------------
!     Assigning wall position and velocity
!---------------------------------------------------------------------
!     

      A = ZERO
      OMEGA_W = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA_W = 2.0d0*Pi*DES_F
         OOMEGAW2 = ONE/(OMEGA_W**2)
         A = DES_GAMMA*GRAV(2)*OOMEGAW2
         SINOMEGAT = SIN(OMEGA_W*S_TIME)
         COSOMEGAT = COS(OMEGA_W*S_TIME)
      END IF

      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = ZERO
      END IF

      DES_WALL_VEL(I,1) = ZERO
      DES_WALL_VEL(I,2) = ZERO
      IF(DIMN.EQ.3) DES_WALL_VEL(I,3) = ZERO

!     west (X)
      IF(I.EQ.1) THEN
         DES_WALL_POS(I,1) = WX1 - DES_R	
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         WALL_NORMAL(1,1) = -ONE
         WALL_NORMAL(1,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(1,3) = ZERO

!     east (X)
      ELSE IF(I.EQ.2) THEN
         DES_WALL_POS(I,1) = EX2 + DES_R
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         WALL_NORMAL(2,1) = ONE
         WALL_NORMAL(2,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(2,3) = ZERO

!     bottom (Y)
      ELSE IF(I.EQ.3) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = BY1 - DES_R + (A*SINOMEGAT)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         DES_WALL_VEL(I,1) = lid_vel
         DES_WALL_VEL(I,2) =  A*OMEGA_W*COSOMEGAT
         WALL_NORMAL(3,1) = ZERO
         WALL_NORMAL(3,2) = -ONE
         IF(DIMN.EQ.3) WALL_NORMAL(3,3) = ZERO

!     top (Y)
      ELSE IF(I.EQ.4) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = TY2 + DES_R
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = DES_POS_NEW(L,3)
         DES_WALL_VEL(I,1) = -lid_vel
         WALL_NORMAL(4,1) = ZERO
         WALL_NORMAL(4,2) = ONE
         IF(DIMN.EQ.3) WALL_NORMAL(4,3) = ZERO

!     south (Z)
      ELSE IF(I.EQ.5) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = SZ1 - DES_R
         WALL_NORMAL(5,1) = ZERO
         WALL_NORMAL(5,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(5,3) = -ONE

!     north (Z)
      ELSE IF(I.EQ.6) THEN
         DES_WALL_POS(I,1) = DES_POS_NEW(L,1)
         DES_WALL_POS(I,2) = DES_POS_NEW(L,2)
         IF(DIMN.EQ.3) DES_WALL_POS(I,3) = NZ2 + DES_R
         WALL_NORMAL(6,1) = ZERO
         WALL_NORMAL(6,2) = ZERO
         IF(DIMN.EQ.3) WALL_NORMAL(6,3) = ONE
      END IF
      
      RETURN
      END SUBROUTINE CFWALLPOSVEL



