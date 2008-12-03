!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLCONTACT(WALL, L, WALLCONTACTI)            C
!>  Purpose: DES - Checking for contact with walls
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLCONTACT(WALL, L, WALLCONTACTI)

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
      
      INTEGER L, I, K, WALL, WALLCONTACTI	
      DOUBLE PRECISION A, OMEGA, OOMEGA2, ASINOMEGAT 
      
!     
!---------------------------------------------------------------------
!     Checking if a particle is in contact with any of the walls
!---------------------------------------------------------------------
!     

      A = ZERO
      OMEGA = ZERO
      ASINOMEGAT = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA = 2.0D0*PI*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         A = DES_GAMMA*GRAV(2)*OOMEGA2
         ASINOMEGAT = A*SIN(OMEGA*S_TIME)
      END IF

      WALLCONTACTI = 0

! west wall (X)
      IF(WALL.EQ.1.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((DES_POS_NEW(L,1)-WX1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

! east wall (X)
      ELSE IF(WALL.EQ.2.AND.(.NOT.DES_PERIODIC_WALLS_X)) THEN
         IF((EX2-DES_POS_NEW(L,1)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

! bottom wall (Y)
      ELSE IF(WALL.EQ.3.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((DES_POS_NEW(L,2)-(BY1+ASINOMEGAT)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

! top wall (Y)
      ELSE IF(WALL.EQ.4.AND.(.NOT.DES_PERIODIC_WALLS_Y)) THEN
         IF((TY2-DES_POS_NEW(L,2)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

! south wall (Z)
      ELSE IF(WALL.EQ.5.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((DES_POS_NEW(L,3)-SZ1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

! north wall (Z)
      ELSE IF(WALL.EQ.6.AND.(.NOT.DES_PERIODIC_WALLS_Z)) THEN
         IF((NZ2-DES_POS_NEW(L,3)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF
      END IF

      RETURN
      END SUBROUTINE CFWALLCONTACT


