!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFWALLCONTACT(WALL, L, STIME, WALLCONTACTI)            C
!  Purpose: DES - Checking for contact with walls                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFWALLCONTACT(WALL, L, STIME, WALLCONTACTI)

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
      DOUBLE PRECISION STIME, A, OMEGA, ASINOMEGAT 
      
!     
!---------------------------------------------------------------------
!     Checking if a particle is in contact with any of the walls
!---------------------------------------------------------------------
!     

      A = 0.0
      OMEGA = 0.0
      IF(DES_F.NE.0.0) THEN
         OMEGA = 2*22*DES_F/7
         IF(UNITS == "CGS") THEN
            A = DES_GAMMA*980/(OMEGA*OMEGA)
         ELSE 
            A = DES_GAMMA*9.81/(OMEGA*OMEGA) 
         END IF
      END IF
      ASINOMEGAT = A*SIN(OMEGA*STIME)

      WALLCONTACTI = 0

      IF(WALL.EQ.1) THEN
         IF((DES_POS_NEW(1,L)-WX1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.2) THEN
         IF((EX2-DES_POS_NEW(1,L)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.3) THEN
         IF((DES_POS_NEW(2,L)-(BY1+ASINOMEGAT)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.4) THEN
         IF((TY2-DES_POS_NEW(2,L)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.5) THEN
         IF((DES_POS_NEW(3,L)-SZ1).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF

      ELSE IF(WALL.EQ.6) THEN
         IF((NZ2-DES_POS_NEW(3,L)).LE.DES_RADIUS(L)) THEN
            WALLCONTACTI = 1
         END IF
      END IF

      RETURN
      END SUBROUTINE CFWALLCONTACT


