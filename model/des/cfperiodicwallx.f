!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_X()                                   C
!  Purpose: DES - Identify particles next to periodix x-walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CF_PERIODIC_WALL_X()
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(L,1).LE.(WX1 + 2D0*RADIUS_EQ)) THEN
            WWALL(1) = WWALL(1) + 1
            WI = WWALL(1) + 1
            IF(WI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits WW', WI
               STOP
            ELSE
               WWALL(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(L,1).GE.(EX2 - 2D0*RADIUS_EQ)) THEN
            EWALL(1) = EWALL(1) + 1 
            EI = EWALL(1) + 1
            IF(EI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits EW', EI
               STOP
            ELSE
               EWALL(EI) = L
            END IF
         END IF

      END DO

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_X


