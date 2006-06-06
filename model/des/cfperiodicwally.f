!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_Y()                                   C
!  Purpose: DES - Identify particles next to periodix y-walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_Y()
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(L,2).LE.(BY1 + 2D0*RADIUS_EQ)) THEN
            BWALL(1) = BWALL(1) + 1
            WI = BWALL(1) + 1
            IF(WI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits BWALL', WI
               STOP
            ELSE
               BWALL(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(L,2).GE.(TY2 - 2D0*RADIUS_EQ)) THEN
            TWALL(1) = TWALL(1) + 1 
            EI = TWALL(1) + 1
            IF(EI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits TWALL', EI
               STOP
            ELSE
               TWALL(EI) = L
            END IF
         END IF

      END DO
      
      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_Y


