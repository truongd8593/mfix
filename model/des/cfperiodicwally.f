!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_Y(BW, TW)                             C
!  Purpose: DES - Identify particles next to periodix y-walls          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_Y(BW, TW)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
      INTEGER BW(1000), TW(1000)
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      BW(1) = 0
      TW(1) = 0
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(L,2).LE.(BY1 + 2*RADIUS_EQ)) THEN
            BW(1) = BW(1) + 1
            WI = BW(1) + 1
            IF(WI.GT.999) THEN
               PRINT *,'Exceeding array limits BW', WI
               STOP
            ELSE
               BW(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(L,2).GE.(TY2 - 2*RADIUS_EQ)) THEN
            TW(1) = TW(1) + 1 
            EI = TW(1) + 1
            IF(EI.GT.999) THEN
               PRINT *,'Exceeding array limits TW', EI
               STOP
            ELSE
               TW(EI) = L
            END IF
         END IF

      END DO
      
      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_Y


