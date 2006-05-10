!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_Z(SW, NW)                             C
!  Purpose: DES - Identify particles next to periodix z-walls          C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_Z(SW, NW)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
      INTEGER SW(1000), NW(1000)
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      SW(1) = 0
      NW(1) = 0
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(L,3).LE.(SZ1 + 2*RADIUS_EQ)) THEN
            SW(1) = SW(1) + 1
            WI = SW(1) + 1
            IF(WI.GT.999) THEN
               PRINT *,'Exceeding array limits SW', WI
               STOP
            ELSE
               SW(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(L,3).GE.(NZ2 - 2*RADIUS_EQ)) THEN
            NW(1) = NW(1) + 1 
            EI = NW(1) + 1
            IF(EI.GT.999) THEN
               PRINT *,'Exceeding array limits NW', EI
               STOP
            ELSE
               NW(EI) = L
            END IF
         END IF

      END DO

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_Z


