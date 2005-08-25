!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_Z(WW, EW)                             C
!  Purpose: DES - Identify particles next to periodix z-walls          C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_Z(WW, EW)
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
      INTEGER WW(1000), EW(1000)
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      WW(1) = 0
      EW(1) = 0
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(3,L).LE.(SZ1 + 2*RADIUS_EQ)) THEN
            WW(1) = WW(1) + 1
            WI = WW(1) + 1
            IF(WI.GT.999) THEN
               PRINT *,'WI', WI
               STOP
            ELSE
               WW(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(3,L).GE.(NZ2 - 2*RADIUS_EQ)) THEN
            EW(1) = EW(1) + 1 
            EI = EW(1) + 1
            IF(EI.GT.999) THEN
               PRINT *,'EI', EI
               STOP
            ELSE
               EW(EI) = L
            END IF
         END IF

      END DO

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_Z


