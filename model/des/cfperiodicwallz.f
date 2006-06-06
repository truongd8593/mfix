!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CF_PERIODIC_WALL_Z()                                   C
!  Purpose: DES - Identify particles next to periodix z-walls          C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  REVISED Aug 24 2005                                                 C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CF_PERIODIC_WALL_Z()
      
      USE discretelement
      IMPLICIT NONE

      INTEGER K, L, WI, EI
!     
!-------------------------------------------------------------------
!     

      WI = 0
      EI = 0	
      
      DO L = 1, PARTICLES

         IF(DES_POS_NEW(L,3).LE.(SZ1 + 2D0*RADIUS_EQ)) THEN
            SWALL(1) = SWALL(1) + 1
            WI = SWALL(1) + 1
            IF(WI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits SWALL', WI
               STOP
            ELSE
               SWALL(WI) = L
            END IF

         ELSE IF(DES_POS_NEW(L,3).GE.(NZ2 - 2D0*RADIUS_EQ)) THEN
            NWALL(1) = NWALL(1) + 1 
            EI = NWALL(1) + 1
            IF(EI.GT.PBP-1) THEN
               PRINT *,'Exceeding array limits NWALL', EI
               STOP
            ELSE
               NWALL(EI) = L
            END IF
         END IF

      END DO

      RETURN
      END SUBROUTINE CF_PERIODIC_WALL_Z


