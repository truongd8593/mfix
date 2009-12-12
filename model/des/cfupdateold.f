!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFUPDATEOLD                                            C
!
!  Purpose: DES - Update arrays                                        
!
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFUPDATEOLD
      
      USE discretelement
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters (no. particles)
      INTEGER LL
! Accounted for particles
      INTEGER PC

!-----------------------------------------------

      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         DES_VEL_OOLD(LL,:) = DES_VEL_OLD(LL,:)
         DES_POS_OLD(LL,:)  = DES_POS_NEW(LL,:)
         DES_VEL_OLD(LL,:)  = DES_VEL_NEW(LL,:)
         OMEGA_OLD(LL,:)    = OMEGA_NEW(LL,:)

         PC = PC + 1
      ENDDO
 

      RETURN

      END SUBROUTINE CFUPDATEOLD



