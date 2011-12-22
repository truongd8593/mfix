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

!!      PC = 0
!$omp parallel do if(max_pip .ge. 10000) default(shared)         &
!$omp private(ll)                    &
!$omp schedule (guided,50)                      	  
      DO LL = 1, MAX_PIP
!!         if (pea(ll,1)) pc = pc + 1 
         IF(.NOT.PEA(LL,1) .or. pea(ll,4)) CYCLE

         DES_VEL_OOLD(LL,:) = DES_VEL_OLD(LL,:)
         DES_POS_OLD(LL,:)  = DES_POS_NEW(LL,:)
         DES_VEL_OLD(LL,:)  = DES_VEL_NEW(LL,:)
         OMEGA_OLD(LL,:)    = OMEGA_NEW(LL,:)

!!         IF(PC .eq. PIP) EXIT
      ENDDO
!$omp end parallel do 

      RETURN

      END SUBROUTINE CFUPDATEOLD



