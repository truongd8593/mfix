!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_RESTART(PARTS)                               C
!  Purpose: Writing DES data for restart                               C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer:Sreekanth Pannala                         Date: 31-Oct-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_RESTART(PARTS)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      INTEGER LL, PARTS, K 

!---------------------------------------------------------------------------

      OPEN (UNIT=901, FILE='des_restart_pos.out', STATUS='unknown')
      OPEN (UNIT=902, FILE='des_restart_vel.out', STATUS='unknown')
      OPEN (UNIT=903, FILE='des_restart_neighbor.out', STATUS='unknown')
      OPEN (UNIT=904, FILE='des_restart_force.out', STATUS='unknown')

      REWIND(901)
      REWIND(902)
      REWIND(903)
      REWIND(904)

      DO LL = LL, PARTICLES
         WRITE (901,*) (DES_POS_NEW(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         WRITE (902,*) (DES_VEL_NEW(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         WRITE (903,*) (NEIGHBOURS(LL,K),K=1,MAXNEIGHBORS)
      END DO

      DO LL = 1, PARTICLES
         WRITE (904,*) (FC(LL,K),K=1,DIMN)
      END DO

      END SUBROUTINE WRITE_DES_RESTART 


