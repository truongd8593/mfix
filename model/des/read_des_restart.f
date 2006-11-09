!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_DES_RESTART                                       C
!  Purpose: Reading DES data for restart                               C
!                                                                      C
!                                                                      C
!  Author: Sreekanth Pannala                          Date: 09-Nov-06  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DES_RESTART

      USE param1      
      USE discretelement
      IMPLICIT NONE

      INTEGER LL, K 

!---------------------------------------------------------------------------

      OPEN (UNIT=901, FILE='des_restart_pos.out', STATUS='unknown')
      OPEN (UNIT=902, FILE='des_restart_vel.out', STATUS='unknown')
      OPEN (UNIT=903, FILE='des_restart_neighbor.out', STATUS='unknown')
      OPEN (UNIT=904, FILE='des_restart_force.out', STATUS='unknown')
      OPEN (UNIT=905, FILE='des_restart_omega.out', STATUS='unknown')
      OPEN (UNIT=904, FILE='des_restart_tow.out', STATUS='unknown')

      REWIND(901)
      REWIND(902)
      REWIND(903)
      REWIND(904)
      REWIND(905)
      REWIND(906)

      DO LL = LL, PARTICLES
         READ (901,*) (DES_POS_OLD(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         READ (902,*) (DES_VEL_OLD(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         READ (903,*) (NEIGHBOURS(LL,K),K=1,MAXNEIGHBORS)
      END DO

      DO LL = 1, PARTICLES
         READ (904,*) (FC(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         READ (905,*) (OMEGA_OLD(LL,K),K=1,DIMN)
      END DO

      DO LL = 1, PARTICLES
         READ (906,*) (TOW(LL,K),K=1,DIMN)
      END DO

      END SUBROUTINE READ_DES_RESTART 


