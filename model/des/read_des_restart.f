!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_DES_RESTART                                       C
!>  Purpose: Reading DES data for restart                               
!                                                                      C
!                                                                      C
!  Author: Sreekanth Pannala                          Date: 09-Nov-06  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DES_RESTART

      USE param1      
      USE run
      USE discretelement
      IMPLICIT NONE

      INTEGER LL, K

!---------------------------------------------------------------------------

      OPEN (UNIT=901, FILE=TRIM(RUN_NAME)//'_DES.RES', FORM='Unformatted', STATUS='unknown')

      REWIND (901)

      READ (901) PARTICLES
      READ (901) IFI
      READ (901) DTSOLID
      READ (901) DES_POS_OLD
      READ (901) DES_VEL_OLD
      READ (901) DES_RADIUS
      READ (901) RO_Sol
      READ (901) NEIGHBOURS
      READ (901) FC
      READ (901) OMEGA_OLD
      READ (901) TOW
      READ (901) DES_U_s
      READ (901) DES_V_s
      READ (901) DES_W_s

      END SUBROUTINE READ_DES_RESTART 


