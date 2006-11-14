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
      USE run
      IMPLICIT NONE

      INTEGER LL, PARTS, K

!---------------------------------------------------------------------------

      OPEN (UNIT=901, FILE=TRIM(RUN_NAME)//'_DES.RES', FORM='Unformatted', STATUS='unknown')

      REWIND (901)

      WRITE (901) PARTICLES
      WRITE (901) IFI
      WRITE (901) DTSOLID
      WRITE (901) DES_POS_NEW
      WRITE (901) DES_VEL_NEW
      WRITE (901) DES_RADIUS
      WRITE (901) RO_Sol
      WRITE (901) NEIGHBOURS
      WRITE (901) FC
      WRITE (901) OMEGA_NEW
      WRITE (901) TOW

      END SUBROUTINE WRITE_DES_RESTART 


