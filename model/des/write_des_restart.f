!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_RESTART                                      C
!  Purpose: Writing DES data for restart
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer:Sreekanth Pannala                         Date: 31-Oct-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_RESTART

      USE param1      
      USE discretelement
      USE run

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------


      OPEN (UNIT=901,FILE=TRIM(RUN_NAME)//'_DES.RES',FORM='Unformatted',STATUS='unknown')

      REWIND (901)

      WRITE (901) PARTICLES
      WRITE (901) IFI
      WRITE (901) DTSOLID

      WRITE (901) DES_POS_NEW
      WRITE (901) DES_VEL_NEW
      WRITE (901) OMEGA_NEW

      WRITE (901) DES_RADIUS
      WRITE (901) RO_Sol

      WRITE (901) NEIGHBOURS

! Needed for particle contact history in tangential direction
      WRITE (901) PFT
      WRITE (901) PN
      WRITE (901) PV


! J. Musser DES boundary condition data
      WRITE (901) PIS
      WRITE (901) PEA

      END SUBROUTINE WRITE_DES_RESTART 

