!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is
!           user-definable.  The user may insert code in this routine
!           or call appropriate user defined subroutines.
!           This routine is not called from an IJK loop, hence
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR3
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      Use usr
      USE fldvar
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!
call abort
      OPEN(6,FILE='POST_Aflame.dat')
      write(6,'(A,G10.3)') 'Adiabatic Flame Temperature = ', T_g(5)
      write(6,'(A,G10.3)') 'P_g = ', P_g(5)
      write(6,'(A,G10.3,A,G10.3)') 'CH4=', X_g(5, 1), ' O2=', X_g(5, 2)
      write(6,'(A,G10.3,A,G10.3)') 'CO2=', X_g(5, 3), ' H2O=', X_g(5, 4)
      write(6,'(A,G10.3)') 'N2 =', X_g(5, 5)
      CLOSE(6)
      RETURN
      END SUBROUTINE USR3
