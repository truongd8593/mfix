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
      USE param 
      USE param1 
      USE fldvar
      USE geometry
      USE indices
      USE compar        !//d
      Use usr
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
      DOUBLE PRECISION num, den
      INTEGER IJK, IJK2, IJK1, I
!
!  Include files defining statement functions here
!
      INCLUDE 'function.inc'
!
!  Insert user-defined code here
!
      num = zero
      den = zero
      do i = imin1, imax1
        IJK1 = FUNIJK(i, JMIN2, 1)
        IJK2 = FUNIJK(i, JMAX1, 1)
        num = num + (X_g(IJK2, 1) * ROP_g(IJK2) * V_g(IJK2))
	den = den + (X_g(IJK1, 1) * ROP_g(IJK1) * V_g(IJK1))
      end do

      write(*,'(//A,G12.5//)')' Conversion = ', (ONE -  (num/den))
      RETURN  
      END SUBROUTINE USR3 
