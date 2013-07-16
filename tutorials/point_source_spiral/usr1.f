!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR1                                                   C
!  Purpose: This routine is called from the time loop and is           C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting or checking errors in quantities   C
!           that vary with time.  This routine is not called from an   C
!           IJK loop, hence all indices are undefined.                 C               C
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
      SUBROUTINE USR1 

      use ps
      use constant
      use geometry
      use physprop
      use run
      use usr

      IMPLICIT NONE

      DOUBLE PRECISION :: lRad,   lTime
      INTEGER :: M

      lTime = time
      do while (lTime > 3.0d0)
         lTime = lTime - 3.0d0
      enddo

      lRad = 2.0d0*Pi*(lTime/3.0d0)

! Calculate the normalized velocity vectors for all phases.
      PS_V_g(1) = 0.12d0/0.15d0
      PS_U_g(1) = cos(lRad)/0.15d0
      PS_W_g(1) = sin(lRad)/0.15d0

      do M=1,MMAX
         PS_V_s(1,M) = 0.12d0/0.15d0
         PS_U_s(1,M) = 0.09d0*cos(lRad)/0.15d0
         PS_W_s(1,M) = 0.09d0*sin(lRad)/0.15d0
      enddo

      RETURN  
      END SUBROUTINE USR1
