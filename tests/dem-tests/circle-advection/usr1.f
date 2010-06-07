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
!
      SUBROUTINE USR1 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!
!  Include modules
!
      Use usr
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE constant
      USE toleranc
      USE compar
      USE run
      USE turb
      USE sendrecv  
      USE discretelement

      IMPLICIT NONE
!
!  Define local variables here
!

!                      Indices
      INTEGER          I, J, K, IJK
      DOUBLE PRECISION XX, YY, ZZ, XM, YM, ZM, LN

      DOUBLE PRECISION ERROR(PARTICLES)
!
!  Include files defining statement functions here
!
      INCLUDE 'function.inc'
!
!     Insert user-defined code here
!     
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         XX = XE(I)
         YY = YN(J)
         XM = XE(I) - 0.5d0*DX(I)
         YM = YN(J) - 0.5d0*DY(J)

         u_g(ijk) = (sin(PI*XX))**2*sin(2.0d0*PI*YM)*cos(PI*time/T_per)
         v_g(ijk) = -sin(2.0d0*PI*XM)*(sin(PI*YM))**2*cos(PI*time/T_per)

!        write(101,*) ijk, i, j, k, time, u_g(ijk), v_g(ijk), w_g(ijk)
!        write(102,*) ijk, i, j, k, time, pi, xx, yy, zz, zm, ym, zm, cos(pi*time/T_per)
!        write(103,*) ijk, i, j, k, time, pi, sin(pi*xx), sin(2*pi*ym), sin(2*pi*zm)
 
!        write(104,*) i, j, k, xx, yy, zz, xm, ym, zm

      END DO

      IF(MOD(time,t_per).lt.1.0e-5) then
      write(201,*) 'Cycle = ', time/T_per 
      DO LN = 1, PARTICLES
         error(ln) = sqrt ((DES_POS_NEW(LN,1)-x_store(ln,1))**2 &
                 + (DES_POS_NEW(LN,2)-x_store(ln,2))**2)
      END DO
      write(201,*) 'Max L1 norm', maxval(error(:))
      ENDIF


      RETURN  
      END SUBROUTINE USR1 
