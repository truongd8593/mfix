!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
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
      SUBROUTINE USR0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      use usr
      use discretelement

      IMPLICIT NONE

      double precision :: t_r

      double precision :: lRad
      double precision :: lGrav

      double precision :: ly, ldydt

      integer lc1
      integer, parameter :: lc_max = 1000

      logical converged

      if(particles /= 1) then
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      endif

      lRad  = des_radius(1)
      h0    = des_pos_new(2,1)
      lGrav = -grav(2)

! Calculate the start time of particle/wall collision.
      time_c =  dsqrt(2.0d0*(h0 - lRad)/lGrav)
! Calculate the particle's velocity at the start of the collision.
      vel_c  = -dsqrt(2.0d0*lGrav*h0)

      b_r  = DES_ETAN_WALL(1)/ (2.0d0 * dsqrt(KN_W * PMASS(1)))
      w0_r = dsqrt(KN_W / PMASS(1))

! Initial loop parameters.
      time_r = time_c + 0.05d0
      lc1 = 0
      converged = .false.

! Calculate the start time for the rebound stage.
      do while (.not.converged)
         lc1 = lc1 + 1

         ly    = y_s2(h0, lRad, b_r, w0_r, lGrav, time_r) - lRad
         ldydt = dydt_s2(h0, lRad, b_r, w0_r, lGrav, time_r)

         if(ldydt /= 0.0d0) then
            t_r = time_r - ly/ldydt
         else
            write(*,"(3x, 'Fatal Error: ldydt == 0')")
            call mfix_exit(0)
         endif

         if(abs(t_r - time_r) < 10.0d-16) then
            converged = .true.
         else
            converged = .false.

            write(*,"(3x,I4,2(3x,F14.8))") lc1, t_r, time_r

         endif

         time_r = t_r

         if(lc1 > lc_max) then
            write(*,"(3x, 'Fatal Error: lc1 > lc_max')")
            call mfix_exit(0)
         endif

      enddo
! Calculate the velocity at the end of the rebound stage.
      vel_r = dydt_s2(h0, lRad, b_r, w0_r, lGrav, time_r)


      write(*,"(//3x,'Collision:')")
      write(*,"(5x,'Time:     ',F14.8)") time_c
      write(*,"(5x,'Velocity: ',F14.8)") vel_c

      write(*,"(/3x,'Rebound:')")
      write(*,"(5x,'Time: ',F14.8)") time_r
      write(*,"(5x,'Velocity: ',F14.8)") vel_r

      write(*,"(//' ')")

      return

      END SUBROUTINE USR0
