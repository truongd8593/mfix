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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use run
      use usr
      use discretelement

      IMPLICIT NONE

      if(particles /= 1) then
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      endif

! Store the intial translational velocity.
      v0 = des_vel_new(1,1)

! Calculate the end of the slip
      end_slip = -(2.0d0*v0)/(7.0d0*MEW_w*Grav(2))


      write(*,"(//3x,'Initial translational velocity:',g18.8)") v0
      write(*,"(3x,'Gravity:  ',g18.8)") GRAV(2)
      write(*,"(3x,'End of slip: ',g18.8)") end_slip

      if(tstop < end_slip) then
         write(*,"(/3x,'Simulation too short. It will not catch')")
         write(*,"( 3x,'the end of particle/wall slip.')")
         call mfix_exit(0)
      endif

      return
      END SUBROUTINE USR0
