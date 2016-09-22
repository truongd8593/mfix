!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR0

      CALL GET_PARTICLES

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GET_PARTICLES()

      use constant
      use desgrid
      use des_allocate
      use discretelement
      use functions
      use physprop
      use mpi_utility

      IMPLICIT NONE

      double precision :: lOffset(3)
      double precision :: lPos(3), lVel(3), lRad, lROs
      integer :: ios, lc1

      lOffset(1) = dg_xStart
      lOffset(2) = dg_yStart
      lOffset(3) = dg_zStart

      open(UNIT=6587, FILE='particle_input.dat', FORM="FORMATTED")

      pip = 0
      read_lp: do
         read(6587,*,iostat=ios) (lPos(lc1),lc1=1,3),&
            lRad, lROs, (lVel(lc1),lc1=1,3)

         if(ios > 0) stop 6587
         if(ios < 0) exit read_lp

         pip = pip + 1

         call particle_grow(pip)
         call set_normal(pip)

         des_pos_new(pip,:) = lPos*1.0d-2 + lOffset
         des_vel_new(pip,:) = lVel*1.0d-2

         des_radius(pip) = lRad*1.0d-2
         ro_sol(pip) = lROs*1.0d3

         omega_new(pip,:) = 0.0d0
      enddo read_lp

      close(6587)

      particles = pip
      call global_all_sum(particles)

      RETURN
      END SUBROUTINE GET_PARTICLES

      END SUBROUTINE WRITE_USR0
