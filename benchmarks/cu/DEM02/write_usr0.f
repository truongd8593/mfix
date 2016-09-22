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

      CALL GEN_PARTICLES

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GEN_PARTICLES()

      use constant
      use desgrid
      use des_allocate
      use discretelement
      use functions
      use physprop
      use mpi_utility

      IMPLICIT NONE

      integer, parameter :: pCount = 1611

      double precision :: ru1, ru2, ru3, ru4, rand3(3)

      double precision :: x0, y0, z0, xLen, yLen, zLen

      double precision :: lRad, lDp
      double precision :: lPos(3), lVel(3), meanVel(3), vScale

      integer :: lc1, fail

      double precision :: gTemp, ldist(3), ldmag

! Initial granular energy (m^2/sec^2)
      double precision, parameter :: T0 = 1.0d-2
      double precision, parameter :: one_third = 1.0d0/3.0d0


      lRad = 0.5d0*D_p0(1)
      lDp = D_p0(1)

      x0 = dg_xStart + lRad
      y0 = dg_yStart + lRad
      z0 = dg_zStart + lRad

      xLen = dg_xEnd - dg_xStart - D_p0(1)
      yLen = dg_yEnd - dg_yStart - D_p0(1)
      zLen = dg_zEnd - dg_zStart - D_p0(1)

      pip = 0
      fail = 0
      seed_lp: do while(pip < pCount)

         if(fail > 50000) then
            write(*,*) 'failed',pip
            stop 23
         endif

         call random_number(rand3)

         lPos(1) = x0 + xLen*rand3(1)
         lPos(2) = y0 + yLen*rand3(2)
         lPos(3) = z0 + zLen*rand3(3)

         do lc1=1, pip
            ldist = des_pos_new(lc1,:) - lPos
            ldmag = sqrt(dot_product(ldist, ldist))

            if(ldmag - ldp < 100.0d-8) then
              fail = fail + 1
              cycle seed_lp
            endif
         enddo

         pip = pip+1
         call particle_grow(pip)
         call set_normal(pip)

         call random_number(ru1)
         call random_number(ru2)
         call random_number(ru3)
         call random_number(ru4)

         lVel(1) = DSQRT(-2.0d0*DLOG(DBLE(ru1)))*COS(2.0d0*PI*ru2)
         lVel(2) = DSQRT(-2.0d0*DLOG(DBLE(ru1)))*SIN(2.0d0*PI*ru2)
         lVel(3) = DSQRT(-2.0d0*DLOG(DBLE(ru3)))*COS(2.0d0*PI*ru4)

         des_pos_new(pip,:) = lPos
         des_vel_new(pip,:) = lVel

         omega_new(pip,:) = 0.0d0

         des_radius(pip) = lRad
         ro_sol(pip) = RO_s0(1)

         meanVel = meanVel + lVel

      enddo seed_lp

! Calc the average mean velocity in each direction
      meanVel = meanVel/dble(pip)

! Subtract mean velocity from the random velocities to get a zero
! mean velocity. Also, calculate the mean granular temperature.
      gTemp = 0.0d0
      do lc1 = 1, pip
         des_vel_new(lc1,:) = des_vel_new(lc1,:) - meanVel
         gTemp = gTemp + dot_product &
            (des_vel_new(lc1,:), des_vel_new(lc1,:))
      enddo
      gTemp = gTemp/(3.0d0*DBLE(pip))

! Scale velocities so the mean granular temperature is equal to
! the targeted valued.
      gTemp = dsqrt(T0/gTemp)

      des_vel_new(:pip,:) = des_vel_new(:pip,:)*gTemp

      particles = pip
      call global_all_sum(particles)


! This write statement keeps the gfortran v6.1.0 version of the code from
! inserting a NAN. This only appears with O3 opt flag.
      if(myPE == PE_IO) write(*,*) 'debug:',des_vel_new(924,:)

      RETURN
      END SUBROUTINE GEN_PARTICLES

      END SUBROUTINE WRITE_USR0
