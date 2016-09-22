!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Author:                                            Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: Write user-defined output                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR1(L)

      use run, only: TIME

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L


      SELECT CASE(L)
      CASE(1); CALL WRITE_TEST_DATA
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_TEST_DATA

      use compar, only: myPE, PE_IO
      use run, only: time
      use discretelement, only: max_pip, particles
      use discretelement, only: des_vel_new
      use functions, only: is_normal
      use mpi_utility, only: global_sum

      implicit none

      integer :: lc1, lc2
      double precision :: ltime, gTemp, gTemp0

      integer, parameter :: fUnit=2030


! This is a test of the setup
      lc2 = 0
      gTemp = 0.0d0
      do lc1=1, max_pip
         if(is_normal(lc1)) then
            gTemp = gTemp + dot_product &
               (des_vel_new(lc1,:),des_vel_new(lc1,:))
            lc2 = lc2 + 1
         endif
      enddo

      call global_sum(gTemp, gTemp0)

      if(myPE == PE_IO) then
         gTemp = gTemp0/(3.0d0*DBLE(particles))
         ltime = time*sqrt(0.1d0)/100.0e-6

         OPEN(UNIT=fUnit,FILE='POST_GRAN_TEMP.dat',&
            POSITION="APPEND",STATUS='OLD')
         write(fUnit,"(3(2x,es13.6))") ltime, ToT0(time), gTemp/0.1d0

         close(fUnit)
      endif

      contains

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      double precision function ToT0(ptime)

      double precision, intent(in) :: pTime

      double precision, parameter :: A = 222.511014d0
      double precision, parameter :: B =  91.382708d0

      double precision :: expBto2

      expBto2 = exp(B*pTime/2.0d0)

      ToT0 = 1.0d0/(expBto2 + (A/B)*sqrt(0.10d0)*(expBto2-1.0d0))
      ToT0 = ToT0**2

      end function ToT0

      END SUBROUTINE WRITE_TEST_DATA
