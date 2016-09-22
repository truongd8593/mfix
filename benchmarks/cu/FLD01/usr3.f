!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Purpose: This routine is called after the time loop ends and is     C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.              C
!           This routine is not called from an IJK loop, hence         C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR3

      use fldvar, only: U_G, V_G

      IMPLICIT NONE

      double precision, parameter :: X(16) = (/&
          0.96880, 0.96090, 0.95310, 0.94530, &
          0.90630, 0.85940, 0.80470, 0.50000, &
          0.23440, 0.22660, 0.15630, 0.09380, &
          0.07810, 0.07030, 0.06250, 0.00000/)

      double precision, parameter :: Y(16) = (/&
          0.97660, 0.96880, 0.96090, 0.95310, &
          0.85160, 0.73440, 0.61720, 0.50000, &
          0.45310, 0.28130, 0.17190, 0.10160, &
          0.07030, 0.06250, 0.05470, 0.00000/)

! U-Velocities from Ghia et al. (1982), Table I at X=0.5
      double precision, parameter :: U_Re0100(16) = (/&
          0.84123, 0.78871, 0.73722, 0.68717, &
          0.23151, 0.00332,-0.13641,-0.20581, &
         -0.21090,-0.15662,-0.10150,-0.06434, &
         -0.04775,-0.04192,-0.03717, 0.00000/)

! V-Velocities from Ghia et al. (1982), Table II at Y=0.5
      double precision, parameter :: V_Re0100(16) = (/&
         -0.05906,-0.07391,-0.08864,-0.10313, &
         -0.16914,-0.22445,-0.24533, 0.05454, &
          0.17527, 0.17507, 0.16077, 0.12317, &
          0.10890, 0.10091, 0.09233, 0.00000/)

      double precision, parameter :: lHalf(16) = 0.5

      call MFIX_to_GHIA(lHalf, Y, U_g, U_Re0100, 'U')
      call MFIX_to_GHIA(X, lHalf, V_g, v_Re0100, 'V')

      RETURN
      contains

!......................................................................!
      SUBROUTINE MFIX_to_GHIA(pX, pY, pVEL, pGhia, pVar)

      use geometry, only: dx, dy
      use functions, only: IS_ON_MYPE_OWNS
      use functions, only: FUNIJK
      use mpi_utility, only: GLOBAL_ALL_SUM
      use compar, only: myPE, PE_IO

      double precision, intent(in) :: pX(:), pY(:)
      double precision, intent(in) :: pVel(:)
      double precision, intent(in) :: pGhia(:)
      character(len=1), intent(in) :: pVar

      integer :: lc1, i, j, k
      double precision :: lOoDx, lOoDy
      double precision :: Sx(0:1), Sy(0:1)

      double precision :: lMFIX(16)

      integer, parameter :: fUnit = 2030

      lOoDx = 1.0d0/dx(1)
      lOoDy = 1.0d0/dy(1)

      lMFIX = 0.0d0

      k=1
      do lc1=1,15

         i = int(lOoDx*pX(lc1))+1
         Sx(0) = (dble(i)*dx(1) - pX(lc1))*lOoDx
         Sx(1) = 1.0d0 - Sx(0)

         j = int(lOoDy*pY(lc1))+1
         Sy(0) = (dble(j)*dy(1) - pY(lc1))*lOoDy
         Sy(1) = 1.0d0 - Sy(0)

         if(is_On_myPE_OWNS(i,j,k)) then
            lMFIX(lc1) = &
               pVel(funijk(i-1,j-1,k))*Sx(0)*Sy(0) + &
               pVel(funijk(i  ,j-1,k))*Sx(1)*Sy(0) + &
               pVel(funijk(i-1,j  ,k))*Sx(0)*Sy(1) + &
               pVel(funijk(i  ,j  ,k))*Sx(1)*Sy(1)
         endif
      enddo

      call global_all_sum(lMFIX)
      if(myPE /= PE_IO) RETURN

! Open file for output
      open(unit=fUnit, file='POST_'//pVar//'G.dat',&
         status='unknown')
! Write data.
      write(fUnit,1000)
      do lc1=1,15
         write(fUnit,1100) X(lc1), pGhia(lc1), lMFIX(lc1)
      enddo
      close(fUnit)

 1000 Format(2/4x,'Position',8x,'Ghia',11x,'MFIX')
 1100 Format(3(2x,es13.6))

      END SUBROUTINE MFIX_to_GHIA

      END SUBROUTINE USR3
