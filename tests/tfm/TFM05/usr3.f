!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Author: Avinash Vaidheeswaran                      Date: July, 2016 C
!  Reviewer: J.Musser                                                  C
!                                                                      C
!  Purpose: Calculates teh exact solution to the Couette flow case and C
!  compares with the MFIX solution.                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR3

      use fldvar, only: U_g
      Use geometry, only: dy
      Use geometry, only: HEIGHT => YLENGTH
      Use geometry, only: XLENGTH

      Use bc, only: Uw => BC_Uw_g
      Use bc, only: delP_x
      use physprop, only: Mu_g0

      Use geometry, only: imin1, jmin1, kmin1
      Use geometry, only: imax1, jmax1, kmax1

      use functions, only: funijk

      Use param1, only: zero, small_number, half

      IMPLICIT NONE

! looping indices
      integer :: ijk, i, j, k
! temporary variable
      double precision  :: TMPdp
! Calculated height of cell center
      double precision  :: yt
! Exact and numerical solutions
      double precision  :: Ug, Ug_MFIX
! Absolute and absolute relative errors
      double precision :: absErr, relErr
! file unit for output data
      integer, parameter :: fUNIT = 2030

! Open file for output.
      OPEN(UNIT=fUNIT, FILE='POST_VEL.dat', &
         POSITION="APPEND",STATUS='OLD')

! Calculate temp value and initial height
      TMPdp = (HALF/MU_g0)*(-delP_x/(xLength))
      yt = -HALF*dy(jmin1)

! generate grid locations for exact solution calculation
      k = kmin1 + (kmax1-kmin1)/2
      i = imin1 + (imax1-imin1)/2
      do j = jmin1, jmax1
! Calculate cell height (assumed uniform spacing)
         yt = yt + dy(j)
! Calculate exact solution
!         Ug = (TMPdp*yt + ((Uw(1)/height) - TMPdp*height))*yt

         Ug = Uw(1)*(yt/height) + &
            (HALF/Mu_g0)*(-delP_x/xLength)*(yt**2 - height*yt)

! Get the MFIX solution
         Ug_MFIX = U_G(funijk(i,j,k))

         absErr = abs(Ug - Ug_MFIX)
         relErr = abs(absErr/max(abs(Ug), SMALL_NUMBER))

         write(fUnit, 1200) yt, Ug, Ug_MFIX, absErr, relErr

      end do

      write(fUnit,"(5/)")
      close(fUnit)


 1200 FORMAT(5(3x,es13.6))

      RETURN
      END SUBROUTINE USR3
