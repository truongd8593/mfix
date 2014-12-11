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
!  Author: S. Venkatesan, M. Syamlal                  Date: 21-JUN-93  C
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
      SUBROUTINE USR0

      use compar
      use cutcell
      use geometry
      use mpi_utility
      use param
      use param1
      use usr
      use output

! User-specified constant
      use constant, only: C
! Total number of solids phases.
      use physprop, only: MMAX
! Initial composition of solids.
      use physprop, only: X_s0
! Flag for variable solids density.
      use run, only: SOLVE_ROs

      implicit none

      INTEGER :: L

! External Function for comparing two numbers.
      LOGICAL, EXTERNAL :: COMPARE

      include 'species.inc'

      Allocate( N_sh (DIMENSION_3, DIMENSION_M) )

! Set the coal proximate and ultimate analysis:
      IF(SOLVE_ROs(1)) THEN
         PAM = X_s0(1,Moisture)  ! Moisture
         PAV = X_s0(1,Volatiles) ! Volatiles
         PAC = X_s0(1,Char)      ! Char
         PAA = X_s0(1,Ash)       ! Ash
      ELSE
          PAC = C(Char)
          PAV = C(Volatiles)
          PAM = C(Moisture)
          PAA = C(Ash)
      ENDIF

! Verify that the proximate analysis sums to one.
      IF(.NOT.compare(ONE,(PAC+PAV+PAM+PAA))) THEN
         IF(myPE == PE_IO) THEN
            write(*,"(3/,3x,'Proximate analysis does NOT sum to one.')")
            write(*,"(3x,'Check values in mfix.dat')")
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

      END SUBROUTINE USR0
