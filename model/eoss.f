!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: EOSS (IJK,M)                                           !
!                                                                      !
!  Author: J.Musser                                   Date: 09-Oct-13  !
!  Reviewer:                                                           !
!                                                                      !
!  Purpose: Equation of state for gas                                  !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION EOSS(IJK, M)

! Global Variables:
!---------------------------------------------------------------------//
! Baseline/initial solids density
      use physprop, only: RO_s0
! Baseline/initial solids mass fractions.
      use physprop, only: X_s0
! Index of inert species.
      use physprop, only: INERT_SPECIES
! Solids mass fraction
      USE fldvar, only: X_s

      implicit none

! Passed Arguments:
!---------------------------------------------------------------------/
! Fluid cell index.
      INTEGER, intent(in) :: IJK
! Solids phase index.
      INTEGER, intent(in) :: M

! Local Variables:
!---------------------------------------------------------------------/
! Alias for inert species index.
      INTEGER :: INERT 


! Set the inert species alias.
      INERT = INERT_SPECIES(M)

! Evaluate the solids EOS.
      EOSS = RO_s0(M) * X_s0(M,INERT) / max(X_s(IJK,M,INERT), 1.0d-8)

      RETURN  
      END FUNCTION EOSS
