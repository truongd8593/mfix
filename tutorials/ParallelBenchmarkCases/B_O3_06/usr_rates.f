!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE sendrecv

      USE toleranc
      USE usr

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

!-----------------------------------------------
      INCLUDE 'species.inc'

      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'

      INCLUDE 'function.inc'

      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'

      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
      DOUBLE PRECISION c_O3   ! Ozone concentration mol/cm^3

! PhseCng:  O3 --> 1.5O2         (mol/cm^3.s)
!---------------------------------------------------------------------//
      c_O3 = (RO_g(IJK)*X_g(IJK,O3)/MW_g(O3))

      RATES(Ozone_Decomp) = (ONE - EP_g(IJK)) * C(1) * c_O3

      RETURN

      END SUBROUTINE USR_RATES
