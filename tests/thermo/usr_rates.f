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
      INTEGER :: L

      L = IJK - 49
      RATES(L) = ONE

      RETURN  
      END SUBROUTINE USR_RATES
