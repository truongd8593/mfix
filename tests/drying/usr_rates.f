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


! Drying:  Liquid --> Vapor         (mol/cm^3.s)
!---------------------------------------------------------------------//
! Note: The rate is artificial and used for the drying test case.


      IF(X_s(IJK,1,Liquid) > ZERO) THEN
         RATES(Drying) = C(1) * ROP_s(IJK, 1) * &
            max(ZERO, (0.5d0 - X_g(IJK,Vapor))) / MW_s(1,Liquid)
      ELSE
         RATES(Drying) = ZERO
      ENDIF

      RETURN  

      END SUBROUTINE USR_RATES
