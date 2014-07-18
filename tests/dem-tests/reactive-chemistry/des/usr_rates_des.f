!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: USR_RATES_DES                                           !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Comments: Write reaction rates in units of moles/ses.               !
!                                                                      !
!  WARNING: Only discrete phase reactions should be specifed here.     !
!  Homogeneous gas phase reactions in DEM simulations must be given    !
!  in the continuum reaction hook (usr_rates.f).                       !
!                                                                      !
!  The call to usr_rates_des is made from inside a particle loop which !
!  is nested inside an IJK loop. Fluid grid calculations independent   !
!  of particle properties can be carried out in des/usr4_des.f to      !
!  reduce redundant calculations.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES_DES(NP, pM, IJK, DES_RATES)

      use des_rxns, only: NO_OF_DES_RXNS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP  ! Index of particle
      INTEGER, INTENT(IN) :: pM  ! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP

! Calculated reaction rates. (reacted moles per sec)
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)


      INCLUDE '../species.inc'


! EX_RXN:    A(g) + 2B(s) --> C(g) + D(s)
!`````````````````````````````````````````````````````````````````````\\
      DES_RATES(EX_RXN) = 3.927d-5 ! (moles/sec)


      RETURN  
      END SUBROUTINE USR_RATES_DES
