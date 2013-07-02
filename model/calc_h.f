!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_H (IJK, M, N)                                       C
!  Purpose: Calculate specific enthalpy of species N in phase M        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-DEC-2007C
!  Reviewer:                                         Date:   C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:None                                           C
!  Variables modified:None                                             C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION CALC_H(refT, M, N) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE fldvar

      USE constant, only: RGAS => GAS_CONST_cal

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      cell, phase and species indices

      DOUBLE PRECISION, INTENT(IN) :: refT   ! Temperature

      INTEGER, INTENT(IN) :: M ! Phase index
      INTEGER, INTENT(IN) :: N ! Species index

      DOUBLE PRECISION ICpoR
      DOUBLE PRECISION lMW

      INTEGER :: IER

      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
!-----------------------------------------------
!
      IER = 0

      if(M == 0)then
         lMW = MW_g(N)
      else
         lMW = MW_s(M,N)
      endif

! Integrate the specific heat from zero to refT
      ICpoR = calc_ICpoR(refT, M, N, IER)

! Evaluate the enthalpy of speices N at refT
      CALC_H = (HfrefoR(M,N)  + ICpoR) * (RGAS / lMW)

      RETURN
      END FUNCTION CALC_H


!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
! Function: DES_CALC_H0                                                !
!                                                                      !
! Purpose: Calculate the enthalpy for des_rrates subroutine.           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION DES_CALC_H0(NP, N)

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
      use discretelement, only: PIJK
      use des_thermo, only: DES_T_s_NEW
      use des_rxns, only: DES_HfrefoR
      use des_rxns, only: DES_MW_s

      implicit none

! Particle index
      INTEGER, INTENT(IN) :: NP
! Species index
      INTEGER, INTENT(IN) :: N
! Index for the discrete solids phase of the particle
      INTEGER M
! Reference temperature
      DOUBLE PRECISION ICpoR
      DOUBLE PRECISION, EXTERNAL :: DES_CALC_ICpoR

      INTEGER :: IER

      IER = 0
      M = PIJK(NP,5)

! Sensible heat (refT to T)
      ICpoR = DES_CALC_ICpoR(DES_T_s_NEW(NP), M,N, IER)

      DES_CALC_H0 = (DES_HfrefoR(M,N)+ ICpoR) * (RGAS/DES_MW_s(M,N))

      END FUNCTION DES_CALC_H0
