!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_H (Tl, M, N)                                       C
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
!  Local variables:                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION CALC_H (Tl, M, N) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE constant
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      cell, phase and species indices
      INTEGER          M, N

!                      Temperature of origination phase
      DOUBLE PRECISION Tl
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
!-----------------------------------------------
!
      IF( M == 0) THEN
        CALC_H =  ( HfrefoR_g(N)  + &
	            (calc_ICpoR(Tl, Thigh_g(N), Tlow_g(N), &
		      Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N)) -  &
			      IC_PGrefoR(N)) &
		   ) * GAS_CONST_cal / MW_g(N)
      ELSE
         CALC_H = ( HfrefoR_s(M, N)  + &
	            (calc_ICpoR(Tl, Thigh_s(M,N), Tlow_s(M,N), &
			     Tcom_s(M,N), Ahigh_s(1,M,N), Alow_s(1,M,N)) -  &
			      IC_PsrefoR(M,N)) &
		   ) * GAS_CONST_cal / MW_s(M,N)
      ENDIF

      RETURN  
      END FUNCTION CALC_H


!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
! Function: DES_CALC_H                                                 !
!                                                                      !
! Purpose: Calculate the enthalpy for des_rrates subroutine.           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION DES_CALC_H(NP, N)

      USE constant
      USE des_rxns
      USE des_thermo
      USE discretelement

      IMPLICIT NONE

! Particle index
      INTEGER, INTENT(IN) :: NP
! Species index
      INTEGER, INTENT(IN) :: N
! Index for the discrete solids phase of the particle
      INTEGER M
! Reference temperature
      DOUBLE PRECISION ICpoR
      DOUBLE PRECISION, EXTERNAL :: CALC_ICpoR

      M = PIJK(NP,5)

!     Enthalpy
      ICpoR = CALC_ICpoR(DES_T_s_NEW(NP), DES_Thigh_s(M,N), &
         DES_Tlow_s(M,N), DES_Tcom_s(M,N), DES_Ahigh_s(1,M,N), &
         DES_Alow_s(1,M,N))

       DES_CALC_H = (DES_HfrefoR_s(M,N)+(ICpoR-DES_IC_PsrefoR(M,N))) * &
		        GAS_CONST_cal / DES_MW_s(M,N)

      END FUNCTION DES_CALC_H
