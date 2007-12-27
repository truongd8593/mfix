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
!  Local variables:                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION CALC_H (IJK, M, N) 
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
      INTEGER          IJK, M, N
      DOUBLE PRECISION, EXTERNAL ::calc_ICpoR
!-----------------------------------------------
!
      IF( M == 0) THEN
        CALC_H =  ( HfrefoR_g(N)  + &
	            (calc_ICpoR(T_G(IJK), Thigh_g(N), Tlow_g(N), &
		      Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N)) -  &
			      IC_PGrefoR(N)) &
		   ) * GAS_CONST_cal / MW_g(N)
      ELSE
         CALC_H = ( HfrefoR_s(M, N)  + &
	            (calc_ICpoR(T_s(IJK, M), Thigh_s(M,N), Tlow_s(M,N), &
			     Tcom_s(M,N), Ahigh_s(1,M,N), Alow_s(1,M,N)) -  &
			      IC_PsrefoR(M,N)) &
		   ) * GAS_CONST_cal / MW_s(M,N)
      ENDIF

      RETURN  
      END FUNCTION CALC_H
