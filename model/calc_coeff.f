!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COEFF(DENSITY, SIZE, SP_HEAT, VISC, COND, DIFF,   C
!       RRATE, DRAG, HEAT_TR, WALL_TR, IER)                            C
!  Purpose: Calculate physical and transport properties, reaction ratesC
!           and exchange rates.                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_COEFF(DENSITY, SIZE, SP_HEAT, VISC, COND, DIFF, RRATE, &
         DRAG, HEAT_TR, WALL_TR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE rxns
      USE compar !//AIKEPARDBG
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Loop indices
      INTEGER          L, M
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), SIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!
!                      Flags to tell whether to calculate or not
      LOGICAL          VISC(0:DIMENSION_M), COND(0:DIMENSION_M),&
                       DIFF(0:DIMENSION_M)
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!
!                      Flag for exchange functions
      LOGICAL          DRAG(0:DIMENSION_M, 0:DIMENSION_M),&
                       HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M),&
                       WALL_TR
!-----------------------------------------------
!
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): begin CALC_COEFF')") myPE    !//AIKEPARDBG
!    call mfix_exit(myPE)   !//AIKEPARDBGSTOP

!     Calculate physical properties
!
      CALL PHYSICAL_PROP (DENSITY, SIZE, SP_HEAT, IER) 

!//AIKEPARDBG
    write(*,"('(PE ',I2,'): aft PHYSICAL_PROP in CALC_COEFF, IER=',I4)") myPE,IER    !//AIKEPARDBG


!
!     Calculate Transport properties
!
      CALL TRANSPORT_PROP (VISC, COND, DIFF, IER) 

!//AIKEPARDBG
    write(*,"('(PE ',I2,'): aft TRANSPORT_PROP in CALC_COEFF, IER=',I4)") myPE,IER    !//AIKEPARDBG
!    call mfix_exit(myPE)   !//AIKEPARDBGSTOP

!
!     Calculate reaction rates and interphase mass transfer
!
      IF (RRATE) THEN 
         IF (NO_OF_RXNS > 0) THEN 
            CALL RRATES0 (IER)                   !rxns defined in mfix.dat and rrates0.f 
         ELSE 
            CALL RRATES (IER)                    !rxns defined in rrates.f 
         ENDIF 
      ENDIF 
!//AIKEPARDBG
    write(*,"('(PE ',I2,'): aft RRATES calls in CALC_COEFF, IER=',I4)") myPE,IER    !//AIKEPARDBG
!    call mfix_exit(myPE)   !//AIKEPARDBGSTOP
      
!
!     Calculate interphase momentum, and energy transfers
!
      CALL EXCHANGE (DRAG, HEAT_TR, WALL_TR, IER) 
!
!     Reset all flags.  The flags need to be set every time this routine is
!     called.
!
      RRATE = .FALSE. 
      WALL_TR = .FALSE. 
      M = 0 
      IF (MMAX + 1 > 0) THEN 
         DENSITY(:MMAX) = .FALSE. 
         SIZE(:MMAX) = .FALSE. 
         SP_HEAT(:MMAX) = .FALSE. 
         VISC(:MMAX) = .FALSE. 
         COND(:MMAX) = .FALSE. 
         DIFF(:MMAX) = .FALSE. 
         L = 0 
         DRAG(:MMAX,:MMAX) = .FALSE. 
         HEAT_TR(:MMAX,:MMAX) = .FALSE. 
         L = MMAX + 1 
         M = MMAX + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE CALC_COEFF 
