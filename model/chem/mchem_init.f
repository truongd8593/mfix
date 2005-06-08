!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM_INIT                                             C
!     Purpose: INITIAL values for rxns calcs                              C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!      
      SUBROUTINE MCHEM_INIT
!
!-------------------
!  M o d u l e s
!-------------------
!
      USE param1
      USE run
      USE physprop
      USE mchem

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------      
!
!             counter
      INTEGER NL
!
!     Time in isat calcualtion
!
      IF (CALL_CHEM .AND. (ISATdt .EQ. UNDEFINED)) THEN
         TIME_isat = TIME + DT
      ELSE
         TIME_isat = TIME + ISATdt
      END IF
!
!     Dimension of ODEs solved in ISAT
!
      NSpec = 0
!
!     For X_G and X_S
! 
      DO NL = 0, MMAX
         NSpec = NMAX(NL) + NSpec
      END DO
!
!     For RO_G
!
      NSpec = NSpec + 1
!
!     For EP_s
      NSpec = NSpec + MMAX
!
!     For T_G and T_S
!
!      IF (ENERGY_EQ) THEN
         NSpec = NSpec + 1 + MMAX
!      END IF
!
!     For D_p
!       
      NSpec = NSpec + MMAX
!
      RXN_source_g = ZERO
      RXN_source_s = ZERO

      RETURN
      END SUBROUTINE MCHEM_INIT
