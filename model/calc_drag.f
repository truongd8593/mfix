!//NOMOD 1117 No modifications necessary
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DRAG(DRAG, IER)                                   C
!  Purpose: Calculate gas-solids and solids-solids drag terms          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for variable grid size capability,           C
!           logic for volume-weighted averaging                        C
!  Author: W. Rogers                                  Date: 20-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: MFIX 2.0 mods                                              C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, U_g, V_g, W_g, MMAX      C
!                        U_s, V_s, W_s, DT                             C
!  Variables modified: I, J, K, IJK, M, DTxF_gs, DTxF_ss               C
!                                                                      C
!  Local variables: UGC, USCM, USCL, VGC, VSCM, VSCL, WGC, WSCM, WSCL, C
!                   VREL, L, LM                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_DRAG(DRAGD, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE drag
      USE compar     !//AIKEPARDBG
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
!                      Solids phase
      INTEGER          M
!
!                      Flag for exchange functions
      LOGICAL          DRAGD(0:DIMENSION_M, 0:DIMENSION_M)
!
!                      Local index for solids phase l
      INTEGER          L
!-----------------------------------------------
!
!
!
!
!  Function subroutines
!
!
!  Local variables
!
!
!
      DO M = 1, MMAX 
         IF (DRAGD(0,M) .AND. RO_G0/=ZERO) CALL DRAG_GS (M, F_GS, IER) 
         DO L = 1, M - 1 
            IF (DRAGD(L,M)) CALL DRAG_SS (L, M, F_SS, IER) 
         END DO 
      END DO 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): eof CALC_DRAG')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP
      
      RETURN  
      END SUBROUTINE CALC_DRAG 
