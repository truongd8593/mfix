!//NOMOD 1117 No modifications necessary

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EXCHANGE(DRAG, HEAT_TR, WALL_TR, IER)                  C
!  Purpose: Calculate interphase mass, momentum, and energy exchange   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
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
      SUBROUTINE EXCHANGE(DRAG, HEAT_TR, WALL_TR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE compar  !//AIKEPARDBG
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IER 
      LOGICAL WALL_TR 
!                      Flag for exchange functions
      LOGICAL, DIMENSION(0:DIMENSION_M,0:DIMENSION_M) :: DRAG, HEAT_TR 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!
!
!     Calculate drag coefficients
!
      CALL CALC_DRAG (DRAG, IER) 
!
!     Calculate interphase heat transfer coefficients
!
      CALL CALC_GAMA (HEAT_TR, IER) 
!
!
!//AIKEPARDBG
      write(*,"('(PE ',I2,'): eof EXCHANGE')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP

      RETURN  
      END SUBROUTINE EXCHANGE 
