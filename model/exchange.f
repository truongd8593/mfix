!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: EXCHANGE                                                C
!  Purpose: Calls routines to drive calculations of the interphase     C
!           mass, momentum, and energy exchange coefficients/terms     C
!           if directed to do so by the corresponding flags            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE EXCHANGE(DRAG, HEAT_TR, IER) 

!-----------------------------------------------
! Module
!-----------------------------------------------
      USE param 
      USE param1 
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Flag for exchange functions
      LOGICAL, INTENT(IN) :: DRAG(0:DIMENSION_M,0:DIMENSION_M),&
                             HEAT_TR(0:DIMENSION_M,0:DIMENSION_M) 
! Error index
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Currently unused flag
      LOGICAL WALL_TR
!-----------------------------------------------


! Calculate drag coefficients
      CALL CALC_DRAG (DRAG, IER) 

! Calculate interphase heat transfer coefficients
      CALL CALC_GAMA (HEAT_TR, IER) 

      RETURN  
      END SUBROUTINE EXCHANGE 


