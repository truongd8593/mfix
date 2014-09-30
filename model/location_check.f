!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: LOCATION_CHECK                                          C
!  Purpose: check calculated and given cell locations for consistency  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE LOCATION_CHECK(CELL_SPECIFIED, CELL_CALCULATED, &
         COUNTER, MESSAGE)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE funits
      USE geometry
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! cell index specified in the input file
      INTEGER, INTENT(IN) :: CELL_SPECIFIED
! cell index calculated for location coordinate specified in the
! data file
      INTEGER, INTENT(IN) :: CELL_CALCULATED
! index for BC or IC
      INTEGER, INTENT(IN) :: COUNTER
! error message to print out
      CHARACTER MESSAGE*(*)
!-----------------------------------------------


! check that the cell_specified in the data input equals to the
! cell calculated.  If not equal, print error message and stop

      IF (CELL_SPECIFIED == CELL_CALCULATED) RETURN

      IF (NO_K) THEN
         IF (MESSAGE(6:6)=='b' .OR. MESSAGE(6:6)=='t') RETURN
      ENDIF
      IF (NO_J) THEN
         IF (MESSAGE(6:6)=='s' .OR. MESSAGE(6:6)=='n') RETURN
      ENDIF
      IF (NO_I) THEN
         IF (MESSAGE(6:6)=='w' .OR. MESSAGE(6:6)=='e') RETURN
      ENDIF

      CALL ERROR_ROUTINE ('location_check', 'consistency error', 0, 2)
      IF(DMP_LOG)WRITE (UNIT_LOG, 1000) MESSAGE, COUNTER, &
         CELL_SPECIFIED, CELL_CALCULATED
      CALL ERROR_ROUTINE (' ', ' ', 1, 3)

 1000 FORMAT(1X,'IC, BC, or IS error for : ',A,/,1X,'IC/BC/IS No  = ',&
      I6,/, 1X,'Cell specified  = ',I6,/,1X,'Cell calculated = ',I6)

      RETURN
      END SUBROUTINE LOCATION_CHECK
