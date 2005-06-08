!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_CHEM                                        C
!  Purpose: check the chemical rxns namelist variables for             C
!           CALL_CHEM or CALL_ISAT                                     C
!                                                                      C
!  Author: NAN XIE                              Date: 02-Aug-04        C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Variables referenced:  CALL_CHEM , CALL_ISAT , ISATdt               C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_CHEM
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param1
      USE run
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!-----------------------------------------------
!
      IF (CALL_CHEM .AND. CALL_ISAT) THEN
         WRITE(*,*) 'CALL_CHEM and CALl_ISAT can not be set true together.'
         CALL MFIX_EXIT(myPE)
      END IF

      IF (CALL_ISAT) THEN
         IF (ISATdt .EQ. UNDEFINED) THEN
            WRITE(*,*) 'ISATdt should be provided if CALL_ISAT is true.'
            CALL MFIX_EXIT(myPE)
         END IF
      END IF


      RETURN
      END SUBROUTINE CHECK_DATA_CHEM
