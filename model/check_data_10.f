!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_09                                          C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_10

      USE ps

      IMPLICIT NONE


!----------------------------------------------------------------------!
!        This is the initial check-in of point souce routines.         !
!----------------------------------------------------------------------!


! Restrict access to point souce routines.
      POINT_SOURCE = .FALSE.

! To-do: Migrate check routines from usr0.f to this routine and include
!        additional checks on user data.

      RETURN
      END SUBROUTINE CHECK_DATA_10
