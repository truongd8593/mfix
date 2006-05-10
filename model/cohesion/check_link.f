!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FUNCTION CHECK_LINK(I,J)                               C
!  Purpose: Integer function to check if two particles are             C
!      linked (ie, their center-to-center distance                     C
!      is less that the sum of their well widths)                      C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      INTEGER FUNCTION CHECK_LINK(I,J)

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

      INTEGER K,I,J
         
      CHECK_LINK = 0
      IF(LINKS(I,1).gt.0)THEN
        K = 2
        DO WHILE(LINKS(I,K).gt.0)
          IF(LINKS(I,K).eq.J)THEN
            CHECK_LINK=1
            GO TO 10
          END IF
          K = K+1
        END DO
      END IF
10    RETURN
      END FUNCTION CHECK_LINK
