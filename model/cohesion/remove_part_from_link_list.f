!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REMOVE_PART_FROM_LINK_LIST(I, J)                       C
!  Purpose: This function will remove two particles form their         C
!           respective linked lists                                    C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE REMOVE_PART_FROM_LINK_LIST(I, J)

!-----MODULES USED
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER K, I, J, JJ
      INTEGER LINK_PARTNER
      INTEGER REMOVED
      LOGICAL ALREADY_EXISTS


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START REMOVE LINK'
      END IF 




!-----Remove particle J from list of particle I
      IF(I.le.PARTICLES)THEN
         REMOVED=0
         K=2
         DO WHILE(LINKS(I,K).gt.0)
           IF(REMOVED.eq.1)THEN
             LINKS(I,K-1)=LINKS(I,K)
             LINKS(I,K)=-1
           ELSE IF(LINKS(I,K).eq.J)THEN
             LINKS(I,1)=LINKS(I,1)-1
             REMOVED=1
             IF(LINKS(I,K+1).lt.0)THEN
               LINKS(I,K)=-1
             END IF
           END IF
         K=K+1
         END DO

         IF(LINKS(I,1).eq.0)THEN
            IS_LINKED(I)=0
         END IF      
      END IF

!-----Remove particle I from the list of particle J
      IF(J.le.PARTICLES)THEN
         REMOVED=0
         K=2
         DO WHILE(LINKS(J,K).gt.0)
           IF(REMOVED.eq.1)THEN
             LINKS(J,K-1)=LINKS(J,K)
             LINKS(J,K)=-1
           ELSE IF(LINKS(J,K).eq.I)THEN
             REMOVED=1
             LINKS(J,1)=LINKS(J,1)-1
             IF(LINKS(I,K+1).lt.0)THEN
               LINKS(J,K)=-1
             END IF
           END IF
         K=K+1
         END DO      

         IF(LINKS(J,1).eq.0)THEN
            IS_LINKED(J)=0
         END IF      

         IF(COHESION_DEBUG.gt.0)THEN
            PRINT *,'**END REMOVE LINK'
         END IF
      END IF

      END SUBROUTINE REMOVE_PART_FROM_LINK_LIST
