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
         DO WHILE(LINKS(K,I).gt.0)
           IF(REMOVED.eq.1)THEN
             LINKS(K-1,I)=LINKS(K,I)
             LINKS(K,I)=-1
           ELSE IF(LINKS(K,I).eq.J)THEN
             LINKS(1,I)=LINKS(1,I)-1
             REMOVED=1
             IF(LINKS(K+1,I).lt.0)THEN
               LINKS(K,I)=-1
             END IF
           END IF
         K=K+1
         END DO

         IF(LINKS(1,I).eq.0)THEN
            IS_LINKED(I)=0
         END IF      
      END IF

!-----Remove particle I from the list of particle J
      IF(J.le.PARTICLES)THEN
         REMOVED=0
         K=2
         DO WHILE(LINKS(K,J).gt.0)
           IF(REMOVED.eq.1)THEN
             LINKS(K-1,J)=LINKS(K,J)
             LINKS(K,J)=-1
           ELSE IF(LINKS(K,J).eq.I)THEN
             REMOVED=1
             LINKS(1,J)=LINKS(1,J)-1
             IF(LINKS(K+1,I).lt.0)THEN
               LINKS(K,J)=-1
             END IF
           END IF
         K=K+1
         END DO      

         IF(LINKS(1,J).eq.0)THEN
            IS_LINKED(J)=0
         END IF      

         IF(COHESION_DEBUG.gt.0)THEN
            PRINT *,'**END REMOVE LINK'
         END IF
      END IF

      END SUBROUTINE REMOVE_PART_FROM_LINK_LIST
