!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: DES_LINKED_LIST_FUNCS_MOD                              !
!                                                                      !
!                                                                      !
!  Reviewer: R. Garg                                  Date: 19-Mar-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!




      MODULE DES_LINKED_LIST_FUNCS


      CONTAINS

      SUBROUTINE merge_part_lists(newlist, baselist)
      USE des_linked_list_data, only : particle

      IMPLICIT NONE
      TYPE(PARTICLE), POINTER  :: newlist, baselist

      TYPE(particle), pointer :: part

      if(.not.associated(newlist)) then
         !do nothing
         return
      endif

      If(.not.associated(baselist)) then
         !just point the baselist to newlist
         baselist => newlist
      else
         nullify(part)
         part => baselist
         DO WHILE (ASSOCIATED(part%next))
            part => part%next
         ENDDO

         part%next => newlist
         newlist%prev => part
      endif

      END SUBROUTINE merge_part_lists

      SUBROUTINE DEALLOC_PART_LIST(LIST_NAME)
      USE des_linked_list_data, only : particle

      IMPLICIT NONE

      TYPE(PARTICLE), POINTER  :: LIST_NAME, part => NULL(), part_old => NULL()

      IF(.not.associated(list_name)) return !do nothing

      part =>  list_name

      DO WHILE (ASSOCIATED(part))
         part_old => part
         part => part%next
         DEALLOCATE(part_old)
      enddo
      RETURN
      END SUBROUTINE DEALLOC_PART_LIST
      SUBROUTINE GEN_AND_ADD_TO_PART_LIST(LIST_NAME, PHASE, POS, VEL,RAD, DENS, STATWT)
      USE discretelement, only : dimn
      USE des_linked_list_data, only : particle

      IMPLICIT NONE

      TYPE(PARTICLE), POINTER  :: LIST_NAME
      INTEGER, INTENT(IN) ::  PHASE
      double precision, INTENT(IN), DIMENSION(DIMN) :: POS, VEL
      double precision, INTENT(IN) :: RAD, DENS, STATWT

      TYPE(PARTICLE), POINTER  :: NEW_PART => NULL()

      ALLOCATE(NEW_PART)

      NEW_PART%M = PHASE
      NEW_PART%INDOMAIN = .true.

      NEW_PART%POSITION(1:DIMN) = POS(1:DIMN)
      NEW_PART%VELOCITY(1:DIMN) = VEL(1:DIMN)
      NEW_PART%STATWT = STATWT

      NEW_PART%RAD = RAD
      NEW_PART%DENS = DENS
      NEW_PART%STATWT = STATWT

      IF(.NOT.ASSOCIATED(LIST_NAME)) THEN !FIRST ENTRY
         LIST_NAME => NEW_PART  !make the new entry as head
      ELSE
         NEW_PART%NEXT => LIST_NAME
         LIST_NAME%PREV => NEW_PART
         LIST_NAME =>  NEW_PART
      END IF
      END SUBROUTINE GEN_AND_ADD_TO_PART_LIST


      SUBROUTINE ADD_TO_PART_LIST(NEW, BASE)
      USE discretelement, only : dimn
      USE des_linked_list_data, only : particle

      IMPLICIT NONE

      !New entry in the linked list
      TYPE(PARTICLE), POINTER  :: NEW

      !Baseline linked list
      TYPE(PARTICLE), POINTER  :: BASE
      IF(.NOT.ASSOCIATED(BASE)) THEN !FIRST ENTRY
         BASE => NEW
      ELSE
         NEW%NEXT => BASE
         BASE%PREV => NEW
         BASE =>  NEW
      END IF

      END SUBROUTINE ADD_TO_PART_LIST

      SUBROUTINE REMOVE_PART_LLIST(PART)
      USE discretelement, only : dimn
      USE des_linked_list_data, only : particle
      IMPLICIT NONE
      LOGICAL :: DELETION
      TYPE(PARTICLE), POINTER  :: PART
      TYPE(PARTICLE), POINTER  :: RIGHT, LEFT

      !CALL INIT_ERR_MSG("REMOVE_PART_LLIST IN MOD DES_LINKED_LIST_DATA")

      RIGHT => PART%NEXT
      LEFT => PART%PREV
      DELETION = .false.

      IF(.NOT.ASSOCIATED(LEFT).AND.(.NOT.ASSOCIATED(RIGHT))) THEN
         !this is the only particle in the list
         DELETION  = .true.
      ELSEIF(.NOT.ASSOCIATED(LEFT)) THEN !LEFT BORDERING POINT, PREVIOUS NULL PTR
         !ONLY CORRECT THE PREV POINTER FOR RIGHT ENTRY
         RIGHT%PREV => NULL()
         DELETION = .TRUE.
      ELSEIF(.NOT.ASSOCIATED(RIGHT)) THEN !RIGHT BORDERING POINT, PREVIOUS NULL PTR
         !ONLY CORRECT THE NEXT POINTER FOR LEFT ENTRY
         LEFT%NEXT => NULL()
         DELETION = .TRUE.
      ELSE
         RIGHT%PREV => LEFT
         LEFT%NEXT => RIGHT
         DELETION = .TRUE.
      END IF
      IF (DELETION) THEN
         DEALLOCATE(PART)
      ENDIF
      END SUBROUTINE REMOVE_PART_LLIST

      END MODULE DES_LINKED_LIST_FUNCS
