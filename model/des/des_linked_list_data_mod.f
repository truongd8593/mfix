!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: DES_LINKED_LIST_DATA_MOD                               !
!                                                                      !
!                                                                      !
!  Reviewer: R. Garg                                  Date: 19-Mar-14  !
!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_LINKED_LIST_DATA
      use error_manager

      IMPLICIT NONE

      TYPE PARTICLE
      INTEGER :: CELL(4)
      LOGICAL :: INDOMAIN
! Solid phase
      INTEGER :: M
      double precision :: RAD, DENS, STATWT
      double precision :: VELOCITY(3), POSITION(3)

      !could be made allocatable later to reduce memory usage for 2-D runs
      TYPE (PARTICLE), POINTER :: NEXT=>NULL() ! NEXT PARTICLE ADDRESS
      TYPE (PARTICLE), POINTER :: PREV=>NULL() ! PREVIOUS PARTICLE ADDRESS

      END TYPE PARTICLE

! This is the linked list of first set of particles
      TYPE (PARTICLE), POINTER :: ORIG_PART_LIST => NULL()
! This is the linked list of deleted set of particles
      TYPE (PARTICLE), POINTER :: DEL_PART_LIST => NULL()
! This is the linked list of remaining set of particles
      TYPE (PARTICLE), POINTER :: REM_PART_LIST => NULL()

      END MODULE DES_LINKED_LIST_DATA


