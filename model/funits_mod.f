      MODULE FUNITS

! Whether this processor should write the log file in DMP mode.
! Usually this flag is true only for PE_IO.  All the PEs may be forced
! to write a log file by setting ENABLE_DMP_LOG to .true. in output_mod.f.
      LOGICAL :: DMP_LOG

! Flag indicated that the log was opened globally.
      LOGICAL :: LOG_WAS_CLOSED = .FALSE.

! RRATES debug file unit number
      INTEGER, PARAMETER :: UNIT_RRATES = 43

! mfix.dat file unit number
      INTEGER, PARAMETER :: UNIT_DAT = 51

! RUN_NAME.OUT file unit number
      INTEGER, PARAMETER :: UNIT_OUT = 52

! RUN_NAME.LOG file unit number. (DEFAULT/Serial 53)
      INTEGER, PARAMETER :: UNIT_LOG = 53

! Temporary (scratch) file unit number
      INTEGER, PARAMETER :: UNIT_TMP = 54

! RUN_NAME.RES file unit number
      INTEGER, PARAMETER :: UNIT_RES = 55

! RUN_NAME.SPx file unit offset number
      INTEGER, PARAMETER :: UNIT_SPX = 60

      CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: NEWUNIT                                                !
!  Purpose: Finds an open i/o unit number; Usage:                      !
!   integer myunit                                                     !
!   open(convert='big_endian',unit=newunit(myunit),file='filename')    !
!                                                                      !
!  Author: A. Choudhary                               Date: 01/21/2015 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      INTEGER FUNCTION newunit(unit)
      IMPLICIT NONE

! optional variable to hold unit number
      INTEGER, INTENT(OUT), OPTIONAL  :: unit

! lower and upper limits to search for available units
      INTEGER, PARAMETER  :: lun_min = 100, lun_max= 999

! check to see if the unit is open
      LOGICAL             :: is_open

! looping variable
      INTEGER             :: lun


      newunit = -1

      DO lun = lun_min, lun_max
        INQUIRE(UNIT=lun, OPENED=is_open)
        IF(.NOT.is_open) THEN
          newunit = lun
          EXIT
        END IF
      END DO

      IF(present(unit)) unit=newunit

      RETURN

      END FUNCTION NEWUNIT


      END MODULE FUNITS
