      MODULE FUNITS

! Whether this processor should write the log file in DMP mode.
! Usually this flag is true only for PE_IO.  All the PEs may be forced
! to write a log file by setting ENABLE_DMP_LOG to .true. in output_mod.f.
      LOGICAL :: DMP_LOG

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

      END MODULE FUNITS
