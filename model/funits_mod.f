      MODULE funits

!           Whether this processor should write the log file in DMP mode.
!           Usually this flag is true only for PE_IO.  All the PEs may be forced
!           to write a log file by setting ENABLE_DMP_LOG to .true. in output_mod.f.
      LOGICAL :: DMP_LOG
!
!              .DAT unit number
      INTEGER  UNIT_DAT
      PARAMETER (UNIT_DAT = 51)
!
!              .OUT unit number
      INTEGER  UNIT_OUT
      PARAMETER (UNIT_OUT = 52)
!
!              .LOG unit number
      INTEGER  UNIT_LOG
      PARAMETER (UNIT_LOG = 53)
!
!              Temporary (scratch) file unit number
      INTEGER  UNIT_TMP
      PARAMETER (UNIT_TMP = 54)
!
!              .RES unit number
      INTEGER  UNIT_RES
      PARAMETER (UNIT_RES = 55)
!
!              .SPx unit offset number
      INTEGER  UNIT_SPX
      PARAMETER (UNIT_SPX = 60)


      END MODULE funits                                                                          
