!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: OUTPUT_MANAGER                                          !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Purpose: Relocate calls to write output files (RES, SPx, VTP). This !
!  was done to simplify the time_march code.                           !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE OUTPUT_MANAGER(BATCHQ_END, FINISHED)

! Global Variables:
!---------------------------------------------------------------------//

      use output, only: RES_TIME, RES_DT
      use output, only: SPX_TIME, SPX_DT
      use output, only: OUT_TIME, OUT_DT
      use output, only: USR_TIME, USR_DT
      use vtk, only:    VTK_TIME, VTK_DT

      use output, only: DISK, DISK_TOT

      use param1, only: N_SPX
      use param, only: DIMENSION_USR
      use vtk, only: DIMENSION_VTK

      use run, only: TIME, DT, TSTOP

      use time_cpu, only: CPU_IO

      use compar, only: myPE, PE_IO

      use discretelement, only: DISCRETE_ELEMENT

      use vtk, only: WRITE_VTK_FILES
      use qmom_kinetic_equation, only: QMOMK

      use param1, only: UNDEFINED

      use vtp, only: write_vtp_file


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag that the the user specified batch time (plus buffer) is met.
      LOGICAL, INTENT(IN) :: BATCHQ_END
! Flag that a steady state case is completed.
      LOGICAL, INTENT(IN) :: FINISHED

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter and counter
      INTEGER :: LC, IDX
! Flag to write NetCDF output
      LOGICAL :: bWRITE_NETCDF_FILES
! Flag that the header (time) has not be written.
      LOGICAL :: HDR_MSG
! SPX file extensions.
      CHARACTER(LEN=35) ::  EXT_END
! Wall time at the start of IO operations.
      DOUBLE PRECISION :: WALL_START

! External function:
!---------------------------------------------------------------------//
! Returns the current wall time.
      DOUBLE PRECISION :: WALL_TIME
!......................................................................!


! Initialize the SPx file extension array.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
! Initial the header flag.
      HDR_MSG = .TRUE.

! Get the current time before any IO operations begin
      WALL_START = WALL_TIME()

! Write restart file, if needed
      IF(CHECK_TIME(RES_TIME) .OR. BATCHQ_END) THEN

         RES_TIME = NEXT_TIME(RES_DT)
         CALL WRITE_RES1
         CALL NOTIFY_USER('.RES;')

         IF(DISCRETE_ELEMENT) THEN
            CALL WRITE_RES0_DES
            CALL NOTIFY_USER('DES.RES;')
         ENDIF

         IF(QMOMK) THEN
            CALL QMOMK_WRITE_RESTART
            CALL NOTIFY_USER('QMOMK.RES;')
         ENDIF

      ENDIF

! Write SPx files, if needed
      IDX = 0
      bWRITE_NETCDF_FILES = .FALSE.

      DO LC=1, N_SPX
         IF(CHECK_TIME(SPX_TIME(LC))) THEN
            SPX_TIME(LC) = NEXT_TIME(SPX_DT(LC))

            CALL WRITE_SPX1(LC, 0)
            CALL NOTIFY_USER('SPx:',EXT_END(LC:LC))

            DISK_TOT = DISK_TOT + DISK(LC)
            IDX = IDX + 1

            bWRITE_NETCDF_FILES = .TRUE.
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST


! Write standard output, if needed
      IF(CHECK_TIME(OUT_TIME)) THEN
         OUT_TIME = NEXT_TIME(OUT_DT)
         CALL WRITE_OUT1
         CALL NOTIFY_USER('.OUT;')
      ENDIF

! Write special output, if needed
      IDX = 0
      DO LC = 1, DIMENSION_USR
         IF(CHECK_TIME(USR_TIME(LC))) THEN
            USR_TIME(LC) = NEXT_TIME(USR_DT(LC))
            CALL WRITE_USR1 (LC)
            CALL NOTIFY_USER('.USR:',EXT_END(LC:LC))
            IDX = IDX + 1
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST

      CALL FLUSH_NOTIFY_USER

! Write vtk file, if needed
      IF(WRITE_VTK_FILES) THEN
         DO LC = 1, DIMENSION_VTK
            IF(CHECK_TIME(VTK_TIME(LC))) THEN
               VTK_TIME(LC) = NEXT_TIME(VTK_DT(LC))
               CALL WRITE_VTU_FILE(LC)
               IF(DISCRETE_ELEMENT) CALL WRITE_VTP_FILE(LC)
            ENDIF
         ENDDO
      ENDIF

! Write NetCDF files.
      IF(bWRITE_NETCDF_FILES) CALL WRITE_NETCDF(0,0,TIME)

! Add the amount of time needed for all IO operations to total.
      CPU_IO = CPU_IO + (WALL_TIME() - WALL_START)

      RETURN

      contains


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION CHECK_TIME(lTIME)

      DOUBLE PRECISION, INTENT(IN) :: lTIME

      IF(DT == UNDEFINED) THEN
         CHECK_TIME = FINISHED
      ELSE
         CHECK_TIME = (TIME+0.1d0*DT>=lTIME).OR.(TIME+0.1d0*DT>=TSTOP)
      ENDIF

      RETURN
      END FUNCTION CHECK_TIME


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      DOUBLE PRECISION FUNCTION NEXT_TIME(lWRITE_DT)

      DOUBLE PRECISION, INTENT(IN) :: lWRITE_DT

      IF (DT /= UNDEFINED) THEN
         NEXT_TIME = (INT((TIME + 0.1d0*DT)/lWRITE_DT)+1)*lWRITE_DT
      ELSE
         NEXT_TIME = lWRITE_DT
      ENDIF

      RETURN
      END FUNCTION NEXT_TIME


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE NOTIFY_USER(MSG, EXT)

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      CHARACTER(len=*), INTENT(IN) :: MSG
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: EXT


      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1000, ADVANCE='NO') TIME
         IF(SCR_LOG) WRITE(*, 1000, ADVANCE='NO') TIME
         HDR_MSG = .FALSE.
      ENDIF

 1000 FORMAT(' ',/' t=',F12.6,' Wrote')

      IF(.NOT.present(EXT)) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1100, ADVANCE='NO') MSG
         IF(SCR_LOG) WRITE(*, 1100, ADVANCE='NO') MSG
      ELSE
         IF(IDX == 0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1110, ADVANCE='NO') MSG, EXT
            IF(SCR_LOG) WRITE(*, 1110, ADVANCE='NO') MSG, EXT
         ELSE
            IF(DMP_LOG) WRITE(UNIT_LOG, 1120, ADVANCE='NO') EXT
            IF(SCR_LOG) WRITE(*, 1120, ADVANCE='NO') EXT
         ENDIF
      ENDIF

 1100 FORMAT(1X,A)
 1110 FORMAT(1X,A,1x,A)
 1120 FORMAT(',',A)

      RETURN
      END SUBROUTINE NOTIFY_USER

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_LIST

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(DMP_LOG) WRITE(UNIT_LOG,1000, ADVANCE='NO')
      IF(SCR_LOG) WRITE(*,1000, ADVANCE='NO')

 1000 FORMAT(';')

      RETURN
      END SUBROUTINE FLUSH_LIST


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_NOTIFY_USER

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START

      use output, only: NLOG

      use run, only: NSTEP

      use error_manager


      DOUBLE PRECISION :: WALL_ELAP, WALL_LEFT, WALL_NOW
      CHARACTER(LEN=9) :: CHAR_ELAP, CHAR_LEFT
      CHARACTER(LEN=4) :: UNIT_ELAP, UNIT_LEFT

      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(.NOT.HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         IF(SCR_LOG) WRITE(*,1000)
      ENDIF

 1000 FORMAT(' ',/' ')

! Write the elapsed time and estimated remaining time
      IF(MOD(NSTEP,NLOG) == 0) THEN
         WALL_NOW = WALL_TIME()
! Calculate the elapsed wall time.
         WALL_ELAP = WALL_NOW - WALL_START
         CALL GET_TUNIT(WALL_ELAP, UNIT_ELAP)
         CHAR_ELAP=''; WRITE(CHAR_ELAP,"(F9.2)") WALL_ELAP
         CHAR_ELAP = trim(adjustl(CHAR_ELAP))
! Estimate the remaining wall time.
         WALL_LEFT = (WALL_NOW-WALL_START)*(TSTOP-TIME)/               &
            max(TIME-TIME_START,1.0d-6)
         CALL GET_TUNIT(WALL_LEFT, UNIT_LEFT)

         IF (DT /= UNDEFINED) THEN
            CHAR_LEFT=''; WRITE(CHAR_LEFT,"(F9.2)") WALL_LEFT
            CHAR_LEFT = trim(adjustl(CHAR_LEFT))
         ELSE
            CHAR_LEFT = '0.0'
            UNIT_LEFT = 's'
         ENDIF

! Notify the user of usage/remaining wall times.
         WRITE(ERR_MSG,2000)                                           &
            'Elapsed:', trim(CHAR_ELAP), trim(UNIT_ELAP),              &
            'Est. Remaining:',trim(CHAR_LEFT), trim(UNIT_LEFT)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 2000 FORMAT('Wall Time - ',2(A,1X,A,A,4X))



      RETURN
      END SUBROUTINE FLUSH_NOTIFY_USER

      END SUBROUTINE OUTPUT_MANAGER
!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INIT_OUTPUT_VARS

      use output, only: RES_TIME, RES_DT
      use output, only: SPX_TIME, SPX_DT
      use output, only: OUT_TIME, OUT_DT
      use output, only: USR_TIME, USR_DT
      use vtk, only:    VTK_TIME, VTK_DT

      use output, only: DISK, DISK_TOT

      use param1, only: N_SPX
      use param, only: DIMENSION_USR
      use vtk, only: DIMENSION_VTK


      use physprop, only: MMAX, NMAX
      use run, only: RUN_TYPE
      use run, only: K_EPSILON
      use param1, only: ZERO
      use geometry, only: IJKMAX2

      use vtk, only: DIMENSION_VTK
      use vtk, only: WRITE_VTK_FILES
      use vtk, only: VTK_TIME, VTK_DT

      use param1, only: UNDEFINED
      use run, only: TIME, DT

      use time_cpu, only: CPU_IO
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START
      use rxns, only: nRR
      use scalars, only: NScalar
      use output, only: ONEMEG

      IMPLICIT NONE

! Disk space needed for one variable and each SPX file
      DOUBLE PRECISION :: DISK_ONE

! External function for geting job wall time
      DOUBLE PRECISION :: WALL_TIME

! Loop counter
      INTEGER :: LC

! Initialize times for writing outputs
      OUT_TIME = TIME
! Initialize the amount of time spent on IO
      CPU_IO = 0.0d0

! Initialize disk space calculations
      DISK_TOT = ZERO
      DISK_ONE = 4.0*IJKMAX2/ONEMEG

      DISK(1) = 1.0*DISK_ONE                           ! EPg
      DISK(2) = 2.0*DISK_ONE                           ! Pg, Ps
      DISK(3) = 3.0*DISK_ONE                           ! Ug, Vg, Wg
      DISK(4) = 3.0*DISK_ONE*MMAX                      ! Us, Vs, Ws
      DISK(5) = 1.0*DISK_ONE*MMAX                      ! ROPs
      DISK(6) = 1.0*DISK_ONE*(MMAX+1)                  ! Tg, Ts
      DISK(7) = 1.0*DISK_ONE*(sum(NMAX(0:MMAX)))       ! Xg, Xs
      DISK(8) = 1.0*DISK_ONE*MMAX                      ! Theta
      DISK(9) = 1.0*DISK_ONE*NScalar                   ! User Scalars
      DISK(10) = nRR*DISK_ONE                          ! ReactionRates
      DISK(11) = merge(2.0*DISK_ONE, ZERO, K_EPSILON)  ! K-Epsilon


      IF (RUN_TYPE == 'NEW') THEN
         RES_TIME = TIME
         SPX_TIME(:N_SPX) = TIME
      ELSE
         IF (DT /= UNDEFINED) THEN
            RES_TIME = RES_DT *                                        &
               (INT((TIME + 0.1d0*DT)/RES_DT) + 1)
            SPX_TIME(:N_SPX) = SPX_DT(:N_SPX) *                        &
               (INT((TIME + 0.1d0*DT)/SPX_DT(:N_SPX)) + 1)
         ENDIF
      ENDIF

      DO LC = 1, DIMENSION_USR
         USR_TIME(LC) = UNDEFINED
         IF (USR_DT(LC) /= UNDEFINED) THEN
            IF (RUN_TYPE == 'NEW') THEN
               USR_TIME(LC) = TIME
            ELSE
               USR_TIME(LC) = USR_DT(LC) *                             &
                  (INT((TIME+0.1d0*DT)/USR_DT(LC))+1)
            ENDIF
         ENDIF
      ENDDO

! Initialize VTK_TIME
      IF(WRITE_VTK_FILES) THEN
         DO LC = 1, DIMENSION_VTK
            VTK_TIME(LC) = UNDEFINED
            IF (VTK_DT(LC) /= UNDEFINED) THEN
               IF (RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2') THEN
                  VTK_TIME(LC) = TIME
               ELSE
                  VTK_TIME(LC) = VTK_DT(LC) *                          &
                     (INT((TIME + 0.1d0*DT)/VTK_DT(LC))+1)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      WALL_START = WALL_TIME()
      TIME_START = TIME

      RETURN
      END SUBROUTINE INIT_OUTPUT_VARS

