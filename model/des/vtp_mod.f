      MODULE vtp

      use mpi_utility
      use cdist

      use desmpi
      use mpi_comm_des
      use error_manager


      IMPLICIT NONE

      INTEGER, PRIVATE :: GLOBAL_CNT
      INTEGER, PRIVATE :: LOCAL_CNT

      INTEGER :: DES_UNIT = 2000

! file unit for ParaView *.pvd data
      INTEGER, PARAMETER :: PVD_UNIT = 2050

! formatted file name
      CHARACTER(LEN=64) :: FNAME_VTP

      INTERFACE VTP_WRITE_DATA
         MODULE PROCEDURE VTP_WRITE_DP1
         MODULE PROCEDURE VTP_WRITE_DP2
         MODULE PROCEDURE VTP_WRITE_I1
      END INTERFACE


      CONTAINS


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP1                                            !
!                                                                      !
! Purpose: Collect and write 1D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(.NOT.PEA(LC,1)) CYCLE
            PC = PC+1
            IF(PEA(LC,4)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") real(DATA(LC))
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (dProcBuf(LOCAL_CNT) )
         allocate (dRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") real(drootbuf(LC))
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(dProcBuf, dRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP1


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP2                                            !
!                                                                      !
! Purpose: Collect and write 2D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP2(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)

      CHARACTER(len=16) :: NOC
      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2

      LB = LBOUND(DATA,1)
      UB = UBOUND(DATA,1)
      NOC=''; WRITE(NOC,*) (UB-LB)+1

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))

         PC = 1
         DO LC1 = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(.NOT.PEA(LC1,1)) CYCLE
            PC = PC+1
            IF(PEA(LC1,4)) CYCLE
            DO LC2=LB, UB
               WRITE(DES_UNIT,1001,ADVANCE="NO") real(DATA(LC1,LC2))
            ENDDO
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (dProcBuf(LOCAL_CNT) )
         allocate (dRootBuf(GLOBAL_CNT))
         allocate (ltemp_array((UB-LB)+1,GLOBAL_CNT))

         DO LC1 = LB, UB
            CALL DES_GATHER(DATA(LC1,:))
            ltemp_array(LC1,:) = drootbuf(:)
         ENDDO

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
            DO LC1=1, GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(DES_UNIT,1001,ADVANCE="NO") &
                     real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array)

      ENDIF


 1000 FORMAT('<DataArray type="Float32" Name="',A,'" NumberOf',        &
         'Components="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP2



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_I1                                             !
!                                                                      !
! Purpose: Collect and write 1D integer arrays to the VTP file. This   !
! routine is designed to collect the data for parallel and serial      !
! runs. This routine also manages the distribted IO case.              !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_I1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      INTEGER, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(.NOT.PEA(LC,1)) CYCLE
            PC = PC+1
            IF(PEA(LC,4)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") DATA(LC)
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (iProcBuf(LOCAL_CNT) )
         allocate (iRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") irootbuf(LC)
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(iProcBuf, iRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(I10,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_I1


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_ELEMENT                                        !
!                                                                      !
! Purpose: Write a string to the VTP file. It masks the need to check  !
! the logical before flushing.                                         !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_ELEMENT(ELEMENT)

      CHARACTER(len=*), INTENT(in) :: ELEMENT

      IF(bDist_IO .OR. myPE == PE_IO) &
         WRITE(DES_UNIT,"(A)") ELEMENT

      RETURN
      END SUBROUTINE VTP_WRITE_ELEMENT



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_OPEN_FILE                                            !
!                                                                      !
! Purpose: This routine opens the VTP file and calcualtes the offsets  !
! for dmp data collection.                                             !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_OPEN_FILE(NoPc)

! Modules
!-----------------------------------------------
      USE run, only: RUN_TYPE, RUN_NAME

      IMPLICIT NONE

      CHARACTER(len=*) :: NoPc

      INTEGER :: NumberOfPoints

! Variables related to gather
      integer lgathercnts(0:numpes-1), lproc

! check whether an error occurs in opening a file
      INTEGER :: IOS
! Integer error flag.
      INTEGER :: IER

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_VTP
! status of the vtp file to be written
      CHARACTER(LEN=8) :: STATUS_VTP


! Initial the global count.
      GLOBAL_CNT = 10
! Calculate the number of 'real' particles on the local process.
      LOCAL_CNT = PIP - iGHOST_CNT

! Distributed IO
      IF(bDIST_IO) THEN
         NumberOfPoints = LOCAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

         WRITE(fname_vtp,'(A,"_DES",I4.4,"_",I5.5,".vtp")') &
            trim(run_name), vtp_findex, mype

! Serial IO
      ELSE

! Calculate the total number of particles system-wide.
         call global_sum(LOCAL_CNT, GLOBAL_CNT)
         NumberOfPoints = GLOBAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

! Set the send count from the local process.
         igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
         lgathercnts = 0
         lgathercnts(myPE) = LOCAL_CNT
         call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
         idispls(0) = 0
         DO lPROC = 1,NUMPEs-1
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
         ENDDO

! set the file name and unit number and open file
         WRITE(fname_vtp,'(A,"_DES_",I5.5,".vtp")') &
            trim(run_name), vtp_findex
      ENDIF

      IER = 0
      IF(bDIST_IO .OR. myPE == PE_IO) THEN

! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_VTP,EXIST=EXISTS_VTP)
! The given file should not exist if the run type is NEW.
         IF(RUN_TYPE == 'NEW' .AND. EXISTS_VTP)THEN
            IER = 1

! A RESTART_1 case may over-write vtp files created during past run
         ELSEIF(RUN_TYPE == 'RESTART_1' .AND. EXISTS_VTP) THEN
            STATUS_VTP = 'REPLACE'

! Set the status to NEW so that any problems in opening the file
! will be accurately reported.
         ELSE
            STATUS_VTP = 'NEW'
         ENDIF

! Open the file and record any erros.
         OPEN(UNIT=DES_UNIT, FILE=FNAME_VTP, STATUS=STATUS_VTP, IOSTAT=IOS)
         IF(IOS /= 0) IER = 1
      ENDIF

      CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         CALL INIT_ERR_MSG("VTP_MOD --> OPEN_VTP")
         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Unable to open VTP file. This could be ',    &
         'causes by a VTP',/'with the same file name already ',        &
         'existing. or an error code returned',/'from the OPEN ',      &
         'function'/'Aborting.')


      END SUBROUTINE VTP_OPEN_FILE



!......................................................................!
! SUBROUTINE: VTP_CLOSE_FILE                                           !
!                                                                      !
! Purpose: This routine closes the vtp file.                           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_CLOSE_FILE


      VTP_FINDEX=VTP_FINDEX+1

      IF(bDist_io .OR. (myPE .eq.pe_IO)) CLOSE(des_unit)


      END SUBROUTINE VTP_CLOSE_FILE


!......................................................................!
! SUBROUTINE: ADD_VTP_TO_PVD                                           !
!                                                                      !
! Purpose: This routine opens the pvd file.                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE ADD_VTP_TO_PVD

      use run, only: RUN_TYPE, RUN_NAME

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index position of desired character
      INTEGER IDX_f, IDX_b
! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_PVD
! Generic input limited to 256 characters
      CHARACTER(LEN=256) INPUT

! formatted file name
      CHARACTER(LEN=64) :: FNAME_PVD = ''

! formatted solids time
      CHARACTER(LEN=12) :: S_TIME_CHAR = ''

      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! IO Status flag
      INTEGER :: IOS

! Variables related to gather
      integer :: IER

!-----------------------------------------------


      CALL INIT_ERR_MSG('VTP_MOD --> ADD_VTP_TO_PVD')

! Initialize the error flag.
      IER = 0

! Obtain the file name and open the pvd file
      FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'

! The PVD file is only written by PE_IO with serial IO.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_PVD,EXIST=EXISTS_PVD)

         IF(FIRST_PASS) THEN

! Open the "NEW" file and write the necessary header information.
            IF(RUN_TYPE /= 'RESTART_1')THEN

! The file exists but first_pass is also true so most likely an existing
! file from an earlier/other run is present in the directory. Exit to
! prevent accidently overwriting the existing file.
               IF(EXISTS_PVD) THEN
                  IER = 1
               ELSE
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,STATUS='NEW')
                  WRITE(PVD_UNIT,"(A)")'<?xml version="1.0"?>'
                  WRITE(PVD_UNIT,"(A)")'<VTKFile type="Collection" &
                     &version="0.1" byte_order="LittleEndian">'
                  WRITE(PVD_UNIT,"(3X,'<Collection>')")
! write two generic lines that will be removed later
                  WRITE(PVD_UNIT,"('SCRAP LINE 1')")
                  WRITE(PVD_UNIT,"('SCRAP LINE 2')")
               ENDIF

! This is the first pass of a restart run. Extra care is needed to make
! sure that the pvd file is ready to accept new data.
            ELSE ! a restart run
               IF(EXISTS_PVD) THEN
! Open the file at the beginning.
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
                     POSITION="REWIND",STATUS='OLD',IOSTAT=IOS)
                  IF(IOS /= 0) IER = 2
               ELSE ! a pvd file does not exist
                  IER = 3
               ENDIF

               IF(IER == 0) THEN
! Loop over the entries in the PVD file, looking for a match to the
! file that is being written. If no match is found, the data will be
! appended to the end of the pvd file, otherwise, the old data will
! be over-written.
                  DO
! Read in the entires of the PVD file.
                     READ(PVD_UNIT,"(A)",IOSTAT=IOS)INPUT
                     IF(IOS > 0) THEN
                        IER = 4
                        EXIT
                     ELSEIF(IOS<0)THEN
! The end of the pvd file has been reached without finding an entry
! matching the current record. Exit the loop.
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
! Find the first instances of file=" and "/> in the read data.
                     IDX_f = INDEX(INPUT,'file="')
                     IDX_b = INDEX(INPUT,'"/>')
! Skip rows that do not contain file data
                     IF(IDX_f == 0 .AND. IDX_b == 0) CYCLE
! Truncate the file name from the read data
                     WRITE (INPUT,"(A)") INPUT(IDX_f+6:IDX_b-1)
! If the file name matches the current VTP record, break the loop to
! over-write this record.
                     IF(TRIM(FNAME_VTP) == TRIM(INPUT)) THEN
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF ! No errors
            ENDIF ! run_type new or restart
! Identify that the files has been created and opened for next pass
            FIRST_PASS = .FALSE.

         ELSE ! not FIRST_PASS
            OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
               POSITION="APPEND",STATUS='OLD',IOSTAT=IOS)
            IF (IOS /= 0) IER = 2
         ENDIF

      ENDIF ! if myPE == PE_IO and not distributed IO


      CAlL GLOBAL_ALL_SUM(IER)
      IF(IER /= 0) THEN
         SELECT CASE(IER)
         CASE(1); WRITE(ERR_MSG,1101) trim(FNAME_PVD)
         CASE(2); WRITE(ERR_MSG,1102) trim(FNAME_PVD)
         CASE(3); WRITE(ERR_MSG,1103) trim(FNAME_PVD)
         CASE(4); WRITE(ERR_MSG,1104) trim(FNAME_PVD)
         CASE DEFAULT; WRITE(ERR_MSG,1105) trim(FNAME_PVD)
         END SELECT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: A PVD file was detected in the run ',        &
         'directory which should',/'not exist for a NEW run.',/        &
         'File: ',A)

 1102 FORMAT('Error 1102: Fatal error status returned while OPENING ', &
         'PVD file.',/'File: ', A)

 1103 FORMAT('Error 1103: PVD file MISSING from run directory.',/      &
         'File: ',A)

 1104 FORMAT('Error 1104: Fatal error status returned while READING ', &
         'PVD file.',/'File: ', A)

 1105 FORMAT('Error 1105:: Fatal unclassified error when processing ', &
         'PVD file.',/'File: ', A)


! If there were no errors, updated the file.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Remove the last two lines written so that additional data can be added
         BACKSPACE(PVD_UNIT)
         BACKSPACE(PVD_UNIT)

! Write the data to the file
         WRITE(PVD_UNIT,"(6X,A,A,A,A,A,A,A)")&
         '<DataSet timestep="',TRIM(iVal(S_TIME)),'" ',&
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,"(3X,A)")'</Collection>'
         WRITE(PVD_UNIT,"(A)")'</VTKFile>'

         CLOSE(PVD_UNIT)
      ENDIF

      CALL FINL_ERR_MSG

! Return to the calling routine
      RETURN

      END SUBROUTINE ADD_VTP_TO_PVD

      END MODULE VTP
