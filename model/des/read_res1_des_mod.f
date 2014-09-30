      MODULE READ_RES1_DES

      use desmpi
      use compar, only: myPE
      use compar, only: PE_IO

      use cdist, only: bDist_IO

      use error_manager

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_READ_RES_DES
      PUBLIC :: FINL_READ_RES_DES

      PUBLIC :: READ_PAR_POS
      PUBLIC :: READ_PAR_COL

      PUBLIC :: READ_RES_DES
      PUBLIC :: READ_RES_pARRAY
      PUBLIC :: READ_RES_cARRAY


      INTERFACE READ_RES_DES
         MODULE PROCEDURE READ_RES_DES_0I
         MODULE PROCEDURE READ_RES_DES_1I
         MODULE PROCEDURE READ_RES_DES_0D
         MODULE PROCEDURE READ_RES_DES_1D
         MODULE PROCEDURE READ_RES_DES_0L
         MODULE PROCEDURE READ_RES_DES_1L
      END INTERFACE

      INTERFACE READ_RES_pARRAY
         MODULE PROCEDURE READ_RES_pARRAY_1I
         MODULE PROCEDURE READ_RES_pARRAY_1D
         MODULE PROCEDURE READ_RES_pARRAY_1L
      END INTERFACE

      INTERFACE READ_RES_cARRAY
         MODULE PROCEDURE READ_RES_cARRAY_1I
         MODULE PROCEDURE READ_RES_cARRAY_1D
         MODULE PROCEDURE READ_RES_cARRAY_1L
      END INTERFACE


      INTEGER, PARAMETER :: RDES_UNIT = 901

      INTEGER :: pIN_COUNT
      INTEGER :: cIN_COUNT

! Send/Recv parameters for Particle arrays:
      INTEGER :: pROOTCNT, pPROCCNT
      INTEGER :: pRECV
      INTEGER, allocatable :: pSCATTER(:)
      INTEGER, allocatable :: pDISPLS(:)

! Variables used for reading restart file
      INTEGER, ALLOCATABLE :: pRestartMap(:)
      INTEGER, ALLOCATABLE :: cRestartMap(:)

! Send/Recv parameters for Particle arrays:
      INTEGER :: cROOTCNT, cPROCCNT
      INTEGER :: cRECV
      INTEGER, allocatable :: cSCATTER(:)
      INTEGER, allocatable :: cDISPLS(:)

      INTEGER, ALLOCATABLE :: iPAR_COL(:,:)

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_READ_RES_DES                                        !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_READ_RES_DES(BASE, lVERSION, lNEXT_REC)

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      DOUBLE PRECISION, INTENT(OUT) :: lVERSION
      INTEGER, INTENT(OUT) :: lNEXT_REC

      CHARACTER(len=32) :: lFNAME

      INTEGER :: lDIMN


      lDIMN = merge(2,3,NO_K)

      allocate(pSCATTER(0:numPEs-1))
      allocate(pDISPLS(0:numPEs-1))

      allocate(cSCATTER(0:numPEs-1))
      allocate(cDISPLS(0:numPEs-1))


      IF(bDIST_IO) THEN

         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

         READ(RDES_UNIT, REC=1) lVERSION
         READ(RDES_UNIT, REC=2) pIN_COUNT
         READ(RDES_UNIT, REC=3) iGHOST_CNT
         READ(RDES_UNIT, REC=4) cIN_COUNT

         IF(PIP > MAX_PIP) THEN
            write(*,*) "From des_read_restart:"
            write(*,*) "Error: The pip is grater than current max_pip"
            write(*,*) "pip=" ,pip,"; max_pip =", max_pip

         ENDIF

         PIP = pIN_COUNT
         COLLISION_NUM = cIN_COUNT

         DO WHILE(COLLISION_NUM > COLLISION_MAX)
            COLLISION_MAX = 2*COLLISION_MAX
         ENDDO
         CALL COLLISION_GROW

      ELSE

         IF(myPE == PE_IO) THEN
            WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
            OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',      &
               STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

            READ(RDES_UNIT, REC=1) pIN_COUNT

            READ(RDES_UNIT, REC=1) lVERSION
            READ(RDES_UNIT, REC=2) pIN_COUNT
!           READ(RDES_UNIT, REC=3) -NOTHING-
            READ(RDES_UNIT, REC=4) cIN_COUNT

         ELSE
            pIN_COUNT = 10
            cIN_COUNT = 10
         ENDIF

         allocate( pRestartMap(pIN_COUNT))
         allocate( cRestartMap(cIN_COUNT))
      ENDIF

      lNEXT_REC = 5

      RETURN
      END SUBROUTINE INIT_READ_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_READ_RES_DES

      IF(bDIST_IO .OR. myPE == PE_IO) close(RDES_UNIT)

      IF(allocated(dPROCBUF)) deallocate(dPROCBUF)
      IF(allocated(dROOTBUF)) deallocate(dROOTBUF)
      IF(allocated(iPROCBUF)) deallocate(iPROCBUF)
      IF(allocated(iROOTBUF)) deallocate(iROOTBUF)

      IF(allocated(pRestartMap)) deallocate(pRestartMap)
      IF(allocated(cRestartMap)) deallocate(cRestartMap)

      IF(allocated(pSCATTER)) deallocate(pSCATTER)
      IF(allocated(pDISPLS)) deallocate(pDISPLS)

      IF(allocated(cSCATTER)) deallocate(cSCATTER)
      IF(allocated(cDISPLS)) deallocate(cDISPLS)



      IF(bDIST_IO) CALL DES_RESTART_NEIGH

      RETURN
      END SUBROUTINE FINL_READ_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_PAR_POS                                             !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_PAR_POS(lNEXT_REC)

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC

      INTEGER :: lDIMN
      INTEGER :: LC1, lPROC
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The number of particles on each process.
      INTEGER :: PAR_CNT(0:NUMPEs-1)

!-----------------------------------------------

      CALL INIT_ERR_MSG("READ_PAR_POS")

      lDIMN = merge(2,3,NO_K)


! All process read positions for distributed IO restarts.
      IF(bDIST_IO) THEN
         DO LC1 = 1, lDIMN
            CALL READ_RES_DES(lNEXT_REC, DES_POS_NEW(LC1,:))
         ENDDO
         RETURN
      ENDIF

      allocate( dPAR_POS(lDIMN, pIN_COUNT))

! Only the IO proccess reads positions.
      IF(myPE == PE_IO) THEN
         DO LC1=1, merge(2,3,NO_K)
            CALL IN_BIN_512(RDES_UNIT, dPAR_POS(LC1,:),                &
               pIN_COUNT, lNEXT_REC)
         ENDDO
      ENDIF

! Use the particle postions and the domain coverage of each process
! to determine which processor each particle belongs.
      CALL MAP_pARRAY_TO_PROC(PAR_CNT)

! Send the particle position data to the individual ranks.
      CALL SCATTER_PAR_POS(PAR_CNT)

! Set up the read/scatter arrary information.
      pPROCCNT = PIP
      pROOTCNT = pIN_COUNT

! Set the recv count for this process.
      pRECV = PIP

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      lScatterCnts(:) = 0; lScatterCnts(mype) = PIP
      CALL GLOBAL_SUM(lScatterCnts,pSCATTER)

! Calculate the displacements for each process in the global array.
      pDispls(0) = 0
      DO lPROC = 1, NUMPEs-1
         pDispls(lPROC) = pDispls(lPROC-1) + pSCATTER(lPROC-1)
      ENDDO

      IF(allocated(dPAR_POS)) deallocate(dPAR_POS)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE READ_PAR_POS


!``````````````````````````````````````````````````````````````````````!
! Subroutine: MAP_pARRAY_TO_PROC                                       !
!                                                                      !
! Purpose: Use the particle positions to determine which processor     !
! they live on and count the number of particles on each process.      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE MAP_pARRAY_TO_PROC(lPAR_CNT)

      implicit none

      INTEGER, INTENT(OUT) :: lPAR_CNT(0:numPEs-1)

! Data dimensionality flag.
      INTEGER :: lDIMN
! Loop counters.
      INTEGER :: LC1, lPROC
! Error flag.
      INTEGER :: IER(0:numPEs-1)
! Local scatter count.
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The X/Y/Z bounds of the physical space "owned" by each process.
      DOUBLE PRECISION :: lxmin(0:NUMPEs-1), lxmax(0:NUMPEs-1)
      DOUBLE PRECISION :: lymin(0:NUMPEs-1), lymax(0:NUMPEs-1)
      DOUBLE PRECISION :: lzmin(0:NUMPEs-1), lzmax(0:NUMPEs-1)
!-----------------------------------------------

      CALL INIT_ERR_MSG("MAP_pARRAY_TO_PROC")

! Initialize the error flag.
      IER = 0

      lDIMN = merge(2, 3, NO_K)

! set the domain range for each processor
      DO lPROC= 0, NUMPEs-1
         lxmin(lproc) = xe(istart1_all(lproc)-1)
         lxmax(lproc) = xe(iend1_all(lproc))
         lymin(lproc) = yn(jstart1_all(lproc)-1)
         lymax(lproc) = yn(jend1_all(lproc))
         lzmin(lproc) = zt(kstart1_all(lproc)-1)
         lzmax(lproc) = zt(kend1_all(lproc))

! modify the range for mass inlet and outlet, as particles injected
! can lie outside the domain and not ghost particles
         IF(istart1_all(lproc).eq.imin1) &
            lxmin(lproc) = xe(istart1_all(lproc)-2)
         IF(iend1_all(lproc).eq.imax1) &
            lxmax(lproc) = xe(iend1_all(lproc)+1)
         IF(jstart1_all(lproc).eq.jmin1) &
            lymin(lproc) = yn(jstart1_all(lproc)-2)
         IF(jend1_all(lproc).eq.jmax1)  &
            lymax(lproc) = yn(jend1_all(lproc)+1)
         IF(kstart1_all(lproc).eq.kmin1 .AND. DO_K) &
            lzmin(lproc) = zt(kstart1_all(lproc)-2)
         IF(kend1_all(lproc).eq.kmax1 .AND. DO_K) &
            lzmax(lproc) = zt(kend1_all(lproc)+1)
      ENDDO

! build the send buffer in PE_IO proc
! first pass to get the count of particles
      IER = 0
      pRestartMap(:) = -1
      lPAR_CNT(:) = 0
      IF(myPE == PE_IO) THEN
         DO LC1 = 1, pIN_COUNT
            DO lPROC=0, NUMPEs-1
               IF(dPAR_POS(1,LC1) >= lxmin(lproc) .AND. &
                  dPAR_POS(1,LC1) <  lxmax(lproc) .AND. &
                  dPAR_POS(2,LC1) >= lymin(lproc) .AND. &
                  dPAR_POS(2,LC1) <  lymax(lproc)) THEN
                  IF(NO_K)THEN
                     lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                     pRestartMap(LC1) = lproc
                     EXIT
                  ELSE
                     IF(dPAR_POS(3,LC1) >= lzmin(lproc) .AND. &
                        dPAR_POS(3,LC1) <  lzmax(lproc)) THEN
                        lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                        pRestartMap(LC1) = lproc
                        EXIT
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO  ! Loop over processes
            IF (pRestartMap(LC1) == -1) then
               IER(myPE) = -1
               WRITE(ERR_MSG,1000) trim(iVal(LC1))
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               IF(NO_K) THEN
                  WRITE(ERR_MSG,1001) dPAR_POS(1:2,LC1)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ELSE
                  WRITE(ERR_MSG,1002) dPAR_POS(1:3,LC1)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ENDIF
            ENDIF
         ENDDO  ! Loop over particles
      ENDIF

 1000 FORMAT('Error 1000: Unable to locat paritcle inside domain:',/&
         3x,'Particle Number:',A)
 1001 FORMAT(3x,'X POS: ',g11.5,/3x,'Y POS: ',g11.5)
 1002 FORMAT(3x,'X POS: ',g11.5,/3x,'Y POS: ',g11.5,/3x,'Z POS: ',g11.5)

! Send out the error flag and exit if needed.
      CALL BCAST(IER, PE_IO)
      IF(IER(PE_IO) /= 0) CALL MFIX_EXIT(myPE)

! PE_IO sends out the number of particles for each process.
      CALL BCAST(lPAR_CNT(0:NUMPES-1), PE_IO)

! Each process stores the number of particles-on-its-process. The error
! flag is set if that number exceeds the maximum.
      PIP = lPAR_CNT(myPE)
      IF(PIP > MAX_PIP) IER(myPE) = MAX_PIP+1

! Global collection of error flags to abort it the max was exceeded.
      CALL GLOBAL_ALL_SUM(IER)
      IF(sum(IER) /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         DO LC1=0, numPEs-1
            IF(IER(LC1) /= 0) THEN
               WRITE(ERR_MSG,"(3(2x,I10))")LC1,IER(LC1)-1,lPAR_CNT(LC1)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF
         ENDDO
         WRITE(ERR_MSG,"('Aborting.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE.,ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Maximum number of particles exceeded.',2/    &
         5x,'Process',5x,'Maximum',7x,'Count')


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAP_pARRAY_TO_PROC



!``````````````````````````````````````````````````````````````````````!
! Subroutine: DES_RESTART_MAP                                          !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE SCATTER_PAR_POS(lPAR_CNT)

      implicit none

! Number of particles on each process.
      INTEGER, INTENT(INOUT) :: lPAR_CNT(0:numPEs-1)
! Dimensionality flag.
      INTEGER :: lDIMN
! Loop counters.
      INTEGER :: LC1, lPROC, lBuf, IER
! Scatter counts.
      INTEGER :: lScatterCNTS(0:NUMPEs-1)


      lDIMN = merge(2,3,NO_K)

! Set up the recv count and allocate the local process buffer.
      iSCR_RECVCNT = PIP*lDIMN
      allocate (dProcBuf(iscr_recvcnt))

! Allocate the buffer for the root.
      IF (myPE == PE_IO) THEN
         allocate (dRootBuf(pIN_COUNT*lDIMN))
      ELSE
         allocate (dRootBuf(10))
      ENDIF

! The IO processor builds drootbuffer and iDISLS
      IF(myPE == PE_IO) THEN
! Determine the offsets for each process and the amount of data that
! is to be scattered to each.
         iDISPLS(0) = 0
         iScatterCnts(0) = lPAR_CNT(0)*lDIMN
         DO lProc = 1, NUMPES-1
            iDispls(lproc) = iDispls(lproc-1) + iScatterCnts(lproc-1)
            iScatterCnts(lproc) = lPAR_CNT(lProc)*lDIMN
         ENDDO
! Copy the position data into the root buffer, mapped to the owner
! process.
         lPAR_CNT(:) = 0
         DO LC1 = 1,pIN_COUNT
            lPROC = pRestartMap(LC1)
            lbuf = iDispls(lProc) + lPAR_CNT(lProc)*lDIMN+1
            dRootBuf(lBuf:lBuf+lDIMN-1) = dPAR_POS(1:lDIMN,LC1)
            lBuf = lBuf + lDIMN
            lPAR_CNT(lProc) = lPAR_CNT(lProc) + 1
         ENDDO
      ENDIF
      CALL DESMPI_SCATTERV(pTYPE=2)

! Unpack the particle data.
      DO LC1 = 1, PIP
         lBuf = (LC1-1)*lDIMN+1
         DES_POS_NEW(1:lDIMN,LC1) = dProcBuf(lBuf:lBuf+lDIMN-1)
         lBuf = lBuf + lDIMN
         PEA(LC1,1) = .TRUE.
      ENDDO

      IF(allocated(dRootBuf)) deallocate(dRootBuf)
      IF(allocated(dProcBuf)) deallocate(dProcBuf)
      IF(allocated(dPAR_POS)) deallocate(dPAR_POS)

      RETURN
      END SUBROUTINE SCATTER_PAR_POS


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_PAR_COL                                             !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_PAR_COL(lNEXT_REC)

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC

      INTEGER :: LC1, lPROC
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The number of particles on each process.
      INTEGER :: COL_CNT(0:NUMPEs-1)

!-----------------------------------------------

      CALL INIT_ERR_MSG("READ_PAR_COL")

! All process read positions for distributed IO restarts.
      IF(bDIST_IO) THEN
         DO LC1 = 1, 2
            CALL READ_RES_DES(lNEXT_REC, COLLISIONS(LC1,:))
         ENDDO
         RETURN
      ENDIF

      allocate(iPAR_COL(2, cIN_COUNT))

! Only the IO proccess reads positions.
      IF(myPE == PE_IO) THEN
         DO LC1=1, 2
            CALL IN_BIN_512i(RDES_UNIT, iPAR_COL(LC1,:),               &
               cIN_COUNT, lNEXT_REC)
         ENDDO
      ENDIF

! Use the particle postions and the domain coverage of each process
! to determine which processor each particle belongs.
      CALL MAP_cARRAY_TO_PROC(COL_CNT)

! Send the particle position data to the individual ranks.
      CALL SCATTER_PAR_COL(COL_CNT)

! Set up the read/scatter arrary information.
      cPROCCNT = COLLISION_NUM
      cROOTCNT = cIN_COUNT

! Set the recv count for this process.
      cRECV = COLLISION_NUM

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      lScatterCnts(:) = 0; lScatterCnts(mype) = COLLISION_NUM
      CALL GLOBAL_SUM(lScatterCnts,cSCATTER)

! Calculate the displacements for each process in the global array.
      cDispls(0) = 0
      DO lPROC = 1, NUMPEs-1
         cDispls(lPROC) = cDispls(lPROC-1) + cSCATTER(lPROC-1)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE READ_PAR_COL


!``````````````````````````````````````````````````````````````````````!
! Subroutine: MAP_cARRAY_TO_PROC                                       !
!                                                                      !
! Purpose: Use the particle positions to determine which processor     !
! they live on and count the number of particles on each process.      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE MAP_cARRAY_TO_PROC(lCOL_CNT)

      implicit none

      INTEGER, INTENT(OUT) :: lCOL_CNT(0:numPEs-1)

! Loop counters.
      INTEGER :: LC1, LC2, lPROC, lBUF
! Error flag.
      INTEGER :: IER(0:numPEs-1)
! Local scatter count.
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
      INTEGER :: lGatherCNTS(0:NUMPEs-1)

      INTEGER :: lROOTCNT, lPROCCNT

!-----------------------------------------------

      CALL INIT_ERR_MSG("MAP_cARRAY_TO_PROC")

! Initialize the error flag.
      IER = 0

! Setup data for particle array data collection:
      lROOTCNT = 10
      lPROCCNT = PIP

! Rank 0 gets the total number of gloabl particles.
      CALL GLOBAL_SUM(lPROCCNT, lROOTCNT)

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      iGath_SendCnt = lPROCCNT

      lGatherCnts = 0
      lGatherCnts(myPE) = lPROCCNT

      CALL GLOBAL_SUM(lGatherCnts, iGatherCnts)

! Calculate the displacements for each process in the global array.
      iDISPLS(0) = 0
      DO lPROC = 1,NUMPES-1
         iDISPLS(lPROC) = iDISPLS(lPROC-1) + iGatherCnts(lPROC-1)
      ENDDO

      allocate(iPROCBUF(lPROCCNT))
      allocate(iROOTBUF(lROOTCNT))

      iProcBuf(1:PIP) = iGLOBAL_ID(1:PIP)
      CALL DESMPI_GATHERV(pTYPE=1)

! Loop over the neighbor pair list and match the read global ID to
! one of the global IDs.
      lCOL_CNT = 0
      IF(myPE == PE_IO) THEN
         cRestartMap = -1
         LC1_LP: DO LC1=1, cIN_COUNT
            DO LC2=1, lROOTCNT
               IF(iPAR_COL(1,LC1) == iROOTBUF(LC2)) THEN
                  DO lPROC = 0, numPEs-1
                     IF(LC2 >= iDISPLS(lPROC)) THEN
                        cRestartMap(LC1) = lPROC
                        lCOL_CNT(lPROC) = lCOL_CNT(lPROC) + 1
                        CYCLE LC1_LP
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO LC1_LP
         DO LC1 = 1, cIN_COUNT
            IF (cRestartMap(LC1) == -1) THEN
               IER(myPE) = -1
               WRITE(ERR_MSG,1000) trim(iVal(LC1))
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               IF(NO_K) THEN
                  WRITE(ERR_MSG,1001) dPAR_POS(1:2,LC1)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ELSE
                  WRITE(ERR_MSG,1002) dPAR_POS(1:3,LC1)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

 1000 FORMAT('Error 1000: Unable to locat paritcle inside domain:',/&
         3x,'Particle Number:',A)
 1001 FORMAT(3x,'X POS: ',g11.5,/3x,'Y POS: ',g11.5)
 1002 FORMAT(3x,'X POS: ',g11.5,/3x,'Y POS: ',g11.5,/3x,'Z POS: ',g11.5)

! Send out the error flag and exit if needed.
      CALL BCAST(IER, PE_IO)
      IF(IER(PE_IO) /= 0) CALL MFIX_EXIT(myPE)

! PE_IO sends out the number of particles for each process.
      CALL BCAST(lCOL_CNT(0:NUMPES-1), PE_IO)

! Each process stores the number of particles-on-its-process. The error
! flag is set if that number exceeds the maximum.
      COLLISION_NUM = lCOL_CNT(myPE)

      DO WHILE(COLLISION_NUM > COLLISION_MAX)
         COLLISION_MAX = 2*COLLISION_MAX
      ENDDO
      CALL COLLISION_GROW

 1100 FORMAT('Error 1100: Maximum number of particles exceeded.',2/    &
         5x,'Process',5x,'Maximum',7x,'Count')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAP_cARRAY_TO_PROC


!``````````````````````````````````````````````````````````````````````!
! Subroutine: SCATTER_PAR_COL                                          !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE SCATTER_PAR_COL(lCOL_CNT)

      implicit none

! Number of particles on each process.
      INTEGER, INTENT(INOUT) :: lCOL_CNT(0:numPEs-1)
! Loop counters.
      INTEGER :: LC1, lPROC, lBuf, IER
! Scatter counts.
      INTEGER :: lScatterCNTS(0:NUMPEs-1)

! Set up the recv count and allocate the local process buffer.
      iSCR_RECVCNT = COLLISION_NUM*2
      allocate (iProcBuf(iscr_recvcnt))

! Allocate the buffer for the root.
      IF (myPE == PE_IO) THEN
         allocate (iRootBuf(cIN_COUNT*2))
      ELSE
         allocate (iRootBuf(10))
      ENDIF

! The IO processor builds drootbuffer and iDISLS
      IF(myPE == PE_IO) THEN
! is to be scattered to each.
         iDISPLS(0) = 0
         iScatterCnts(0) = lCOL_CNT(0)*2
         DO lProc = 1, NUMPES-1
            iDispls(lproc) = iDispls(lproc-1) + iScatterCnts(lproc-1)
            iScatterCnts(lproc) = lCOL_CNT(lProc)*2
         ENDDO
! Copy the position data into the root buffer, mapped to the owner
! process.
         lCOL_CNT(:) = 0
         DO LC1 = 1,cIN_COUNT
            lPROC = cRestartMap(LC1)
            lbuf = iDispls(lProc) + lCOL_CNT(lProc)*2+1
            iRootBuf(lBuf:lBuf+2-1) = iPAR_COL(1:2,LC1)
            lBuf = lBuf + 2
            lCOL_CNT(lProc) = lCOL_CNT(lProc) + 1
         ENDDO
      ENDIF

! The IO processor builds drootbuffer and iDISLS
      IF(myPE == PE_IO) THEN
! Determine the offsets for each process and the amount of data that
! is to be scattered to each.
         iDISPLS(0) = 0
         iScatterCnts(0) = lCOL_CNT(0)*2
         DO lProc = 1, NUMPES-1
            iDispls(lproc) = iDispls(lproc-1) + iScatterCnts(lproc-1)
            iScatterCnts(lproc) = lCOL_CNT(lProc)*2
         ENDDO
! Copy the position data into the root buffer, mapped to the owner
! process.
         lCOL_CNT(:) = 0
         DO LC1 = 1,cIN_COUNT
            lPROC = cRestartMap(LC1)
            lbuf = iDispls(lProc) + lCOL_CNT(lProc)*2+1
            iRootBuf(lBuf:lBuf+2-1) = iPAR_COL(1:2,LC1)
            lBuf = lBuf + 2
            lCOL_CNT(lProc) = lCOL_CNT(lProc) + 1
         ENDDO
      ENDIF
      CALL DESMPI_SCATTERV(pTYPE=1)

! Unpack the particle data.
      DO LC1 = 1, COLLISION_NUM
         lBuf = (LC1-1)*2+1
         COLLISIONS(1:2,LC1) = iProcBuf(lBuf:lBuf+2-1)
         lBuf = lBuf + 2
      ENDDO

      IF(allocated(iRootBuf)) deallocate(iRootBuf)
      IF(allocated(iProcBuf)) deallocate(iProcBuf)

      RETURN
      END SUBROUTINE SCATTER_PAR_COL



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0I(lNEXT_REC, INPUT_I)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
!      INTEGER, INTENT(OUT) :: INPUT_I
      INTEGER :: INPUT_I

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
         CALL BCAST(INPUT_I)
      ENDIF

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_1I                                              !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1I(lNEXT_REC, INPUT_I)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: INPUT_I(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_I)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_I, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0D                                          !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
         CALL BCAST(INPUT_D, PE_IO)
      ENDIF
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1D                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1D(lNEXT_REC, INPUT_D)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_D)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_D, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0L                                          !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0L(lNEXT_REC, OUTPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L

      INTEGER :: OUTPUT_I

      OUTPUT_L = .TRUE.

      IF(bDIST_IO)THEN
         READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
         CALL BCAST(OUTPUT_I, PE_IO)
      ENDIF

      IF(OUTPUT_I == 1) OUTPUT_L = .TRUE.
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1L                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1L(lNEXT_REC, INPUT_L)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: INPUT_L(:)

      INTEGER, ALLOCATABLE :: INPUT_I(:)

      INTEGER :: lSIZE, LC1

      lSIZE = size(INPUT_I)
      ALLOCATE( INPUT_I(lSIZE))

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_I, PE_IO)
      ENDIF

      DO LC1=1, LSIZE
         IF(INPUT_I(LC1) == 1) THEN
            INPUT_L(LC1) = .TRUE.
         ELSE
            INPUT_L(LC1) = .FALSE.
         ENDIF
      ENDDO

      IF(allocated(INPUT_I)) deallocate(INPUT_I)

      RETURN
      END SUBROUTINE READ_RES_DES_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1I(lNEXT_REC, OUTPUT_I)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, pIN_COUNT, lNEXT_REC)
      ELSE

         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, PIP
            OUTPUT_I(LC1) = iProcBuf(LC1)
         ENDDO

      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1I



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1D(lNEXT_REC, OUTPUT_D)

      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(pPROCCNT))
      allocate(dROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, pIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_D(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512(RDES_UNIT, lBUF_D, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               dRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_D(LC1)
            ENDDO

            deallocate(lBUF_D)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=2)
         DO LC1=1, PIP
            OUTPUT_D(LC1) = dProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1L(lNEXT_REC, OUTPUT_L)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         allocate(lBUF_I(pIN_COUNT))
         CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)
         DO LC1=1,pIN_COUNT
            IF(lBUF_I(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
         deallocate(lBUF_I)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, PIP
            IF(iProcBuf(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1I(lNEXT_REC, OUTPUT_I)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, cIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, COLLISION_NUM
            OUTPUT_I(LC1) = iProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_cARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1D(lNEXT_REC, OUTPUT_D)

      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(cPROCCNT))
      allocate(dROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER


      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, cIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_D(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512(RDES_UNIT, lBUF_D, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               dRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_D(LC1)
            ENDDO

            deallocate(lBUF_D)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=2)
         DO LC1=1, COLLISION_NUM
            OUTPUT_D(LC1) = dProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1L(lNEXT_REC, OUTPUT_L)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      IF(bDIST_IO) THEN
         allocate(lBUF_I(cIN_COUNT))
         CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)
         DO LC1=1,cIN_COUNT
            IF(lBUF_I(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
         deallocate(lBUF_I)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, COLLISION_NUM
            IF(iProcBuf(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1L


      END MODULE READ_RES1_DES
