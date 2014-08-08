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

      PUBLIC :: DES_RESTART_MAP

      PUBLIC :: READ_RES_DES
      PUBLIC :: READ_RES_DES_NPA


      INTERFACE READ_RES_DES
         MODULE PROCEDURE READ_RES_DES_0I
         MODULE PROCEDURE READ_RES_DES_1I
         MODULE PROCEDURE READ_RES_DES_0D
         MODULE PROCEDURE READ_RES_DES_1D
         MODULE PROCEDURE READ_RES_DES_0L
         MODULE PROCEDURE READ_RES_DES_1L
      END INTERFACE

      INTERFACE READ_RES_DES_NPA
         MODULE PROCEDURE READ_RES_DES_NPA_1I
         MODULE PROCEDURE READ_RES_DES_NPA_1D
         MODULE PROCEDURE READ_RES_DES_NPA_1L
      END INTERFACE


      INTEGER, PARAMETER :: RDES_UNIT = 901

      INTEGER :: IN_COUNT

      CONTAINS


!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_READ_RES_DES                                        !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_READ_RES_DES(BASE, lNEXT_REC)

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      INTEGER, INTENT(OUT) :: lNEXT_REC

      CHARACTER(len=32) :: lFNAME

      INTEGER :: lDIMN


      lDIMN = merge(2,3,NO_K)

      IF(bDIST_IO) THEN

         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

         READ(RDES_UNIT, REC=1) PIP
         READ(RDES_UNIT, REC=2) iGHOST_CNT

         IF(PIP > MAX_PIP) THEN
            write(*,*) "From des_read_restart:"
            write(*,*) "Error: The pip is grater than current max_pip"
            write(*,*) "pip=" ,pip,"; max_pip =", max_pip

         ENDIF

         IN_COUNT = PIP

      ELSE

         IF(myPE == PE_IO) THEN
            WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
            OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',      &
               STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

            READ(RDES_UNIT, REC=1) IN_COUNT
         ELSE
            IN_COUNT = 10
         ENDIF

         allocate( dPAR_POS( lDIMN, IN_COUNT))
         allocate( iRestartMap(IN_COUNT))

      ENDIF

      lNEXT_REC = 3

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

      IF(bDIST_IO) CALL DES_RESTART_NEIGH

      RETURN
      END SUBROUTINE FINL_READ_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: DES_RESTART_MAP                                          !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE DES_RESTART_MAP(lNEXT_REC)

!      use desmpi, only:

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC


      INTEGER :: lDIMN

      INTEGER :: LC1, lPROC, lBuf, IER

      INTEGER :: lScatterCNTS(0:NUMPEs-1)

      INTEGER :: lPAR_CNT(0:NUMPEs-1)

      DOUBLE PRECISION :: TMP_D

      DOUBLE PRECISION :: lxmin(0:NUMPEs-1), lxmax(0:NUMPEs-1)
      DOUBLE PRECISION :: lymin(0:NUMPEs-1), lymax(0:NUMPEs-1)
      DOUBLE PRECISION :: lzmin(0:NUMPEs-1), lzmax(0:NUMPEs-1)
!-----------------------------------------------

      CALL INIT_ERR_MSG("DES_RESTART_MAP")

      lDIMN = merge(2, 3, NO_K)

      IF(myPE == PE_IO) THEN
         DO LC1=1,lDIMN
            CALL IN_BIN_512(RDES_UNIT, dPAR_POS(LC1,:), &
               IN_COUNT, lNEXT_REC)
         ENDDO
      ENDIF

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
      iRestartMap(:) = -1
      lPAR_CNT(:) = 0
      IF(myPE == PE_IO) THEN
         DO LC1 = 1, IN_COUNT
            DO lPROC=0, NUMPEs-1
               IF(dPAR_POS(1,LC1) >= lxmin(lproc) .AND. &
                  dPAR_POS(1,LC1) <  lxmax(lproc) .AND. &
                  dPAR_POS(2,LC1) >= lymin(lproc) .AND. &
                  dPAR_POS(2,LC1) <  lymax(lproc)) THEN
                  IF(NO_K)THEN
                     lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                     iRestartMap(LC1) = lproc
                     EXIT
                  ELSE
                     IF(dPAR_POS(3,LC1) >= lzmin(lproc) .AND. &
                        dPAR_POS(3,LC1) <  lzmax(lproc)) THEN
                        lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                        iRestartMap(LC1) = lproc
                        EXIT
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO  ! Loop over processes
            IF (iRestartMap(LC1) == -1) then
               IER = 1
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

      CALL BCAST(IER, PE_IO)
      IF(IER /= 0) CALL MFIX_EXIT(myPE)

      CALL BCAST(lPAR_CNT(0:NUMPES-1), PE_IO)

      PIP = lPAR_CNT(myPE)
      IF(PIP > MAX_PIP) THEN
         WRITE(*,602) PIP, MAX_PIP
         call des_mpi_stop
      ENDIF

      iSCR_RECVCNT = PIP*lDIMN
      allocate (dProcBuf(iscr_recvcnt))

      IF (myPE == PE_IO) THEN
         allocate (dRootBuf(IN_COUNT*lDIMN))
      ELSE
         allocate (dRootBuf(10))
      ENDIF

! The IO processor builds drootbuffer and iDISLS
      IF(myPE == PE_IO) THEN

         iDISPLS(0) = 0
         iScatterCnts(0) = lPAR_CNT(0)*lDIMN

         DO lProc = 1, NUMPES-1
            iDispls(lproc) = iDispls(lproc-1) + iScatterCnts(lproc-1)
            iScatterCnts(lproc) = lPAR_CNT(lProc)*lDIMN
         ENDDO

         lPAR_CNT(:) = 0
         DO LC1 = 1,IN_COUNT
            lproc = irestartmap(LC1)
            lbuf = iDispls(lProc) + lPAR_CNT(lProc)*lDIMN+1
            dRootBuf(lBuf:lBuf+lDIMN-1) = dPAR_POS(1:lDIMN,LC1)
            lBuf = lBuf + lDIMN
            lPAR_CNT(lProc) = lPAR_CNT(lProc) + 1
         ENDDO
      ENDIF
      call desmpi_scatterv(ptype=2)

! Unpack the particles on each processor.
      DO LC1 = 1, PIP
         lBuf = (LC1-1)*lDIMN+1
         DES_POS_NEW(1:lDIMN,LC1) = dProcBuf(lBuf:lBuf+lDIMN-1)
         lBuf = lBuf + lDIMN
         PEA(LC1,1) = .TRUE.
      ENDDO

      IF(allocated(dRootBuf)) deallocate(dRootBuf)
      IF(allocated(dProcBuf)) deallocate(dProcBuf)
      IF(allocated(dPAR_POS)) deallocate(dPAR_POS)

! Set up the scatter arrays for future read/scatter operations.
      allocate(iProcBuf(PIP))
      allocate(iRootBuf(IN_COUNT))

      allocate(dProcBuf(PIP))
      allocate(dRootBuf(IN_COUNT))

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      lScatterCnts(:) = 0
      lScatterCnts(mype) = PIP
      iScr_RecvCnt = PIP

      CALL GLOBAL_SUM(lScatterCnts,iScatterCnts)

! Calculate the displacements for each process in the global array.
      iDispls(0) = 0
      DO lPROC = 1, NUMPEs-1
         iDispls(lPROC) = iDispls(lPROC-1) + iScatterCnts(lPROC-1)
      ENDDO


      CALL FINL_ERR_MSG


 600  FORMAT(/2X,'From: DES_RESTART_MAP: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 601  FORMAT(/2X,'From: DES_RESTART_MAP: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

 611  FORMAT(/2X,'From: DES_RESTART_MAP: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain'/,2x,'... DELETE IT!')



 602  FORMAT(/2X,'From: DES_RESTART_MAP: ',/2X,&
         'ERROR: Particles in the processor ',I10,&
         'exceeds MAX_PIP', I10)

      RETURN
      END SUBROUTINE DES_RESTART_MAP



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
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1I(lNEXT_REC, OUTPUT_I)

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


      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, IN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(IN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, IN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, IN_COUNT
               lPROC = iRestartMap(LC1)
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
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1D(lNEXT_REC, OUTPUT_D)

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


      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, IN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_D(IN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512(RDES_UNIT, lBUF_D, IN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, IN_COUNT
               lPROC = iRestartMap(LC1)
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
! Subroutine: READ_RES_DES_1D                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1L(lNEXT_REC, OUTPUT_L)

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


      IF(bDIST_IO) THEN
         allocate(lBUF_I(IN_COUNT))
         CALL IN_BIN_512i(RDES_UNIT, lBUF_I, IN_COUNT, lNEXT_REC)
         DO LC1=1,IN_COUNT
            IF(lBUF_I(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
         deallocate(lBUF_I)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(IN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, IN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, IN_COUNT
               lPROC = iRestartMap(LC1)
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

      RETURN
      END SUBROUTINE READ_RES_DES_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_NPA_1I(lNEXT_REC, INPUT_I)

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
      END SUBROUTINE READ_RES_DES_NPA_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_NPA_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_NPA_1D(lNEXT_REC, INPUT_D)

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
      END SUBROUTINE READ_RES_DES_NPA_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_NPA_1L(lNEXT_REC, INPUT_L)

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
      END SUBROUTINE READ_RES_DES_NPA_1L



      END MODULE READ_RES1_DES
