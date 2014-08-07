      MODULE RES_DES

      use desmpi
      use compar, only: myPE
      use compar, only: PE_IO

      use cdist, only: bDist_IO

      PRIVATE

      PUBLIC :: INIT_WRITE_RES_DES
      PUBLIC :: FINL_WRITE_RES_DES

      PUBLIC :: WRITE_RES_DES

      INTERFACE WRITE_RES_DES
         MODULE PROCEDURE WRITE_RES_DES_0I
         MODULE PROCEDURE WRITE_RES_DES_1I
         MODULE PROCEDURE WRITE_RES_DES_0D
         MODULE PROCEDURE WRITE_RES_DES_1D
         MODULE PROCEDURE WRITE_RES_DES_0L
         MODULE PROCEDURE WRITE_RES_DES_1L
      END INTERFACE



      INTEGER, PARAMETER :: RDES_UNIT = 901

      INTEGER :: OUT_COUNT

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: OPEN_RES_DES                                             !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE OPEN_RES_DES(BASE)

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      CHARACTER(len=32) :: lFNAME

      IF(bDIST_IO) THEN
         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

      ELSEIF(myPE == PE_IO) THEN
         WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)
      ENDIF

      END SUBROUTINE OPEN_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_WRITE_RES_DES                                       !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_WRITE_RES_DES(BASE, lNEXT_REC)


      use mpi_utility, only: GLOBAL_SUM

      use desmpi, only: iProcBuf, dProcBuf
      use desmpi, only: iRootBuf, dRootBuf

      use desmpi, only: iGath_SendCnt

      use desmpi, only: iDisPls

      use discretelement, only: PIP, iGHOST_CNT

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      INTEGER, INTENT(OUT) :: lNEXT_REC

! Number of real particles on local rank
      INTEGER :: lParCnt
! Total number of real particles.
      INTEGER :: lGloCnt

      INTEGER :: lGatherCnts(0:NUMPEs-1)

      CHARACTER(len=32) :: lFNAME


      lNEXT_REC = 1

      IF(bDIST_IO) THEN

         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

         allocate (dPROCBUF(PIP)) ! local process particle count
         allocate (iPROCBUF(PIP)) ! local process particle count

         OUT_COUNT = PIP

         CALL WRITE_RES_DES(lNEXT_REC, PIP)
         CALL WRITE_RES_DES(lNEXT_REC, iGHOST_CNT)

      ELSE

! Root process is informed of the total number of (real) particles.
         lGLOCNT = 10
         lPARCNT = PIP - iGHOST_CNT

! Rank 0 gets the total number of gloabl particles.
         CALL GLOBAL_SUM(lPARCNT, lGLOCNT)

! Allocate variables for gather.
! double precision arryas:
         allocate (dPROCBUF(lparcnt)) ! local process particle count
         allocate (dROOTBUF(lglocnt)) ! global particle count (root)
! integer arrays:
         allocate (iPROCBUF(lPARCNT)) ! local process particle count
         allocate (iROOTBUF(lGLOCNT)) ! global particle count (root)

! Construct an array for the Root process that states the number of
! (real) particles on each process.
         iGath_SendCnt = lParCnt

         lGatherCnts = 0
         lGatherCnts(myPE) = lParCnt

         CALL GLOBAL_SUM(lGatherCnts, iGatherCnts)

! Calculate the displacements for each process in the global array.
         iDispls(0) = 0
         DO lPROC = 1,NUMPES-1
            iDISPLS(lPROC) = iDISPLS(lPROC-1) + iGatherCnts(lPROC-1)
         ENDDO

         IF(myPE == PE_IO) THEN
            WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
            OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',      &
               STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)
         ENDIF

         OUT_COUNT = lGLOCNT

         CALL WRITE_RES_DES(lNEXT_REC, lGLOCNT)
         CALL WRITE_RES_DES(lNEXT_REC, 0)

      ENDIF

      RETURN
      END SUBROUTINE INIT_WRITE_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_WRITE_RES_DES

      IF(bDIST_IO .OR. myPE == PE_IO) close(RDES_UNIT)

      IF(allocated(dPROCBUF)) deallocate(dPROCBUF)
      IF(allocated(dROOTBUF)) deallocate(dROOTBUF)
      IF(allocated(iPROCBUF)) deallocate(iPROCBUF)
      IF(allocated(iROOTBUF)) deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE FINL_WRITE_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0I                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0I(lNEXT_REC, INPUT_I)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I

      IF(bDIST_IO .OR. myPE == PE_IO) &
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1I                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1I(lNEXT_REC, INPUT_I, pLOC2GLB)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)
      LOGICAL, INTENT(IN), OPTIONAL :: pLOC2GLB

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1, LC2

      lLOC2GLB = .FALSE.
      IF(present(pLOC2GLB)) lLOC2GLB = pLOC2GLB

      IF(bDIST_IO) THEN
         LC1 = 1

         IF(lLOC2GLB) THEN
            DO LC2 = 1, MAX_PIP
               IF(LC1 > PIP) EXIT
               IF(.NOT. PEA(LC1,1)) CYCLE
               iProcBuf(LC1) = iGLOBAL_ID(INPUT_I(LC2))
               LC1 = LC1 + 1
            ENDDO
         ELSE
            DO LC2 = 1, MAX_PIP
               IF(LC1 > PIP) EXIT
               IF(.NOT. PEA(LC1,1)) CYCLE
               iProcBuf(LC1) = INPUT_I(LC2)
               LC1 = LC1 + 1
            ENDDO
         ENDIF
         CALL OUT_BIN_512i(RDES_UNIT, iProcBuf, OUT_COUNT, lNEXT_REC)

      ELSE
         CALL DES_GATHER(INPUT_I, lLOC2GLB)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512i(RDES_UNIT,iROOTBUF, OUT_COUNT, lNEXT_REC)
      ENDIF

      RETURN
      END SUBROUTINE WRITE_RES_DES_1I




!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0D                                         !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D

      IF(bDIST_IO .OR. myPE == PE_IO) & 
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_D

      lNEXT_REC = lNEXT_REC + 1
      
      RETURN
      END SUBROUTINE WRITE_RES_DES_0D




!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1D                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1D(lNEXT_REC, INPUT_D)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

! Loop counters
      INTEGER :: LC1, LC2

      IF(bDIST_IO) THEN
         LC1 = 1
         DO LC2 = 1, MAX_PIP
            IF(LC1 > PIP) EXIT
            IF(.NOT. PEA(LC1,1)) CYCLE
            dProcBuf(LC1) = INPUT_D(LC2)
            LC1 = LC1 + 1
         ENDDO
         CALL OUT_BIN_512(RDES_UNIT, dProcBuf, OUT_COUNT, lNEXT_REC)
      ELSE
         CALL DES_GATHER(INPUT_D)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512(RDES_UNIT,dRootBuf, OUT_COUNT, lNEXT_REC)
      ENDIF

      RETURN
      END SUBROUTINE WRITE_RES_DES_1D



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0L                                         !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0L(lNEXT_REC, INPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L

      INTEGER :: INPUT_I

      INPUT_I = merge(1,0,INPUT_L)

      IF(bDIST_IO .OR. myPE == PE_IO) &
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0L

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1D                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1L(lNEXT_REC, INPUT_L)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA
 
      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

! Loop counters
      INTEGER :: LC1, LC2

      IF(bDIST_IO) THEN
         LC1 = 1
         DO LC2 = 1, MAX_PIP
            IF(LC1 > PIP) EXIT
            IF(.NOT. PEA(LC1,1)) CYCLE
            iProcBuf(LC1) = merge(1,0,INPUT_L(LC2))
            LC1 = LC1 + 1
         ENDDO
         CALL OUT_BIN_512(RDES_UNIT, iProcBuf, OUT_COUNT, lNEXT_REC)
      ELSE
         CALL DES_GATHER(INPUT_L)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512(RDES_UNIT,iRootBuf, OUT_COUNT, lNEXT_REC)
      ENDIF

      RETURN
      END SUBROUTINE WRITE_RES_DES_1L


      END MODULE RES_DES
