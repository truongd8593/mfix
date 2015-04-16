!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! None

! Local variables
!---------------------------------------------------------------------//
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
      DOUBLE PRECISION, SAVE :: OUT_TIME
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 1.0d-5

      IF(FIRST_PASS) THEN
         CALL WRITE_DES_OUT(S_TIME)
         OUT_TIME = OUT_dT
         FIRST_PASS = .FALSE.
      ELSE
         IF(S_TIME >= OUT_TIME) THEN
            CALL WRITE_DES_Out(S_TIME)
            OUT_TIME = OUT_TIME + OUT_dT
         ENDIF
      ENDIF

      RETURN

      contains

!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_DES_Out(lTime)

      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION lTime


! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: uPos1 = 2030
      INTEGER, PARAMETER :: uPos2 = 2031


      double precision :: v1, v2
      double precision :: w1, w2

      double precision :: y1_K1, y1_K2, y1_K3, y1_K4
      double precision :: y2_K1, y2_K2, y2_K3, y2_K4
      double precision :: x1_K1, x1_K2, x1_K3, x1_K4
      double precision :: x2_K1, x2_K2, x2_K3, x2_K4

      DOUBLE PRECISION, SAVE :: RK4_TIME = 0.0d0
      DOUBLE PRECISION :: RK4_DT, RK4_DT_LAST
      DOUBLE PRECISION TIME_INTERVAL
      DOUBLE PRECISION, PARAMETER :: RK4_DT_DEFAULT = 1.0d-6
      INTEGER :: I, RK4_STEPS
      LOGICAL :: NOISY

! Position and velocity
      double precision, save :: y1 = 0.045d0
      double precision, save :: y2 = 0.135d0
      double precision, save :: x1 = 0.000d0
      double precision, save :: x2 = 0.000d0

! Open the files.
      FNAME1 = 'POST_POS1.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=uPos1,FILE=FNAME1,STATUS='NEW')
         WRITE(uPos1,"(8X,'Time',9X,'Stage', &
            8X,'Pos',13X,'Pos_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uPos1,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      FNAME2 = 'POST_POS2.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=uPos2,FILE=FNAME2,STATUS='NEW')
         WRITE(uPos2,"(8X,'Time',9X,'Stage', &
            8X,'Pos',13X,'Pos_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uPos2,FILE=FNAME2,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF


! Calculate the value for the RK4 solutions.
      TIME_INTERVAL = lTime - RK4_TIME
      IF(TIME_INTERVAL .LE. RK4_DT) THEN
         RK4_STEPS = 1
         RK4_DT = TIME_INTERVAL
         RK4_DT_LAST = UNDEFINED
      ELSE
         RK4_STEPS = floor(real(TIME_INTERVAL/RK4_DT_DEFAULT))
         RK4_DT = RK4_DT_DEFAULT
         RK4_DT_LAST = lTime - (RK4_TIME + RK4_STEPS*RK4_DT)
      ENDIF



      NOISY = .FALSE.
      IF(NOISY) THEN
         write(*,*) ''
         write(*,*) ''
         write(*,*) ' Entering RK4 calculations:'
         write(*,*) '    lTime: ', lTime
         write(*,*) '    RK4_TIME: ', RK4_TIME
         write(*,*) '    TIME_INTERVAL: ', TIME_INTERVAL
         write(*,*) '    RK4_STEPS: ', RK4_STEPS
         write(*,*) '    RK4_DT :', RK4_DT
         write(*,*) '    RK4_DT_LAST: ',RK4_DT_LAST
      ENDIF



      DO I=1, RK4_STEPS
         CALL RK4_V4(RK4_DT, y1, x1, y2, x2)
         RK4_TIME = RK4_TIME + RK4_DT
      ENDDO


      IF(RK4_DT_LAST .NE. UNDEFINED) THEN
         CALL RK4_V4(RK4_DT_LAST, y1, x1, y2, x2)
         RK4_TIME = RK4_TIME + RK4_DT_LAST
      ENDIF



      IF(NOISY) THEN
         write(*,*) ' Leaving RK4 calculations:'
         write(*,*) '    RK4_TIME: ', RK4_TIME
         write(*,*) ' Last calculation values:'
         write(*,*)'    y1: ',y1
         write(*,*)'    y2: ',y2
         write(*,*)'    x1: ',x1
         write(*,*)'    x2: ',x2
      ENDIF

      if(y1 .NE. y1 .OR. y2 .NE. y2 .OR. x1 .NE. x1 .OR. x2 .NE. x2) then
         write(*,*)' NAN found... exiting'
         stop
      endif


! Write the results to a file.
      WRITE(uPos1,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, y1,        &
         DES_POS_new(2,1), (ABS(y1 - DES_POS_new(2,1))/ABS(y1))*100

      WRITE(uPos2,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, y2,        &
         DES_POS_new(2,2), (ABS(y2 - DES_POS_new(2,2))/ABS(y1))*100

      CLOSE(uPos1)
      CLOSE(uPos2)


      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
