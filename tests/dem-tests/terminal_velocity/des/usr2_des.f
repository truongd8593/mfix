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
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 1.0d-3

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
      INTEGER, PARAMETER :: uPOS = 2030
      INTEGER, PARAMETER :: uVEL = 2031

      DOUBLE PRECISION, SAVE :: RK4_TIME = 0.0d0
      DOUBLE PRECISION :: RK4_DT, RK4_DT_LAST
      DOUBLE PRECISION TIME_INTERVAL
      DOUBLE PRECISION, PARAMETER :: RK4_DT_DEFAULT = 1.0d-6
      INTEGER :: I, RK4_STEPS

! Position and velocity
      double precision, save :: RK4_POS(3) = (/0.5d0, 9.25d0, 0.5d0/)
      double precision, save :: RK4_VEL(3) = (/0.0d0, 0.00d0, 0.0d0/)

      double precision :: POS_ERR, VEL_ERR

! Open the files.
      FNAME1 = 'POST_POS.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=uPOS,FILE=FNAME1,STATUS='NEW')
         WRITE(uPOS,"(8X,'Time',17X, &
            'Pos',13X,'Pos_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uPOS,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

      FNAME2 = 'POST_VEL.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=uVEL,FILE=FNAME2,STATUS='NEW')
         WRITE(uVEL,"(8X,'Time',17x, &
            'Vel',13X,'Vel_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uVEL,FILE=FNAME2,&
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

! Take the total number of RK4 steps.
      DO I=1, RK4_STEPS
         CALL RK4_V2b3(RK4_DT, RK4_POS, RK4_VEL)
         RK4_TIME = RK4_TIME + RK4_DT
      ENDDO

! Take the 'last' RK4 step with accounts for the difference in times.
      IF(RK4_DT_LAST .NE. UNDEFINED) THEN
         CALL RK4_V2b3(RK4_DT, RK4_POS, RK4_VEL)
         RK4_TIME = RK4_TIME + RK4_DT_LAST
      ENDIF

! Calculate erors.
      POS_ERR = (ABS(RK4_POS(2) - DES_POS_NEW(1,2))/ABS(RK4_POS(2)))*100
      VEL_ERR = (ABS(RK4_VEL(2) - DES_VEL_NEW(1,2))/ABS(RK4_VEL(2)))*100

! Write the results to file.
      WRITE(uPOS,"(3x,F15.8,5X,F15.8,2(3x,F15.8))")                     &
         lTime, RK4_POS(2), DES_POS_NEW(1,2), POS_ERR

      WRITE(uVEL,"(3x,F15.8,5X,F15.8,2(3x,F15.8))")                     &
         lTime, RK4_VEL(2), DES_VEL_NEW(1,2), VEL_ERR

      CLOSE(uPOS)
      CLOSE(uVEL)


      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
