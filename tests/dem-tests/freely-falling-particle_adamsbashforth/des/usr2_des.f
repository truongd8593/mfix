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

      INCLUDE 'usrnlst.inc'

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

      Use des_rxns
      Use des_thermo
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
      INTEGER, PARAMETER :: uPos = 2030
      INTEGER, PARAMETER :: uVel = 2031

! Analytic position and velocity
      double precision :: lPos_Y, lVel_Y
! Absolute relative error between MFIX solution and analytic solution.
      double precision :: Pos_rErr, Vel_rErr

      double precision :: lGrav
      double precision :: lRad
      integer :: lStage

! Open the files.
      FNAME1 = 'POST_POS.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=uPos,FILE=FNAME1,STATUS='NEW')
         WRITE(uPos,"(8X,'Time',9X,'Stage', &
            8X,'Pos',13X,'Pos_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uPos,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      FNAME2 = 'POST_VEL.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=uVel,FILE=FNAME2,STATUS='NEW')
         WRITE(uVel,"(8X,'Time',9X,'Stage', &
            8X,'Vel',13X,'Vel_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uVel,FILE=FNAME2,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Set local variables.
      lGrav = -grav(2)
      lRad  = des_radius(1)
      lStage = 0

! Calculate the position and velocity of the particle

! Stage 1: Free fall
      if(lTime < time_c) then
         lStage = 1
         lPos_Y = y_s1(h0, lGrav, lTime)
         lVel_Y = dydt_s1(lGrav, lTime)
! Stage 2: Contact
      elseif( lTime < time_r) then
         lStage = 2
         lPos_Y = y_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
         lVel_Y = dydt_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
! Stage 3: Rebound
      else
         lStage = 3
         lPos_Y = y_s3(lRad, lGrav, lTime)
         lVel_Y = dydt_s3(lGrav, lTime)
      endif

! Calculate the absolute relative error.
       Pos_rErr = (ABS(lPos_Y - DES_POS_new(1,2))/ABS(lPos_Y))*100
       Vel_rErr = (ABS(lVel_Y - DES_VEL_new(1,2))/ABS(lVel_Y))*100

! Write the results to a file.
      WRITE(uPos,"(3x,F15.8,5x,I1,5X,F15.8,2(3x,F15.8))") lTime, &
         lStage, lPos_Y,DES_POS_new(1,2),Pos_rErr
      CLOSE(uPos)


      WRITE(uVel,"(3x,F15.8,5x,I1,5X,F15.8,2(3x,F15.8))")lTime, &
         lStage, lVel_Y, DES_VEL_new(1,2),Vel_rErr
      CLOSE(uVel)

      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
