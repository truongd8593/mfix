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
      DOUBLE PRECISION, PARAMETER :: OUT_dT = 5.0d-4

      IF(FIRST_PASS) THEN
         CALL WRITE_DES_OUT(S_TIME)
         OUT_TIME = OUT_dT
         FIRST_PASS = .FALSE.
      ELSE
         IF(S_TIME >= OUT_TIME .AND. S_TIME <= END_SLIP) THEN
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
      INTEGER, PARAMETER :: uPos = 2030
      INTEGER, PARAMETER :: uVel = 2031

! Analytic position and velocity
      double precision :: xVel, aVel
! Absolute relative error between MFIX solution and analytic solution.
      double precision :: xVel_rErr, aVel_rErr

      double precision :: lGrav
      double precision :: lRad

! Open the files.
      FNAME1 = 'POST_XVEL.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=uPos,FILE=FNAME1,STATUS='NEW')
         WRITE(uPos,"(8X,'Time', 17X, &
            'xVel',13X,'aVEL_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uPos,FILE=FNAME1,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      FNAME2 = 'POST_AVEL.dat'
      INQUIRE(FILE=FNAME2,EXIST=F_EXISTS2)
      IF (.NOT.F_EXISTS2) THEN
         OPEN(UNIT=uVel,FILE=FNAME2,STATUS='NEW')
         WRITE(uVel,"(8X,'Time', 17X, &
            'aVel',13X,'aVEL_MFIX',10X,'REL ERR')")
      ELSE
         OPEN(UNIT=uVel,FILE=FNAME2,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF

! Set local variables.
      lGrav = -grav(2)
      lRad  = des_radius(1)

! Calculate the position and velocity of the particle
       xVel = v0 - MEW_W*lGrav*lTime
       aVel = -(5.0d0*MEW_W*lGrav*lTime) / (2.0d0*lRad)

! Calculate the absolute relative error.
       xVel_rErr = (ABS(xVel - DES_VEL_NEW(1,1))/ABS(xVel))*100
       aVel_rErr = (ABS(aVel - OMEGA_NEW(3,1))/ABS(aVel))*100

! Write the results to a file.
      WRITE(uPos,"(3x,F15.8,5X,F15.8,2(3x,F15.8))") lTime, &
         xVel, DES_VEL_NEW(1,1), xVel_rErr
      CLOSE(uPos)


      WRITE(uVel,"(3x,F15.8,5X,F15.8,2(3x,F15.8))")lTime, &
          aVel, OMEGA_NEW(3,1), aVel_rErr
      CLOSE(uVel)

      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
