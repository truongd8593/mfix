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

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lTime

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER(LEN=64) :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: outfile = 2030
      INTEGER :: angle_ii,particle_index
      DOUBLE PRECISION, DIMENSION(PARTICLES) :: angles,rebounds

      character(LEN=10) :: angle_s

! Analytic position and velocity
      double precision :: lPos_Y, lVel_Y
! Absolute relative error between MFIX solution and analytic solution.
      double precision :: Pos_rErr
      double precision, DIMENSION(9), save :: Pos_x_old = 1
      double precision, DIMENSION(9), save :: Pos_y_old = 1

      double precision :: lGrav
      double precision :: lRad
      double precision :: rebound_angle
      integer :: lStage

angles = (/0.239726, 2.294521, 4.292237, 9.942922, 20.55936, 30.37671, 39.62329, 49.78311, 60.05708/)
rebounds = (/0.119419, 1.763049, 2.957663, 6.691719, 12.36445, 21.78089, 32.54549, 45.10513, 57.3652/)

do angle_ii = 1, size(angles)

   write(angle_s,'(I2)') angle_ii

! Open the files.
      FNAME1 = 'POST_ALPHA_'//angle_s//'.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=outfile,FILE=FNAME1,STATUS='NEW')
         WRITE(outfile,"(8X,'Time',9X,'Stage', &
            8X,'Pos',13X,'X_Pos_MFIX',6X,'Y_Pos_MFIX',6X,'REL ERR',6X,'calc_rebound',1X,'exp_rebound')")
      ELSE
         OPEN(UNIT=outfile,FILE=FNAME1,&
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

       particle_index = angle_ii

! Calculate the absolute relative error.
       Pos_rErr = (ABS(lPos_Y - DES_POS_new(2,particle_index))/ABS(lPos_Y))*100

       rebound_angle = 45*atan2(DES_POS_new(1,particle_index)-Pos_x_old(particle_index),DES_POS_new(2,particle_index)-Pos_y_old(particle_index))/ATAN(1.D0)
       pos_x_old(particle_index)  = DES_POS_new(1,particle_index)
       pos_y_old(particle_index)  = DES_POS_new(2,particle_index)

! Write the results to a file.
      WRITE(outfile,"(3x,F15.8,5x,I1,5X,F15.8,7(3x,F15.8))") lTime, &
         lStage, lPos_Y,DES_POS_new(1,particle_index),DES_POS_new(2,particle_index),Pos_rErr,rebound_angle,rebounds(particle_index)
      CLOSE(outfile)
   enddo

      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
