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
      double precision, DIMENSION(PARTICLES) :: coefs

      double precision :: lGrav
      double precision :: lRad
      double precision :: coefficient
      integer :: lStage

angles = (/2.075099,3.913043,9.664032,20.09881,30.23715,39.3083,49.56522,59.88142/)
coefs = (/0.490244, 0.343902, 0.191463, 0.37561, 0.542683, 0.692683, 0.802439/)

do angle_ii = 1, size(angles)

   write(angle_s,'(I2)') angle_ii

! Open the files.
      FNAME1 = 'POST_COEF_'//angle_s//'.dat'
      INQUIRE(FILE=FNAME1,EXIST=F_EXISTS1)
      IF (.NOT.F_EXISTS1) THEN
         OPEN(UNIT=outfile,FILE=FNAME1,STATUS='NEW')
         WRITE(outfile,"(8X,'Time',9X,'Stage', &
            8X,'Pos',13X,'X_Pos_MFIX',6X,'Y_Pos_MFIX',6X,'REL ERR',6X,'calc_coef',1X,'exp_coef')")
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

       coefficient = 390/sqrt(dot_product(DES_VEL_new(:,particle_index),DES_VEL_new(:,particle_index)))

! Write the results to a file.
      WRITE(outfile,"(3x,F15.8,5x,I1,5X,F15.8,7(3x,F15.8))") lTime, &
         lStage, lPos_Y,DES_POS_new(1,particle_index),DES_POS_new(2,particle_index),Pos_rErr,coefficient,coefs(particle_index)
      CLOSE(outfile)
   enddo
   print *,""

      END SUBROUTINE WRITE_DES_Out

      END SUBROUTINE USR2_DES
