!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_INLET                                        !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !                                
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CHECK_DES_INLET

      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE param1
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables 
!-----------------------------------------------
      DOUBLE PRECISION TMP

!-----------------------------------------------


      WRITE(*,*) '---------- START CHECK_DES_INLET ---------->'

!  Check to see if there is sufficient boundary condition information
      IF(DES_BC_X_w == UNDEFINED .OR. &
         DES_BC_X_e == UNDEFINED .OR. &
         DES_BC_Y_s == UNDEFINED .OR. &
         DES_BC_Y_n == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF 
      IF(DIMN == 3)THEN
        IF(DES_BC_Z_b == UNDEFINED .OR. &
           DES_BC_Z_t == UNDEFINED) THEN
         WRITE (UNIT_LOG, 1000)
         CALL MFIX_EXIT(myPE)
        ENDIF
      ENDIF

! Check if the rate of particles entering the system is given in 
! volumetric flow rate and if so convert to particles per time step
      

      IF(DES_BC_VOLFLOW_s /= UNDEFINED)THEN
        IF(UNITS == 'CGS')THEN
          TMP = (DES_BC_VOLFLOW_s * DTSOLID) / (PI * DES_RADIUS(1)**3)
          IF(TMP .LT. 1)THEN
            PI_FACTOR = FLOOR(real(1.0/TMP))
            PI_COUNT = 1
          ELSE
            PI_FACTOR = 1
            PI_COUNT = CEILING(real(PI_COUNT))
          ENDIF
        ENDIF
        IF(UNITS == 'SI')THEN
          WRITE(*,*) '     UNITS = SI NOT COMPLETE -- HARD STOP'
          STOP
        ENDIF
      ENDIF

      WRITE(*,*) '     PI_FACTOR: ', PI_FACTOR !<--------------------------------   REMOVE
      WRITE(*,*) '     PI_COUNT: ', PI_COUNT !<----------------------------------   REMOVE

! Check if the rate of particles entering the system is given in 
! mass flow rate and if so convert to particles per time step


      IF(DES_BC_MASSFLOW_s /= UNDEFINED)THEN
        WRITE(*,*) '      DES_BC_MASSFLOW_s NOT COMPLETE - HARD STOP'
        STOP
      ENDIF



!      PRINT*,'WORKING -- HARD STOP FROM CHECK_DES_INLET' !<------------   REMOVE
!      STOP !<----------------------------------------------------------   REMOVE

 1000 FORMAT(/1X,70('*')//&
         ' From: CHECK_DES_INLET -',/&
         ' Message: Insufficient boundary condition infomation',//&
         1X,70('*')/)

      WRITE(*,*) '<---------- END CHECK_DES_INLET ----------'      
      RETURN
      END SUBROUTINE CHECK_DES_INLET
























































