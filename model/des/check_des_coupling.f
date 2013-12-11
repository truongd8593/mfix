!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_COUPLING                                      !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input on coupling settings and write some       !
!  information to the screen and log file.                             !
!                                                                      !
!  Comments: Moved org code from CHECK_DES_DATA.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_COUPLING

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Utilize cutcell geometry.
      USE cutcell, only: CARTESIAN_GRID
! Runtime Flag: Interpolate DEM field quanties.
      USE discretelement, only: DES_INTERP_MEAN_FIELDS
! Runtime Flag: Invoke gas/solids coupled simulation.
      USE discretelement, only: DES_CONTINUUM_COUPLED
! Runtime Flag: Interpolate DEM field quanties.
      USE discretelement, only: DES_INTERP_ON
! File unit for LOG messages.
      USE funits, only: UNIT_LOG

      USE mpi_utility

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
!  NONE


! Verify that the mean fields are interpolated where required.
      IF(.NOT.DES_INTERP_MEAN_FIELDS) THEN
! MPPIC and DES/CUTCELL simulations require that the mean feilds be
! interpolated without regard to user specifications.
         IF(MPPIC.OR.CARTESIAN_GRID) THEN
            DES_INTERP_MEAN_FIELDS = .TRUE.
            IF(myPe.EQ.pe_IO) WRITE(*,1000) 
            IF(DMP_LOG) WRITE(UNIT_LOG,1000) 

! DES_INTERP_MEAN_FIELDS is invoked if des_interp_on is true to remain
! consistent with previous implementations.
         ELSEIF(DES_CONTINUUM_COUPLED.and.DES_INTERP_ON) THEN
            IF(myPe.eq.pe_IO) WRITE(*,1001)
            IF(DMP_LOG) WRITE(UNIT_LOG,1001)
            DES_INTERP_MEAN_FIELDS= .TRUE.

         ELSE
            IF(myPe.eq.pe_IO) WRITE(*,1002)  
            IF(DMP_LOG) WRITE(UNIT_LOG,1002)
         ENDIF

      ELSE
         IF(myPe.EQ.pe_IO) WRITE(*,1003)
         IF(DMP_LOG) WRITE(UNIT_LOG,1003)
      ENDIF
            
      RETURN

 1000 FORMAT(' MPPIC or CG run detected. Using backward interpolation',&
         ' for DES fields.')

 1001 FORMAT(' DES_INTERP_ON specified:  Using backward interpolation',&
         ' for DES fields.')

 1002 FORMAT(' No mean field interpolation specified.',/' Mean DES',   &
         ' fields obained by cell averaging.')

 1003 FORMAT(' DES_INTERP_MEAN_FIELDS specified: Using backward',      &
         ' interpolation for DES fields.')

      END SUBROUTINE CHECK_DES_COUPLING
