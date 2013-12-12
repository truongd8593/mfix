!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_HYBRID                                        !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for TFM/DEM hybrid model.            !
!                                                                      !
!  Comments: Moved org code from CHECK_DES_DATA.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_HYBRID

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke gas/solids coupled simulation.
      USE discretelement, only: DES_CONTINUUM_COUPLED
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Invoke QMOMK model.
      USE qmom_kinetic_equation, only: QMOMK
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! Number of TFM solids phases.
      USE physprop, only: MMAX
! User specified kinetic theory model.
      USE run, only: KT_TYPE
! Slope limiter parameter (0 < C _FAC <= 1.0)
      USE constant, only: C_f
! File unit for LOG messages.
      USE funits, only: UNIT_LOG
! Parameter to determine an uninitialized value.
      USE param1, only: UNDEFINED

      USE mpi_utility


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Integer error flag.
      INTEGER :: IER


! Initialize the local error flag.
      IER = 0

! des_continuum_coupled must be true if des_continuum_hybrid              
      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
         IF(myPE == PE_IO) WRITE(*, 1000)
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         IER = 1
      ENDIF

! Hybrid model does not work with GHD or QMOMK.
      IF (TRIM(KT_TYPE)=='GHD') THEN
         IF(myPE == PE_IO) WRITE(*, 1001)
         IF(DMP_LOG) WRITE(UNIT_LOG,1001)
         IER = 1
      ENDIF
      IF (QMOMK) THEN
         IF(myPE == PE_IO) WRITE(*, 1002)
         IF(DMP_LOG) WRITE(UNIT_LOG,1002)
         IER = 1
      ENDIF

      IF (C_F==UNDEFINED .AND. MMAX>0 .AND. DES_MMAX >0) THEN
         IF(myPE == PE_IO) WRITE(*, 1003)
         IF(DMP_LOG) WRITE(UNIT_LOG, 1003)
         IER = 1
      ENDIF                 

! MPPIC and Hybrid models don't work together.
      IF (MPPIC) THEN
         IF(myPE == PE_IO) WRITE(*, 1004)
         IF(DMP_LOG) WRITE(UNIT_LOG, 1004)
         IER = 1
      ENDIF         
      
      IF(IER /= 0) CALL MFIX_EXIT(myPE)

      RETURN


 1000 FORMAT(/1X,70('*')/' From: CHECK_DES_HYBRID',/' Error 1000:',     &
         ' DES_CONTINUUM_COUPLED must be to true when using ',/         &
         ' DES_CONTINUUM_HYBRID.', /1X,70('*')/)

 1001 FORMAT(/1X,70('*')/' From: CHECK_DES_HYBRID',/' Error 1001:',     &
         ' KT_TYPE cannot be set to GHD when using',                    &
         ' DES_CONTINUUM_HYBRID.',/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: CHECK_DES_HYBRID',/' Error 1002:',    &
         ' QMOMK cannot be TRUE when using DES_CONTINUUM_HYBRID.',/1X, &
         70('*')/)

 1003 FORMAT(/1X,70('*')/' From: CHECK_DES_HYBRID',/' Error 1003:',    &
         ' C_F must be defined for hybrid simulations containing',/    &
         ' continuum and discrete solids. (MMAX>=1 and DES_MMAX>=1)',/ &
         1X,70('*')/)

 1004 FORMAT(/1X,70('*')/' From: CHECK_DES_HYBRID',/' Error 1004:',    &
         'MPPIC and DES_CONTINUUM_HYBRID models are incompatible.',/   &
         1X,70('*')/)


         END SUBROUTINE CHECK_DES_HYBRID
