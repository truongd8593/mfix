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

!-----------------------------------------------
! Modules 
!-----------------------------------------------      

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

! des_continuum_coupled must be true if des_continuum_hybrid              
      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
         WRITE(UNIT_LOG, 1094)
         CALL MFIX_EXIT(myPE)
      ENDIF
! DES_CONTINUUM_HYBRID does not work with GHD or QMOMK
      IF (TRIM(KT_TYPE)=='GHD') THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1091)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF (QMOMK) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1092)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF (C_F==UNDEFINED .AND. MMAX>0 .AND. DES_MMAX >0) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1095)
         CALL MFIX_EXIT(myPE)
      ENDIF                 

      IF (MPPIC) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1093)
         CALL MFIX_EXIT(myPE)
      ENDIF         
      

      RETURN


 1090 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'DES_MMAX not defined in mfix.dat. It must be defined',/10X,&
         'when using DES_CONTINUUM_HYBRID',/1X,70('*')/)

 1091 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'KT_TYPE cannot be set to GHD when using ',&
         'DES_CONTINUUM_HYBRID.',/1X,70('*')/)
 1092 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'QMOMK cannot be TRUE when using ',&
         'DES_CONTINUUM_HYBRID.',/1X,70('*')/)
 1093 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'MPPIC and DES_CONTINUUM_HYBRID cannot both be TRUE.',&
         /1X,70('*')/)
 1094 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'DES_CONTINUUM_COUPLED must be to true when using ',/10X,&
         'DES_CONTINUUM_HYBRID.',&
         /1X,70('*')/)
 1095 FORMAT(/1X,70('*')//' From: CHECK_DES_HYBRID',/' Message: ',&
         'C_F must be defined when DES_CONTINUUM_HYBRID TRUE ',/10X,&
         'and both continuum and discrete solids phases are ',&
         'present'/10X,'(MAX>=1 and DES_MMAX>=1).',/1X,70('*')/)

         END SUBROUTINE CHECK_DES_HYBRID
