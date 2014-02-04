!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_CONTROL

! Global Variables:
!---------------------------------------------------------------------//
! Time intervalue between updating the RES and SPx files.
      use output, only: RES_DT, SPX_DT
! Time-step intervalue between updating the .LOG file.
      use output, only: NLOG
! Flag: Use the K-Epsilon model
      use run, only: K_EPSILON
! Number of arrays to store in SPA
      use rxns, only: nRR

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, LARGE_NUMBER
! Number of SPx files.
      USE param1, only: N_SPX

! Global Module proceedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: LC


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_OUTPUT_CONTROL")


! Check the values specified for the RES file.
      IF (RES_DT==UNDEFINED)THEN
         WRITE(ERR_MSG,1000) 'RES_DT'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(RES_DT <= ZERO) THEN
         WRITE(ERR_MSG,1002) 'RES_DT', RES_DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Check the SPx Files
      SPx_LP: DO LC = 1, N_SPX

! Disable writing the .SPA file if nRR is unspecified.
         IF(LC == 10) THEN
            IF(nRR == 0) THEN
               IF (SPX_DT(LC) == UNDEFINED) SPX_DT(LC) = LARGE_NUMBER
               CYCLE SPx_LP
            ENDIF

! Disable writing the .SPB file if K-Epsilon is unspecified.
         ELSEIF(LC == 11) THEN
            IF(.NOT.K_Epsilon) THEN
               IF (SPX_DT(LC)==UNDEFINED)SPX_DT(LC) = LARGE_NUMBER
               CYCLE SPx_LP
            ENDIF

! Verify the remaining SPx files.
         ELSE
            IF(SPX_DT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) iVar('SPX_DT',LC)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

            ELSEIF(SPX_DT(LC) <= ZERO) THEN
               WRITE(ERR_MSG,1001) iVar('SPX_DT',LC), SPX_DT(LC)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO SPx_LP

! Verify that the LOG frequency is valid.
      IF(NLOG <= 0) THEN
         WRITE(ERR_MSG,1003) 'NLOG', NLOG
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',E14.6,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_OUTPUT_CONTROL
