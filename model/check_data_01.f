!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_01                                          C
!  Purpose: check the run control namelist section                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: allow Si unit                                              C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:  RUN_NAME, DESCRIPTION, UNITS, TIME, TSTOP    C
!                         DT, RUN_TYPE                                 C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_01 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, N, L 
      Character*80  Line(1)
!-----------------------------------------------
!
!
!
!
!                      Solids phase
!
      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' ' 
!
      IF ( (UNITS /= 'CGS').AND. (UNITS /= 'SI') ) CALL ERROR_ROUTINE ('check_data_01', &
         'UNITS not specified or illegal in mfix.dat', 1, 1) 
!
      IF (DT /= UNDEFINED) THEN 
         IF (TIME==UNDEFINED .OR. TIME<ZERO) CALL ERROR_ROUTINE (&
            'check_data_01', 'TIME not specified OR < 0.0 in mfix.dat', 1, 1) 
      ELSE 
         TIME = ZERO 
      ENDIF 
!
      IF (DT /= UNDEFINED) THEN 
         IF (TSTOP==UNDEFINED .OR. TSTOP<TIME) CALL ERROR_ROUTINE (&
            'check_data_01', 'TSTOP not specified OR TSTOP < TIME in mfix.dat'&
            , 1, 1) 
      ENDIF 
!
!      IF (DT.EQ.UNDEFINED .OR. DT.LT.ZERO)
!     &          CALL ERROR_ROUTINE ('check_data_01',
!     &          'DT not specified OR DT < 0.0 in mfix.dat',1,1)
      IF (DT==UNDEFINED .OR. DT<ZERO) THEN 
         ODT = ZERO 
      ELSE 
         ODT = ONE/DT 
      ENDIF 
      IF (.NOT.(RUN_TYPE=='NEW' .OR. RUN_TYPE=='RESTART_1' .OR. RUN_TYPE==&
         'RESTART_2')) CALL ERROR_ROUTINE ('check_data_01', &
         'RUN_TYPE not specified or illegal in mfix.dat', 1, 1) 
      IF (FRICTION) THEN 
         IF (.NOT.GRANULAR_ENERGY) CALL ERROR_ROUTINE ('check_data_01', &
            'For FRICTION=.T., GRANULAR_ENERGY must be turned on', 1, 1) 
         IF (SAVAGE>2 .OR. SAVAGE<0) CALL ERROR_ROUTINE ('check_data_01', &
            'Value of SAVAGE should be 0, 1, or 2', 1, 1) 
      ENDIF  
!
! sof: cannot use both Ahmadi and Simonin models at the same time.
      IF (AHMADI .AND. SIMONIN) CALL ERROR_ROUTINE ('check_data_01', &
            'Cannot set both AHMADI = .T. and SIMONIN = .T.', 1, 1)  
!
! sof: cannot use both L_scale0 and K-epsilon models at the same time.
      IF (K_Epsilon .AND. L_SCALE0/=ZERO) CALL ERROR_ROUTINE ('check_data_01', &
            'Cannot set both K_Epsilon = .T. and L_SCALE0 /= ZERO', 1, 1)
!
!  Set variable ANY_SPECIES_EQ
!
      ANY_SPECIES_EQ = .FALSE. 
      M = 0 
      IF (MMAX + 1 > 0) THEN 
         ANY_SPECIES_EQ = ANY_SPECIES_EQ .OR. ANY(SPECIES_EQ(:MMAX)) 
         M = MMAX + 1 
      ENDIF 
      
!
!  Check phase specification for Scalars
!

      DO N = 1, NScalar
        IF(Phase4Scalar(N) < 0 .OR. Phase4Scalar(N) > MMAX) THEN
	  WRITE (Line,'(A, I3, A, I4, A)')&
	  'Phase4Scalar( ', N, ') = ', Phase4Scalar(N), ' is invalid'

	  CALL ERROR_ROUTINE ('check_data_01', Line, 1, 1) 
	END IF
      END DO
      
!
! CHECK THE DISCRETIZATION METHDS
!

!      IF (DISCRETIZE(1) > 1 .OR. DISCRETIZE(2) > 1) THEN 
!         DISCRETIZE(1) = 0
!	 DISCRETIZE(2) = 0
!	 IF(DMP_LOG)WRITE (UNIT_LOG, 1500) 
!      ENDIF 
!      IF (.not.DEF_COR) THEN 
!        DEF_COR= .true.
!	IF(DMP_LOG)WRITE (UNIT_LOG, 1501) 
!      ENDIF


      IF (FPFOI) THEN
         DO L = 1,8
            IF(DISCRETIZE(L).LE.1) DISCRETIZE(L) = 2
         END DO
	 IF(DMP_LOG)WRITE (UNIT_LOG, 1502) 
      ENDIF
      
      IF(Chi_scheme)THEN
        IF(DISCRETIZE(7).NE.3 .and. DISCRETIZE(7).NE.5) THEN
	   IF(DMP_LOG)WRITE (UNIT_LOG, 1503)
	   CALL Mfix_exit(0) 
	ENDIF

      ENDIF


      RETURN  
 1500 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Continuity equations can be discretized only with FOUP.',/&
         'DISCRETIZE (1 or 2 or both) reset to 0.',G12.5,/1X,70('*')/) 
 1501 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Only deferred correction method can be used.',/&
         'DEF_COR reset to .true.',G12.5,/1X,70('*')/) 
 1502 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Discretize>=2 has to be used along with',/&
         'fourth order scheme',G12.5,/1X,70('*')/) 
 1503 FORMAT(/1X,70('*')//' From: CHECK_DATA_01',/' Message: ',&
         'Only schemes 3 and 5 are Chi-scheme enabled for',/&
	 'species equations [DISCRETIZE(7)].',/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_01 
