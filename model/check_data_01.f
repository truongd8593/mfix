!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_01                                          C
!  Purpose: check the run control namelist section                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
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
      USE run
      USE physprop
      USE indices
      USE scalars
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
      INTEGER :: M, N 
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
      IF (UNITS /= 'CGS') CALL ERROR_ROUTINE ('check_data_01', &
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

      RETURN  
      END SUBROUTINE CHECK_DATA_01 
