!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:                                                        C
!    CHECK_ONE_AXIS(NA,DIMEN,ALENGTH,DA,AXIS,AXIS_INDEX,NO_IJK, SHIFT) C
!  Purpose: check geometry data for one axis                           C
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
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: PERCENT_ERROR, N_SPECIFIED, TEMP_STOR, LC          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_ONE_AXIS(NA, DIMEN, ALENGTH, DA, AXIS, AXIS_INDEX, &
         NO_IJK, SHIFT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE funits 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER NA, DIMEN 
      DOUBLE PRECISION ALENGTH 
      LOGICAL NO_IJK, SHIFT 
      CHARACTER AXIS, AXIS_INDEX 
!//EFD use explicit dimension for DA
!     DA should be dimensioned DA(DIMEN) rather than DA(0:DIMEN+1) to be able to use
!     the logic from previous versions that assumed DA(1) as the first element.  An error
!     check has been added to ensure that DX, DY and DZ definitions in mfix.dat starts
!     with the zeroth element; i.e. DA(1).
      DOUBLE PRECISION, DIMENSION(DIMEN) :: DA 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: PERCENT_ERROR = 1.0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: N_SPECIFIED, LC 
      DOUBLE PRECISION :: TEMP_STOR 
!-----------------------------------------------
!
! passed arguments:
!    NA = number of axis cells (IMAX,JMAX,KMAX)
!    ALENGTH = axis length (XLENGTH,YLENGTH,ZLENGTH)
!    DA = cell sizes (DX,DY,DZ)
!    AXIS = axis checked ('X','Y','Z')
!    AXIS_INDEX = index associated with AXIS ('I','J','K')
!    NO_IJK = Flag that specifies whether variation along that axis is
!             considered (passed variable for NO_I, NO_J, or NO_K)
!
! local variables:
!    PERCENT_ERROR - PERCENT ERROR ALLOWED IN AXIS LENGTH CHECKS
!    N_SPECIFIED   - NUMBER OF ITEMS SPECIFIED FROM NA,ALENGTH,DA
!    TEMP_STOR     - TEMPORARY STORAGE
!    LC            - LOOP COUNTER
!
!
! 0) Ensure that if DA is defined then it starts with DA(1); i.e. DX(0), DY(0) or DZ(0)
!
      if(.not.no_ijk)then
        IF( DA(2) /= UNDEFINED .AND. DA(1) == UNDEFINED) THEN
          CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
          IF(DMP_LOG)WRITE (UNIT_LOG, 1001) AXIS
          CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
        ENDIF
      endif
!
! 1) MAKE SURE AT LEAST TWO OF NA, ALENGTH, DA ARE SPECIFIED
!
      N_SPECIFIED = 0 
      IF (NA /= UNDEFINED_I) N_SPECIFIED = N_SPECIFIED + 1 
      IF (ALENGTH /= UNDEFINED) N_SPECIFIED = N_SPECIFIED + 1 
      IF (DA(1) /= UNDEFINED) N_SPECIFIED = N_SPECIFIED + 1 
      IF (N_SPECIFIED < 2) THEN 
         CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) AXIS, AXIS, AXIS, AXIS_INDEX 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
!
! 2) NUMBER OF CELLS NOT SPECIFIED - calculate NA based on
!    input that was specified
!
      IF (NA == UNDEFINED_I) THEN
        if(no_ijk)then
	  na = 1
	else 
          IF (DA(2) == UNDEFINED) THEN 
            TEMP_STOR = ALENGTH/DA(1) 
            NA = NINT(TEMP_STOR) 
            IF (NA - 1 > 0) THEN 
               DA(2:NA) = DA(1) 
            ENDIF 
          ELSE 
            NA = DIMEN 
            DO LC = 2, DIMEN 
               IF (DA(LC) == UNDEFINED) THEN 
                  NA = LC - 1 
                  EXIT 
               ENDIF 
            END DO 
          ENDIF 
        endif 
        GO TO 700 
      ENDIF 
!
      IF (NA>=0 .AND. NA<=DIMEN) THEN 
!
! 3) AXIS LENGTH NOT SPECIFIED - calculate ALENGTH based on
!    input that was specified
!
         IF (ALENGTH == UNDEFINED) THEN 
	   if(no_ijk)then
             ALENGTH = DA(1) 
	   else
             IF (DA(2) == UNDEFINED) THEN 
               IF (NA - 1 > 0) THEN 
                  DA(2:NA) = DA(1) 
               ENDIF 
             ENDIF 
             ALENGTH = 0.0 
             IF (NA > 0) THEN 
               ALENGTH = SUM(DA(:NA)) 
             ENDIF
	   endif 
         ENDIF 
!
! 4) CELL SIZE NOT SPECIFIED - calculate NON_VARIABLE DA based on
!    input that was specified
!
         IF (DA(1) == UNDEFINED) THEN 
            TEMP_STOR = ALENGTH/DBLE(NA) 
            IF (NA > 0) THEN 
               DA(:NA) = TEMP_STOR 
            ENDIF 
         ENDIF 
!
! 5) ALL 3 SPECIFIED
!
         if(.not.no_ijk)then
           IF (DA(2) == UNDEFINED) THEN 
             IF (NA - 1 > 0) THEN 
               DA(2:NA) = DA(1) 
             ENDIF 
           ENDIF
	 endif 
!
! 6) CHECK CONSISTENCY OF AXIS INPUT
!
      ENDIF 
  700 CONTINUE 
      IF (NA<0 .OR.  .NOT.NO_IJK .AND. NA>DIMEN-2) THEN 
         CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) AXIS_INDEX, NA, AXIS_INDEX, DIMEN, AXIS_INDEX, &
            AXIS_INDEX 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
!
      TEMP_STOR = 0.0 
      DO LC = 1, NA 
         IF (DA(LC)<=0.0 .OR. DA(LC)==UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1200) AXIS, LC 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
         TEMP_STOR = TEMP_STOR + DA(LC) 
      END DO 
      IF (ALENGTH <= 0.0) THEN 
         CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1300) AXIS 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      TEMP_STOR = 100.0*ABS(TEMP_STOR - ALENGTH)/ALENGTH 
      IF (TEMP_STOR > PERCENT_ERROR) THEN 
         CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1400) AXIS, AXIS 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' %(AXIS_LENGTH - SUM(DAs))/AXIS_LENGTH = ', &
            TEMP_STOR 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' AXIS LENGTH = ', ALENGTH 
         DO LC = 1, NA 
            IF(DMP_LOG)WRITE (UNIT_LOG, *) ' next DA = ', DA(LC) 
         END DO 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
!
      DO LC = NA + 1, DIMEN 
         IF (SHIFT .AND. DA(LC)/=UNDEFINED) THEN 
            CALL ERROR_ROUTINE ('check_one_axis', 'AXIS error', 0, 2) 
            IF(DMP_LOG)WRITE (UNIT_LOG, 1500) AXIS, LC, AXIS_INDEX, NA 
            CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
         ENDIF 
      END DO 
      RETURN  
!
 1000 FORMAT(1X,'not enough info supplied for ',A1,'-axis',/,1X,&
         'AT LEAST TWO of      ',A1,'LENTGH , D',A1,' , ',A1,&
         'MAX     must be specified') 
 1001 FORMAT(1X, 'The grid specification must start with D',A1,'(0)') 
 1100 FORMAT(1X,'BAD VALUE FOR ',A1,'MAX = ',I6,/,1X,'DIMENSION_',A1,&
         ' IN param.inc = ',I6,/,1X,A1,&
         'MAX+2 must be less than or equal to DIMENSION_',A1) 
 1200 FORMAT(1X,'D',A1,'(',I6,') is not positive or not specified') 
 1300 FORMAT(1X,A1,'LENGTH is not postive') 
 1400 FORMAT(1X,A1,'LENGTH and D',A1,'s not consistent') 
 1500 FORMAT(1X,'D',A1,'(',I6,') has been set.',/,1X,'Only ',A1,'MAX (',I6,&
         ' ) of these are needed.') 
!
      END SUBROUTINE CHECK_ONE_AXIS 
