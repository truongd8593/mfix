CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name:                                                        C
C    CHECK_ONE_AXIS(NA,DIMEN,ALENGTH,DA,AXIS,AXIS_INDEX,NO_IJK, SHIFT) C
C  Purpose: check geometry data for one axis                           C
C                                                                      C
C  Author: P. Nicoletti                               Date: 27-NOV-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: PERCENT_ERROR, N_SPECIFIED, TEMP_STOR, LC          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CHECK_ONE_AXIS 
     &   (NA,DIMEN,ALENGTH,DA,AXIS,AXIS_INDEX,NO_IJK, SHIFT)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'funits.inc'
      DOUBLE PRECISION  PERCENT_ERROR
      PARAMETER (PERCENT_ERROR = 1.0)
C
      DOUBLE PRECISION    ALENGTH , DA(*) , TEMP_STOR
      INTEGER             NA , N_SPECIFIED , LC, DIMEN
      CHARACTER           AXIS*1 , AXIS_INDEX*1
      LOGICAL             NO_IJK, SHIFT 
C
C passed arguments: 
C    NA = number of axis cells (IMAX,JMAX,KMAX)
C    ALENGTH = axis length (XLENGTH,YLENGTH,ZLENGTH)
C    DA = cell sizes (DX,DY,DZ)
C    AXIS = axis checked ('X','Y','Z')
C    AXIS_INDEX = index associated with AXIS ('I','J','K')
C    NO_IJK = Flag that specifies whether variation along that axis is
C             considered (passed variable for NO_I, NO_J, or NO_K)
C
C local variables:
C    PERCENT_ERROR - PERCENT ERROR ALLOWED IN AXIS LENGTH CHECKS
C    N_SPECIFIED   - NUMBER OF ITEMS SPECIFIED FROM NA,ALENGTH,DA
C    TEMP_STOR     - TEMPORARY STORAGE
C    LC            - LOOP COUNTER
C
C
C 1) MAKE SURE AT LEAST TWO OF NA, ALENGTH, DA ARE SPECIFIED
C
      N_SPECIFIED = 0
      IF (NA.NE.UNDEFINED_I) N_SPECIFIED = N_SPECIFIED + 1
      IF (ALENGTH.NE.UNDEFINED) N_SPECIFIED = N_SPECIFIED + 1
      IF (DA(1).NE.UNDEFINED) N_SPECIFIED = N_SPECIFIED + 1
      IF (N_SPECIFIED.LT.2) THEN
         CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
         WRITE (UNIT_LOG,1000) AXIS,AXIS,AXIS,AXIS_INDEX
         CALL ERROR_ROUTINE (' ',' ',1,3)
      END IF
C
C 2) NUMBER OF CELLS NOT SPECIFIED - calculate NA based on
C    input that was specified
C
      IF (NA.EQ.UNDEFINED_I) THEN
         IF (DA(2).EQ.UNDEFINED) THEN
            TEMP_STOR = ALENGTH / DA(1)
            NA        = NINT(TEMP_STOR)
            DO 100 LC = 2,NA
               DA(LC) = DA(1)
100         CONTINUE
         ELSE
            NA = DIMEN
            DO 200 LC = 2,DIMEN
               IF (DA(LC).EQ.UNDEFINED) THEN
                  NA = LC - 1
                  GOTO 300
               END IF
200         CONTINUE
         END IF
300      GOTO 700
      END IF
C
      IF (NA.LT.0 .OR. NA.GT.DIMEN) GOTO 700
C
C 3) AXIS LENGTH NOT SPECIFIED - calculate ALENGTH based on
C    input that was specified
C
      IF (ALENGTH.EQ.UNDEFINED) THEN
         IF (DA(2).EQ.UNDEFINED) THEN
            DO 400 LC = 2,NA
               DA(LC) = DA(1)
400         CONTINUE
         END IF
         ALENGTH = 0.0
         DO 500 LC = 1,NA
            ALENGTH = ALENGTH + DA(LC)
500      CONTINUE
      END IF
C
C 4) CELL SIZE NOT SPECIFIED - calculate NON_VARIABLE DA based on
C    input that was specified
C
      IF (DA(1).EQ.UNDEFINED) THEN
         TEMP_STOR = ALENGTH / DBLE(NA)
         DO 600 LC = 1,NA
            DA(LC) = TEMP_STOR
600      CONTINUE
      END IF
C
C 5) ALL 3 SPECIFIED
C
      IF (DA(2).EQ.UNDEFINED) THEN
         DO 650 LC = 2,NA
            DA(LC) = DA(1)
650      CONTINUE
      END IF
C
C 6) CHECK CONSISTENCY OF AXIS INPUT
C
700   CONTINUE
C
      IF (NA.LT.0 .OR. (.NOT.NO_IJK .AND. NA.GT.DIMEN-2)) THEN
         CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
         WRITE(UNIT_LOG,1100) AXIS_INDEX , NA , AXIS_INDEX, DIMEN ,
     &                        AXIS_INDEX, AXIS_INDEX
         CALL ERROR_ROUTINE (' ',' ',1,3)
      END IF
C
      TEMP_STOR = 0.0
      DO 710 LC = 1,NA
         IF (DA(LC).LE.0.0 .OR.DA(LC).EQ.UNDEFINED) THEN
            CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
            WRITE(UNIT_LOG,1200) AXIS,LC
            CALL ERROR_ROUTINE (' ',' ',1,3)
         END IF
         TEMP_STOR = TEMP_STOR + DA(LC)
710   CONTINUE
C
      IF (ALENGTH.LE.0.0) THEN
          CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
          WRITE(UNIT_LOG,1300) AXIS
          CALL ERROR_ROUTINE (' ',' ',1,3)
      END IF
      TEMP_STOR = 100.0 * ABS(TEMP_STOR - ALENGTH) / ALENGTH
      IF (TEMP_STOR.GT.PERCENT_ERROR) THEN
          CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
          WRITE(UNIT_LOG,1400) AXIS , AXIS
          WRITE (UNIT_LOG,*)
     &      ' %(AXIS_LENGTH - SUM(DAs))/AXIS_LENGTH = ' , temp_stor
          WRITE (UNIT_LOG,*) ' AXIS LENGTH = ' , ALENGTH
          DO 720 LC = 1,NA
             WRITE (UNIT_LOG,*) ' next DA = ' , DA(LC)
720       CONTINUE
          CALL ERROR_ROUTINE (' ',' ',1,3)
      END IF
C
      DO 730 LC = NA+1,DIMEN
         IF (SHIFT .AND. DA(LC).NE.UNDEFINED) THEN
            CALL ERROR_ROUTINE ('check_one_axis','AXIS error',0,2)
            WRITE(UNIT_LOG,1500) AXIS , LC , AXIS_INDEX , NA
            CALL ERROR_ROUTINE (' ',' ',1,3)
         END IF
730   CONTINUE
C
      RETURN
C
1000  FORMAT(1X,'not enough info supplied for ',A1,'-axis',/,
     &       1X,'AT LEAST TWO of      ',A1,'LENTGH , D',A1,' , ' ,
     &       A1,'MAX     must be specified')
1100  FORMAT(1X,'BAD VALUE FOR ',A1,'MAX = ',I6,/,
     &       1X,'DIMENSION_',A1,' IN param.inc = ',I6,/,
     &       1X,A1,'MAX+2 must be less than or equal to DIMENSION_',A1)
1200  FORMAT(1X,'D',A1,'(',I6,') is not positive or not specified')
1300  FORMAT(1X,A1,'LENGTH is not postive')
1400  FORMAT(1X,A1,'LENGTH and D',A1,'s not consistent')
1500  FORMAT(1X,'D',A1,'(',I6,') has been set.',/,
     &       1X,'Only ',A1,'MAX (', I6, ' ) of these are needed.')
C
      END
