CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CHECK_DATA_03 (SHIFT)                                  C
C  Purpose: check the geometry and discretization namelist section     C
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
C  Variables referenced: IJKMAX2                                       C
C  Variables modified: IMAX, XLENGTH, DX, JMAX, YLENGTH, DY, KMAX      C
C                      ZLENGTH, DZ                                     C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CHECK_DATA_03 (SHIFT)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'funits.inc'
C
C  Local variables
C
C                    Shift or not shift DX, DY and DZ values
      LOGICAL        SHIFT
C
C  Check XMIN specification
C
      IF(XMIN .LT. ZERO) THEN
        WRITE(UNIT_LOG, 990)
      ENDIF
C
C If no variation in a direction is considered, the number of cells in
C that direction should be 1
C
      IF(NO_I) THEN
        IF(IMAX .EQ. UNDEFINED_I) IMAX = 1
        IF(DX(1) .EQ. UNDEFINED .AND. XLENGTH .EQ. UNDEFINED) THEN
          DX(1) = ONE
          XLENGTH = ONE
        ENDIF
      ENDIF
      IF(NO_J) THEN
        IF(JMAX .EQ. UNDEFINED_I) JMAX = 1
        IF(DY(1) .EQ. UNDEFINED .AND. YLENGTH .EQ. UNDEFINED) THEN
          DY(1) = ONE
          YLENGTH = ONE
        ENDIF
      ENDIF
      IF(NO_K) THEN
        IF(KMAX .EQ. UNDEFINED_I) KMAX = 1
        IF(DZ(1) .EQ. UNDEFINED .AND. ZLENGTH .EQ. UNDEFINED) THEN
          IF(COORDINATES .EQ. 'CYLINDRICAL') THEN
            DZ(1) = 8. * ATAN(ONE)
            ZLENGTH = 8. * ATAN(ONE)
          ELSE
            DZ(1) = ONE
            ZLENGTH = ONE
          ENDIF
        ENDIF
      ENDIF
      IF(NO_I .AND. IMAX .GT. 1)THEN
        WRITE(UNIT_LOG, 1000)
        STOP
      ENDIF
      IF(NO_J .AND. JMAX .GT. 1)THEN
        WRITE(UNIT_LOG, 1100)
        STOP
      ENDIF
      IF(NO_K .AND. KMAX .GT. 1)THEN
        WRITE(UNIT_LOG, 1200)
        STOP
      ENDIF
C
C CHECK THE DATA FOR THE INDIVIDUAL AXES
C this must be changed if something other than 'NEW' or 'RESTART_1'
C
      CALL CHECK_ONE_AXIS 
     &  (IMAX,DIMENSION_I,XLENGTH,DX,'X','I', NO_I, SHIFT)
      CALL CHECK_ONE_AXIS
     &  (JMAX,DIMENSION_J,YLENGTH,DY,'Y','J', NO_J, SHIFT)
      CALL CHECK_ONE_AXIS 
     &  (KMAX,DIMENSION_K,ZLENGTH,DZ,'Z','K', NO_K, SHIFT)
C
      DO_I = .NOT.NO_I
      DO_J = .NOT.NO_J
      DO_K = .NOT.NO_K
C
      IF(COORDINATES .EQ. 'CYLINDRICAL')THEN
        CYLINDRICAL = .TRUE.
        IF(CYCLIC_X .OR. CYCLIC_X_PD)THEN
          WRITE(UNIT_LOG, 1250)
          STOP
        ENDIF
      ELSEIF(COORDINATES .EQ. 'CARTESIAN')THEN
        CYLINDRICAL = .FALSE.
      ELSE
        WRITE(UNIT_LOG, 1300)
        STOP
      ENDIF
C
C calculate IMAX1, IMAX2, etc. and shift the DX,DY,DZ arrays to take into
C account the fictitious cells
C
      CALL SET_MAX2
      IF(SHIFT) CALL SHIFT_DXYZ  !only for new and restart_1 runs
C
C  Ensure that the cell sizes across cyclic boundaries are comparable
C
      IF(CYCLIC_X .OR. CYCLIC_X_PD)THEN
        IF(DX(IMIN1) .NE. DX(IMAX1)) THEN
          WRITE(UNIT_LOG, 1400)DX(IMIN1), DX(IMAX1)
          STOP
        ENDIF
      ENDIF
C
      IF(CYCLIC_Y .OR. CYCLIC_Y_PD)THEN
        IF(DY(JMIN1) .NE. DY(JMAX1)) THEN
          WRITE(UNIT_LOG, 1410)DY(JMIN1), DY(JMAX1)
          STOP
        ENDIF
      ENDIF
C
      IF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYLINDRICAL)THEN
        IF(DZ(KMIN1) .NE. DZ(KMAX1)) THEN
          WRITE(UNIT_LOG, 1420)DZ(KMIN1), DZ(KMAX1)
          STOP
        ENDIF
      ENDIF
C
C CHECK THE TOTAL DIMENSION
C
      IF (IJKMAX2 .GT. DIMENSION_3) THEN
         CALL ERROR_ROUTINE ('check_data_03','dimension error',0,2)
         WRITE(UNIT_LOG,*) '(IMAX+2)*(JMAX+2)*(KMAX+2) = ' , IJKMAX2
         WRITE(UNIT_LOG,*) 'DIMENSION_3 in param.inc = ' , DIMENSION_3
         WRITE(UNIT_LOG,*) 'modify param.inc and recompile .... or'
         WRITE(UNIT_LOG,*) 'change dimensions in mfix.dat' , 
     &                          ' ... whichever is appropriate'
         CALL ERROR_ROUTINE (' ',' ',1,3)
      END IF
         
      RETURN
990   FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'XMIN should not be less than zero', 70('*')/)
1000  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'IMAX should be 1, since NO_I is true',/1X, 70('*')/)
1100  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'JMAX should be 1, since NO_J is true',/1X, 70('*')/)
1200  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'KMAX should be 1, since NO_K is true',/1X, 70('*')/)
1250  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'Cyclic bc for X not allowed in cylindrical coordinates',/1X,
     &  70('*')/)
1300  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'COORDINATES specified is illegal',/1X, 70('*')/)
1400  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'Cells adjacent to cyclic boundaries must be of same size:',/
     & 'DX(IMIN1) = ', G12.5, '     DX(IMAX1) = ', G12.5,
     & /1X, 70('*')/)
1410  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'Cells adjacent to cyclic boundaries must be of same size:',/
     & 'DY(JMIN1) = ', G12.5, '     DY(JMAX1) = ', G12.5,
     & /1X, 70('*')/)
1420  FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',
     & 'Cells adjacent to cyclic boundaries must be of same size:',/
     & 'DZ(KMIN1) = ', G12.5, '     DZ(KMAX1) = ', G12.5,
     & /1X, 70('*')/)
      END
