!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_03 (SHIFT)                                  C
!  Purpose: check the geometry and discretization namelist section     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:  Replaced STOP with mfix_exit() to abort all processors    C
!                                                                      C
!  Author:   Aeolus Res. Inc.                         Date: 07-AUG-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2                                       C
!  Variables modified: IMAX, XLENGTH, DX, JMAX, YLENGTH, DY, KMAX      C
!                      ZLENGTH, DZ                                     C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_03(SHIFT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE funits 
      USE compar
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL SHIFT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL COMPARE 
!-----------------------------------------------
      IF (XMIN < ZERO) WRITE (UNIT_LOG, 990) 
!
! If no variation in a direction is considered, the number of cells in
! that direction should be 1
!
      IF (NO_I) THEN 
!
         WRITE (UNIT_LOG, 995)                   !disabled 
         call mfix_exit(myPE)
!
!        IF(IMAX .EQ. UNDEFINED_I) IMAX = 1
!        IF(DX(1) .EQ. UNDEFINED .AND. XLENGTH .EQ. UNDEFINED) THEN
!          DX(1) = ONE
!          XLENGTH = ONE
!        ENDIF
      ENDIF 
      IF (NO_J) THEN 
!
         WRITE (UNIT_LOG, 996)                   !disabled 
         call mfix_exit(myPE) 

!
!        IF(JMAX .EQ. UNDEFINED_I) JMAX = 1
!        IF(DY(1) .EQ. UNDEFINED .AND. YLENGTH .EQ. UNDEFINED) THEN
!          DY(1) = ONE
!          YLENGTH = ONE
!        ENDIF
      ENDIF 


      IF (NO_K) THEN 
         IF (KMAX == UNDEFINED_I) KMAX = 1 
         IF (DZ(1)==UNDEFINED)THEN
	    IF(ZLENGTH==UNDEFINED) THEN 
              IF (COORDINATES == 'CYLINDRICAL') THEN 
                DZ(1) = 8.*ATAN(ONE) 
                ZLENGTH = 8.*ATAN(ONE) 
              ELSE 
                DZ(1) = ONE 
                ZLENGTH = ONE 
              ENDIF
	    ELSE
	      DZ(1) = ZLENGTH 
	    ENDIF
	 ELSE
	    IF(ZLENGTH==UNDEFINED) THEN
	      ZLENGTH = DZ(1)
	    ELSE
	      IF(.NOT.COMPARE(ZLENGTH,DZ(1)))THEN
                WRITE (UNIT_LOG, 997) 
                call mfix_exit(myPE) 
	      ENDIF
	    ENDIF 
         ENDIF 
      ENDIF 
      IF (NO_I .AND. IMAX>1) THEN 
         WRITE (UNIT_LOG, 1000) 
         call mfix_exit(myPE) 
      ENDIF 
      IF (NO_J .AND. JMAX>1) THEN 
         WRITE (UNIT_LOG, 1100) 
         call mfix_exit(myPE) 
      ENDIF 
      IF (NO_K .AND. KMAX>1) THEN 
         WRITE (UNIT_LOG, 1200) 
         call mfix_exit(myPE) 
      ENDIF 
!
! CHECK THE DATA FOR THE INDIVIDUAL AXES
! this must be changed if something other than 'NEW' or 'RESTART_1'
!

      CALL CHECK_ONE_AXIS (IMAX, DIMENSION_I, XLENGTH, DX, 'X', 'I', NO_I, &
         SHIFT)
      CALL CHECK_ONE_AXIS (JMAX, DIMENSION_J, YLENGTH, DY, 'Y', 'J', NO_J, &
         SHIFT)
      CALL CHECK_ONE_AXIS (KMAX, DIMENSION_K, ZLENGTH, DZ, 'Z', 'K', NO_K, &
         SHIFT)
!
      DO_I = .NOT.NO_I
      DO_J = .NOT.NO_J
      DO_K = .NOT.NO_K
!
      IF (COORDINATES == 'CYLINDRICAL') THEN 
         CYLINDRICAL = .TRUE. 
         IF (CYCLIC_X .OR. CYCLIC_X_PD) THEN 
            WRITE (UNIT_LOG, 1250) 
            call mfix_exit(myPE) 
         ENDIF 
      ELSE IF (COORDINATES == 'CARTESIAN') THEN 
         CYLINDRICAL = .FALSE. 
      ELSE 
         WRITE (UNIT_LOG, 1300) 
         call mfix_exit(myPE) 
      ENDIF 
!
! calculate IMAX1, IMAX2, etc. and shift the DX,DY,DZ arrays to take into
! account the fictitious cells
!
!      CALL SET_MAX2    ... moved to allocate_arrays
      IF (SHIFT) CALL SHIFT_DXYZ                 !only for new and restart_1 runs 
!
!  Ensure that the cell sizes across cyclic boundaries are comparable
!
      IF (CYCLIC_X .OR. CYCLIC_X_PD) THEN 
         IF (DX(IMIN1) /= DX(IMAX1)) THEN 
            WRITE (UNIT_LOG, 1400) DX(IMIN1), DX(IMAX1) 
            call mfix_exit(myPE) 
         ENDIF 
      ENDIF 
!
      IF (CYCLIC_Y .OR. CYCLIC_Y_PD) THEN 
         IF (DY(JMIN1) /= DY(JMAX1)) THEN 
            WRITE (UNIT_LOG, 1410) DY(JMIN1), DY(JMAX1) 
            call mfix_exit(myPE) 
         ENDIF 
      ENDIF 
!
      IF (CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYLINDRICAL) THEN 
         IF (DZ(KMIN1) /= DZ(KMAX1)) THEN 
            WRITE (UNIT_LOG, 1420) DZ(KMIN1), DZ(KMAX1) 
            call mfix_exit(myPE) 
         ENDIF 
      ENDIF 
!
! CHECK THE TOTAL DIMENSION
!


      IF (IJKMAX3 > DIMENSION_3G) THEN 
         CALL ERROR_ROUTINE ('check_data_03', 'global dimension error', 0, 2) 
         WRITE (UNIT_LOG, *) '(IMAX+2+1)*(JMAX+2+1)*(KMAX+2+1) = ', IJKMAX3 
         WRITE (UNIT_LOG, *) 'DIMENSION_3 in allocate_arrays = ', DIMENSION_3G 
         WRITE (UNIT_LOG, *) 'modify allocate_arrays and recompile .... or' 
         WRITE (UNIT_LOG, *) 'change dimensions in mfix.dat', &
            ' ... whichever is appropriate' 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      IF (IJKsize3_all(myPE) > DIMENSION_3) THEN 
         CALL ERROR_ROUTINE ('check_data_03', 'subdomain dimension error', 0, 2) 
         WRITE (UNIT_LOG, *) &
	 '(Iend3-Istart3+1)*(Jend3-Jstart3+1)*(Kend3-Kstart3+1) = ', IJKsize3 
         WRITE (UNIT_LOG, *) 'DIMENSION_3 in allocate_arrays = ', DIMENSION_3 
         WRITE (UNIT_LOG, *) 'modify allocated_arrays and recompile .... or' 
         WRITE (UNIT_LOG, *) 'change dimensions in mfix.dat', &
            ' ... whichever is appropriate' 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      
!
      RETURN  
  990 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'XMIN should not be less than zero',70('*')/) 
  995 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Error: ',&
         'This has been disabled.  Use one cell in I direction',/&
         ' and make the East and West walls free-slip',70('*')/) 
  996 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Error: ',&
         'This has been disabled.  Use one cell in J direction',/&
         ' and make the North and South walls free-slip',70('*')/) 
  997 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Error: ',&
         'DZ(1) and ZLENGTH are not equal!  (Recall NO_K=.true.)' ,70('*')/) 
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'IMAX should be 1, since NO_I is true',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'JMAX should be 1, since NO_J is true',/1X,70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'KMAX should be 1, since NO_K is true',/1X,70('*')/) 
 1250 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'Cyclic bc for X not allowed in cylindrical coordinates',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'COORDINATES specified is illegal',/1X,70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'Cells adjacent to cyclic boundaries must be of same size:',/&
         'DX(IMIN1) = ',G12.5,'     DX(IMAX1) = ',G12.5,/1X,70('*')/) 
 1410 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'Cells adjacent to cyclic boundaries must be of same size:',/&
         'DY(JMIN1) = ',G12.5,'     DY(JMAX1) = ',G12.5,/1X,70('*')/) 
 1420 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'Cells adjacent to cyclic boundaries must be of same size:',/&
         'DZ(KMIN1) = ',G12.5,'     DZ(KMAX1) = ',G12.5,/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_03 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 375 Added new check for IJKsize3_all() > DIMENSION_3=DIMENSION_3L
!// 990 Replace STOP with exitMPI
