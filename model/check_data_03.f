!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_03                                           C
!  Purpose: check the geometry and discretization namelist section     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose: Replaced STOP with mfix_exit() to abort all processors     C
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
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_03(SHIFT) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE bc
      USE funits 
      USE compar
      USE mpi_utility
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! shift dx, dy and dz values (true only for new and restart_1 runs)
      LOGICAL, INTENT(IN) :: SHIFT 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      LOGICAL :: COMPARE 
!-----------------------------------------------






! CHECK THE TOTAL DIMENSION
      IF (IJKMAX3 > DIMENSION_3G) THEN 
         CALL ERROR_ROUTINE ('check_data_03', &
            'global dimension error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            '(IMAX+2+1)*(JMAX+2+1)*(KMAX+2+1) = ', IJKMAX3 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'DIMENSION_3 in allocate_arrays = ', DIMENSION_3G 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'modify allocate_arrays and recompile .... or' 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'change dimensions in mfix.dat', &
            ' ... whichever is appropriate' 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 

      IF (IJKsize3_all(myPE) > DIMENSION_3) THEN 
         CALL ERROR_ROUTINE ('check_data_03', &
            'subdomain dimension error', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            '(Iend3-Istart3+1)*(Jend3-Jstart3+1)*',&
            '(Kend3-Kstart3+1) = ',IJKsize3 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'DIMENSION_3 in allocate_arrays = ', DIMENSION_3 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'modify allocated_arrays and recompile .... or' 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) &
            'change dimensions in mfix.dat', &
            ' ... whichever is appropriate' 
         CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
      ENDIF 
      
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
         'DZ(1) and ZLENGTH are not equal!  ',&
         '(Recall NO_K=.true.)' ,70('*')/) 
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'IMAX should be 1, since NO_I is true',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'JMAX should be 1, since NO_J is true',/1X,70('*')/) 

 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'KMAX should be 1, since NO_K is true',/1X,70('*')/) 

 1250 FORMAT(/1X,70('*')//' From: CHECK_DATA_03',/' Message: ',&
         'Cyclic bc for X not allowed in cylindrical coordinates',&
         /1X,70('*')/) 
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

