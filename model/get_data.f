!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_DATA                                               C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Put version.release in LOG file                            C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add CALL to SET_L_scale                                    C
!  Author: W. Sams                                    Date: 10-MAY-94  C
!  Review:                                            Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_NAME, RUN_TYPE, ID_VERSION, ID_NODE       C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_DATA 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
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
 
! 
!                    Shift or not shift DX, DY and DZ values 
      LOGICAL        SHIFT 
!-----------------------------------------------
!
!
!
!
! This module call routines to initialize the namelist variables,
! read in the namelist variables from the ascii input file,
! checks that the input is valid, and opens the files.
!
      CALL INIT_NAMELIST 
      CALL READ_NAMELIST (0) 
!
!

      CALL ALLOCATE_ARRAYS
      
!
!
      IF (RUN_NAME == UNDEFINED_C) THEN 
         WRITE (*, 1000) 
         CALL MFIX_EXIT 
      ENDIF 
!
! open files
!
      IF (RUN_TYPE == 'RESTART_3') THEN 
         RUN_TYPE = 'RESTART_1' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0 
         CALL READ_NAMELIST (0)                  ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_1' 
         SHIFT = .FALSE. 
      ELSE IF (RUN_TYPE == 'RESTART_4') THEN 
         RUN_TYPE = 'RESTART_2' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0 
         CALL READ_NAMELIST (0)                  ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_2' 
         SHIFT = .FALSE. 
      ELSE 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         SHIFT = .TRUE. 
      ENDIF 
!
! write header in the .LOG file
!
      CALL WRITE_HEADER 
!
!  Check data and do some preliminary computations
!
      CALL START_LOG 
!
      CALL CHECK_DATA_01                         ! run_control input 
      CALL CHECK_DATA_02                         ! output_control input 
      CALL CHECK_DATA_03 (SHIFT)                 ! geometry input 
!
!  Set X, X_E, oX, oX_E ... etc.
!
      CALL SET_GEOMETRY 
!
      CALL SET_L_SCALE 
!
      CALL CHECK_DATA_04                         ! solid phase section 
      CALL CHECK_DATA_05                         ! gas phase section 
!
!  Set constants
!
      CALL SET_CONSTANTS 
!
      CALL CHECK_DATA_06                         ! initial condition section 
      CALL CHECK_DATA_07                         ! boundary condition section 
      CALL CHECK_DATA_08                         ! Internal surfaces section 
      CALL CHECK_DATA_09                         ! Chemical reactions section 
!
! close .LOG file
!
      CALL END_LOG 
!
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: GET_DATA.',/' Message: ',&
         'RUN_NAME not specified in mfix.dat',/1X,70('*')/) 
      END SUBROUTINE GET_DATA 


