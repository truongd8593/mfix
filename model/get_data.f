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
!
!  Revision Number:                                                    C
!  Purpose:  call GRIDMAP_INIT to handle the domain decomposition and  C
!            arrangement of all appropriate indices.                   C
!            Introduced MPI_Barrier for RESTART file read situation    C
!            (see comment !//S 0815 )
!
!  Author:   Aeolus Res. Inc.                         Date: 04-SEP-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      USE compar      
      USE gridmap     !// 300 Domain decomposition partioner

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

      CALL CHECK_DATA_00

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): reached end of read_namelist')") myPE	!//AIKEPARDBG
!      call mfix_exit(myPE)	!//AIKEPARDBG


!// 300 0905 Partition the domain and set indices
      call SET_MAX2
      call GRIDMAP_INIT
!//? 0906 GRIDMAP_INIT gives following error:
!//?          Assertion error: ** sendrecv_init: invalid kk  9
!
!//? before allocate, must do something with KMAX,DZ(),ZLENGTH for each PEs and also Global values
      CALL ALLOCATE_ARRAYS


!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft ALLOCATE_ARRAYS in get_data')") myPE	!//AIKEPARDBG
!      call mfix_exit(myPE)	!//AIKEPARDBG
      
!
!
      IF (RUN_NAME == UNDEFINED_C) THEN 
         WRITE (*, 1000) 
         CALL MFIX_EXIT 
      ENDIF 

!
! open files
!

!//S 0815 May be worth to replace following MPI_Barrier calls with a wrapper
!//S      equivalent to avoid direct inclusion of MPI library thru out the source code

      IF (RUN_TYPE == 'RESTART_3') THEN 
         RUN_TYPE = 'RESTART_1' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0

!         write(*,"('(PE ',I2,'): aft READ_RES0')") myPE	!//AIKEPARDBG
         call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

         CALL READ_NAMELIST (0)                  ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_1' 
         SHIFT = .FALSE. 
      ELSE IF (RUN_TYPE == 'RESTART_4') THEN 
         RUN_TYPE = 'RESTART_2' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 

         CALL READ_RES0  
         call MPI_Barrier(MPI_COMM_WORLD,mpierr)  !//PAR_I/O enforce barrier here

         CALL READ_NAMELIST (0)              ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_2' 
         SHIFT = .FALSE. 
      ELSE 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         SHIFT = .TRUE. 
      ENDIF 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): end of RUN_TYPE branches in get_data')") myPE	!//AIKEPARDBG
!      call mfix_exit(myPE)	!//AIKEPARDBG

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
!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft call chk_data_01 in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG


      CALL CHECK_DATA_02                         ! output_control input 
!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): aft call chk_data_02 in get_data')") myPE  !//AIKEPARDBG
!      call mfix_exit(myPE)  !//AIKEPARDBG

      CALL CHECK_DATA_03 (SHIFT)                 ! geometry input 
!//AIKEPARDBGSTOP 0907
!      write(*,"('(PE ',I2,'): aft call chk_data_03 in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

!
!  Set X, X_E, oX, oX_E ... etc.
!
      CALL SET_GEOMETRY 
!
      CALL SET_L_SCALE 
!
!//AIKEPARDBGSTOP 0920
!      write(*,"('(PE ',I2,'): aft call set_L_scale in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

      CALL CHECK_DATA_04                         ! solid phase section 
      CALL CHECK_DATA_05                         ! gas phase section 
!
!  Set constants
!
      CALL SET_CONSTANTS 
!
!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft call set_constants in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG


      CALL CHECK_DATA_06                         ! initial condition section 


!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft call check_data_06 in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

      CALL CHECK_DATA_07                         ! boundary condition section 
!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft call check_data_07 in get_data')") myPE !//AIKEPARDBG
      CALL CHECK_DATA_08                         ! Internal surfaces section 
      CALL CHECK_DATA_09                         ! Chemical reactions section 
!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft call check_data_09 in get_data')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG
      

!
! close .LOG file
!
      CALL END_LOG 
!
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: GET_DATA.',/' Message: ',&
         'RUN_NAME not specified in mfix.dat',/1X,70('*')/) 
      END SUBROUTINE GET_DATA 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_00                                          C
!  Purpose: check the distributed parallel namelist variables          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 14-DEC-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:  NODESI , NODESJ , NODESK                     C
!                         DT, RUN_TYPE                                 C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_00
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param1 
      USE compar
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
!-----------------------------------------------
!
      IF (NODESI .EQ. UNDEFINED_I) THEN
          WRITE (*,*) ' NODESI not found in MFIX.DAT'
          CALL MFIX_EXIT
      END IF
!
      IF (NODESJ .EQ. UNDEFINED_I) THEN
          WRITE (*,*) ' NODESJ not found in MFIX.DAT'
          CALL MFIX_EXIT
      END IF
!
      IF (NODESK .EQ. UNDEFINED_I) THEN
          WRITE (*,*) ' NODESK not found in MFIX.DAT'
          CALL MFIX_EXIT
      END IF
!
!
      RETURN  
      END SUBROUTINE CHECK_DATA_00
 
