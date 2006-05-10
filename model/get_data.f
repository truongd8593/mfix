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
!            (see comment !//S  )
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
      USE gridmap
      USE discretelement
      
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
      LOGICAL        CYCLIC_X_BAK, CYCLIC_Y_BAK, CYCLIC_Z_BAK, CYLINDRICAL_BAK
!-----------------------------------------------
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE
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

      CYCLIC_X_BAK = CYCLIC_X
      CYCLIC_Y_BAK = CYCLIC_Y
      CYCLIC_Z_BAK = CYCLIC_Z
      CYLINDRICAL_BAK = CYLINDRICAL

      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE.
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE.
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE.
      
      DO_K = .NOT.NO_K
!
      IF (COORDINATES == 'CYLINDRICAL') CYLINDRICAL = .TRUE.

      IF (CYLINDRICAL .AND. COMPARE(ZLENGTH,8.D0*ATAN(ONE)) .AND. DO_K) &
         CYCLIC_Z = .TRUE.

!// Partition the domain and set indices
      call SET_MAX2
      call GRIDMAP_INIT

!\\SP Copying back logical variables to original values to retain the structure of the code
!     For cylindrical case making CYCLIC_Z = .TRUE., makes CYCLIC TRUE in set_geometry and
!     which in turn leads to fixing the pressure in source_pp_g.

      CYCLIC_X = CYCLIC_X_BAK
      CYCLIC_Y = CYCLIC_Y_BAK
      CYCLIC_Z = CYCLIC_Z_BAK
      CYLINDRICAL = CYLINDRICAL_BAK




!     write(*,*) 'ISTART-IEND'
!     write(*,*) ISTART3, ISTART2, ISTART, ISTART1, IEND1, IEND, IEND2, IEND3
!     write(*,*) 'JSTART-JEND'
!     write(*,*) JSTART3, JSTART2, JSTART, JSTART1, JEND1, JEND, JEND2, JEND3
!     write(*,*) 'KSTART-KEND'
!     write(*,*) KSTART3, KSTART2, KSTART, KSTART1, KEND1, KEND, KEND2, KEND3
!
! MUST use k-epsilon model and granular temperature PDE with Simonin
! and Ahmadi models (must be done before allocate_arrays). Sof --> 02/01/05
      IF(SIMONIN .OR. AHMADI) THEN
        K_Epsilon = .TRUE.
	GRANULAR_ENERGY = .TRUE.
      ENDIF

      CALL ALLOCATE_ARRAYS    

!     DES
      IF(NO_I.OR.NO_J.OR.NO_K) DIMN = 2       
      CALL DES_ALLOCATE_ARRAYS
      CALL DES_INIT_ARRAYS 
!
!
      IF (RUN_NAME == UNDEFINED_C) THEN 
         WRITE (*, 1000) 
         CALL MFIX_EXIT(myPE) 
      ENDIF 

!
! open files
!

      IF (RUN_TYPE == 'RESTART_3') THEN 
         RUN_TYPE = 'RESTART_1' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  

         CALL READ_NAMELIST (0)                  ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_1' 
         SHIFT = .FALSE. 
      ELSE IF (RUN_TYPE == 'RESTART_4') THEN 
         RUN_TYPE = 'RESTART_2' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 

         CALL READ_RES0  
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  

         CALL READ_NAMELIST (0)              ! to modify the .RES data with .DAT data 
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
!     CHEM & ISAT begin (nan xie)
!  check rxns
      CALL CHECK_DATA_CHEM      
!     CHEM & ISAT end (nan xie)
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
      if( numPEs > 1 ) then
        IF (NODESI .EQ. UNDEFINED_I .AND. NODESJ .EQ. UNDEFINED_I &
	    .AND. NODESK .EQ. UNDEFINED_I) THEN
          WRITE (*,*) ' No grid partitioning data (NODESI, NODESJ, or NODESK) in mfix.dat'
          CALL MFIX_EXIT(myPE)
        END IF
      endif
!
      IF (NODESI .EQ. UNDEFINED_I) NODESI = 1
      IF (NODESJ .EQ. UNDEFINED_I) NODESJ = 1
      IF (NODESK .EQ. UNDEFINED_I) NODESK = 1
!
!
      RETURN  
      END SUBROUTINE CHECK_DATA_00
 
!// Comments on the modifications for DMP version implementation      
!// 300 Partition the domain and set indices
!//PAR_I/O enforce barrier here
