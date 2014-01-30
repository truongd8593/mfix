!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GET_DATA                                                C
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
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:  call GRIDMAP_INIT to handle the domain decomposition and  C
!            arrangement of all appropriate indices.                   C
!            Introduced MPI_Barrier for RESTART file read situation    C
!  Author:   Aeolus Res. Inc.                         Date: 04-SEP-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!  kinetic equation                                                    C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!								       C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_NAME, RUN_TYPE, ID_VERSION, ID_NODE       C
!  Variables modified: None                                            C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GET_DATA 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE funits 
      USE compar      
      USE gridmap
      USE discretelement
      USE des_thermo
      USE des_rxns
      USE leqsol
      USE parallel
      USE qmom_kinetic_equation
      USE mfix_pic
      USE cutcell
      USE dashboard

      USE error_manager


      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! shift DX, DY and DZ values 
      LOGICAL :: SHIFT 
      LOGICAL :: CYCLIC_X_BAK, CYCLIC_Y_BAK, CYCLIC_Z_BAK, &
                 CYLINDRICAL_BAK
!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE
!-----------------------------------------------


! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST 
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0) 

! Initialize the error manager. This call occurs after the mfix.dat
! is read so that message verbosity can be set and the .LOG file 
! can be opened.
      CALL INIT_ERROR_MANAGER

! Check the minimum geometry requirements.
      CALL CHECK_GEOMETRY_PREREQS

! Check the minimum dmp requirements.
      CALL CHECK_DMP_PREREQS

! Set up the physical domain decomposition.
      CALL SET_MAX2

! Set constants
      CALL SET_CONSTANTS 

! Partition the domain and set indices
      CALL GRIDMAP_INIT

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS


!
      CALL ALLOCATE_ARRAYS    

! Alberto Passalacqua - QMOMK
      IF (QMOMK) THEN
         CALL QMOMK_ALLOCATE_ARRAYS
      ENDIF



! open files
      IF (RUN_TYPE == 'RESTART_3') THEN 
         RUN_TYPE = 'RESTART_1' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  

         CALL READ_NAMELIST (0)     ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_1' 
         SHIFT = .FALSE. 
      ELSE IF (RUN_TYPE == 'RESTART_4') THEN 
         RUN_TYPE = 'RESTART_2' 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         CALL READ_RES0  
!        call MPI_Barrier(MPI_COMM_WORLD,mpierr)  

         CALL READ_NAMELIST (0)     ! to modify the .RES data with .DAT data 
         RUN_TYPE = 'RESTART_2' 
         SHIFT = .FALSE. 
      ELSE 
         CALL OPEN_FILES (RUN_NAME, RUN_TYPE, N_SPX) 
         SHIFT = .TRUE. 
      ENDIF 


! write header in the .LOG file
      CALL WRITE_HEADER 


! Check data and do some preliminary computations
      CALL START_LOG 

      CALL CHECK_DATA_01                         ! run_control input 
      CALL CHECK_DATA_02                         ! output_control input 
      CALL CHECK_DATA_03 (SHIFT)                 ! geometry input 

! Set X, X_E, oX, oX_E ... etc.
      CALL SET_GEOMETRY 

      CALL SET_L_SCALE 

      CALL CHECK_DATA_04                         ! solid phase section 
      CALL CHECK_DATA_05                         ! gas phase section 
      CALL CHECK_DATA_06                         ! initial condition section 
!      IF(.NOT.CARTESIAN_GRID) CALL CHECK_DATA_07                         ! boundary condition section 
      CALL CHECK_DATA_07                         ! boundary condition section 
      CALL CHECK_DATA_08                         ! Internal surfaces section 
      CALL CHECK_DATA_09                         ! Chemical reactions section       
      CALL CHECK_DATA_10

      CALL CHECK_DATA_ODEPACK                    ! Stiff Chemistry Solver

! close .LOG file
      CALL END_LOG 

!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER
            
      DTMIN =  LARGE_NUMBER                                                                                                                        
      DTMAX = -LARGE_NUMBER
             
      NIT_MIN = MAX_NIT
      NIT_MAX = 0
             
      N_DASHBOARD = 0


      IF (OPT_PARALLEL) THEN
         IS_SERIAL = .FALSE.
         DO_TRANSPOSE = .FALSE.
         MINIMIZE_DOTPRODUCTS = .TRUE.
         SOLVER_STATISTICS = .TRUE.
         DEBUG_RESID = .FALSE.
         LEQ_SWEEP(1:2) = 'ASAS'
         LEQ_METHOD(1:2) = 2
         LEQ_METHOD(3:9) = 1
      ENDIF


      RETURN  

 1000 FORMAT(/1X,70('*')//' From: GET_DATA.',/' Message: ',&
         'RUN_NAME not specified in mfix.dat',/1X,70('*')/)  

      END SUBROUTINE GET_DATA
