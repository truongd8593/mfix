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


! This module call routines to initialize the namelist variables,
! read in the namelist variables from the ascii input file,
! checks that the input is valid, and opens the files.
      CALL INIT_NAMELIST 
      CALL READ_NAMELIST (0) 

! Basic geometry checks. Needed before DMP decomp.
      CALL CHECK_DATA_00

! Set constants
      CALL SET_CONSTANTS 

      CYCLIC_X_BAK = CYCLIC_X
      CYCLIC_Y_BAK = CYCLIC_Y
      CYCLIC_Z_BAK = CYCLIC_Z
      CYLINDRICAL_BAK = CYLINDRICAL

      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE.
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE.
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE.
      
      DO_K = .NOT.NO_K

      IF (COORDINATES == 'CYLINDRICAL') CYLINDRICAL = .TRUE.

      IF (CYLINDRICAL .AND. COMPARE(ZLENGTH,8.D0*ATAN(ONE)) .AND. DO_K) &
         CYCLIC_Z = .TRUE.

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

      IF(DISCRETE_ELEMENT) THEN
! make these checks here before the other continuum check routines are called.      
         IF(DES_CONTINUUM_HYBRID.AND.MPPIC) then 
            if(myPE.eq.pe_IO) WRITE(*,'(3(/,2x,A))') 'DISCRETE ELEMENT IS TRUE', &
            'BUT BOTH DES_CONTINUUM_HYBRID AND MPPIC CANNOT BE SET TO TRUE', &
            'Problem in input file; TERMINAL ERROR: exiting the simulation'
            
            if(DMP_LOG) WRITE(UNIT_LOG,'(3(/,2x,A))') 'DISCRETE ELEMENT IS TRUE', &
            'BUT BOTH DES_CONTINUUM_HYBRID AND MPPIC CANNOT BE SET TO TRUE', &
            'Problem in input file; TERMINAL ERROR: exiting the simulation'
            
            CALL MFIX_EXIT(myPE) 
         ENDIF

      ELSE
         
! If discrete_element is .false. then overwrite the following user 
! DES logicals which may be set to true in the input file.  Only need
! to set those that may impact continuum aspects
! This is placed here (after reading of namelists) to ensure that
! correct dem related flags are set when other continuum check_data
! routines (such as check_data_04, ....06, ...07, ...20) are called          
         DES_CONTINUUM_COUPLED = .FALSE.
         DES_INTERP_ON = .FALSE.
         DES_CONTINUUM_HYBRID = .FALSE.
         TSUJI_DRAG = .FALSE.
         PRINT_DES_DATA = .FALSE.
         MPPIC = .FALSE. 
         DES_ONEWAY_COUPLED = .false. 
         DES_CONV_EQ = .FALSE. ! No convection - ref by SOLVE_ENERGY_EQ
         ANY_DES_SPECIES_EQ = .FALSE.
      ENDIF

! Partition the domain and set indices
      call SET_MAX2
      call GRIDMAP_INIT


! Copying back logical variables to original values to retain the structure of the code
! For cylindrical case making CYCLIC_Z = .TRUE., makes CYCLIC TRUE in set_geometry and
! which in turn leads to fixing the pressure in source_pp_g.  (SP)
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


      CALL ALLOCATE_ARRAYS    

! Alberto Passalacqua - QMOMK
      IF (QMOMK) THEN
         CALL QMOMK_ALLOCATE_ARRAYS
      ENDIF


      IF (RUN_NAME == UNDEFINED_C) THEN 
         WRITE (*, 1000) 
         CALL MFIX_EXIT(myPE) 
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


      RETURN  

 1000 FORMAT(/1X,70('*')//' From: GET_DATA.',/' Message: ',&
         'RUN_NAME not specified in mfix.dat',/1X,70('*')/)  

      END SUBROUTINE GET_DATA 
      

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_DATA_00                                           !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 20-Sep-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DATA_00


! Global Variables:
!---------------------------------------------------------------------//
! This processes rank.
      use compar, only: myPE
! Number of ranks.
      use compar, only: numPEs
! DMP grid partitioning data:
      use compar, only: NODESI  ! Partitions along x-axis
      use compar, only: NODESJ  ! Partitions along y-axis
      use compar, only: NODESK  ! Partitions along z-axis
! Flag: All ranks report errors.
      use output, only: ENABLE_DMP_LOG
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! Domain partitions in various directions.
      use geometry, only: IMAX 
      use geometry, only: JMAX
      use geometry, only: KMAX
! Runtime flag specifying 2D simulations
      use geometry, only: NO_K

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED_I

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! The approximate number of domain partitions that will be assigned
! to each process for DMP runs. (eg IMAX/NODESI)
      INTEGER :: LAYERS

! Local Parameters:
!---------------------------------------------------------------------//
! The minimum number of computational cell layers requried.
      INTEGER, PARAMETER :: DMP_MIN = 4
      

! This turns on error messaging from all processes.
      IF(ENABLE_DMP_LOG) DMP_LOG = .TRUE.

! Verify that the domain decomposition was specifed.
      IF(IMAX == UNDEFINED_I .OR. JMAX == UNDEFINED_I .OR. &
         (.NOT.NO_K .AND. KMAX == UNDEFINED_I) ) THEN
         IF(DMP_LOG) WRITE(*,1002)
         CALL MFIX_EXIT(myPE) 
      ENDIF

! Verify that DMP partitioning information is provided given if
! there is more than one rank.
      IF( numPEs > 1 ) then
         IF(NODESI .EQ. UNDEFINED_I .AND.                              &
            NODESJ .EQ. UNDEFINED_I .AND.                              &
            NODESK .EQ. UNDEFINED_I) THEN
            IF(DMP_LOG) THEN
! This message is only passed to standard out because the unit log file
! has not been opened at this point.
               WRITE(*,1000); WRITE(*,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

! Initialize NODE values if undefined. If this is a DMP run, then a 
! warning message is passed to the user.
      IF (NODESI .EQ. UNDEFINED_I) THEN
         IF((numPEs > 1) .AND. DMP_LOG) WRITE(*,1001)'I','I'
         NODESI = 1

! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESI > 1) THEN
         LAYERS=int(IMAX/NODESI)
         IF(LAYERS < DMP_MIN) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1003) 'X', DMP_MIN, 'I', 'I', LAYERS
               WRITE(*,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF (NODESJ .EQ. UNDEFINED_I) THEN
         IF((numPEs > 1) .AND. DMP_LOG) WRITE(*,1001)'J','J'
         NODESJ = 1

! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESJ > 1) THEN
         LAYERS=int(JMAX/NODESJ)
         IF(LAYERS < DMP_MIN) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1003) 'Y', DMP_MIN, 'J', 'J', LAYERS
               WRITE(*,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      IF (NODESK .EQ. UNDEFINED_I) THEN
         IF((numPEs > 1) .AND. DMP_LOG) WRITE(*,1001)'K','K'
         NODESK = 1
! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESK > 1) THEN
         LAYERS=int(KMAX/NODESK)
         IF(LAYERS < DMP_MIN) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1003) 'Z', DMP_MIN, 'K', 'K', LAYERS
               WRITE(*,9999)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN  

 1000 FORMAT(//1X,70('*'),/' From: CHECK_DATA_00',/' Error 1000:',     &
         '  No DMP grid partitioning data provided in mfix.dat. ',     &
         '  NODESI, NODESJ, and NODESK are all undefined.')

 1001 FORMAT(/' From: CHECK_DATA_00',/' Warning 1001:',                &
         ' Setting NODES',A1,' to default: NODES',A1,'=1.',/)

 1002 FORMAT(//1X,70('*'),/' From: CHECK_DATA_00',/' Error 1002:',     &
         'IMAX or JMAX or KMAX not specified in mfix.dat',/1X,70('*')/)

 1003 FORMAT(//1X,70('*'),/' From: CHECK_DATA_00',/' Error 1003:',     &
         ' Too many DMP partitions partitions specified for ',A1,      &
         ' axis.',/' There must be at least ',I2,' computational',     &
         ' cells per DMP parition.',/' >>> Computational Cells/DMP',   &
         ' Partition = int(',A1,'MAX/NODES',A1,' = ',I2)


 9999 FORMAT(/' Please refer to the Readme file on the required input',&
         ' and make',/' the necessary corrections to the data file.',  &
         /1X,70('*')//)

      END SUBROUTINE CHECK_DATA_00
