!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: run_control.inc                                        C
!  Purpose: Common block containing run control data                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision: 1                                                         C
!  Purpose: add node name,version                                      C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      MODULE run


      Use param
      Use param1


!
!     Run Control
!
!                      Main filename to be used for output files  Name must
!                      still be legal after extensions are added to it.
      CHARACTER*60     RUN_NAME
!
!                      Brief description of the problem.
      CHARACTER*60     DESCRIPTION
!
!                      Units for data input and output: CGS.
      CHARACTER*16     UNITS
!
!                      Type of run: NEW, RESTART
      CHARACTER*16     RUN_TYPE
!
!                      Variable which triggers automatic restart
      LOGICAL          AUTOMATIC_RESTART
!
!                      counter to keep track of how many auto_retart were performed
      INTEGER          ITER_RESTART
!
!                      computer node name/id
      CHARACTER*64     ID_NODE
!
!                      version.release of software
      CHARACTER*10     ID_VERSION
!
!                      Start-time of the run.
      DOUBLE PRECISION TIME
!
!                      Stop-time of the run.
      DOUBLE PRECISION TSTOP
!
!                      Time step.
      DOUBLE PRECISION DT
!
!                      1./Time step.
      DOUBLE PRECISION oDT
!
!                      Number of times steps completed.
      INTEGER          NSTEP
!
!AE TIME 091501   Declare a new variable to use on CN with RESTART cases
!                      Number of time steps when restart file was read
      INTEGER          NSTEPRST      
!
!                      Discretization scheme for different equations
      INTEGER          DISCRETIZE(9)
      
!
!                      Use Chi scheme for discretizing certain equation sets
!                       (species mass fractions)
      LOGICAL          Chi_scheme
      
!
!                      RUN ID info
!
      INTEGER          ID_MONTH
      INTEGER          ID_DAY
      INTEGER          ID_YEAR
      INTEGER          ID_HOUR
      INTEGER          ID_MINUTE
      INTEGER          ID_SECOND
!
!                      If .TRUE. solve X momentum equations
      LOGICAL          MOMENTUM_X_EQ(0:DIM_M)
!
!                      If .TRUE. solve Y momentum equations
      LOGICAL          MOMENTUM_Y_EQ(0:DIM_M)
!
!                      If .TRUE. solve Z momentum equations
      LOGICAL          MOMENTUM_Z_EQ(0:DIM_M)
!
!                      If .TRUE. solve energy equations
      LOGICAL          ENERGY_EQ
!
!                      If .TRUE. use the deferred correction method
      LOGICAL          DEF_COR
!
!                      If .TRUE. use the fourth order interpolation
      LOGICAL          FPFOI
!
!                      If .TRUE. solve granular energy equations
      LOGICAL          GRANULAR_ENERGY
!
!                      If .TRUE. solve K_Epsilon turbulence eq.
      LOGICAL          K_Epsilon
!
!                      If .TRUE. use Simonin model (k_epsilon will
!                      automatically be set to true in check_data_02.
      LOGICAL          SIMONIN
!
!                      If .TRUE. use Ahmadi model (k_epsilon will
!                      automatically be set to true in check_data_02.
      LOGICAL          AHMADI
!
!                      If .TRUE. use Jenkins small friction BC
      LOGICAL          JENKINS
!
!                      If .TRUE. use Yu and Standish correlation to
!                      compute ep_star
      LOGICAL          YU_STANDISH
!
!                      If .TRUE. use Fedors and Landel correlation to
!                      compute ep_star
      LOGICAL          FEDORS_LANDEL
!
!                      If .TRUE. solve species balance equations
      LOGICAL          SPECIES_EQ(0:DIM_M)
!
!                      If .TRUE. one of the species equations is being solved
      LOGICAL          ANY_SPECIES_EQ
!
!                      If .TRUE. call user-defined subroutines
      LOGICAL          CALL_USR
!                       If .TRUE. solve population balance  equations
      LOGICAL          Call_DQMOM
!
!
!                      If .TRUE. use Model-B momentum equations
      LOGICAL          Model_B
!                      Options are syam_obrien (default), Umf corrected
!                      (if C(2) and C(3) are defined), gidaspow, Wen_Yu
!                      koch_hill
      CHARACTER(64)    DRAG_TYPE
      
!	Parameter used to calculate lubrication interactions between different
!	particles in HYS drag model
      LOGICAL		use_def_lam_HYS

      DOUBLE PRECISION	lam_HYS 
!
! start anuj 4/20
!	               If. TRUE. calculate frictional stress terms
      LOGICAL	       FRICTION
!
 
!
! start loezos
!	               If. TRUE. calculate frictional stress terms
      LOGICAL	       SHEAR

!AE TIME 041601        If .TRUE. activate 2nd order accurate time implementation
      LOGICAL          CN_ON

!AEOLUS STOP Trigger mechanism to terminate MFIX normally before batch queue terminates
      LOGICAL          CHK_BATCHQ_END
      DOUBLE PRECISION BATCH_WALLCLOCK
      DOUBLE PRECISION TERM_BUFFER



!		       If 0: use S:S
!	               If 1: use the form of Savage to compute S:S
!		       If 2: use combination of both
!	               for frictional stress terms
      INTEGER	       SAVAGE		
! end anuj 4/20
!
! sof: added SCHAEFFER keyword (02/16/2005)
!	               default set to .TRUE. use Scheffer frictional stress
      LOGICAL	       SCHAEFFER
! sof: end
! sp: added BLENDING_STRESS keyword (02/8/2006)
!	               default set to .FALSE. do not blend
!	               frictional/kinetic stresses
      LOGICAL	       BLENDING_STRESS
! sp: added TANH_BLEND & SIGM_BLEND keyword (10/24/2006)
      LOGICAL	       TANH_BLEND ! default set to true
      LOGICAL	       SIGM_BLEND ! default set to false
! sp: end
! sof: added in case a user wants the code to automatically restart for DT < DT_MIN
      LOGICAL	       AUTO_RESTART
! sof: end 5/24/2005
 
!      parameters for dynamically adjusting time step
 
!
!                      +1 -> increase dt; -1 decrease dt
      INTEGER          DT_dir
!
!                      number of steps since last adjustment
      INTEGER          STEPS_tot
!
!                      number of iterations since last adjustment
      INTEGER          NIT_tot
!
!                      iterations per second for last dt
      DOUBLE PRECISION NITos
!
!                      Maximum Time step.
      DOUBLE PRECISION DT_MAX
!
!                      Minimum Time step.
      DOUBLE PRECISION DT_MIN
!
!                      Time step adjustment factor (<1.0)
      DOUBLE PRECISION DT_FAC
!
!                      The previous time step used in iterate (before it is
!                      changed by adjust_dt)
      DOUBLE PRECISION DT_prev
!
!                      Slope limiter parameter (0 < C _FAC <= 1.0)
      DOUBLE PRECISION C_FAC
!
!                      If .TRUE. reduce time step when residuals do not decrease
      LOGICAL          DETECT_STALL
!
! String which controls reduction of global sums for residual calculations
      LOGICAL          DEBUG_RESID        
!
!

! loezos                     Shear Vel
      DOUBLE PRECISION V_sh
! loezos
         common /run_dp/ time      !for Linux
!
!     CHEM & ISAT begin (nan xie)
!                      If .TRUE. call ODEPACK
      LOGICAL          CALL_DI
!                      If .TRUE. calculate the change of diameter due to grow
      LOGICAL          CALL_GROW
!                      If .TRUE. call isat subroutines
      LOGICAL          CALL_ISAT
!                      time step for isat calculation
      DOUBLE PRECISION ISATdt
!     CHEM & ISAT end (nan xie)
!
!
!     JEG added
!                       for m > 1 option is IA_nonep.  
!                       for m = 1 option is GD_99
      CHARACTER(64)     KT_TYPE
!
!                       for m > 1 options are lebowitz, modified_lebowitz,
!                       mansoori, modified_mansoori.  default = lebowitz
!                       for m = 1 then carnahan and starling rdf used
      CHARACTER(64)     RDF_TYPE 
!
!     JEG end      

      END MODULE run                                                                             
