!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: run                                                    C
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
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE run

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use param
      Use param1
!-----------------------------------------------


! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER*60     RUN_NAME

! Brief description of the problem.
      CHARACTER*60     DESCRIPTION

! Units for data input and output: CGS.
      CHARACTER*16     UNITS

! Type of run: NEW, RESTART
      CHARACTER*16     RUN_TYPE

! Variable which triggers automatic restart
      LOGICAL          AUTOMATIC_RESTART

! counter to keep track of how many auto_retart were performed
      INTEGER          ITER_RESTART

! computer node name/id
      CHARACTER*64     ID_NODE

! version.release of software
      CHARACTER*10     ID_VERSION

! Start-time of the run.
      DOUBLE PRECISION TIME

! Stop-time of the run.
      DOUBLE PRECISION TSTOP

! Time step.
      DOUBLE PRECISION DT

! 1./Time step.
      DOUBLE PRECISION oDT

! Number of times steps completed.
      INTEGER          NSTEP

! AE 091501: Declare a new variable to use on CN with RESTART cases
!            Number of time steps when restart file was read
      INTEGER          NSTEPRST      

! Discretization scheme for different equations
      INTEGER          DISCRETIZE(9)
      
! Use Chi scheme for discretizing certain equation sets
!  (species mass fractions)
      LOGICAL          Chi_scheme

! RUN ID info
      INTEGER          ID_MONTH
      INTEGER          ID_DAY
      INTEGER          ID_YEAR
      INTEGER          ID_HOUR
      INTEGER          ID_MINUTE
      INTEGER          ID_SECOND

! If .TRUE. solve X momentum equations
      LOGICAL          MOMENTUM_X_EQ(0:DIM_M)

! If .TRUE. solve Y momentum equations
      LOGICAL          MOMENTUM_Y_EQ(0:DIM_M)

! If .TRUE. solve Z momentum equations
      LOGICAL          MOMENTUM_Z_EQ(0:DIM_M)

! If .TRUE. solve energy equations
      LOGICAL          ENERGY_EQ

! If .TRUE. use the deferred correction method
      LOGICAL          DEF_COR

! If .TRUE. use the fourth order interpolation
      LOGICAL          FPFOI

! If .TRUE. solve granular energy equations
      LOGICAL          GRANULAR_ENERGY

! If .TRUE. solve K_Epsilon turbulence eq.
      LOGICAL          K_Epsilon

! If .TRUE. include added (virtual) mass in momentum eq.
      LOGICAL          Added_Mass

! phase number where added mass is applied.
      INTEGER          M_AM

! If .TRUE. use Simonin model (k_epsilon will
! automatically be set to true in check_data_02.
      LOGICAL          SIMONIN

! If .TRUE. use Ahmadi model (k_epsilon will
! automatically be set to true in check_data_02.
      LOGICAL          AHMADI

! If .TRUE. use Jenkins small friction BC
      LOGICAL          JENKINS
!                      If .TRUE. use revised phip for JJ BC 
      LOGICAL		BC_JJ_M
!                      If .TRUE. output PHIP to JJ_PHIP.dat
      LOGICAL		PHIP_OUT_JJ
!  			used to write specularity
      INTEGER		PHIP_OUT_ITER

! If .TRUE. use Yu and Standish correlation to
! compute ep_star
      LOGICAL          YU_STANDISH

! If .TRUE. use Fedors and Landel correlation to
! compute ep_star
      LOGICAL          FEDORS_LANDEL

! If .TRUE. solve species balance equations
      LOGICAL          SPECIES_EQ(0:DIM_M)

! If .TRUE. one of the species equations is being solved
      LOGICAL          ANY_SPECIES_EQ

! If .TRUE. call user-defined subroutines
      LOGICAL          CALL_USR

!  If .TRUE. solve population balance  equations
      LOGICAL          Call_DQMOM

! If .TRUE. use Model-B momentum equations
      LOGICAL          Model_B

! Drag model: see drag_gs for full options 
! syam_obrien (default & Umf corrected if C(2) and C(3) are defined), 
! gidaspow, Wen_Yu, koch_hill, BVK, HYS
      CHARACTER(64)    DRAG_TYPE
      
! Parameter used to calculate lubrication interactions between 
! different particles in HYS drag model
      DOUBLE PRECISION LAM_HYS

! anuj 4/20:  If .TRUE. calculate frictional stress terms
      LOGICAL          FRICTION
! Form of friction model: 
!             If 0: use S:S
!             If 1: use the form of Savage to compute S:S
!             If 2: use combination of both for frictional stress terms
      INTEGER          SAVAGE

! sof (02/16/2005):
! use Scheffer frictional stress (default set to .TRUE.)
      LOGICAL          SCHAEFFER

! sp (02/8/2006): blending frictional/kinetic stresses 	
! (default set to .FALSE. do not blend)
      LOGICAL          BLENDING_STRESS
      LOGICAL          TANH_BLEND ! default set to true
      LOGICAL          SIGM_BLEND ! default set to false

! loezos:  If .TRUE. treat system as if shearing
      LOGICAL          SHEAR
! Shear Vel
      DOUBLE PRECISION V_sh

! AE 04/16/01: If .TRUE. activate 2nd order accurate time implementation
      LOGICAL          CN_ON

! AEOLUS:
! AE: STOP Trigger mechanism to terminate MFIX normally before batch 
! queue terminates flag variable to check for end of batch queue when 
! set to TRUE check performed at the beginning of each time step and
! termination of mfix triggered after saving all files if condition
! is met
      LOGICAL          CHK_BATCHQ_END
! variable to store the total wall clock duration of the batch queue 
! session wall clock time specified in seconds
! for jaguarcnl@NCCS max wall clock limit is 2.5 hr limit up to 512
! processors
      DOUBLE PRECISION BATCH_WALLCLOCK
! variable to set a buffer time before the batch queue session ends to
! make sure once MFIX is triggered to shutdown, there is sufficient 
! time to save files, make copies to HPSS storage before batch queue
! time runs out. Current logic in MFIX checks for:
!    if CPU_TIME > (BATCH_WALLCLOCK - TERM_BUFFER) then 
!    save all .RES .SP files and trigger shutdown      
      DOUBLE PRECISION TERM_BUFFER

! sof (5/24/2005): added in case a user wants the code to automatically 
! restart for DT < DT_MIN
      LOGICAL          AUTO_RESTART

! parameters for dynamically adjusting time step
! +1 -> increase dt; -1 decrease dt
      INTEGER          DT_dir

! number of steps since last adjustment
      INTEGER          STEPS_tot

! number of iterations since last adjustment
      INTEGER          NIT_tot

! iterations per second for last dt
      DOUBLE PRECISION NITos

! Maximum Time step.
      DOUBLE PRECISION DT_MAX

! Minimum Time step.
      DOUBLE PRECISION DT_MIN

! Time step adjustment factor (<1.0)
      DOUBLE PRECISION DT_FAC

! The previous time step used in iterate (before it is
! changed by adjust_dt)
      DOUBLE PRECISION DT_prev

! Slope limiter parameter (0 < C _FAC <= 1.0)
      DOUBLE PRECISION C_FAC

! If .TRUE. reduce time step when residuals do not decrease
      LOGICAL          DETECT_STALL

! String which controls reduction of global sums for residual
! calculations
      LOGICAL          DEBUG_RESID        

       common /run_dp/ time      !for Linux

! kinetic theory model: see calc_mu_s for details
! for m > 1 option is IA_nonep, GHD,
! for m = 1 option is GD_99
      CHARACTER(64)    KT_TYPE

! radial distribution function options: see g_0 for details
! for m > 1 options are lebowitz, modified_lebowitz,
! mansoori, modified_mansoori.  default = lebowitz
! for m = 1 then carnahan and starling rdf used
      CHARACTER(64)    RDF_TYPE 

! Flag to invoke the variable solids density model.
      LOGICAL          SOLID_RO_V
!end

      END MODULE RUN
