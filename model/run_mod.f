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
!                      Discretization scheme for different equations
      INTEGER          DISCRETIZE(9)
      
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
!                      If .TRUE. solve species balance equations
      LOGICAL          SPECIES_EQ(0:DIM_M)
!
!                      If .TRUE. one of the species equations is being solved
      LOGICAL          ANY_SPECIES_EQ
!
!                      If .TRUE. call user-defined subroutines
      LOGICAL          CALL_USR
!
!                      If .TRUE. use Model-B momentum equations
      LOGICAL          Model_B
!
! start anuj 4/20
!	               If. TRUE. calculate frictional stress terms
      LOGICAL	       FRICTION
!
 
!
! start loezos
!	               If. TRUE. calculate frictional stress terms
      LOGICAL	       SHEAR


!		       If 0: use S:S
!	               If 1: use the form of Savage to compute S:S
!		       If 2: use combination of both
!	               for frictional stress terms
      INTEGER	       SAVAGE		
! end anuj 4/20
 
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
!

! loezos                     Shear Vel
      DOUBLE PRECISION V_sh
! loezos
         common /run_dp/ time      !for Linux

      END MODULE run                                                                             
