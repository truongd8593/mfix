CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: run_control.inc                                        C
C  Purpose: Common block containing run control data                   C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision: 1                                                         C
C  Purpose: add node name,version                                      C
C  Author: P.Nicoletti                                Date: 07-FEB-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References: None                                C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: None                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C     Run Control
C
C                      Main filename to be used for output files  Name must
C                      still be legal after extensions are added to it.
      CHARACTER*60     RUN_NAME
C
C                      Brief description of the problem.
      CHARACTER*60     DESCRIPTION
C
C                      Units for data input and output: CGS.
      CHARACTER*16     UNITS
C
C                      Type of run: NEW, RESTART
      CHARACTER*16     RUN_TYPE
C
C                      computer node name/id
      CHARACTER*64     ID_NODE
C
C                      version.release of software
      CHARACTER*10     ID_VERSION
C
C                      Start-time of the run.
      DOUBLE PRECISION TIME
C
C                      Stop-time of the run.
      DOUBLE PRECISION TSTOP
C
C                      Time step.
      DOUBLE PRECISION DT
C
C                      1./Time step.
      DOUBLE PRECISION oDT
C
C                      Number of times steps completed.
      INTEGER          NSTEP
C
C                      Discretization scheme for different equations
      INTEGER          DISCRETIZE(8)
C
C                      RUN ID info
C
      INTEGER          ID_MONTH
      INTEGER          ID_DAY
      INTEGER          ID_YEAR
      INTEGER          ID_HOUR
      INTEGER          ID_MINUTE
      INTEGER          ID_SECOND
C
C                      If .TRUE. solve X momentum equations
      LOGICAL          MOMENTUM_X_EQ(0:DIMENSION_M)
C
C                      If .TRUE. solve Y momentum equations
      LOGICAL          MOMENTUM_Y_EQ(0:DIMENSION_M)
C
C                      If .TRUE. solve Z momentum equations
      LOGICAL          MOMENTUM_Z_EQ(0:DIMENSION_M)
C
C                      If .TRUE. solve energy equations
      LOGICAL          ENERGY_EQ
C
C                      If .TRUE. solve granular energy equations
      LOGICAL          GRANULAR_ENERGY
C
C                      If .TRUE. solve species balance equations
      LOGICAL          SPECIES_EQ(0:DIMENSION_M)
C
C                      If .TRUE. one of the species equations is being solved
      LOGICAL          ANY_SPECIES_EQ
C
C                      If .TRUE. call user-defined subroutines
      LOGICAL          CALL_USR
C
C                      If .TRUE. use Model-B momentum equations
      LOGICAL          Model_B
C
C      parameters for dynamically adjusting time step

C
C                      +1 -> increase dt; -1 decrease dt
      INTEGER          DT_dir     
C
C                      number of steps since last adjustment
      INTEGER          STEPS_tot 
C
C                      number of iterations since last adjustment
      INTEGER          NIT_tot   
C
C                      iterations per second for last dt
      DOUBLE PRECISION NITos      
C
C                      Maximum Time step.
      DOUBLE PRECISION DT_MAX
C
C                      Minimum Time step.
      DOUBLE PRECISION DT_MIN
C
C                      Time step adjustment factor (<1.0)
      DOUBLE PRECISION DT_FAC
C
C                      If .TRUE. reduce time step when residuals do not decrease
      LOGICAL          DETECT_STALL
C
C
      COMMON / RUN_CONTROL_DP /
     &       TIME, TSTOP, DT, oDT, NITos, DT_MAX, DT_MIN, DT_FAC
C
      COMMON / RUN_CONTROL_I /
     &       NSTEP, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR,
     &       ID_MINUTE , ID_SECOND, DT_dir, STEPS_tot, NIT_tot,
     &       DISCRETIZE
C
      COMMON / RUN_CONTROL_C /
     &       RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, ID_NODE, ID_VERSION
C
      COMMON / RUN_CONTROL_L /
     &       ENERGY_EQ, SPECIES_EQ, ANY_SPECIES_EQ, CALL_USR, Model_B,
     &       MOMENTUM_X_EQ, MOMENTUM_Y_EQ, MOMENTUM_Z_EQ, DETECT_STALL,
     &       GRANULAR_ENERGY
