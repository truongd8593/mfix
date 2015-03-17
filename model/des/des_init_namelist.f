!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Added some interpolation based inputs                     C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_INIT_NAMELIST

      USE param1
      USE discretelement
      USE mfix_pic
      USE des_bc
      USE des_thermo
      USE des_rxns
      USE pic_bc
      USE particle_filter

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      INCLUDE 'desnamelist.inc'



!#####################################################################!
!                             Run Control                             !
!#####################################################################!




!#####################################################################!
!                           Physical Parameters                       !
!#####################################################################!



!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required="false"
!  dem="true" pic="true">
!  <description>
!    Reports mass based on Lagrangian particles and continuum
!    representation. Useful to ensure mass conservation between
!    Lagrangian and continuum representations. Recommended use for
!    debugging purposes.
!  </description>
!  <dependent keyword="DES_INTERP_MEAN_FIELDS" value=".TRUE."/>
      DES_REPORT_MASS_INTERP = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false"
!  dem="true" pic="true">
!  <description>
!    Allows writing of discrete particle data to output files. Relevant
!    to both granular and coupled simulations.
!  </description>
      PRINT_DES_DATA = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    The frequency at which particle data is written. This only applies
!    to pure granular simulations. For coupled simulation, the output
!    frequency is controlled by SPX_DT(1).
!  </description>
!  <dependent keyword="PRINT_DES_DATA" value=".True."/>
!  <dependent keyword="DES_CONTINUUM_COUPLED" value=".False."/>
!  <conflict keyword="DES_CONTINUUM_COUPLED" value=".True."/>
!  <conflict keyword="MPPIC" value=".True."/>
      DES_SPX_DT = LARGE_NUMBER
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    The frequency at which _DES.RES file is written. This only applies
!    to pure granular simulations as the restart frequency is governed
!    by RES_DT for coupled simulations.
!  </description>
!  <dependent keyword="DES_CONTINUUM_COUPLED" value=".False."/>
!  <conflict keyword="DES_CONTINUUM_COUPLED" value=".True."/>
!  <conflict keyword="MPPIC" value=".True."/>
      DES_RES_DT = LARGE_NUMBER
!</keyword>

!<keyword category="Output Control" required="false"
!   dem="true" pic="true">
!  <description> The output file format for DES data.</description>
!  <valid value="PARAVIEW" note="ParaView formatted files (.vtp)"/>
!  <valid value="TECPLOT" note="Tecplot formatted files (.dat)"/>
      DES_OUTPUT_TYPE = "PARAVIEW"
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Runtime flag to generate debugging information. Additional data for
!    FOCUS_PARTICLE is saved.
!  </description>
      DEBUG_DES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>
!    Specify particle number for particle level debugging details.
!  </description>
!  <dependent keyword="DEBUG_DES" value=".TRUE."/>
      FOCUS_PARTICLE = 0
!</keyword>

!<keyword category="Output Control" required="false" pic="true">
!  <description>
!    Flag to print processor level parcel seeding statistics for inflow
!    BC with PIC model.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_SEEDING_STATS = .false.
!</keyword>

!<keyword category="Output Control" required="false" pic="true">
!  <description>
!     Flag to print processor level parcel deletion statistics for
!     outflow BC with PIC model. Not recommended for production runs.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_DELETION_STATS = .false.
!</keyword>




!#####################################################################!
! DEM/PIC COMMON:      Discrete Element Simulation                    !
!#####################################################################!


!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    Number of particles to be read in from the particle_input.dat file.
!    This value is overwritten when using automatic particle generation.
!    A simulation with a mass inflow BC can start without solids by
!    setting PARTICLES = 0.
!  </description>
!  <range min="0" max="+Inf" />
      PARTICLES = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    Automatically generate the initial particle position and velocity
!    data based on the parameters specified for each initial condition
!    (IC) region.
!  </description>
!  <valid value=".TRUE." note="Generate particle configuration based
!    on the initial condition parameters. Data provided in the
!    particle_input.dat file, if present, is ignored. "/>
!  <valid value=".FALSE." note="Particle position and velocity data are
!    provided in the particle_input.dat file. A runtime error occurs if
!    this file is not provided."/>
      GENER_PART_CONFIG = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    To switch between pure granular or coupled simulations of carried
!    and dispersed phase flows.
!  </description>
!  <valid value=".true." note="Performs coupled simulations. "/>
      DES_CONTINUUM_COUPLED = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false">
!  <description>Run one-way coupled simulations. The fluid does not
! see the particles in terms of drag force. The effect of particle volume
! is still felt by the fluid through non-unity voidage values.
! </description>
      DES_ONEWAY_COUPLED = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    Use interpolate to calculate field variables at the particle's
!    position when evaluating fuild-particle interactions.
!  </description>
      DES_INTERP_ON = .FALSE.
!</keyword>

!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    Use interpolation to compute dispersed phase average fields, such
!    as, solids volume fraction. This is enabled automatically when
!    using particle interpolation (DES_INTERP_ON), when using the
!    MPPIC solids model, or the DES cut-cell implementation.
!  </description>
!  <valid value=".TRUE."
!    note="Use interpolation to calculate field quantities."/>
!  <valid value=".FALSE."
!    note="Use cell-based averaging to calculate field quantities."/>
      DES_INTERP_MEAN_FIELDS = .FALSE.
!</keyword>


!<keyword category="Discrete Element Simulation" required="false">
!  <description>
!    Solve a diffusion equation to smooth Lagrangian data mapped to
!    the Eulerian grid (e.g., volume fraction).
!  </description>
      DES_DIFFUSE_MEAN_FIELDS = .FALSE.
!</keyword>


      DES_INTERP_SCHEME = 'NONE'
      DES_INTERP_WIDTH = UNDEFINED

      DES_DIFFUSE_WIDTH = UNDEFINED

      DES_EXPLICITLY_COUPLED = .FALSE.

!<keyword category="Discrete Element Simulation" required="false" dem="true">
!  <description>
!    Time stepping scheme.
!  </description>
!  <valid value="EULER"
!    note="First-Order Euler Scheme."/>
!  <valid value="ADAMS BASHFORTH"
!    note="Second order ADAMS BASHFORTH scheme (DEM only)"/>
      DES_INTG_METHOD = 'EULER'
!</keyword>


      DES_USR_VAR_SIZE = 0


!#####################################################################!
! DEM ONLY:            Discrete Element Model                         !
!#####################################################################!

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    The number of iterations of a pure granular simulation to let
!    the initial particle configuration settle before a coupled
!    gas-solid is started.
!  </description>
!  <range min="0" max="+Inf" />
      NFACTOR = 10
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Maximum number of steps through a DEM loop before a neighbor
!    search will be performed. The search may be called earlier
!    based on other logic.
!  </description>
!  <range min="0.0" max="+Inf" />
      NEIGHBOR_SEARCH_N = 25
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to set the neighbor search algorithm.
!  </description>
!  <valid value="1" note="N-Square search algorithm (most expensive)"/>
!  <valid value="2-4" note="Grid-Based Neighbor Search (Recommended)"/>
      DES_NEIGHBOR_SEARCH = 4
!</keyword>

!<keyword category="Discrete Element Model" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the I-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times the
!    maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element Model" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the J-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times
!    the maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_JMAX = UNDEFINED_I
!</keyword>


!<keyword category="Discrete Element Model" required="false"
!  dem="true" pic="true">
!  <description>
!    Number of des grid cells in the K-direction. If left undefined,
!    then it is set by MFIX such that its size equals three times
!    the maximum particle diameter with a minimum of 1 cell.
!  </description>
      DESGRIDSEARCH_KMAX = UNDEFINED_I
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Ratio of the distance (imaginary sphere radius) to particle radius
!    that is allowed before a neighbor search is performed. This works
!    in conjunction with the logic imposed by NEIGHBOR_SEARCH_N in
!    deciding calls to the neighbor search algorithm.
!  </description>
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Effectively increase the radius of a particle (multiple of the sum
!    of particle radii) during the building of particle neighbor list.
!  </description>
      FACTOR_RLM = 1.2
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to use van der Hoef et al. (2006) model for adjusting the
!    rotation of the contact plane. See the MFIX-DEM documentation.
!  </description>
      USE_VDH_DEM_MODEL = .FALSE.
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Collision model for the soft-sphere approach used in DEM model.
!    All models require specifying the following parameters: DES_EN_INPUT,
!    DES_EN_WALL_INPUT, MEW, and MEW_W.
!  </description>
!  <valid value="LSD" note="The linear spring-dashpot model.
!    Requires: KN, KN_W, KT_FAC, KT_W_FAC, DES_ETAT_FAC, DES_ETAT_W_FAC."/>
!  <valid value="HERTZIAN" note="The Hertzian model.
!    Requires: DES_ET_INPUT, DES_ET_WALL_INPUT, E_YOUNG, EW_YOUNG
!    V_POISSON, VW_POISSON."/>
      DES_COLL_MODEL = 'LSD'
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Normal spring constant for inter-particle collisions needed when
!    using the linear spring-dashpot collision model.
!  </description>
      KN = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential spring constant to normal spring constant
!    for inter-particle collisions. Use it to specify the tangential
!    spring constant for particle-particle collisions as KT_FAC*KN.
!    Required when using the linear spring-dashpot collision model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_FAC = 2.d0/7.d0
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem=.true.>
!  <description>
!    Normal spring constant for particle-wall collisions. Needed when
!    using the linear spring-dashpot collision model.
!  </description>
      KN_W = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential spring constant to normal spring constant
!    for particle-wall collisions. Use it to specify the tangential
!    spring constant for particle-wall collisions as KT_W_FAC*KN_W.
!    Needed when using the linear spring-dashpot collision model.
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_W_FAC = 2.d0/7.d0
!</keyword>

!<keyword category="Discrete Element Model" required="false" dem="true"
!  <description>
!    Inter-particle Coulomb friction coefficient.
!  </description>
! <range min="0.0" max="1.0" />
      MEW = UNDEFINED
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Particle-wall Coulomb friction coefficient.
!  </description>
! <range min="0.0" max="1.0" />
      MEW_W = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    The normal restitution coefficient for inter-particle collisions
!    used to determine the inter-particle normal damping factor.
!
!    Values should be defined for a single dimensional array. For
!    example, a simulation with three solids phases (MMAX=3) needs
!    six values: en11, en12, en13; en22 en 23; en33.
!  </description>
!  <range min="0.0" max="1.0" />
      DES_EN_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    The normal restitution coefficient for particle-wall collisions
!    used to determine the particle-wall normal damping factor.
!
!    Values should be defined in a single dimensional array. For
!    example, a simulation with three solids phases (MMAX=3) needs
!    three values: enw1, enw2, enw3.
!  </description>
!  <range min="0.0" max="1.0" />
      DES_EN_WALL_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Tangential restitution coefficient for inter-particle collisions.
!    Values are defined in a one dimensional array. This is required
!    input when using the Hertzian collision model.
! </description>
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
! <range min="0.0" max="1.0" />
      DES_ET_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Tangential restitution coefficient for particle wall collisions.
!    Values are defined in a one dimensional array. This is required
!    input when using the Hertzian collision model.
!  </description>
! <range min="0.0" max="1.0" />
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      DES_ET_WALL_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false" dem="true">
!  <description>
!    Ratio of the tangential damping factor to the normal damping factor
!    for inter-particle collisions.  Required for the linear spring-
!    dashpot model collision model
!  </description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
!  <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
!   reverts to default value of 0.5" />
      DES_ETAT_FAC = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
! <description>
!    Ratio of the tangential damping factor to the normal damping
!    factor for particle-wall collisions. Required for the linear
!    spring-dashpot model for soft-spring collision modelling under
!    DEM. For the Hertzian model, the tangential damping coefficients
!    have to be explicitly specified and specification of this
!    variable is not required.
! </description>
! <dependent keyword="DES_COLL_MODEL" value="LSD"/>
! <range min="0.0" max="1.0" />
! <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
! will revert to default value of 0.5" />
      DES_ETAT_W_FAC = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>Youngs modulus for the wall. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      EW_YOUNG = UNDEFINED
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>Poisson ratio for the wall. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      VW_POISSON = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>Youngs modulus for the particle. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <arg index="1" id="Phase" min="1" max="DES_MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      E_YOUNG(:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>Poissons ratio for the particle. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <arg index="1" id="Phase" min="1" max="DES_MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      V_POISSON(:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to enable/disable cohesion model.
!  </description>
      USE_COHESION = .FALSE.
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Flag to turn on the use Hamaker van der Waals forces.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VAN_DER_WAALS = .FALSE.
!</keyword>


! for cohesion: van der waals
!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Hamaker constant used in particle-particle cohesive interactions.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      HAMAKER_CONSTANT = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Hamaker constant used in particle-wall cohesive interactions.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_HAMAKER_CONSTANT = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Maximum separation distance above which van der Waals forces are
!    not implemented.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VDW_OUTER_CUTOFF = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Minimum separation distance below which van der Waals forces are
!    calculated using a surface adhesion model.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VDW_INNER_CUTOFF = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Maximum separation distance above which van der Waals forces are
!    not implemented (particle-wall interactions).
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_VDW_OUTER_CUTOFF = ZERO
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Minimum separation distance below which van der Waals forces are
!    calculated using a surface adhesion model (particle-wall
!    interactions).
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_VDW_INNER_CUTOFF = UNDEFINED
!</keyword>


!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Mean radius of surface asperities that influence the cohesive force
!    following a model. See H. Rumpf, Particle Technology, Chapman & Hall,
!    London/New York, 1990.
!  </description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      Asperities = ZERO
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Specify the Nusselt number correlation used for particle-gas
!    convection.
!  </description>
!  <valid value="RANZ_1952" note="Ranz, W.E. and Marshall, W.R. (1952).
!    Chemical Engineering Progress, 48: 141-146 and 173-180"/>
      DES_CONV_CORR = 'RANZ_1952'
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Minimum separation distance between the surfaces of two contacting
!    particles.
!  </description>
      DES_MIN_COND_DIST = UNDEFINED
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>
!    Fluid lens proportion constant used to calculate the radius of
!    the fluid lens that surrounds a particle. This parameter is used
!    in the  particle-fluid-particle conduction model.
!  </description>
      FLPC = 1.0d0/5.0d0
!</keyword>

!<keyword category="Discrete Element Model" required="false">
!  <description>Emissivity of solids phase M.</description>
!  <arg index="1" id="Phase" min="1" max="DES_MMAX"/>
      DES_Em(:DIM_M) = UNDEFINED
!</keyword>



!#####################################################################!
!                          Particle In Cell                           !
!#####################################################################!


!<keyword category="Particle In Cell" required="false">
!  <description>
!    Turn on snider's version of frictional model.
!    Does not run very stably.
!  </description>
      MPPIC_SOLID_STRESS_SNIDER = .false.
!</keyword>


!<keyword category="Particle In Cell" required="false">
!  <description>
!    First coefficient of restitution for the frictional stress model
!    in the MPPIC model. See the MPPIC documentation for more details.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN1 = UNDEFINED
!</keyword>


!<keyword category="Particle In Cell" required="false">
!  <description>
!    Second coefficient of restitution for the frictional stress model
!    in the MPPIC model. See the MPPIC documentation for more details.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN2 = UNDEFINED
!</keyword>


!<keyword category="Particle In Cell" required="false">
!  <description>
!    Normal coefficient of restitution for parcel-wall collisions
!    in the MPPIC model.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN_WALL = UNDEFINED
!</keyword>


!<keyword category="Particle In Cell" required="false">
!  <description> Tangential coefficient of restitution for
! parcel-wall collisions in the MPPIC model.
! Currently not implemented in the code.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_ET_WALL = 1.0
!</keyword>


!<keyword category="Particle In Cell" required="false">
!  <description> Turn on the implicit treatment for interphase drag force.
! Valid only for MPPIC model..
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_PDRAG_IMPLICIT = .false.
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>
!     Variable to decide if special treatment is needed or not in the
!     direction of gravity in the frictional stress tensor. See the
!     MPPIC documentation for details.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_GRAV_TREATMENT = .true.
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>
!    A run time flag to report minimum value and location of gas
!    voidage. This is useful only for debugging and is not
!    recommended for production runs.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_MIN_EPG = .FALSE.
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>
!    P_s term in the frictional stress model of Snider.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PSFAC_FRIC_PIC = 100
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>
!    Beta term in the frictional stress model of Snider.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      FRIC_EXP_PIC = 2.5
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>
!    Non-singularity term (epsilon) in the frictional stress model of
!    Snider.
!  </description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      FRIC_NON_SING_FAC = 1E-07
!</keyword>

!<keyword category="Particle In Cell" required="false">
!  <description>CFL number used to decide maximum time
! step size for parcels evolution equations.
! Relevant to MPPIC model only.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      CFL_PIC = 0.1
!</keyword>


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                            UNSUPPORTED KEYWORDS                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! Logical to force the inlet to operate with an ordered boundary
! condition. This may be useful during long simulations or if the
! inlet appears to be taking a long time to randomly place particles.
      FORCE_ORD_BC = .FALSE.

! Lees-Edwards boundary condition to simulate homogeneous shear
! problem with periodic boundary conditions. Not supported in this
! version.
      DES_LE_BC = .FALSE.

! Relative velocity needed for Lees-Edwards BC.
! Not supported in this version.
      DES_LE_REL_VEL = UNDEFINED

! Direction of shear for Lees-Edwards BC.
! Not supported in this version. </description>
      DES_LE_SHEAR_DIR = UNDEFINED_C

! des wall boundaries: wall velocities. I think they probably
! defined for the Lees-Edwards BC's
      DES_BC_Uw_s(:,:) = ZERO
      DES_BC_Vw_s(:,:) = ZERO
      DES_BC_Ww_s(:,:) = ZERO


! These need to be inialized to 0, but they are not part of the namelist
      VTP_FINDEX = 0
      TECPLOT_FINDEX = 0

! not a well supported feature and not generic either. So removing
! from namelists
      DES_CALC_BEDHEIGHT = .FALSE.
      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
