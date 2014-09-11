!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: INIT_NAMELIST                                           !
!  Purpose: initialize the NAMELIST variables                          !
!                                                                      !
!  Author: P. Nicoletti                               Date: 26-NOV-91  !
!                                                                      !
!  Keyword Documentation Format:                                       !
!                                                                      !
!<keyword category="category name" required=.TRUE./.FALSE.             !
!                                    legacy=.TRUE./.FALSE.>            !
!  <description></description>                                         !
!  <arg index="" id="" max="" min=""/>                                 !
!  <dependent keyword="" value="DEFINED"/>                             !
!  <conflict keyword="" value="DEFINED"/>                              !
!  <valid value="" note="" alias=""/>                                  !
!  <range min="" max="" />                                             !
!  MFIX_KEYWORD=INIT_VALUE                                             !
!</keyword>                                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE INIT_NAMELIST

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE bc
      USE ps
      USE fldvar
      USE constant
      USE indices
      USE is
      USE toleranc
      USE scales
      USE ur_facs
      USE leqsol
      USE residual
      USE rxns
      USE scalars
      USE compar
      USE parallel
      USE cdist
      USE stiff_chem
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counters
      INTEGER :: LC, LCM, M, N
! Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'namelist.inc'




!#####################################################################!
!                             Run Control                             !
!#####################################################################!



!<keyword category="Run Control" required=.TRUE.>
!  <description> Name used to create output files. The name should
!    generate legal file names after appending extensions.
!    Ex: Given the input, RUN_NAME = "bub01", MFIX will generate
!    the output files: BUB01.LOG, BUB01.OUT, BUB01.RES, etcs.
!  </description>
      RUN_NAME = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Problem description. Limited to 60 characters.</description>
      DESCRIPTION = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required=.TRUE.>
!  <description> Simulation input/output units.</description>
!  <valid value="cgs" note="all input and output in cgs units (g, cm, s, cal)."/>
!  <valid value="si" note="all input and output in si units (kg, m, s, j)."/>
      UNITS = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required=.TRUE.>
!  <description>type of run.</description>
!  <valid value="new" note="new run."/>
!  <valid value="RESTART_1" note="Traditional restart."/>
!  <valid value="RESTART_2"
!    note="Start a new run with initial conditions from a .RES file
!      created from another run."/>
      RUN_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Start-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TIME = UNDEFINED
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Stop-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TSTOP = UNDEFINED
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Starting time step. If left undefined, a steady-state
!    calculation is performed.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT = UNDEFINED
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Maximum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MAX = ONE
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Minimum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MIN = 1.0D-6
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Factor for adjusting time step. Must be less than 1.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="1" />
      DT_FAC = 0.9D0
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Reduce the time step if the residuals stop decreasing.
!  </description>
!  <valid value=".FALSE." note="Continue iterating if residuals stall."/>
!  <valid value=".TRUE."  note="Reduce time step if residuals stall."/>
      DETECT_STALL = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Flag to restart the code when DT < DT_MIN.
!  </description>
      AUTO_RESTART = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Shared gas-pressure formulation. See Syamlal, M. and Pannala, S.
!    “Multiphase continuum formulation for gas-solids reacting flows,”
!    chapter in Computational Gas-Solids Flows and Reacting Systems:
!    Theory, Methods and Practice, S. Pannala, M. Syamlal and T.J.
!    O’Brien (editors), IGI Global, Hershey, PA, 2011.
!  </description>
!  <valid value=".FALSE." note="Use Model A"/>
!  <valid value=".TRUE."  note="Use Model B. Bouillard, J.X.,
!    Lyczkowski, R.W., Folga, S., Gidaspow, D., Berry, G.F. (1989).
!    Canadian Journal of Chemical Engineering 67:218–229."/>
      MODEL_B = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Flag to enable/disable solving the X-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve X-momentum equations."/>
!  <valid value=".FALSE." note="The X velocity initial conditions
!   persist throughout the entire simulation."/>
      MOMENTUM_X_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Flag to enable/disable solving the Y-momentum equations.
! </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Y-momentum equations."/>
!  <valid value=".FALSE." note="The Y velocity initial conditions
!   persist throughout the entire simulation."/>
      MOMENTUM_Y_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Flag to enable/disable solving the Z-momentum equations.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Z-momentum equations."/>
!  <valid value=".FALSE." note="The Z velocity initial conditions
!   persist throughout the entire simulation."/>
      MOMENTUM_Z_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Solve energy equations.</description>
!  <valid value=".FALSE." note="Do not solve energy equations."/>
!  <valid value=".TRUE." note="Solve energy equations."/>
      ENERGY_EQ = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Solve species transport equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".FALSE." note="Solve species equations."/>
!  <valid value=".TRUE." note="Do not solve species equations."/>
      SPECIES_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Granular energy formulation selection.</description>
!  <valid value=".FALSE."
!    note="Use algebraic granular energy equation formulation."/>
!  <valid value=".TRUE."
!    note="Use granular energy transport equation (PDE) formulation."/>
      GRANULAR_ENERGY = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Solids phase stress model.</description>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <valid value="AHMADI"
!    note="Cao and Ahmadi (1995). Int. J. Multiphase Flow 21(6), 1203."/>
!  <valid value="GD_99"
!     note="Garzo and Dufty (1999). Phys. Rev. E 59(5), 5895."/>
!  <valid value="GHD"
!    note="Garzo, Hrenya and Dufty (2007). Phys. Rev. E 76(3), 31304"/>
!  <valid value="IA_NONEP"
!     note="Iddir & Arastoopour (2005). AIChE J. 51(6), 1620"/>
!  <valid value="LUN_1984"
!    note="Lun et al (1984). J. Fluid Mech., 140, 223."/>
!  <valid value="SIMONIN"
!    note="Simonin (1996). VKI Lecture Series, 1996-2"/>
      KT_TYPE = "LUN_1984"
!</keyword>

! Retired keyword for specifying Ahmadi KT Theory.
! Use: KT_TYPE = "AHMADI"
      AHMADI = .FALSE.

! Retired keyword for specifying Simonin KT Theory.
! Use: KT_TYPE = "SIMONIN"
      SIMONIN = .FALSE.

!<keyword category="Run Control" required=.FALSE.>
!  <description>Jenkins small frictional boundary condition.</description>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <dependent keyword="PHI_W" value="DEFINED"/>
!  <valid value=".FALSE." note=""/>
!  <valid value=".TRUE."
!    note="Use the Jenkins small frictional boundary condition."/>
      JENKINS = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Solids stress model selection.</description>
!  <valid value=".FALSE." note="Use the Schaeffer solids stress model."/>
!  <valid value=".TRUE."  note="Use the Princeton solids stress model"/>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <dependent keyword="PHI_W" value="DEFINED"/>
      FRICTION = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    For a term appearing in the frictional stress model
!    invoked with friction = .TRUE.
!  </description>
!  <valid value="0" note="Use S:S in the frictional stress model."/>
!  <valid value="1" note="Use an alternate form suggested by Savage."/>
!  <valid value="2" note="An appropriate combination of above."/>
!  <dependent keyword="friction" value=".TRUE."/>
      SAVAGE = 1
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Schaeffer frictional stress tensor formulation. </description>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <valid value=".TRUE." note="Use the Schaeffer model"/>
!  <valid value=".FALSE." note="Do not use the Schaeffer model."/>
      SCHAEFFER = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Blend the Schaeffer stresses with that of kinetic theory around Ep*.
!  </description>
      BLENDING_STRESS = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Hyperbolic tangent function for blending frictional stress models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".TRUE."/>
!  <conflict keyword="SIGM_BLEND" value=".TRUE."/>
      TANH_BLEND = .TRUE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    A scaled and truncated sigmoidal function for blending
!    frictional stress  models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".TRUE."/>
!  <conflict keyword="TANH_BLEND" value=".TRUE."/>
      SIGM_BLEND = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Correlation to compute maximum packing for polydisperse systems.
!  </description>
!  <valid value=".TRUE."
!    note="Use the Yu and Standish correlation."/>
!  <valid value=".FALSE."
!    note="Do not use the Yu and Standish correlation."/>
      YU_STANDISH = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    use fedors and landel correlation to compute maximum
!    packing for a binary (only) mixture of powders.
!  </description>
      FEDORS_LANDEL = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>call user-defined subroutines.</description>
!  <valid value=".TRUE." note="call user-defined subroutines."/>
!  <valid value=".FALSE." note="do not call user-defined subroutines."/>
      CALL_USR = .FALSE.
!</keyword>


!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    The number of user defined chemical reactions stored
!    in the *.SPA file.
!  </description>
      nRR = 0
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    when activated the k-epsilon turbulence model (for single
!    -phase flow) is solved using standard wall functions.
!  </description>
!  <conflict keyword="L_SCALE0" value="DEFINED"/>
      K_Epsilon = .FALSE.
!</keyword>


!<keyword category="Run Control" required=.FALSE.>
!  <description>
!     Available gas-solids drag models.
!     Note: The extension _PCF following the specified drag model
!     indicates that the polydisperse correction factor is available.
!     For PCF details see:
!     * Van der Hoef MA, Beetstra R, Kuipers JAM. (2005)
!       Journal of Fluid Mechanics.528:233-254.
!     * Beetstra, R., van der Hoef, M. A., Kuipers, J.A.M. (2007).
!       AIChE Journal, 53:489-501.
!     * Erratum (2007), AIChE Journal, Volume 53:3020
!  </description>
!
!  <valid value="SYAM_OBRIEN" note="Syamlal M, O'Brien TJ (1988).
!   International Journal of Multiphase Flow 14:473-481.
!   Two additional parameters may be specified: DRAG_C1, DRAG_D1"/>
!
!  <valid value="GIDASPOW" note="Ding J, Gidaspow D (1990).
!   AIChE Journal 36:523-538"/>
!
!  <valid value="GIDASPOW_BLEND" note="Lathouwers D, Bellan J (2000).
!    Proceedings of the 2000 U.S. DOE
!        Hydrogen Program Review NREL/CP-570-28890."/>
!
!  <valid value="WEN_YU" note="Wen CY, Yu YH (1966).
!   Chemical Engineering Progress Symposium Series 62:100-111."/>
!
!  <valid value="KOCH_HILL" note="Hill RJ, Koch DL, Ladd JC (2001).
!   Journal of Fluid Mechanics, 448: 213-241. and 448:243-278."/>
!
!  <valid value="BVK" note="Beetstra, van der Hoef, Kuipers (2007).
!   Chemical Engineering Science 62:246-255"/>
!
!  <valid value="HYS" note="Yin, X, Sundaresan, S. (2009).
!   AIChE Journal 55:1352-1368
!   This model has a lubrication cutoff distance, LAM_HYS, that can be
!   specified."/>
!
!  <valid value="USER_DRAG" note="Invoke user-defined drag law. (usr_drag.f)"/>
!
!  <valid value="GIDASPOW_PCF" note="see GIDASPOW"/>
!  <valid value="GIDASPOW_BLEND_PCF" note="see GIDASPOW_BLEND"/>
!  <valid value="WEN_YU_PCF" note="see WEN_YU"/>
!  <valid value="KOCH_HILL_PCF" note="see KOCH_HILL"/>
!
      DRAG_TYPE = 'SYAM_OBRIEN'
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This are determined using the Umf spreadsheet.
!  </description>
      drag_c1 = 0.8d0
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This are determined using the Umf spreadsheet.
!  </description>
      drag_d1 = 2.65d0
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    If use_def_lam_hys is set to .FALSE. the user is able to specify a
!    value for the lubrication cutoff distance (lam_hys).  in practice
!    this number should be on the order of the mean free path of the
!    gas for smooth particles, or the rms roughness of a particle if
!    they are rough (if particle roughness is larger than the mean
!   free path).
!  </description>
!  <dependents>USE_DEF_LAM_HYS</dependents>
      LAM_HYS = UNDEFINED
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Radial distribution function at contact for polydisperse systems.
!    Do not specify any RDF for monodisperse systems because Carnahan-
!    Starling is the model only available.
!
!    Carnahan, N.F. and Starling K.E., (1969).
!    The Journal of Chemical Physics, Vol. 51(2):635-636.
!  </description>
!
!  <valid value="LEBOWITZ" note="Lebowitz, J.L. (1964)
!   The Physical Review, A133, 895-899"/>
!
!  <valid value="MODIFIED_LEBOWITZ" note=""/>
!
!  <valid value="MANSOORI" note="
!   Mansoori, GA, Carnahan N.F., Starling, K.E. Leland, T.W. (1971).
!    The Journal of Chemical Physics, Vol. 54:1523-1525."/>
!
!  <valid value="MODIFIED_MANSOORI" note="van Wachem, B.G.M., Schouten, J.C.,
!    van den Bleek, C.M., Krishna, R. and Sinclair, J. L. (2001)
!    AIChE Journal 47:1035–1051."/>
      RDF_TYPE = 'LEBOWITZ'
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    when activated the added (or virtual) mass force effectively
!    acts to increase the inertia of the dispersed phase, which
!    tends to stabilize simulations of bubbly gas-liquid flows.
!  </description>
      Added_Mass = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Subgrid models.
!  </description>
!
!  <valid value="Igci" note="
!   Igci, Y., Pannala, ., Benyahia, S., and Sundaresan S. (2012).
!   Industrial & Engineering Chemistry Research, 2012, 51(4):2094-2103"/>
!
!  <valid value="Milioli" note="
!   Milioli, C.C., Milioli, F. E., Holloway, W., Agrawal, K. and
!   Sundaresan, S. (2013). AIChE Journal, 59:3265–3275."/>
!
      SUBGRID_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>
!    Ratio of filter size to computational cell size.
!  </description>
      filter_size_ratio = 2.0D0
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Flag for subgrid wall effects.</description>
!  <valid value=".FALSE." note="Do not include wall effects."/>
!  <valid value=".TRUE." note="Include subgrid wall effects."/>
      SUBGRID_Wall = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Provide detailed logging of negative density errors.</description>
!  <valid value=".FALSE." note="Do not log negative density errors."/>
!  <valid value=".TRUE." note="Log negative density errors."/>
      REPORT_NEG_DENSITY = .FALSE.
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>Number of user-defined scalar transport equations to solve.</description>
!  <range min="0" max="DIM_SCALAR" />
      NScalar = 0
!</keyword>

!<keyword category="Run Control" required=.FALSE.>
!  <description>The phase convecting the indexed scalar transport equation.</description>
!  <arg index="1" id="Scalar Equation" min="0" max="DIM_SCALAR"/>
!  <range min="0" max="DIM_M" />
      Phase4Scalar(:DIM_SCALAR) = UNDEFINED_I
!</keyword>



!#####################################################################!
!                           Physcial Parameters                       !
!#####################################################################!

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>User defined constants.</description>
      C(:DIMENSION_C) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>Name of user-defined constant. (20 character max)</description>
      C_NAME(:DIMENSION_C) = '....................'
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>coefficient of restitution for particle-particle collisions.</description>
      C_E = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    GHD Theory: Coefficient of restitution for particle-particle collisions.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <arg index="2" id="Phase" min="0" max="DIM_M"/>
      r_p(:DIM_M, :DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>coefficient of restitution for particle-wall collisions.</description>
      E_W = 1.D0
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>specularity coefficient associated with particle-wall collisions.</description>
      PHIP = 0.6D0
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    output the variable specularity coefficient when bc_jj_m is
!    .TRUE.. The specularity coefficient will be stored in reactionrates
!    array for post-processing by post-mfix. user needs to set nrr to 1
!    for this purpose. be careful with this setting when reacting flow
!    is simulated.
!  </description>
      PHIP_OUT_JJ=.FALSE.
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    specify the value of specularity coefficient when the
!    normalized slip velocity goes to zero when bc_jj_m is
!    .TRUE.. this variable is calculated internally in the
!    code. do not modify unless an accurate number is known.
!  </description>
!  <dependents>bc_jj_m</dependents>
      phip0 = undefined
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    coefficient of friction between the particles of two solids phases.
!  </description>
      C_F = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!     angle of internal friction (in degrees). set this value
!     to zero to turn off plastic regime stress calculations.
!  </description>
      PHI = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    Angle of internal friction (in degrees) at walls. Set this
!    value to non-zero (phi_w = 11.31 means tan_phi_w = mu = 0.2)
!    when using jenkins or bc_jj_m boundary condition.
!  </description>
      PHI_W = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    minimum solids fraction above which friction sets in.  (when
!    friction = .TRUE.)</description>
!  <dependents>friction</dependents>
      EPS_F_MIN = 0.5D0
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    Maximum solids volume fraction at packing for polydisperse
!    systems (more than one solids phase used). The value of
!    EP_star may change during the computation if solids phases
!    with different particle diameters are specified and
!    Yu_Standish or Fedors_Landel correlations are used.
!  </description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <range min="0" max="1-ep_star" />
      EP_S_MAX(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    used in calculating the initial slope of segregation: see
!    Gera et al. (2004) - recommended value 0.3. increasing this
!    coefficient results in decrease in segregation of particles
!    in binary mixtures.
!  </description>
      SEGREGATION_SLOPE_COEFFICIENT=0.D0
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    value of turbulent length initialized. this may be overwritten
!    in specific regions with the keyword ic_l_scale.
!</description>
!  <conflict keyword="K_EPSILON" value=".TRUE."/>
      L_SCALE0 = ZERO
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    maximum value of the turbulent viscosity of the fluid.
!  </description>
      MU_GMAX = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>excluded volume in boyle-massoudi stress.</description>
!  <valid value="0.0" note="b-m stress is turned off."/>
      V_EX = ZERO
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>reference pressure.</description>
      P_REF = ZERO
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>scale factor for pressure.</description>
      P_SCALE = ONE
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>gravitational acceleration. </description>
      GRAVITY = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    X-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_X = ZERO
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    Y-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Y = ZERO
!</keyword>

!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>
!    Z-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Z = ZERO
!</keyword>


!<keyword category="Physical Parameters" required=.FALSE.>
!  <description>disperse phase number where the added mass applies.</description>
      M_AM = UNDEFINED_I
!</keyword>




!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum number of iterations.</description>
      MAX_NIT = 500
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>factor to normalize the gas continuity equation residual.</description>
      NORM_G = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>factor to normalize the solids continuity equation residual.</description>
      NORM_S = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (continuity+momentum).</description>
      TOL_RESID = 1.0D-3
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (granular energy).</description>
      TOL_RESID_Th = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (energy).</description>
      TOL_RESID_T = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (species balance).</description>
      TOL_RESID_X = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (scalar balances.)</description>
      TOL_RESID_Scalar = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>maximum residual at convergence (K_Epsilon Model)</description>
      TOL_RESID_K_Epsilon = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description> Minimum residual for declaring divergence.
!    This parameter is useful for incompressible fluid simulations because
!    velocity residuals can take large values for the second iteration (e.g., 1e+8)
!    before drop down to smaller values for third  third iteration (e.g., 0.1).
!  </description>
      TOL_DIVERGE = 1.0D+4
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>
!    The code declares divergence if the velocity anywhere in the domain
!    exceeds a maximum value.  this maximum value is automatically
!    determined from the boundary values. the user may scale the maximum
!    value by adjusting this scale factor.
!  </description>
      MAX_INLET_VEL_FAC = ONE
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Number of iterations in the linear equation solver.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
      LEQ_IT(1) = 20
      LEQ_IT(2) = 20
      LEQ_IT(3) = 5
      LEQ_IT(4) = 5
      LEQ_IT(5) = 5
      LEQ_IT(6) = 15
      LEQ_IT(7) = 15
      LEQ_IT(8) = 15
      LEQ_IT(9) = 15
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>LEQ Solver selection.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="1" note="SOR"/>
!  <valid value="2" note="BiCGSTAB"/>
!  <valid value="3" note="GMRES"/>
!  <valid value="5" note="CG"/>
      LEQ_METHOD(:) = 2
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Linear equation sweep direction.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="RSRS" note="(Red/Black Sweep, Send Receive) repeated twice"/>
!  <valid value="ISIS" note="(Sweep in I, Send Receive) repeated twice"/>
!  <valid value="JSJS" note="(Sweep in J, Send Receive) repeated twice"/>
!  <valid value="KSKS" note="(Sweep in K, Send Receive) repeated twice"/>
!  <valid value="ASAS" note="(All Sweep, Send Receive) repeated twice"/>
      LEQ_SWEEP(:) = 'RSRS'
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Linear Equation tolerance.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <dependent keyword="LEQ_METHOD" value="2"/>
!  <dependent keyword="LEQ_METHOD" value="3"/>
      LEQ_TOL(:) = 1.0D-4
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Linear precondition used for LEQ solver sweeps.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="NONE" note="No preconditioner"/>
!  <valid value="LINE" note="Line relaxation"/>
!  <valid value="DIAG" note="Diagonal Scaling"/>
      LEQ_PC(:) = 'LINE'
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Under relaxation factors.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
      UR_FAC(1)  = 0.8D0             !pressure
      UR_FAC(2)  = 0.5D0             !rho, ep
      UR_FAC(3)  = 0.5D0             !U
      UR_FAC(4)  = 0.5D0             !V
      UR_FAC(5)  = 0.5D0             !W
      UR_FAC(6)  = 1.0D0             !T
      UR_FAC(7)  = 1.0D0             !X
      UR_FAC(8)  = 0.5D0             !Th
      UR_FAC(9)  = 0.8D0             !Scalar
!</keyword>

!<keyword category="run control" required=.FALSE.>
!  <description>
!    Use deferred correction method for implementing higher order
!    discretization.
!  </description>
!  <valid value=".FALSE." note="use down-wind factor method (default)."/>
!  <valid value=".TRUE."
!    note="use deferred correction method for implementing higher order discretization."/>
      DEF_COR  =  .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>discretization scheme of equations.</description>
!  <valid value="1" note="first-order upwinding (using down-wind factors)."/>
!  <valid value="3" note="smart."/>
!  <valid value="2" note="superbee (recommended method)."/>
!  <valid value="5" note="quickest (does not work)."/>
!  <valid value="4" note="ultra-quick."/>
!  <valid value="7" note="van leer."/>
!  <valid value="6" note="muscl."/>
!  <valid value="8" note="minmod."/>
!  <valid value="0" note="first-order upwinding."/>
      DISCRETIZE(:) = 0
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>
!    four point fourth order interpolation and is upstream biased. if
!    this scheme is chosen and discretize(*) < 2, discretize(*) is
!    defaulted to 2. if you chose this scheme, set the c_fac value
!    between 0 and 1.
!  </description>
!  <dependent keyword="c_fac" value="DEFINED"/>
      FPFOI = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>
!    factor used in the universal limiter (when fpfoi is set .TRUE.) and
!    can be any value in the set (0,1). the choice of 1 will give
!    (diffusion) first order upwinding and as this value becomes closer
!    to 0 the scheme becomes more compressive.
!  </description>
!  <range min="0.0" max="1.0" />
!  <dependent keyword="fpfoi" value=".TRUE."/>
      C_FAC = UNDEFINED
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>temporal discretization scheme.</description>
!  <valid value=".FALSE."
!    note="Implicit Euler based temporal discretization scheme employed
!      (first order accurate in time)."/>
!  <valid value=".TRUE."
!    note="crank-nicholson based temporal discretization scheme employed
!      (second order accurate in time excluding the restart timestep
!      which is first order)."/>
      CN_ON = .FALSE.
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>
!    chi-scheme, proposed by darwish and moukalled (2003), is activated.
!    this scheme guarantees that the set of differenced species mass
!    balance equations has the property that the sum of mass fractions
!    add up to one. when a flux limiter is used with (higher order)
!    spatial discretization schemes it is not guaranteed that the mass
!    fractions add up to one. this problem may be rectified by
!    activating the chi-scheme.
!  </description>
      Chi_scheme = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>
!    The implicitness calculation of the gas-solids drag coefficient
!    may be underrelaxed by changing ur_f_gs, which takes values
!    between 0 to 1:
!  </description>
!  <range min="0" max="1" />
      UR_F_gs = 1.0D0
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Unknown</description>
      UR_Kth_sml = 1.0D0
!</keyword>


!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Solve transpose of linear system. (BICGSTAB ONLY).</description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      DO_TRANSPOSE = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Frequency to check for convergence. (BICGSTAB ONLY)</description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      icheck_bicgs = 1
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Sets optimal LEQ flags for parallel runs.</description>
      OPT_PARALLEL = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Use do-loop assignment over direct vector assignment.</description>
      USE_DOLOOP = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required=.FALSE.>
!  <description>Calculate dot-products more efficiently (Serial runs only.)</description>
      IS_SERIAL = .TRUE.
!</keyword>


!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!


!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>coordinates used in the simulation.</description>
!  <valid value="cartesian" note="cartesian coordinates."/>
!  <valid value="cylindrical" note="cylindrical coordinates."/>
      COORDINATES = UNDEFINED_C
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>(do not use.)</description>
!  <valid value=".FALSE." note="x (r) direction is considered."/>
!  <valid value=".TRUE." note="x (r) direction is not considered."/>
      NO_I = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>number of cells in the x (r) direction.</description>
      IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    Cell sizes in the x (r) direction. Enter values from DX(0) to
!    DX(IMAX-1). (Use uniform mesh size with higher-order
!    discretization methods.  Also in cylindrical coordinates DX
!    should be kept uniform for strict momentum conservation.)
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_I"/>
      DX(:DIM_I) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    The inner radius in the simulation of an annular cylindrical region.
!  </description>
      XMIN = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>reactor length in the x (r) direction.</description>
      XLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>(do not use.)</description>
!  <valid value=".FALSE. note="y direction is considered."/>
!  <valid value=".TRUE." note="y direction is not considered."/>
      NO_J = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>number of cells in the y direction.</description>
      JMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    Cell sizes in the y direction. Enter values from DY(0) to
!    DY(IMAX-1). (Use uniform mesh size with second-order
!    discretization methods.)
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_J"/>
      DY(:DIM_J) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>reactor length in the y direction.</description>
      YLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description></description>
!  <valid value=".FALSE." note="z(theta) direction is considered."/>
!  <valid value=".TRUE." note="z(theta) direction is not considered."/>
      NO_K = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>number of cells in the z (() direction.</description>
      KMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    Cell sizes in the z (theta) direction. Enter values from DZ(0) to
!    DZ(IMAX-1). (Use uniform mesh size with second-order discretization
!    methods.)
!  </description>
!  <arg index="1" id="Cell" min="0" max="DIM_K"/>
      DZ(:DIM_K) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>reactor length in the z (theta) direction.</description>
      ZLENGTH = UNDEFINED
!</keyword>


!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    Flag for making the x-direction cyclic without pressure drop. no other
!    boundary conditions for the x-direction should be specified.
!</description>
!  <valid value=".FALSE." note="no cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="cyclic condition at x-boundary."/>
      CYCLIC_X = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    flag for making the x-direction cyclic with pressure drop. if the
!    keyword flux_g is given a value this becomes a cyclic boundary
!    condition with specified mass flux. no other boundary conditions
!    for the x-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="no cyclic condition at x-boundary."/>
!  <valid value=".TRUE." note="cyclic condition with pressure drop at x-boundary."/>
      CYCLIC_X_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    fluid pressure drop across xlength when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      DELP_X = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    flag for making the y-direction cyclic without pressure drop. no
!    other boundary conditions for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="no cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="cyclic condition at x-boundary."/>
      CYCLIC_Y = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    flag for making the y-direction cyclic with pressure drop. if the
!    keyword flux_g is given a value this becomes a cyclic boundary
!    condition with specified mass flux. no other boundary conditions
!    for the y-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="no cyclic condition at y-boundary."/>
!  <valid value=".TRUE." note="cyclic condition with pressure drop at y-boundary."/>
      CYCLIC_Y_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    fluid pressure drop across ylength when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      DELP_Y = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    flag for making the z-direction cyclic without pressure drop. no
!    other boundary conditions for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="no cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="cyclic condition at z-boundary."/>
      CYCLIC_Z = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    flag for making the z-direction cyclic with pressure drop. if the
!    keyword flux_g is given a value this becomes a cyclic boundary
!    condition with specified mass flux. no other boundary conditions
!    for the z-direction should be specified.
!  </description>
!  <valid value=".FALSE." note="no cyclic condition at z-boundary."/>
!  <valid value=".TRUE." note="cyclic condition with pressure drop at z-boundary."/>
      CYCLIC_Z_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    fluid pressure drop across zlength when a cyclic boundary condition
!    with pressure drop is imposed in the z-direction.
!  </description>
      DELP_Z = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    if .TRUE. imposes a mean shear on the flow field as a linear
!    function of _x_ coordinate. this feature should only be used when
!    cyclic_x=.TRUE. also, the keyword v-sh needs to be set.
!  </description>
!  <dependents>cyclic_x v_sh</dependents>
      SHEAR = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    Specifies the mean _y_ velocity component at the eastern boundary
!    of the domain (v_sh), and the mean _y_ velocity (-v_sh) at the
!    western boundary of the domain.
!  </description>
      V_sh=0d0
!</keyword>


!<keyword category="Geometry and Discretization" required=.FALSE.>
!  <description>
!    if a value is specified (in units of g/cm2.s), the domain-averaged gas
!    flux is held constant at that value in simulations over a periodic
!    domain.  a pair of boundaries specified as periodic with fixed
!    pressure drop is then treated as periodic with fixed mass flux.
!    Even for this case a pressure drop must also be specified, which
!    is used as the initial guess in the simulations.
!  </description>
      Flux_g = UNDEFINED
!</keyword>




!#####################################################################!
!                               Gas Phase                             !
!#####################################################################!




!<keyword category="Gas Phase" required=.FALSE.>
!  <description>
!    Specified constant gas density. this value may be set to zero to
!    make the drag zero and to simulate granular flow in a vacuum. for
!    this case, users may turn off solving for gas momentum equations
!    to accelerate convergence.
!  </description>
      RO_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>specified constant gas viscosity.</description>
      MU_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>specified constant gas conductivity.</description>
      K_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>specified constant gas specific heat.</description>
      C_PG0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>specified constant gas diffusivity.</description>
      DIF_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>average molecular weight of gas.</description>
      MW_AVG = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>molecular weight of gas species n.</description>
!  <arg index="1" id="Species" min="1" max="DIM_N_G"/>
      MW_G(:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>Number of species comprising the gas phase.</description>
      NMAX_g = UNDEFINED_I
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>Name of gas phase species n as it appears in the materials database.</description>
      SPECIES_g = UNDEFINED_C
!</keyword>

!<keyword category="Gas Phase" required=.FALSE.>
!  <description>User defined name for gas phase species n.</description>
      SPECIES_ALIAS_g = UNDEFINED_C
!</keyword>



!#####################################################################!
!                            Solids Phase                             !
!#####################################################################!


!<keyword category="Solids Phase" required=.FALSE.>
!  <description>number of solids phases.</description>
      MMAX = 1
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Initial particle diameters.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      D_P0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Specified constant solids density.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      RO_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Baseline species mass fraction.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      X_s0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Specified constant solids species density.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      RO_Xs0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Index of inert solids phase species.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      INERT_SPECIES(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>specified constant granular viscosity. if this value is
!    specified, then the kinetic theory calculation is turned off and
!    p_s = 0 and lambda_s = -2/3 mu_s0.
!  </description>
      MU_S0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>specified constant solids conductivity.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      K_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>specified constant solids specific heat.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      C_PS0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>specified constant solids diffusivity.</description>
      DIF_S0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Molecular weight of solids phase-m, species n.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      MW_S(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>Number of species comprising solids phase m.</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      NMAX_s(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>
!    Name of solids phase m, species n as it appears in the materials
!    database.
!</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      SPECIES_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>User defined name for solids phase m, species n</description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_s"/>
      SPECIES_ALIAS_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>packed bed void fraction.</description>
      EP_STAR = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required=.FALSE.>
!  <description>
!    Indicates that the solids phase forms a packed bed with a void
!    fraction ep_star.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
      CLOSE_PACKED(:DIM_M) = .TRUE.
!</keyword>


!<keyword category="Solids Phase" required=.FALSE.>
!  <description>specified constant granular viscosity. if this value is
!    specified, then the kinetic theory calculation is turned off and
!    p_s = 0 and lambda_s = -2/3 mu_s0.
!  </description>
!  <arg index="1" id="Phase" min="1" max="DIM_M"/>
!  <valid value='TFM' note='Two-fluid Model (coninuum).' />
!  <valid value='DEM' note='Discrete Element Model' />
!  <valid value='PIC' note='Multiphase-Particle in Cell' />
      SOLIDS_MODEL(:DIM_M) = 'TFM'
!</keyword>



!#####################################################################!
!                         Initial Conditions                          !
!#####################################################################!


      DO LC = 1, DIMENSION_IC

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>x coordinate of the west face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>x coordinate of the east face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>y coordinate of the south face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>y coordinate of the north face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>z coordinate of the bottom face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>z coordinate of the top face.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>i index of the west-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>i index of the east-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>j index of the south-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>j index of the north-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>k index of the bottom-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>k index of the top-most wall.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Type of initial condition. Mainly used in restart runs to overwrite
!    values read from the .RES file by specifying it as _PATCH_. The
!    user needs to be careful when using the _PATCH_ option, since the
!    values from the .RES file are overwritten and no error checking is
!    done for the patched values.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, MFIX will set up a hydrostatic pressure profile,
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Initial solids pressure in the IC region. Usually, this value is
!    specified as zero.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_P_STAR(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Turbulence length scale in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_L_SCALE(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Initial bulk density (rop_s = ro_s x ep_s) of solids phase-m in the
!    IC region. Users need to specify this IC only for polydisperse flow
!    (MMAX > 1). Users must make sure that summation of ( IC_ROP_s(ic,m)
!    / RO_s(m) ) over all solids phases is equal to ( 1.0 - IC_EP_g(ic)).
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Initial solids volume fraction of solids phase-m in the IC region.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial solids phase-m granular temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Gas phase radiation coefficient in the IC region. Modify file
!    rdtn2.inc to change the source term.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RG(LC) = ZERO
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Gas phase radiation temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RG(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>
!    Solids phase-m radiation coefficient in the IC region. Modify file
!    radtn2.inc to change the source term.
!  </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RS(LC,:DIM_M) = ZERO
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Solids phase-m radiation temperature in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RS(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial mass fraction of gas species n.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         IC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial mass fraction of gas species n.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         IC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial value of Scalar n.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
        IC_SCALAR(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial value of K in K-Epsilon.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Initial value of Epsilon in K-Epsilon.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
         IC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Flag for inflating initial lattice distribution
! to the entire IC region. </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
          IC_DES_FIT_TO_REGION(LC) = .FALSE.
!</keyword>


!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Flag to specify the initial constant number
! of particles per cell for the PIC method initialization.
!Statistical weight of parcels will be calculated by the code.</description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
!  <conflict keyword="IC_PIC_CONST_STATWT" value="DEFINED"/>
          IC_PIC_CONST_NPC(LC, :DIM_M) = 0
!</keyword>


!<keyword category="Initial Condition" required=.FALSE.>
!  <description>Flag to specify the initial constant statistical
! weight for computational particles/parcels. Actual number of
! parcels will be automatically computed. </description>
!  <arg index="1" id="IC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
!  <conflict keyword="IC_PIC_CONST_NPC" value="DEFINED"/>
          IC_PIC_CONST_STATWT(LC, :DIM_M) = ZERO
!</keyword>
      ENDDO




!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!
      DO LC = 1, DIMENSION_BC


!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>x coordinate of the west face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>x coordinate of the east face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>i index of the west-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>i index of the east-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>j index of the south-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>j index of the north-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>k index of the bottom-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>k index of the top-most cell.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!
!  <valid value='DUMMY'
!    note='The specified boundary condition is ignored. This is
!      useful for turning off some boundary conditions without having
!      to delete them from the file.' />
!
!  <valid value='MASS_INFLOW' alias='MI'
!    note='Mass inflow rates for gas and solids phases are
!      specified at the boundary.'/>
!
!  <valid value='MASS_OUTFLOW' alias='MO'
!    note='The specified values of gas and solids mass outflow
!      rates at the boundary are maintained, approximately. This
!      condition should be used sparingly for minor outflows, when
!      the bulk of the outflow is occurring through other constant
!      pressure outflow boundaries.' />
!
!  <valid value='P_INFLOW' alias='PI'
!    note='Inflow from a boundary at a specified constant
!      pressure. To specify as the west, south, or bottom end of
!      the computational region, add a layer of wall cells to the
!      west, south, or bottom of the PI cells. Users need to specify
!      all scalar quantities and velocity components. The specified
!      values of fluid and solids velocities are only used initially
!      as MFIX computes these values at this inlet boundary.' />
!
!  <valid value='P_OUTFLOW' alias='PO'
!    note='Outflow to a boundary at a specified constant pressure.
!      To specify as the west, south, or bottom end of the computational
!      region, add a layer of wall cells to the west, south, or bottom of
!      the PO cells.' />
!
!  <valid value='FREE_SLIP_WALL' alias='FSW'
!    note='Velocity gradients at the wall vanish. If BC_JJ_PS is
!      equal to 1, the Johnson-Jackson boundary condition is used for
!      solids.  A FSW is equivalent to using a PSW with hw=0.' />
!
!  <valid value='NO_SLIP_WALL' alias='NSW'
!    note='All components of the velocity vanish at the wall. If
!      BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is
!      used for solids.  A NSW is equivalent to using a PSW with vw=0
!      and hw undefined.' />
!
!  <valid value='PAR_SLIP_WALL' alias='PSW'
!    note='Partial slip at the wall implemented as
!      dv/dn + hw (v - vw) = 0, where n is the normal pointing from the
!      fluid into the wall. The coefficients hw and vw should be
!      specified. For free slip set hw = 0. For no slip leave hw
!      undefined (hw=+inf) and set vw = 0. To set hw = +inf, leave it
!      unspecified. If BC_JJ_PS is equal to 1, the Johnson-Jackson
!      boundary condition is used for solids.' />
         BC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_UW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_UW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_WW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_WW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Johnson and Jackson partial slip BC.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <valid value='0'
!    note='Do not use Johnson and Jackson partial slip bc. If
!      granular energy transport equation is not solved./>
!  <valid value='1'
!    note='Use Johnson and Jackson partial slip bc. If granular
!      energy transport equation is solved. />
         BC_JJ_PS(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified wall value, THETAw_M, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Transfer coefficient, Hw, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Hw for granular energy bc.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant flux, C, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_HW_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified gas phase wall temperature, Tw_g, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_TW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_C_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Solids phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Solids phase hw for heat transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified solids phase wall temperature, Tw_s, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_TW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant solids phase heat flux, C, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_C_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Gas phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified wall gas species mass fraction, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant gas species mass flux, C, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Solid phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_HW_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified solids species mass fraction at the wall, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_XW_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant solids species mass flux, C, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_C_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Scalar transfer coefficient, Hw, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_HW_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified scalar value at the wall, ScalarW, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_ScalarW(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>
!    Specified constant scalar flux, C, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_C_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Bulk density of solids phase at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids volume fraction at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids phase-m granular temperature at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Mass fraction of gas species n at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Mass fraction of solids phase-m, species n at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>x-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>x-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>The interval at the beginning when the normal
!    velocity at the boundary is equal to BC_Jet_g0. When restarting,
!    run this value and BC_Jet_g0 should be specified such that the
!    transient jet continues correctly. MFIX does not store the jet
!    conditions. For MASS_OUTFLOW boundary conditions, BC_DT_0 is
!    the time period to average and print the outflow rates. The
!    adjustment of velocities to get a specified mass or volumetric
!    flow rate is based on the average outflow rate.
!  </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Value of normal velocity during the initial interval BC_DT_0.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_G0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>The interval when normal velocity is equal to BC_Jet_gh.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_H(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Value of normal velocity during the interval BC_DT_h.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GH(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>The interval when normal velocity is equal to BC_JET_gL.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_DT_L(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Value of normal velocity during the interval BC_DT_L.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_JET_GL(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Boundary value for user-defined scalar equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Boundary value of K for K-Epsilon Equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Boundary value of Epsilon for K-Epsilon Equation.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
         BC_VELMAG_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Magnitude of gas velocity in a specified boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Flag to specify the constant number
! of computational particles per cell for the PIC solids inflow BC.
!Statistical weight of parcels will be calculated by the code.</description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="BC_PIC_CONST_STATWT" value="DEFINED"/>
!  <dependent keyword="SOLIDS_MODEL" value="PIC"/>
          BC_PIC_MI_CONST_NPC(LC, :DIM_M) = 0
!</keyword>


!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Flag to specify the constant statistical
! weight for inflowing computational particles/parcels. Actual number of
! parcels will be automatically computed. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <conflict keyword="IC_PIC_CONST_NPC" value="DEFINED"/>
          BC_PIC_MI_CONST_STATWT(LC, :DIM_M) = ZERO
!</keyword>

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Flag to make the PO BC invisible to discrete solids.
! Set this flag to.FALSE.to remove this BC for discrete solids. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
         BC_PO_APPLY_TO_DES(LC) = .TRUE.
!</keyword>


!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Flag to make the inflow plane invisible to discrete solids.
! Set this flag to.FALSE.to remove to inflow plane. </description>
!  <arg index="1" id="BC" min="1" max="DIMENSION_IC"/>
         BC_MI_AS_WALL_FOR_DES(LC)  = .TRUE.
!</keyword>

         BC_ROP_G(LC) = UNDEFINED
      ENDDO

!<keyword category="Boundary Condition" required=.FALSE.>
!  <description>Use the modified Johnson and Jackson partial slip BC
!    with variable specularity coefficient.
!  </description>
!  <dependent keyword="e_w" value="DEFINED"/>
!  <dependent keyword="phi_w" value="DEFINED"/>
         BC_JJ_M = .FALSE.
!</keyword>




!#####################################################################!
!                         Internal Surfaces                           !
!#####################################################################!
      DO LC = 1, DIMENSION_IS


!<keyword category="Internal Surface" required=.FALSE.>
!  <description>x coordinate of the west face or edge.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>x coordinate of the east face or edge.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>y coordinate of the south face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>y coordinate of the north face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>z coordinate of the bottom face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>z coordinate of the top face or edge</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>i index of the west-most cell.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>i index of the east-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>j index of the south-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>j index of the north-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>k index of the bottom-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>k index of the top-most cell</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>Type of internal surface</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
!  <valid value="IMPERMEABLE"
!    note="No gas or solids flow through the surface." alias="IP"/>
!  <valid value="SEMIPERMEABLE" aliase='SP'
!    note="Gas flows through the surface with an additional resistance.
!      Solids velocity through the surface is set to zero or to a user-
!      specified fixed value (i.e., solids momentum equation for this
!      direction is not solved)." />
         IS_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>permeability</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_PC(LC,1) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required=.FALSE.>
!  <description>Inertial resistance coefficient.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
         IS_PC(LC,2) = ZERO
!</keyword>


!<keyword category="Internal Surface" required=.FALSE.>
!  <description>Value of fixed solids velocity through semipermeable surfaces.</description>
!  <arg index="1" id="IS" min="1" max="DIMENSION_IS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IS_VEL_S(LC,:DIM_M) = ZERO
!</keyword>
      ENDDO


!#####################################################################!
!                     Point Source Mass Inlets                        !
!#####################################################################!
      DO LC = 1, DIMENSION_PS

!<keyword category="Point Source" required=.FALSE.>
!  <description>x coordinate of the west face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>x coordinate of the east face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>y coordinate of the south face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>y coordinate of the north face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>z coordinate of the top face or edge.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>i index of the west-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>i index of the east-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>j index of the south-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>j index of the north-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>k index of the bottom-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>k index of the top-most cell.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_K_T(LC) = UNDEFINED_I
!</keyword>


!<keyword category="Point Source" required=.FALSE.>
!  <description>x-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>y-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>z-component of incoming gas velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Gas mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Temperature of incoming gas.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Gas phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
         PS_X_G(LC,:DIM_N_g) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>x-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>y-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>z-component of incoming solids velocity.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Solids mass flow rate through the point source.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Temperature of incoming solids.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         PS_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required=.FALSE.>
!  <description>Solids phase incoming species n mass fraction.</description>
!  <arg index="1" id="PS" min="1" max="DIMENSION_PS"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         PS_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

      ENDDO


!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required=.FALSE.>
!  <description>interval at which restart (.res) file is updated.</description>
      RES_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Interval at which .SPX files are updated. </description>
!  <valid value="SP1" note="void fraction (EP_G"/>
!  <valid value="SP2" note="Gas pressure (P_G), and Solids pressure (P_star)"/>
!  <valid value="SP3" note="Gas veloicty (U_G, V_G, W_G)"/>
!  <valid value="SP4" note="Solids veloicty (U_S, V_S, W_S)"/>
!  <valid value="SP5" note="Solids builk density (ROP_s)"/>
!  <valid value="SP6" note="Gas and solids temperature (T_G, T_S)"/>
!  <valid value="SP7" note="Gas and solids mass fractions (X_G, X_S)"/>
!  <valid value="SP8" note="Granular temperature (THETA_M)"/>
!  <valid value="SP9" note="User defined scalars. (SCALAR)"/>
!  <valid value="SPA" note="Reaction Rates (ReactionRates)"/>
!  <valid value="SPB" note="Turbulence quantities (K_TURB_G, E_TURB_G)"/>
      SPX_DT(:N_SPX) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description> Interval at which standard output (.OUT) file is updated.
!    Only run configuration information is written if left undefined. Otherwise
!    all field variables for the entire domain are written in ASCII
!    format to the .OUT file at OUT_DT intervals.
!  </description>
      OUT_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Number of time steps between .LOG file updates.</description>
      NLOG = 25
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description> Display the residuals on the screen and provide
!    messages about convergence on the screen and in the .LOG file.
!  </description>
      FULL_LOG = .FALSE.
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Specifies the residuals to display. </description>
!  <arg index="1" id="Residual Index" max="8" min="1"/>
!  <valid value="P0" note="Gas pressure"/>
!  <valid value="PM" note="Solids phase M pressure"/>
!  <valid value="R0" note="Gas density"/>
!  <valid value="RM" note="Solids phase M density"/>
!  <valid value="U0" note="Gas phase U-velocity"/>
!  <valid value="V0" note="Gas phase V-velocity"/>
!  <valid value="W0" note="Gas phase W-velocity"/>
!  <valid value="UM" note="Solids phase M U-velocity"/>
!  <valid value="VM" note="Solids phase M V-velocity"/>
!  <valid value="WM" note="Solids phase M W-velocity"/>
!  <valid value="T0" note="Gas temperature"/>
!  <valid value="TM" note="Solids phase M temperature"/>
!  <valid value="X0NN" note="Gas phase species NN mass fraction"/>
!  <valid value="XMNN" note="Solids phase M species NN mass fraction"/>
!  <valid value="K0" note="K-Epsilon model residuals"/>
      RESID_STRING(:8) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Display residuals by equation.  </description>
      GROUP_RESID = .FALSE.
!</keyword>


      DO LC=1, DIMENSION_USR
!<keyword category="Output Control" required=.FALSE.>
!  <description> Intervals at which subroutine write_usr1 is called. </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_DT(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: x coordinate of the west face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: x coordinate of the east face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: y coordinate of the south face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: y coordinate of the north face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: z coordinate of the top face or edge.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: i index of the west-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: i index of the east-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: j index of the south-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: j index of the north-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: k index of the bottom-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: k index of the top-most cell.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: Type of user-defined ouput: Binary of ASCII.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook:
!    Variables to be written in the user-defined output files.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_VAR(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook:
!    Format for writing user-defined (ASCII) output file.
!  </description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_FORMAT(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>UDF Hook: File extension for the user-defined output.</description>
!  <arg index="1" id="USR" max="DIMENSION_USR" min="1"/>
         USR_EXT(LC) = UNDEFINED_C
!</keyword>
      ENDDO


!<keyword category="Output Control" required=.FALSE.>
!  <description> Frequency to perform an overall species mass balance.
!    Leaving undefined suppresses the mass balance calculations which
!    can slightly extend run time.
!  </description>
      REPORT_MASS_BALANCE_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Use distributed IO :: Each rank generates RES/SPx files.</description>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
      bDist_IO = .FALSE.
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Restart a unified IO run as distributed IO.</description>
!  <dependent keyword="RUN_TYPE" value="RESTART_2"/>
!  <dependent keyword="bDist_IO" value=".TRUE."/>
      bStart_with_one_RES = .FALSE.
!</keyword>

!<keyword category="Output Control" required=.FALSE.>
!  <description>Flag to write variable in NetCDF output file. NetCDF support
!    is not included in MFIX by default. The executable must be compiled and
!    linked with an appropriate NetCDF library to use this functionality.
!
!    Variable Index List:
!     1: void fraction (EP_G)
!     2: Gas pressure (P_G)
!     3: Solids pressure (P_star)
!     4: Gas veloicty (U_G, V_G, W_G)
!     5: Solids veloicty (U_S, V_S, W_S)
!     6: Solids builk density (ROP_s)
!     7: Gas temperature (T_G)
!     8: Gas and solids temperature (T_S)
!     9: Gas mass fractions (X_G)
!    10: Solids mass fractions (X_S)
!    11: Granular temperature (THETA_M)
!    12: User defined scalars. (SCALAR)
!    13: Reaction Rates (ReactionRates)
!    14: Turbulence quantities (K_TURB_G, E_TURB_G)
!  </description>
!  <arg index="1" id="NetCDF Variable Reference" max="20" min="1"/>
!  <valid value=".TRUE." note="Write variable in NetCDF output."/>
!  <valid value=".FALSE." note="Do not include variable in NetCDF output."/>
      bWrite_netCDF(:20) = .FALSE.
!</keyword>


!#####################################################################!
!                        Chemical Reactions                           !
!#####################################################################!


!<keyword category="Chemical Reactions" required=.FALSE.>
!  <description>Flag to use stiff chemistry solver (Direct Integration).</description>
!  <conflict keyword="USE_RRATES" value=".TRUE."/>
      STIFF_CHEMISTRY = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE.>
!  <description>
!    Maximum number of internal steps ODEPACK may use to integrate
!    over the time interval.Leaving this value unspecified permits
!    an unlimited number of steps. Thee stiff solver reports the
!    number of cells that exceed the number of steps as 'incomplete'.
!  </description>
!  <dependent keyword="STIFF_CHEMISTRY" value=".TRUE."/>
!  <conflict keyword="USE_RRATES" value=".TRUE."/>
      STIFF_CHEM_MAX_STEPS = UNDEFINED_I
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE.>
!  <description>Flag to use legacy chemcial reaction UDFs.</description>
      USE_RRATES = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE. legacy=.TRUE.>
!  <description>
!    Names of gas and solids phase species as it appears in the
!    materials database. The first NMAX(0) are the names of gas
!    species. The next NMAX(1) are the names of solids phase-1
!    species, etc.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      SPECIES_NAME(:DIM_N_ALL) = UNDEFINED_C
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE.>
!  <description>
!    Number of species in phase m. Note that the gas phase is indicated
!    as m=0.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      NMAX = UNDEFINED_I
!</keyword>

!<keyword category="Chemical Reactions" legacy=.TRUE.>
!  <description>Flag for previous stiff solver.</description>
      CALL_DI   = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE. legacy=.TRUE.>
!  <description>
!   Flag to specify variable solids diameter in original stiff chem solver.
!   (Non-Functional, Removed in 2013-2 Release)
!  </description>
      CALL_GROW = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE. legacy=.TRUE.>
!  <description>
!   Flag to use ISAT tables with original DI solver.
!   (Non-Functional, Removed in 2013-2 Release)
!  </description>
      CALL_ISAT = .FALSE.      ! Legacy Keyword
!  </description>
!</keyword>

!<keyword category="Chemical Reactions" required=.FALSE. legacy=.TRUE.>
!  <description>
!   Specified constant call to ISAT functions.
!   (Non-Functional, Removed in 2013-2 Release)
!  </description>
      ISATdt = UNDEFINED
!  </description>
!</keyword>



!#####################################################################!
!                    Parallelization Control                          !
!#####################################################################!


!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>number of grid blocks in x-direction.</description>
      NODESI = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>number of grid blocks in y-direction.</description>
      NODESJ = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>number of grid blocks in z-direction.</description>
      NODESK = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>Print out additional statistics for parallel runs</description>
      solver_statistics = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>Group residuals to reduce global collectives.</description>
      DEBUG_RESID = .TRUE.
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>All ranks write error messages.</description>
      ENABLE_DMP_LOG = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required=.FALSE.>
!  <description>Print the index layout for debugging.</description>
      DBGPRN_LAYOUT = .FALSE.
!</keyword>


!#####################################################################!
!                       Batch Queue Environment                       !
!#####################################################################!


!<keyword category="Batch Queue Environment" required=.FALSE.>
!  <description>enables clean termination feature.</description>
      CHK_BATCHQ_END = .FALSE.
!</keyword>

!<keyword category="Batch Queue Environment" required=.FALSE.>
!  <description>total wall-clock duration of the job, in seconds.</description>
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
!</keyword>

!<keyword category="Batch Queue Environment" required=.FALSE.>
!  <description>buffer time when initiating clean termination, in seconds.</description>
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job
!</keyword>



!#####################################################################!
!          Direct Quadrature Method of Moments (DQMOM)                !
!#####################################################################!


!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required=.FALSE.>
!  <description>variable to decide if the population balance equations are solved.</description>
      Call_DQMOM = .FALSE.
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required=.FALSE.>
!  <description>success-factor for aggregation.</description>
      AGGREGATION_EFF=0.D0
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required=.FALSE.>
!  <description>success-factor for breakage.</description>
      BREAKAGE_EFF=0.D0
!</keyword>








! ---------------------------------- questionable namelist entries below








!<keyword category="category name" required=.FALSE.>
!  <description>Variable which triggers an automatic restart.</description>
      AUTOMATIC_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required=.FALSE.>
!  <description>ATUO_RESTART counter.</description>
      ITER_RESTART = 1
!</keyword>



! NO_OF_RXNS is not a keyword. However, it is initialized here so that
! if there are no reactions, this value is assigned.
      NO_OF_RXNS = UNDEFINED_I


      U_G0 = UNDEFINED
      V_G0 = UNDEFINED
      W_G0 = UNDEFINED
      U_S0(:DIM_M) = UNDEFINED
      V_S0(:DIM_M) = UNDEFINED
      W_S0(:DIM_M) = UNDEFINED


      PHIP_OUT_ITER=0





      CALL DES_INIT_NAMELIST

      CALL QMOMK_INIT_NAMELIST

      CALL USR_INIT_NAMELIST

      CALL CARTESIAN_GRID_INIT_NAMELIST

      RETURN
      END SUBROUTINE INIT_NAMELIST
