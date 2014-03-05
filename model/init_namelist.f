!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: INIT_NAMELIST                                           C
!  Purpose: initialize the NAMELIST variables                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 26-NOV-91  C
!                                                                      C
!  Keyword Documentation Format:                                       C
!                                                                      C
!<keyword category="category name" required="true/false"               C
!                                    legacy="true/false">              C
!  <description></description>                                         C
!  <arg index="" id="" max="" min=""/>                                 C
!  <dependent keyword="" value="DEFINED"/>                             C
!  <conflict keyword="" value="DEFINED"/>                              C
!  <valid value="" note="" alias=""/>                                  C
!  <range min="" max="" />                                             C
!  MFIX_KEYWORD=INIT_VALUE                                             C
!</keyword>                                                            C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

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



!<keyword category="Run Control" required="true">
!  <description> Name used to create output files. The name should
!    generate legal file names after appending extensions.
!    Ex: Given the input, RUN_NAME = "bub01", MFIX will generate
!    the output files: BUB01.LOG, BUB01.OUT, BUB01.RES, etcs.
!  </description>
      RUN_NAME = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Problem description. Limited to 60 characters.</description>
      DESCRIPTION = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description> Simulation input/output units.</description>
!  <valid value="cgs" note="all input and output in cgs units (g, cm, s, cal)."/>
!  <valid value="si" note="all input and output in si units (kg, m, s, j)."/>
      UNITS = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="true">
!  <description>type of run.</description>
!  <valid value="new" note="new run."/>
!  <valid value="RESTART_1" note="Traditional restart."/>
!  <valid value="RESTART_2" 
!    note="Start a new run with initial conditions from a .RES file
!      created from another run."/>
      RUN_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Start-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TIME = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Stop-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TSTOP = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Starting time step. If left undefined, a steady-state 
!    calculation is performed.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT = UNDEFINED
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Maximum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MAX = ONE
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Minimum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MIN = 1.0D-6
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Factor for adjusting time step. Must be less than 1.
!  </description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="1" />
      DT_FAC = 0.9D0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description></description>
!  <valid value=".FALSE." 
!    note="Do not reduce time step for stalled iterations."/>
!  <valid value=".true." 
!    note="Reduce time step if the residuals sum does not decrease."/>
      DETECT_STALL = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Flag to restart the code when DT < DT_MIN. 
!  </description>
      AUTO_RESTART = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>momentum equations.</description>
!  <valid value=".false." note="model A"/>
!  <valid value=".true." note="model B"/>
      MODEL_B = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve X-momentum equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".true." note="Solve X-momentum equations."/>
!  <valid value=".false." note="Do not solve X-momentum equations."/>
      MOMENTUM_X_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve Y-momentum equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".TRUE." note="Solve Y-momentum equations."/>
!  <valid value=".FALSE." note="Do not solve Y-momentum equations."/>
      MOMENTUM_Y_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve Z-momentum equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".true." note="Solve Z-momentum equations."/>
!  <valid value=".false." note="Do not solve Z-momentum equations."/>
      MOMENTUM_Z_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve energy equations.</description>
!  <valid value=".false." note="Do not solve energy equations."/>
!  <valid value=".true." note="Solve energy equations."/>
      ENERGY_EQ = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solve species transport equations.</description>
!  <arg index="1" id="Phase" min="0" max="DIM_M"/>
!  <valid value=".false." note="Solve species equations."/>
!  <valid value=".true." note="Do not solve species equations."/>
      SPECIES_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Granular energy formulation selection.</description>
!  <valid value=".FALSE."
!    note="Use algebraic granular energy equation formulation."/>
!  <valid value=".TRUE."
!    note="Use granular energy transport equation (PDE) formulation."/>
      GRANULAR_ENERGY = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
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

!<keyword category="Run Control" required="false">
!  <description>Jenkins small frictional boundary condition.</description>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <dependent keyword="PHI_W" value="DEFINED"/>
!  <valid value=".FALSE." note=""/>
!  <valid value=".true." 
!    note="Use the Jenkins small fricational boundary condition."/>
      JENKINS = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Solids stress model selection.</description>
!  <valid value=".false." note="Use the Schaeffer solids stress model."/>
!  <valid value=".true." note="Use the Princeton solids stress model"/>
!  <dependent keyword="GRANULAR_ENERGY" value=".TRUE."/>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <dependent keyword="PHI_W" value="DEFINED"/>
      FRICTION = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    For a term appearing in the frictional stress model 
!    invoked with friction = .true.
!  </description>
!  <valid value="0" note="use s:s in the frictional stress model."/>
!  <valid value="2" note="an appropriate combination of the above two forms."/>
!  <valid value="1" note="use an alternate form suggested by savage."/>
!  <dependent keyword="friction" value=".true."/>
      SAVAGE = 1
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Schaeffer frictional stress tensor formulation. </description>
!  <dependent keyword="PHI" value="DEFINED"/>
!  <valid value=".TRUE." note="Use the Schaeffer model"/>
!  <valid value=".FALSE." note="Do not use the Schaeffer model."/>
!    If FRICTION=.FALSE. then the model will have no frictional viscosity."/>
      SCHAEFFER = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Blend the Schaeffer stresses with that of kinetic theory around Ep*.
!  </description>
      BLENDING_STRESS = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Hyperbolic tangent function for blending frictional stress models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".true."/>
!  <conflict keyword="SIGM_BLEND" value=".true."/>
      TANH_BLEND = .TRUE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    A scaled and truncated sigmoidal function for blending
!    frictional stress  models.
!  </description>
!  <dependent keyword="BLENDING_STRESS" value=".true."/>
!  <conflict keyword="TANH_BLEND" value=".true."/>
      SIGM_BLEND = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Correlation to compute maximum packing for polydisperse systems.
!  </description>
!  <valid value=".TRUE."
!    note="Use the Yu and Standish correlation."/>
!  <valid value=".FALSE."
!    note="Do not use the Yu and Standish correlation."/>
      YU_STANDISH = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    use fedors and landel correlation to compute maximum
!    packing for a binary (only) mixture of powders.
!  </description>
      FEDORS_LANDEL = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>call user-defined subroutines.</description>
!  <valid value=".true." note="call user-defined subroutines."/>
!  <valid value=".false." note="do not call user-defined subroutines."/>
      CALL_USR = .FALSE.
!</keyword>


!<keyword category="Run Control" required="false">
!  <description>
!    the number of user defined chemical reactions stored
!    in the *.spa file.  see section 4.10 chemical reactions.
!  </description>
      nRR = 0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    when activated the k-epsilon turbulence model (for single
!    -phase flow) is solved using standard wall functions.
!  </description>
!  <conflict keyword="L_SCALE0" value="DEFINED"/>
      K_Epsilon = .FALSE.
!</keyword>


!<keyword category="Run Control" required="false">
!  <description>drag model </description>
!  <valid value="syam_obrien" note=""/>
!  <valid value="gidaspow" note=""/>
!  <valid value="gidaspw_pcf" note=""/>
!  <valid value="gidaspow_blend" note=""/>
!  <valid value="gidaspow_blend_pcf" note=""/>
!  <valid value="wen_yu" note=""/>
!  <valid value="wen_yu_pcf" note=""/>
!  <valid value="koch_hill" note=""/>
!  <valid value="koch_hill_pcf" note=""/>
!  <valid value="bvk" note=""/>
!  <valid value="hys" note=""/>
      DRAG_TYPE = 'SYAM_OBRIEN'
!</keyword>



!<keyword category="Run Control" required="false">
!  <description> Sinlge particle drag correlation. Note that names
!    containing xRE are correlations multiplied by Reynods number.
!    A default for each DRAG_TYPE is set if left UNDEFINED.
!  </description>
!  <valid value="C_DS_SN"    note="Schiller and Naumann (1933)"/>
!  <valid value="C_DSxRE_DV" note="Dalla Valle (1948)"/>
!  <valid value="C_DS_DEL"   note="Dellion et al. (2005)"/>
!  <valid value="C_DSxRE_TL" note="Turton and Levenspiel (1986)"/>
      CD_FUNCTION = UNDEFINED_C
!</keyword>


!<keyword category="Run Control" required="false">
!  <description>
!    radial distribution function at contact for polydisperse systems.
!  </description>
!  <valid value="lebowitz" note=""/>
!  <valid value="modified_lebowitz" note=""/>
!  <valid value="mansoori" note=""/>
!  <valid value="modified_mansoori" note=""/>
      RDF_TYPE = 'LEBOWITZ'
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    when activated the added (or virtual) mass force effectively 
!    acts to increase the inertia of the dispersed phase, which
!    tends to stabilize simulations of bubbly gas-liquid flows.
!  </description>
      Added_Mass = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    subgrid models include: igci and milioli
!  </description>
!  <valid value="igci" note=""/>
!  <valid value="milioli" note=""/>
      SUBGRID_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>
!    Ratio of filter size to computational cell size.
!  </description>
      filter_size_ratio = 2.0D0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Flag for subgrid wall effects.</description>
!  <valid value=".FALSE." note="Do not include wall effects."/>
!  <valid value=".TRUE." note="Include subgrid wall effects."/>
      SUBGRID_Wall = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Provide detailed logging of negative density errors.</description>
!  <valid value=".FALSE." note="Do not log negative density errors."/>
!  <valid value=".TRUE." note="Log negative density errors."/>
      REPORT_NEG_DENSITY = .FALSE.
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>Number of user-defined scalar transport equations to solve.</description>
!  <range min="0" max="DIM_SCALAR" />
      NScalar = 0
!</keyword>

!<keyword category="Run Control" required="false">
!  <description>The phase convecting the indexed scalar transport equation.</description>
!  <arg index="1" id="Scalar Equation" min="0" max="DIM_SCALAR"/>
!  <range min="0" max="DIM_M" />
      Phase4Scalar(:DIM_SCALAR) = UNDEFINED_I
!</keyword>



!#####################################################################!
!                           Physcial Parameters                       !
!#####################################################################!

!<keyword category="Physical Parameters" required="false">
!  <description>User defined constants.</description>
      C(:DIMENSION_C) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>Name of user-defined constant. (20 character max)</description>
      C_NAME(:DIMENSION_C) = '....................'
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>coefficient of restitution for particle-particle collisions.</description>
      C_E = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    GHD Theory: Coefficient of restitution for particle-particle collisions.
!  </description>
      r_p(:DIM_M, :DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>coefficient of restitution for particle-wall collisions.</description>
      E_W = 1.D0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>specularity coefficient associated with particle-wall collisions.</description>
      PHIP = 0.6D0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    output the variable specularity coefficient when bc_jj_m is 
!    .true.. The specularity coefficient will be stored in reactionrates
!    array for post-processing by post-mfix. user needs to set nrr to 1
!    for this purpose. be careful with this setting when reacting flow 
!    is simulated.
!  </description>
      PHIP_OUT_JJ=.false.
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    specify the value of specularity coefficient when the
!    normalized slip velocity goes to zero when bc_jj_m is
!    .true.. this variable is calculated internally in the
!    code. do not modify unless an accurate number is known.
!  </description>
!  <dependents>bc_jj_m</dependents>
      phip0 = undefined
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    coefficient of friction between the particles of two solids phases.
!  </description>
      C_F = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!     angle of internal friction (in degrees). set this value
!     to zero to turn off plastic regime stress calculations.
!  </description>
      PHI = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Angle of internal friction (in degrees) at walls. Set this
!    value to non-zero (phi_w = 11.31 means tan_phi_w = mu = 0.2)
!    when using jenkins or bc_jj_m boundary condition.
!  </description>
      PHI_W = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    minimum solids fraction above which friction sets in.  (when
!    friction = .true.)</description>
!  <dependents>friction</dependents>
      EPS_F_MIN = 0.5D0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Maximum solids volume fraction at packing for polydisperse 
!    systems (more than one solids phase used). The value of 
!    EP_star may change during the computation if solids phases
!    with different particle diameters are specified and 
!    Yu_Standish or Fedors_Landel correlations are used.
!  </description>
!  <range min="0" max="1-ep_star" />
      EP_S_MAX(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    used in calculating the initial slope of segregation: see 
!    Gera et al. (2004) - recommended value 0.3. increasing this
!    coefficient results in decrease in segregation of particles
!    in binary mixtures.
!  </description>
      SEGREGATION_SLOPE_COEFFICIENT=0.D0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    value of turbulent length initialized. this may be overwritten 
!    in specific regions with the keyword ic_l_scale.
!</description>
!  <conflict keyword="K_EPSILON" value=".TRUE."/>
      L_SCALE0 = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    maximum value of the turbulent viscosity of the fluid.
!  </description>
      MU_GMAX = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>excluded volume in boyle-massoudi stress.</description>
!  <valid value="0.0" note="b-m stress is turned off."/>
      V_EX = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>reference pressure.</description>
      P_REF = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>scale factor for pressure.</description>
      P_SCALE = ONE
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>gravitational acceleration. </description>
      GRAVITY = UNDEFINED
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    X-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction. 
!  </description>
      GRAVITY_X = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Y-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Y = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Z-component of gravitational acceleration vector. By default, the
!    gravity force acts in the negative y-direction.
!  </description>
      GRAVITY_Z = ZERO
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This are determined using the Umf spreadsheet.
!  </description>
      drag_c1 = 0.8d0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    Quantity for calibrating Syamlal-O'Brien drag correlation using Umf
!    data.  This are determined using the Umf spreadsheet.
!  </description>
      drag_d1 = 2.65d0
!</keyword>

!<keyword category="Physical Parameters" required="false">
!  <description>
!    If use_def_lam_hys is set to .false. the user is able to specify a
!    value for the lubrication cutoff distance (lam_hys).  in practice
!    this number should be on the order of the mean free path of the
!    gas for smooth particles, or the rms roughness of a particle if
!    they are rough (if particle roughness is larger than the mean 
!   free path).
!  </description>
!  <dependents>USE_DEF_LAM_HYS</dependents>
      LAM_HYS = UNDEFINED
!</keyword>


!<keyword category="Physical Parameters" required="false">
!  <description>disperse phase number where the added mass applies.</description>
      M_AM = UNDEFINED_I
!</keyword>




!#####################################################################!
!                          Numerical Parameters                       !
!#####################################################################!



!<keyword category="Numerical Parameters" required="false">
!  <description>maximum number of iterations.</description>
      MAX_NIT = 500
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>factor to normalize the gas continuity equation residual.</description>
      NORM_G = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>factor to normalize the solids continuity equation residual.</description>
      NORM_S = UNDEFINED
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (continuity+momentum).</description>
      TOL_RESID = 1.0D-3
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (granular energy).</description>
      TOL_RESID_Th = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (energy).</description>
      TOL_RESID_T = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (species balance).</description>
      TOL_RESID_X = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (scalar balances.)</description>
      TOL_RESID_Scalar = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>maximum residual at convergence (K_Epsilon Model)</description>
      TOL_RESID_K_Epsilon = 1.0D-4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description> Minimum residual for declaring divergence.
!    This parameter is useful for incompressible fluid simulations because
!    velocity residuals can take large values for the second iteration (e.g., 1e+8)
!    before drop down to smaller values for third  third iteration (e.g., 0.1). 
!  </description>
      TOL_DIVERGE = 1.0D+4
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    The code declares divergence if the velocity anywhere in the domain
!    exceeds a maximum value.  this maximum value is automatically 
!    determined from the boundary values. the user may scale the maximum
!    value by adjusting this scale factor.
!  </description>
      MAX_INLET_VEL_FAC = ONE
!</keyword>

!<keyword category="Numerical Parameters" required="false">
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

!<keyword category="Numerical Parameters" required="false">
!  <description>LEQ Solver selection.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="1" note="SOR"/>
!  <valid value="2" note="BiCGSTAB"/>
!  <valid value="3" note="GMRES"/>
!  <valid value="5" note="CG"/>
      LEQ_METHOD(:) = 2
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>Linear equation sweep direction.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="RSRS" note="(Red/Black Sweep, Send Receive) repeated twice"/>
!  <valid value="ISIS" note="(Sweep in I, Send Receive) repeated twice"/>
!  <valid value="JSJS" note="(Sweep in J, Send Receive) repeated twice"/>
!  <valid value="KSKS" note="(Sweep in K, Send Receive) repeated twice"/>
!  <valid value="ASAS" note="(All Sweep, Send Receive) repeated twice"/>
      LEQ_SWEEP(:) = 'RSRS'
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>Linear Equation tolerance.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <dependent keyword="LEQ_METHOD" value="2"/>
!  <dependent keyword="LEQ_METHOD" value="3"/>
      LEQ_TOL(:) = 1.0D-4
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>Linear precondition used for LEQ solver sweeps.</description>
!  <arg index="1" id="Equation ID Number" min="1" max="9"/>
!  <valid value="NONE" note="No preconditioner"/>
!  <valid value="LINE" note="Line relaxation"/>
!  <valid value="DIAG" note="Diagonal Scaling"/>
      LEQ_PC(:) = 'LINE'
!</keyword>

!<keyword category="Numerical Parameters" required="false">
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

!<keyword category="run control" required="false">
!  <description>
!    Use deferred correction method for implementing higher order 
!    discretization.
!  </description>
!  <valid value=".false." note="use down-wind factor method (default)."/>
!  <valid value=".true." 
!    note="use deferred correction method for implementing higher order discretization."/>
      DEF_COR  =  .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
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

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    four point fourth order interpolation and is upstream biased. if 
!    this scheme is chosen and discretize(*) < 2, discretize(*) is 
!    defaulted to 2. if you chose this scheme, set the c_fac value
!    between 0 and 1.
!  </description>
!  <dependent keyword="c_fac" value="DEFINED"/>
      FPFOI = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    factor used in the universal limiter (when fpfoi is set .true.) and 
!    can be any value in the set (0,1). the choice of 1 will give 
!    (diffusion) first order upwinding and as this value becomes closer
!    to 0 the scheme becomes more compressive.
!  </description>
!  <range min="0.0" max="1.0" />
!  <dependent keyword="fpfoi" value=".true."/>
      C_FAC = UNDEFINED
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>temporal discretization scheme.</description>
!  <valid value=".false."
!    note="Implicit Euler based temporal discretization scheme employed 
!      (first order accurate in time)."/>
!  <valid value=".true." 
!    note="crank-nicholson based temporal discretization scheme employed
!      (second order accurate in time excluding the restart timestep
!      which is first order)."/>
      CN_ON = .FALSE.
!</keyword>


!<keyword category="Numerical Parameters" required="false">
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

!<keyword category="Numerical Parameters" required="false">
!  <description>
!    The implicitness calculation of the gas-solids drag coefficient
!    may be underrelaxed by changing ur_f_gs, which takes values 
!    between 0 to 1:
!  </description>
!  <range min="0" max="1" />
      UR_F_gs = 1.0D0
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Unknown</description>
      UR_Kth_sml = 1.0D0
!</keyword>


!<keyword category="Numerical Parameters" required="false">
!  <description>Solve transpose of linear system. (BICGSTAB ONLY).</description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      DO_TRANSPOSE = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Frequency to check for convergence. (BICGSTAB ONLY)</description>
!  <dependent keyword="LEQ_METHOD" value="2"/>
      icheck_bicgs = 1
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Sets optimal LEQ flags for parallel runs.</description>
      OPT_PARALLEL = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Use do-loop assignment over direct vector assignment.</description>
      USE_DOLOOP = .FALSE.
!</keyword>

!<keyword category="Numerical Parameters" required="false">
!  <description>Calculate dot-products more efficiently (Serial runs only.)</description>
      IS_SERIAL = .TRUE.
!</keyword>


!#####################################################################!
!                      Geometry and Discretization                    !
!#####################################################################!


!<keyword category="Geometry and Discretization" required="false">
!  <description>coordinates used in the simulation.</description>
!  <valid value="cartesian" note="cartesian coordinates."/>
!  <valid value="cylindrical" note="cylindrical coordinates."/>
      COORDINATES = UNDEFINED_C
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(do not use.)</description>
!  <valid value=".false." note="x (r) direction is considered."/>
!  <valid value=".true." note="x (r) direction is not considered."/>
      NO_I = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>number of cells in the x (r) direction.</description>
      IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the x (r) direction. Enter values from DX(0) to 
!    DX(IMAX-1). (Use uniform mesh size with higher-order 
!    discretization methods.  Also in cylindrical coordinates DX
!    should be kept uniform for strict momentum conservation.)
!  </description>
      DX(:DIM_I) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    The inner radius in the simulation of an annular cylindrical region.
!  </description>
      XMIN = ZERO
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>reactor length in the x (r) direction.</description>
      XLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>(do not use.)</description>
!  <valid value=".false" note="y direction is considered."/>
!  <valid value=".true." note="y direction is not considered."/>
      NO_J = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>number of cells in the y direction.</description>
      JMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the y direction. Enter values from DY(0) to 
!    DY(IMAX-1). (Use uniform mesh size with second-order 
!    discretization methods.)
!  </description>
      DY(:DIM_J) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>reactor length in the y direction.</description>
      YLENGTH = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description></description>
!  <valid value=".false." note="z(theta) direction is considered."/>
!  <valid value=".true." note="z(theta) direction is not considered."/>
      NO_K = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>number of cells in the z (() direction.</description>
      KMAX = UNDEFINED_I
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Cell sizes in the z (theta) direction. Enter values from DZ(0) to 
!    DZ(IMAX-1). (Use uniform mesh size with second-order discretization
!    methods.)
!  </description>
      DZ(:DIM_K) = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>reactor length in the z (theta) direction.</description>
      ZLENGTH = UNDEFINED
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Flag for making the x-direction cyclic without pressure drop. no other 
!    boundary conditions for the x-direction should be specified.
!</description>
!  <valid value=".false." note="no cyclic condition at x-boundary."/>
!  <valid value=".true." note="cyclic condition at x-boundary."/>
      CYCLIC_X = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    flag for making the x-direction cyclic with pressure drop. if the 
!    keyword flux_g is given a value this becomes a cyclic boundary 
!    condition with specified mass flux. no other boundary conditions
!    for the x-direction should be specified.
!  </description>
!  <valid value=".false." note="no cyclic condition at x-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at x-boundary."/>
      CYCLIC_X_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    fluid pressure drop across xlength when a cyclic boundary condition
!    with pressure drop is imposed in the x-direction.
!  </description>
      DELP_X = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    flag for making the y-direction cyclic without pressure drop. no
!    other boundary conditions for the y-direction should be specified.
!  </description>
!  <valid value=".false." note="no cyclic condition at y-boundary."/>
!  <valid value=".true." note="cyclic condition at x-boundary."/>
      CYCLIC_Y = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    flag for making the y-direction cyclic with pressure drop. if the
!    keyword flux_g is given a value this becomes a cyclic boundary 
!    condition with specified mass flux. no other boundary conditions 
!    for the y-direction should be specified.
!  </description>
!  <valid value=".false." note="no cyclic condition at y-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at y-boundary."/>
      CYCLIC_Y_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    fluid pressure drop across ylength when a cyclic boundary condition
!    with pressure drop is imposed in the y-direction.
!  </description>
      DELP_Y = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    flag for making the z-direction cyclic without pressure drop. no
!    other boundary conditions for the z-direction should be specified.
!  </description>
!  <valid value=".false." note="no cyclic condition at z-boundary."/>
!  <valid value=".true." note="cyclic condition at z-boundary."/>
      CYCLIC_Z = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    flag for making the z-direction cyclic with pressure drop. if the 
!    keyword flux_g is given a value this becomes a cyclic boundary 
!    condition with specified mass flux. no other boundary conditions
!    for the z-direction should be specified.
!  </description>
!  <valid value=".false." note="no cyclic condition at z-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at z-boundary."/>
      CYCLIC_Z_PD = .FALSE.
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    fluid pressure drop across zlength when a cyclic boundary condition 
!    with pressure drop is imposed in the z-direction.
!  </description>
      DELP_Z = UNDEFINED
!</keyword>

!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    if .true. imposes a mean shear on the flow field as a linear 
!    function of _x_ coordinate. this feature should only be used when
!    cyclic_x=.true. also, the keyword v-sh needs to be set.
!  </description>
!  <dependents>cyclic_x v_sh</dependents>
      SHEAR = .FALSE.
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
!  <description>
!    Specifies the mean _y_ velocity component at the eastern boundary 
!    of the domain (v_sh), and the mean _y_ velocity (-v_sh) at the 
!    western boundary of the domain.
!  </description>
      V_sh=0d0
!</keyword>


!<keyword category="Geometry and Discretization" required="false">
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




!<keyword category="Gas Phase" required="false">
!  <description>
!    Specified constant gas density. this value may be set to zero to 
!    make the drag zero and to simulate granular flow in a vacuum. for
!    this case, users may turn off solving for gas momentum equations 
!    to accelerate convergence.
!  </description>
      RO_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>specified constant gas viscosity.</description>
      MU_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>specified constant gas conductivity.</description>
      K_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>specified constant gas specific heat.</description>
      C_PG0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>specified constant gas diffusivity.</description>
      DIF_G0 = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>average molecular weight of gas.</description>
      MW_AVG = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>molecular weight of gas species n.</description>
      MW_G(:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>Number of species comprising the gas phase.</description>
      NMAX_g = UNDEFINED_I
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>Name of gas phase species n as it appears in the materials database.</description>
      SPECIES_g = UNDEFINED_C
!</keyword>

!<keyword category="Gas Phase" required="false">
!  <description>User defined name for gas phase species n.</description>
      SPECIES_ALIAS_g = UNDEFINED_C
!</keyword>



!#####################################################################!
!                            Solids Phase                             !
!#####################################################################!


!<keyword category="Solids Phase" required="false">
!  <description>number of solids phases.</description>
      MMAX = 1
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Initial particle diameters.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
      D_P0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Specified constant solids density.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
      RO_S0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Baseline species mass fraction.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      X_s0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Specified constant solids species density.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="INERT_SPECIES" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      RO_Xs0(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Index of inert solids phase species.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
!  <dependent keyword="X_s0" value="DEFINED"/>
!  <dependent keyword="RO_Xs0" value="DEFINED"/>
!  <conflict keyword="RO_s0" value="DEFINED"/>
      INERT_SPECIES(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>specified constant granular viscosity. if this value is
!    specified, then the kinetic theory calculation is turned off and
!    p_s = 0 and lambda_s = -2/3 mu_s0.
!  </description>
      MU_S0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>specified constant solids conductivity.</description>
      K_S0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>specified constant solids specific heat.</description>
      C_PS0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>specified constant solids diffusivity.</description>
      DIF_S0 = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Molecular weight of solids phase-m, species n.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
      MW_S(:DIM_M,:DIM_N_s) = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>Number of species comprising solids phase m.</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
      NMAX_s(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>
!    Name of solids phase m, species n as it appears in the materials 
!    database.
!</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
      SPECIES_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>User defined name for solids phase m, species n</description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <arg index="2" id="Species index" min="1" max="DIM_N_s"/>
      SPECIES_ALIAS_s(:DIM_M,:DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>packed bed void fraction.</description>
      EP_STAR = UNDEFINED
!</keyword>

!<keyword category="Solids Phase" required="false">
!  <description>
!    Indicates that the solids phase forms a packed bed with a void
!    fraction ep_star.
!  </description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
      CLOSE_PACKED(:DIM_M) = .TRUE.
!</keyword>


!<keyword category="Solids Phase" required="false">
!  <description>specified constant granular viscosity. if this value is
!    specified, then the kinetic theory calculation is turned off and
!    p_s = 0 and lambda_s = -2/3 mu_s0.
!  </description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <valid value='TFM' description='Two-fluid Model (coninuum).' />
!  <valid value='DEM' description='Discrete Element Model' />
!  <valid value='MPPIC' description='Multiphase-Particle in Cell' />
      SOLIDS_MODEL(:DIM_M) = 'TFM'
!</keyword>


!<keyword category="Solids Phase" required="false">
!  <description>Partice shape factor. </description>
!  <arg index="1" id="Solids phase index" min="1" max="DIM_M"/>
!  <dependent keyword="CD_FUNCTION" value="C_DS_DEL"/>
      PSI_s(:DIM_M) = UNDEFINED
!</keyword>


!#####################################################################!
!                         Initial Conditions                          !
!#####################################################################!


      DO LC = 1, DIMENSION_IC

!<keyword category="Initial Condition" required="false">
!  <description>x coordinate of the west face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>x coordinate of the east face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>y coordinate of the south face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>y coordinate of the north face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>z coordinate of the bottom face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>z coordinate of the top face.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>i index of the west-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>i index of the east-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>j index of the south-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>j index of the north-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>k index of the bottom-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>k index of the top-most wall.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Type of initial condition. Mainly used in restart runs to overwrite
!    values read from the .RES file by specifying it as _PATCH_. The 
!    user needs to be careful when using the _PATCH_ option, since the
!    values from the .RES file are overwritten and no error checking is
!    done for the patched values.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial gas pressure in the IC region. If this quantity is not
!    specified, MFIX will set up a hydrostatic pressure profile, 
!    which varies only in the y-direction.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids pressure in the IC region. Usually, this value is 
!    specified as zero.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_P_STAR(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Turbulence length scale in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_L_SCALE(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial bulk density (rop_s = ro_s x ep_s) of solids phase-m in the
!    IC region. Users need to specify this IC only for polydisperse flow
!    (MMAX > 1). Users must make sure that summation of ( IC_ROP_s(ic,m)
!    / RO_s(m) ) over all solids phases is equal to ( 1.0 – IC_EP_g(ic)).
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Initial solids volume fraction of solids phase-m in the IC region.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m granular temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Gas phase radiation coefficient in the IC region. Modify file 
!    rdtn2.inc to change the source term.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RG(LC) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Gas phase radiation temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RG(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>
!    Solids phase-m radiation coefficient in the IC region. Modify file
!    radtn2.inc to change the source term.
!  </description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_GAMA_RS(LC,:DIM_M) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Solids phase-m radiation temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_T_RS(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         IC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial mass fraction of gas species n.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="species" min="1" max="DIM_N_g"/>
         IC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial mass fraction of gas species n.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_s"/>
         IC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of Scalar n.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
        IC_SCALAR(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of K in K-Epsilon.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of Epsilon in K-Epsilon.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Flag for inflating IC region to cover full domain.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
          IC_DES_FIT_TO_REGION(LC) = .FALSE. 
!</keyword>

      ENDDO




!#####################################################################!
!                        Boundary Conditions                          !
!#####################################################################!
      DO LC = 1, DIMENSION_BC


!<keyword category="Boundary Condition" required="false">
!  <description>x coordinate of the west face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>x coordinate of the east face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y coordinate of the south face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y coordinate of the north face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z coordinate of the top face or edge.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>i index of the west-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>i index of the east-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>j index of the south-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>j index of the north-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>k index of the bottom-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>k index of the top-most cell.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <valid value='DUMMY' 
!    description='The specified boundary condition is ignored. This is
!      useful for turning off some boundary conditions without having 
!      to delete them from the file.' />
!  <valid value='MASS_INFLOW' alias='MI'
!    description='Mass inflow rates for gas and solids phases are 
!      specified at the boundary.'/>
!  <valid value='MASS_OUTFLOW' alias='MO'
!    description='The specified values of gas and solids mass outflow
!      rates at the boundary are maintained, approximately. This 
!      condition should be used sparingly for minor outflows, when 
!      the bulk of the outflow is occurring through other constant
!      pressure outflow boundaries.' />
!  <valid value='P_INFLOW' alias='PI'
!    description='Inflow from a boundary at a specified constant 
!      pressure. To specify as the west, south, or bottom end of 
!      the computational region, add a layer of wall cells to the
!      west, south, or bottom of the PI cells. Users need to specify
!      all scalar quantities and velocity components. The specified 
!      values of fluid and solids velocities are only used initially
!      as MFIX computes these values at this inlet boundary.' />
!  <valid value='P_OUTFLOW' alias='PO'
!    description='Outflow to a boundary at a specified constant pressure.
!      To specify as the west, south, or bottom end of the computational
!      region, add a layer of wall cells to the west, south, or bottom of
!      the PO cells.' />
!  <valid value='FREE_SLIP_WALL' alias='FSW'
!    description='Velocity gradients at the wall vanish. If BC_JJ_PS is
!      equal to 1, the Johnson-Jackson boundary condition is used for 
!      solids.  A FSW is equivalent to using a PSW with hw=0.' />
!  <valid value='NO_SLIP_WALL' alias='NSW'
!    description='All components of the velocity vanish at the wall. If
!      BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is
!      used for solids.  A NSW is equivalent to using a PSW with vw=0 
!      and hw undefined.' />
!  <valid value='PAR_SLIP_WALL' alias='PSW'
!    description='Partial slip at the wall implemented as 
!      dv/dn + hw (v – vw) = 0, where n is the normal pointing from the
!      fluid into the wall. The coefficients hw and vw should be
!      specified. For free slip set hw = 0. For no slip leave hw
!      undefined (hw=∞) and set vw = 0. To set hw = ∞, leave it
!      unspecified. If BC_JJ_PS is equal to 1, the Johnson-Jackson
!      boundary condition is used for solids.' />
         BC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_UW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_UW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_VW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_WW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_WW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Johnson and Jackson partial slip bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <valid value='0'
!    description='Do not use Johnson and Jackson partial slip bc. If 
!      granular energy transport equation is not solved./>
!  <valid value='1' 
!    description='Use Johnson and Jackson partial slip bc. If granular
!      energy transport equation is solved. />
         BC_JJ_PS(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall value, THETAw_M, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Transfer coefficient, Hw, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Hw for granular energy bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant flux, C, in diffusion boundary condition:
!    d(Theta_M)/dn + Hw (THETA_M - THETAw_M) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified gas phase wall temperature, Tw_g, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_TW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas phase heat flux, C, in diffusion boundary condition:
!    d(T_g)/dn + Hw (T_g - Tw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Solids phase heat transfer coefficient, Hw, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Solids phase hw for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_HW_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified solids phase wall temperature, Tw_s, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
         BC_TW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant solids phase heat flux, C, in diffusion boundary condition:
!    d(T_s)/dn + Hw (T_s - Tw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Gas phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified wall gas species mass fraction, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant gas species mass flux, C, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Solid phase species mass transfer coefficient, Hw, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_HW_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified solids species mass fraction at the wall, Xw, in diffusion boundary condition:
!    d(X_g)/dn + Hw (X_g - Xw_g) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_XW_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant solids species mass flux, C, in diffusion boundary condition:
!    d(X_s)/dn + Hw (X_s - Xw_s) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_S"/>
         BC_C_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Scalar transfer coefficient, Hw, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_HW_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified scalar value at the wall, ScalarW, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_ScalarW(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>
!    Specified constant scalar flux, C, in diffusion boundary condition:
!    d(Scalar)/dn + Hw (Scalar - ScalarW) = C, where n is the fluid-to-wall normal.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_C_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Bulk density of solids phase at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volume fraction at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_EP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m granular temperature at the BC plane.</description>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of gas species n at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Species" min="1" max="DIM_N_G"/>
         BC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of solids phase-m, species n at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Species" min="1" max="DIM_N_s"/>
         BC_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>x-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>x-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Solid phase" min="1" max="DIM_M"/>
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval at the beginning when the normal
!    velocity at the boundary is equal to BC_Jet_g0. When restarting,
!    run this value and BC_Jet_g0 should be specified such that the
!    transient jet continues correctly. MFIX does not store the jet
!    conditions. For MASS_OUTFLOW boundary conditions, BC_DT_0 is 
!    the time period to average and print the outflow rates. The 
!    adjustment of velocities to get a specified mass or volumetric
!    flow rate is based on the average outflow rate.
!  </description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_DT_0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the initial interval BC_DT_0.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_JET_G0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_Jet_gh.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_DT_H(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_h.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_JET_GH(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_JET_gL.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_DT_L(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_L.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_JET_GL(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value for user-defined scalar equation.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Scalar Eq." min="1" max="DIM_SCALAR"/>
         BC_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value of K for K-Epsilon Equation.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_K_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Boundary value of Epsilon for K-Epsilon Equation.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_E_Turb_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specifed boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VELMAG_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Magnitude of gas velocity in a specifed boundary region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <arg index="2" id="Phase" min="1" max="DIMENSION_BC"/>
         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED
!</keyword>

         BC_APPLY_TO_MPPIC(LC) = .true.
         BC_ROP_G(LC) = UNDEFINED
      ENDDO

!<keyword category="Boundary Condition" required="false">
!  <description>Use the modified Johnson and Jackson partial slip BC 
!    with variable specularity coefficient. 
!  </description>
!  <dependent keyword="e_w" value="DEFINED"/>
!  <dependent keyword="phi_w" value="DEFINED"/>
         BC_JJ_M = .false.
!</keyword>




!#####################################################################!
!                         Internal Surfaces                           !
!#####################################################################!
      DO LC = 1, DIMENSION_IS


!<keyword category="Internal Surface" required="false">
!  <description>x coordinate of the west face or edge.</description>
         IS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>x coordinate of the east face or edge.</description>
         IS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>y coordinate of the south face or edge</description>
         IS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>y coordinate of the north face or edge</description>
         IS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>z coordinate of the bottom face or edge</description>
         IS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>z coordinate of the top face or edge</description>
         IS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>i index of the west-most cell.</description>
         IS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>i index of the east-most cell</description>
         IS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>j index of the south-most cell</description>
         IS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>j index of the north-most cell</description>
         IS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>k index of the bottom-most cell</description>
         IS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>k index of the top-most cell</description>
         IS_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Type of internal surface</description>
!  <valid value="IMPERMEABLE" 
!    note="No gas or solids flow through the surface." alias="IP"/>
!  <valid value="SEMIPERMEABLE" aliase='SP'
!    note="Gas flows through the surface with an additional resistance.
!      Solids velocity through the surface is set to zero or to a user-
!      specified fixed value (i.e., solids momentum equation for this
!      direction is not solved)." />
         IS_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>permeability</description>
         IS_PC(LC,1) = LARGE_NUMBER
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Inertial resistance coefficient.</description>
         IS_PC(LC,2) = ZERO
!</keyword>


!<keyword category="Internal Surface" required="false">
!  <description>Value of fixed solids velocity through semipermeable surfaces.</description>
         IS_VEL_S(LC,:DIM_M) = ZERO
!</keyword>
      ENDDO


!#####################################################################!
!                     Point Source Mass Inlets                        !
!#####################################################################!
      DO LC = 1, DIMENSION_PS

!<keyword category="Point Source" required="false">
!  <description>x coordinate of the west face or edge.</description>
         PS_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>x coordinate of the east face or edge.</description>
         PS_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y coordinate of the south face or edge.</description>
         PS_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y coordinate of the north face or edge.</description>
         PS_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z coordinate of the bottom face or edge.</description>
         PS_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z coordinate of the top face or edge.</description>
         PS_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>i index of the west-most cell.</description>
         PS_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>i index of the east-most cell.</description>
         PS_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>j index of the south-most cell.</description>
         PS_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>j index of the north-most cell.</description>
         PS_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>k index of the bottom-most cell.</description>
         PS_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>k index of the top-most cell.</description>
         PS_K_T(LC) = UNDEFINED_I
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>x-component of incoming gas velocity.</description>
         PS_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y-component of incoming gas velocity.</description>
         PS_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z-component of incoming gas velocity.</description>
         PS_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas mass flow rate through the point source.</description>
         PS_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming gas.</description>
         PS_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas phase incoming species n mass fraction.</description>
         PS_X_G(LC,:DIM_N_g) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>x-component of incoming solids velocity.</description>
         PS_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y-component of incoming solids velocity.</description>
         PS_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z-component of incoming solids velocity.</description>
         PS_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids mass flow rate through the point source.</description>
         PS_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming solids.</description>
         PS_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids phase incoming species n mass fraction.</description>
         PS_X_S(LC,:DIM_M,:DIM_N_S) = UNDEFINED
!</keyword>

      ENDDO


!#####################################################################!
!                          Output Control                             !
!#####################################################################!

!<keyword category="Output Control" required="false">
!  <description>interval at which restart (.res) file is updated.</description>
      RES_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>interval at which .spx files are updated.
!    SP1: void fraction (EP_G)
!    SP2: Gas pressure (P_G), and Solids pressure (P_star)
!    SP3: Gas veloicty (U_G, V_G, W_G)
!    SP4: Solids veloicty (U_S, V_S, W_S)
!    SP5: Solids builk density (ROP_s)
!    SP6: Gas and solids temperature (T_G, T_S)
!    SP7: Gas and solids mass fractions (X_G, X_S)
!    SP8: Granular temperature (THETA_M)
!    SP9: User defined scalars. (SCALAR)
!    SPA: Reaction Rates (ReactionRates)
!    SPB: Turbulence quantities (K_TURB_G, E_TURB_G)
!  </description>
      SPX_DT(:N_SPX) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Interval at which standard output (.OUT) file is updated.
!    Only run configuration information is written if left undefined. Otherwise
!    all field variables for the entire domain are written in ASCII 
!    format to the .OUT file at OUT_DT intervals.
!  </description>
      OUT_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Number of time steps between .LOG file updates.</description>
      NLOG = 25
!</keyword>

!<keyword category="Output Control" required="false">
!  <description> Display the residuals on the screen and provide
!    messages about convergence on the screen and in the .LOG file.
!  </description>
      FULL_LOG = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Specifies the residuals to display.
!    P0  : Gas pressure
!    PM  : Solids phase M pressure
!    R0  : Gas density
!    RM  : Solids phase M density
!    U0  : Gas phase U-velocity
!    V0  : Gas phase V-velocity
!    W0  : Gas phase W-velocity
!    Um  : Solids phase M U-velocity
!    Vm  : Solids phase M V-velocity
!    Wm  : Solids phase M W-velocity
!    T0  : Gas temperature
!    Tm  : Solids phase M temperature
!    X0NN: Gas phase species NN mass fraction
!    XMNN: Solids phase M species NN mass fraction
!    K0  : K-Epsilon model residuals
!  </description>
!  <arg index="1" id="Residual Index" max="8" min="1"/>
      RESID_STRING(:8) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Display residuals by equation.  </description>
      GROUP_RESID = .FALSE.
!</keyword>


      DO LC=1, DIMENSION_USR
!<keyword category="Output Control" required="false">
!  <description> Intervals at which subroutine write_usr1 is called. </description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_DT(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: x coordinate of the west face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_X_W(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: x coordinate of the east face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_X_E(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: y coordinate of the south face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_Y_S(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: y coordinate of the north face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_Y_N(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: z coordinate of the bottom face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_Z_B(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: z coordinate of the top face or edge.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_Z_T(LC) = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: i index of the west-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_I_W(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: i index of the east-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_I_E(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: j index of the south-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_J_S(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: j index of the north-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_J_N(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: k index of the bottom-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_K_B(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: k index of the top-most cell.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_K_T(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: Type of user-defined ouput: Binary of ASCII.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: 
!    Variables to be written in the user-defined output files.
!  </description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_VAR(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: 
!    Format for writing user-defined (ASCII) output file.
!  </description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_FORMAT(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>UDF Hook: File extension for the user-defined output.</description>
!  <arg index="1" id="USR Index" max="DIMENSION_USR" min="1"/>
         USR_EXT(LC) = UNDEFINED_C
!</keyword>
      ENDDO


!<keyword category="Output Control" required="false">
!  <description> Frequency to perform an overall species mass balance.
!    Leaving undefined suppresses the mass balance calculations which
!    can slightly extend run time.
!  </description>
      REPORT_MASS_BALANCE_DT = UNDEFINED
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Use distributed IO :: Each rank generates RES/SPx files.</description>
!  <dependent keyword="SPECIES_EQ" value=".TRUE."/>
      bDist_IO = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
!  <description>Restart a unified IO run as distributed IO.</description>
!  <dependent keyword="RUN_TYPE" value="RESTART_2"/>
!  <dependent keyword="bDist_IO" value=".TRUE."/>
      bStart_with_one_RES = .FALSE.
!</keyword>

!<keyword category="Output Control" required="false">
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


!<keyword category="Chemical Reactions" required="false">
!  <description>Flag to use stiff chemistry solver (Direct Integration).</description>
!  <conflict keyword="USE_RRATES" value=".TRUE."/>
      STIFF_CHEMISTRY = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false">
!  <description>Flag to use legacy chemcial reaction UDFs.</description>
      USE_RRATES = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false" legacy="true">
!  <description>
!    Names of gas and solids phase species as it appears in the 
!    materials database. The first NMAX(0) are the names of gas
!    species. The next NMAX(1) are the names of solids phase-1 
!    species, etc.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      SPECIES_NAME(:DIM_N_ALL) = UNDEFINED_C
!</keyword>

!<keyword category="Chemical Reactions" required="false">
!  <description>
!    Number of species in phase m. Note that the gas phase is indicated
!    as m=0.
!  </description>
!  <dependent keyword="USE_RRATES" value=".TRUE."/>
      NMAX = UNDEFINED_I
!</keyword>

!<keyword category="Chemical Reactions" legacy="true">
!  <description>Flag for previous stiff solver.</description>
      CALL_DI   = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false" legacy="true">
!  <description>
!   Flag to specify variable solids diameter in original stiff chem solver.
!   (Non-Functional, Removed in 2013-2 Release)
!  </description>
      CALL_GROW = .FALSE.
!</keyword>

!<keyword category="Chemical Reactions" required="false" legacy="true">
!  <description>
!   Flag to use ISAT tables with original DI solver.
!   (Non-Functional, Removed in 2013-2 Release)
!  </description>
      CALL_ISAT = .FALSE.      ! Legacy Keyword
!  </description>
!</keyword>

!<keyword category="Chemical Reactions" required="false" legacy="true">
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


!<keyword category="Parallelization Control" required="false">
!  <description>number of grid blocks in x-direction.</description>
      NODESI = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>number of grid blocks in y-direction.</description>
      NODESJ = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>number of grid blocks in z-direction.</description>
      NODESK = UNDEFINED_I
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print out additional statistics for parallel runs</description>
      solver_statistics = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Group residuals to reduce global collectives.</description>
      DEBUG_RESID = .TRUE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>All ranks write error messages.</description>
      ENABLE_DMP_LOG = .FALSE.
!</keyword>

!<keyword category="Parallelization Control" required="false">
!  <description>Print the index layout for debugging.</description>
      DBGPRN_LAYOUT = .FALSE.
!</keyword>


!#####################################################################!
!                       Batch Queue Environment                       !
!#####################################################################!


!<keyword category="Batch Queue Environment" required="false">
!  <description>enables clean termination feature.</description>
      CHK_BATCHQ_END = .FALSE.
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>total wall-clock duration of the job, in seconds.</description>
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
!</keyword>

!<keyword category="Batch Queue Environment" required="false">
!  <description>buffer time when initiating clean termination, in seconds.</description>
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job
!</keyword>



!#####################################################################!
!          Direct Quadrature Method of Moments (DQMOM)                !
!#####################################################################!


!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>variable to decide if the population balance equations are solved.</description>
      Call_DQMOM = .FALSE.
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>success-factor for aggregation.</description>
      AGGREGATION_EFF=0.D0
!</keyword>

!<keyword category="Direct Quadrature Method of Moments (DQMOM)" required="false">
!  <description>success-factor for breakage.</description>
      BREAKAGE_EFF=0.D0
!</keyword>








! ---------------------------------- questionable namelist entries below








!<keyword category="category name" required="false">
!  <description>Variable which triggers an automatic restart.</description>
      AUTOMATIC_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
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
