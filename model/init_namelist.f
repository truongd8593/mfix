!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: INIT_NAMELIST                                           C
!  Purpose: initialize the NAMELIST variables                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 26-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Initialize Phi and Phi_w                                   C
!  Author: M. Syamlal                                 Date: 11-FEB-93  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Add L_scale0, L_scale                                      C
!  Author: W. Sams                                    Date: 04-MAY-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To call DES_Init_Namelist                                  C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: To call CARTESIAN_GRID_INIT_NAMELIST                       C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!  Keyword Documentation Format:                                       C
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

!<keyword category="run control" required="true">
!  <description>name used to create output files. the name should be legal after extensions are added to it; e.g., for run name bub01, the output files bub01.log, bub01.out, bub01.res, etc., will be created.</description>
      RUN_NAME = UNDEFINED_C
!</keyword>

!<keyword category="run control" required="false">
!  <description>problem description in 60 characters.</description>
      DESCRIPTION = UNDEFINED_C
!</keyword>

!<keyword category="run control" required="true">
!  <description>units for data input and output.</description>
!  <valid value="cgs" note="all input and output in cgs units (g, cm, s, cal)."/>
!  <valid value="si" note="all input and output in si units (kg, m, s, j)."/>
      UNITS = UNDEFINED_C
!</keyword>

!<keyword category="run control" required="true">
!  <description>type of run.</description>
!  <valid value="new" note="new run."/>
!  <valid value="restart_2" note="start a new run with initial conditions from a .res file created from another run."/>
!  <valid value="restart_1" note="normal restart run. initial conditions from .res file."/>
!  <valid value="restart_3" note="continue old run as in restart_1, but any input data not given in mfix.dat is read from the .res file. (do not use)"/>
!  <valid value="restart_4" note="start a new run as in restart_2, but any input data not given in mfix.dat is read from the .res file. (do not use)"/>
      RUN_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="run control" required="false">
!  <description>start-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TIME = UNDEFINED
!</keyword>

!<keyword category="run control" required="false">
!  <description>stop-time of the run.</description>
!  <range min="0.0" max="+Inf" />
      TSTOP = UNDEFINED
!</keyword>

!<keyword category="run control" required="false">
!  <description>starting time step. if dt is not defined, a steady-state calculation will be performed.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT = UNDEFINED
!</keyword>

!<keyword category="run control" required="false">
!  <description>maximum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MAX = ONE
!</keyword>

!<keyword category="run control" required="false">
!  <description>minimum time step.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="+Inf" />
      DT_MIN = 1.D-6
!</keyword>

!<keyword category="category name" required="false">
!  <description>factor for adjusting time step. should be less than 1.</description>
!  <dependent keyword="TIME" value="DEFINED"/>
!  <dependent keyword="TSTOP" value="DEFINED"/>
!  <range min="0.0" max="1" />
      DT_FAC = 0.9D0
!</keyword>

!<keyword category="category name" required="false">
!  <description></description>
!  <valid value=".false." note="do not reduce time step for stalled iterations."/>
!  <valid value=".true." note="reduce time step if the residuals sum does not decrease."/>
      DETECT_STALL = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>enables clean termination feature.</description>
      CHK_BATCHQ_END = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>total wall-clock duration of the job, in seconds.</description>
      BATCH_WALLCLOCK = 9000.0    ! set to 2.5 hrs for jaguarcnl w/ nproc<=512
!</keyword>

!<keyword category="category name" required="false">
!  <description>buffer time when initiating clean termination, in seconds.</description>
      TERM_BUFFER = 180.0         ! set to 3 minutes prior to end of job
!</keyword>

!<keyword category="category name" required="false">
!  <description>Solve energy equations.</description>
!  <valid value=".false." note="Do not solve energy equations."/>
!  <valid value=".true." note="Solve energy equations."/>
      ENERGY_EQ = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Use deferred correction method for implementing higher order discretization.</description>
!  <valid value=".false." note="use down-wind factor method (default)."/>
!  <valid value=".true." note="use deferred correction method for implementing higher order discretization."/>
      DEF_COR  =  .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>factor used in the universal limiter (when fpfoi is set .true.) and can be any value in the set (0,1). the choice of 1 will give (diffusion) first order upwinding and as this value becomes closer to 0 the scheme becomes more compressive.</description>
!  <range min="0.0" max="1.0" />
!  <dependent keyword="fpfoi" value=".true."/>
      C_FAC = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>four point fourth order interpolation and is upstream biased. if this scheme is chosen and discretize(*) < 2, discretize(*) is defaulted to 2. if you chose this scheme, set the c_fac value between 0 and 1.</description>
!  <dependent keyword="c_fac" value="DEFINED"/>
      FPFOI = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>use the granular energy transport equation (pde) as opposed to the algebraic (alg) equation formulation.</description>
      GRANULAR_ENERGY = .FALSE.
!</keyword>

! do not use revised JJ BC
!<keyword category="Boundary Condition" required="false">
!  <description>Use the modified Johnson and Jackson partial slip BC with variable specularity coefficient. Must set e_w and phi_w with this BC.</description>
!  <dependent keyword="e_w" value="DEFINED"/>
!  <dependent keyword="phi_w" value="DEFINED"/>
      BC_JJ_M = .false.
!</keyword>

!<keyword category="category name" required="false">
!  <description>output the variable specularity coefficient when bc_jj_m is .true.. the specularity coefficient will be stored in reactionrates array for post-processing by post-mfix. user needs to set nrr to 1 for this purpose. be careful with this setting when reacting flow is simulated.</description>
      PHIP_OUT_JJ=.false.
!</keyword>

      PHIP_OUT_ITER=0

! SUBGRID model, Sebastien Dartevelel, LANL, May 2013

!<keyword category="category name" required="false">
!  <description>subgrid models include: igci and milioli</description>
!  <valid value="igci" note=""/>
!  <valid value="milioli" note=""/>
      SUBGRID_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>Include wall effect term in subgrid model.</description>
      SUBGRID_Wall = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Ratio of filter size to computational cell size.</description>
      filter_size_ratio = 2.0D0
!</keyword>

!<keyword category="category name" required="false">
!  <description>when activated the k-epsilon turbulence model (for single-phase flow) is solved using standard wall functions.</description>
      K_Epsilon = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>when activated the added (or virtual) mass force effectively acts to increase the inertia of the dispersed phase, which tends to stabilize simulations of bubbly gas-liquid flows.</description>
      Added_Mass = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>disperse phase number where the added mass applies.</description>
      M_AM = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>use simonin model (see ~mfix/doc/simonin_ahmadi_models.pdf for details).</description>
      SIMONIN = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>use ahmadi model (see ~mfix/doc/simonin_ahmadi_models.pdf for details).</description>
      AHMADI = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>use jenkins small frictional boundary condition (see ~mfix/doc/simonin_ahmadi_models.pdf for details).</description>
      JENKINS = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>use yu and standish correlation to compute maximum packing for polydisperse systems.</description>
      YU_STANDISH = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>use fedors and landel correlation to compute maximum packing for a binary (only) mixture of powders.</description>
      FEDORS_LANDEL = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      AUTOMATIC_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      AUTO_RESTART = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      ITER_RESTART = 1
!</keyword>

! peter: 7/15
!<keyword category="category name" required="false">
!  <description>specifies the mean _y_ velocity component at the eastern boundary of the domain (v_sh), and the mean _y_ velocity (-v_sh) at the western boundary of the domain.</description>
      V_sh=0d0
!</keyword>

! jeg: 4/01/2005
!<keyword category="category name" required="false">
!  <description>solids phase stress model. To specify a kinetic theory granular_energy must be .true.  If kt_type is left undefined (default) the viscous model is based on the theory of lun et al. (1984).</description>
!  <dependent keyword="GRANULAR_ENERGY" value=".true."/>
!  <valid value="ia_nonep" note="iddir & arastoopour (aiche j, 2005)"/>
!  <valid value="gd_99" note="garzo and dufty (pre, 1999)"/>
!  <valid value="ghd" note="garzo, hrenya and dufty (pre, 2007)"/>
!  <valid value="ia_nonep" note="iddir & arastoopour (aiche j, 2005)"/>
      KT_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>radial distribution function at contact for polydisperse systems. </description>
!  <valid value="lebowitz" note=""/>
!  <valid value="modified_lebowitz" note=""/>
!  <valid value="mansoori" note=""/>
!  <valid value="modified_mansoori" note=""/>
      RDF_TYPE = 'LEBOWITZ'
!</keyword>

! anuj: 04/20
!<keyword category="category name" required="false">
!  <description>use the schaeffer model when .false., or use the princeton model when .true.</description>
!  <valid value=".false." note="use the schaeffer model"/>
!  <valid value=".true." note="use the princeton model"/>
      FRICTION = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>for a term appearing in the frictional stress model invoked with friction = .true.</description>
!  <dependent keyword="friction" value=".true."/>
!  <valid value="0" note="use s:s in the frictional stress model."/>
!  <valid value="2" note="an appropriate combination of the above two forms."/>
!  <valid value="1" note="use an alternate form suggested by savage."/>
      SAVAGE = 1
!</keyword>

! sof: 02/16/2005
!<keyword category="category name" required="false">
!  <description>if set to false with friction = .false., then the model will not have any frictional viscosity.</description>
!  <valid value=".true." note=""/>
!  <valid value=".false." note=""/>
      SCHAEFFER = .TRUE.
!</keyword>

! sp: 02/08/2006
!<keyword category="category name" required="false">
!  <description>This will turn on the blending function to blend the Schaeffer stresses with that of kinetic theory around epsilon*. The default is hyperbolic tangent function for blending (TANH_BLEND = .TRUE.) and one could also utilize a scaled and truncated sigmoidal function (SIGM_BLEND=.TRUE.).</description>
      BLENDING_STRESS = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>hyperbolic tangent function for blending</description>
!  <dependent keyword="BLENDING_STRESS" value=".true."/>
!  <conflict keyword="SIGM_BLEND" value=".true."/>
      TANH_BLEND      = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>a scaled and truncated sigmoidal function for blending</description>
!  <dependent keyword="BLENDING_STRESS" value=".true."/>
!  <conflict keyword="TANH_BLEND" value=".true."/>
      SIGM_BLEND      = .FALSE.
!</keyword>

! sp: 06/15/2007
!<keyword category="category name" required="false">
!  <description>if set to false, the residuals are grouped into far fewer global collectives and this should only have impact on parallel runs</description>
      DEBUG_RESID     = .TRUE.
!</keyword>

! DISTIO distributed IO
!<keyword category="category name" required="false">
!  <description>NA</description>
      bDist_IO            = .false.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      bStart_with_one_RES = .false.
!</keyword>

! loezos:
!<keyword category="category name" required="false">
!  <description>if .true. imposes a mean shear on the flow field as a linear function of _x_ coordinate. this feature should only be used when cyclic_x=.true. also, the keyword v-sh needs to be set.</description>
!  <dependents>cyclic_x v_sh</dependents>
      SHEAR = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>drag model (see drag_gs.f for details).</description>
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

!<keyword category="category name" required="false">
!  <description>Quantity for calibrating Syamlal-O'Brien drag correlation using Umf data.  This are determined using the Umf spreadsheet.</description>
      drag_c1 = 0.8d0
!</keyword>

!<keyword category="category name" required="false">
!  <description>Quantity for calibrating Syamlal-O'Brien drag correlation using Umf data.  This are determined using the Umf spreadsheet.</description>
      drag_d1 = 2.65d0
!</keyword>

!<keyword category="category name" required="false">
!  <description>if use_def_lam_hys is set to .false. the user is able to specify a value for the lubrication cutoff distance (lam_hys).  in practice this number should be on the order of the mean free path of the gas for smooth particles, or the rms roughness of a particle if they are rough (if particle roughness is larger than the mean free path).</description>
!  <dependents>USE_DEF_LAM_HYS</dependents>
      LAM_HYS = UNDEFINED
!</keyword>

! AE: 041601 Set the default to 1st order accurate time implementation
!<keyword category="category name" required="false">
!  <description>temporal discretization scheme.</description>
!  <valid value=".false." note="Implicit Euler based temporal discretization scheme employed (first order accurate in time)."/>
!  <valid value=".true." note="crank-nicholson based temporal discretization scheme employed (second order accurate in time excluding the restart timestep which is first order)."/>
      CN_ON = .FALSE.
!</keyword>

      IF (DIM_M + 1 > 0) THEN
!<keyword category="category name" required="false">
!  <description>Solve X-momentum equations of phase m.</description>
!  <valid value=".true." note="Solve X-momentum equations of phase m."/>
!  <valid value=".false." note="Do not solve X-momentum equations of phase m. Beware of inconsistencies when the momentum equations are turned off; e.g., 2-D developing flow with only Y-momentum should not specify no-slip-walls."/>
         MOMENTUM_X_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Solve Y-momentum equations of phase m.</description>
!  <valid value=".true." note="Solve Y-momentum equations of phase m."/>
!  <valid value=".false." note="Do not solve Y-momentum equations of phase m."/>
         MOMENTUM_Y_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Solve Z-momentum equations of phase m.</description>
!  <valid value=".true." note="Solve Z-momentum equations of phase m."/>
!  <valid value=".false." note="Do not solve Z-momentum equations of phase m."/>
         MOMENTUM_Z_EQ(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Solve species equations (m=0 indicates gas phase).</description>
!  <valid value=".false." note="do not solve species equations of phase m."/>
!  <valid value=".true." note="solve species equations of phase m.  to solve species equation with no chemical reactions, copy the file mfix/model/rrates.f into run directory and remove the first two executable lines as explained in the comments.  no other change is needed in that file."/>
         SPECIES_EQ(:DIM_M) = .TRUE.
!</keyword>

      ENDIF

!<keyword category="category name" required="false">
!  <description>call user-defined subroutines.</description>
!  <valid value=".true." note="call user-defined subroutines."/>
!  <valid value=".false." note="do not call user-defined subroutines."/>
      CALL_USR = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>momentum equations.</description>
!  <valid value=".false." note="model A"/>
!  <valid value=".true." note="model B"/>
      MODEL_B = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>discretization scheme for seven types of equations.</description>
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

!<keyword category="category name" required="false">
!  <description>chi-scheme, proposed by darwish and moukalled (2003), is activated. this scheme guarantees that the set of differenced species mass balance equations has the property that the sum of mass fractions add up to one. when a flux limiter is used with (higher order) spatial discretization schemes it is not guaranteed that the mass fractions add up to one. this problem may be rectified by activating the chi-scheme.</description>
      Chi_scheme = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of solid phases to solve the population balance equations.</description>
      NScalar = 0
!</keyword>

      Phase4Scalar(:) = UNDEFINED_I

!<keyword category="category name" required="false">
!  <description>the number of user defined chemical reactions stored in the *.spa file.  see section 4.10 chemical reactions.</description>
      nRR = 0
!</keyword>

!<keyword category="category name" required="false">
!  <description>variable to decide if the population balance equations are solved.</description>
      Call_DQMOM = .FALSE.
!</keyword>


! INITIALIZE THE OUTPUT CONTROL SECTION
!<keyword category="category name" required="false">
!  <description>NA</description>
      report_mass_balance_dt = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>interval in number of time steps at which .log file is written.</description>
!  <valid value="report_mass_balance_dt " note="if a value is defined, say 0.1 s, an overall species mass balance is performed and reported in the log file. the over all mass balance calculations may slightly slow down the run."/>
!  <valid value="resid_string" note="specify residuals to be printed as 4-character strings. "/>
      NLOG = 25
!</keyword>

!<keyword category="category name" required="false">
!  <description>if true, display the residuals on the screen and messages about convergence on the screen and in the .log file.</description>
      FULL_LOG = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>interval at which restart (.res) file is updated.</description>
      RES_DT = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>interval at which .spx files are updated.</description>
!  <valid value=".sp4" note="solids velocity (u_s, v_s, w_s)."/>
!  <valid value=".sp7" note="gas and solids mass fractions (x_g, x-s)."/>
!  <valid value=".sp6" note="gas and solids temperature (t_g, t_s1, t_s2)."/>
!  <valid value=".sp1" note="void fraction (ep_g)."/>
!  <valid value=".sp3" note="gas velocity (u_g, v_g, w_g)."/>
!  <valid value=".sp2" note="gas pressure, solids pressure (p_g, p_star)."/>
!  <valid value=".sp9" note="user defined scalars."/>
!  <valid value=".sp8" note="granular temperature (g)."/>
!  <valid value=".sp5" note="solids density (rop_s)."/>
!  <valid value=".spa" note="reaction rates. (see section 4.11)"/>
!  <valid value=".spb" note=""/>
      SPX_DT(:N_SPX) = UNDEFINED
!</keyword>
      LC = N_SPX + 1

!<keyword category="category name" required="false">
!  <description>interval at which standard output (.out) file is updated.</description>
      OUT_DT = UNDEFINED
!</keyword>

      DO LC = 1, DIMENSION_USR
!<keyword category="category name" required="false">
!  <description>interval at which user-defined outputs are written from the subroutine write_usr1.</description>
         USR_DT(LC) = UNDEFINED
!</keyword>

         USR_TYPE(LC) = UNDEFINED_C
         USR_VAR(LC) = UNDEFINED_C
         USR_FORMAT(LC) = UNDEFINED_C
         USR_EXT(LC) = UNDEFINED_C
      END DO
      DO LC = 1, 8
         RESID_STRING(LC) = UNDEFINED_C
      END DO

!<keyword category="category name" required="false">
!  <description>coordinates used in the simulation.</description>
!  <valid value="cartesian" note="cartesian coordinates."/>
!  <valid value="cylindrical" note="cylindrical coordinates."/>
      COORDINATES = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>(do not use.)</description>
!  <valid value=".false." note="x (r) direction is considered."/>
!  <valid value=".true." note="x (r) direction is not considered."/>
      NO_I = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>(do not use.)</description>
!  <valid value=".false" note="y direction is considered."/>
!  <valid value=".true." note="y direction is not considered."/>
      NO_J = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description></description>
!  <valid value=".false." note="z(theta) direction is considered."/>
!  <valid value=".true." note="z(theta) direction is not considered."/>
      NO_K = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of cells in the x (r) direction.</description>
      IMAX = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of cells in the y direction.</description>
      JMAX = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of cells in the z (() direction.</description>
      KMAX = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of solids phases.</description>
      MMAX = 1
!</keyword>

!<keyword category="category name" required="false">
!  <description>the inner radius in the simulation of an annular cylindrical region.</description>
      XMIN = ZERO
!</keyword>

!<keyword category="category name" required="false">
!  <description>reactor length in the x (r) direction.</description>
      XLENGTH = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>reactor length in the y direction.</description>
      YLENGTH = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>reactor length in the z (theta) direction.</description>
      ZLENGTH = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Cell sizes in the x (r) direction. Enter values from DX(0) to DX(IMAX-1). (Use uniform mesh size with higher-order discretization methods.  Also in cylindrical coordinates DX should be kept uniform for strict momentum conservation.)</description>
      DX(:DIM_I) = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Cell sizes in the y direction. Enter values from DY(0) to DY(IMAX-1). (Use uniform mesh size with second-order discretization methods.)</description>
      DY(:DIM_J) = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Cell sizes in the z (theta) direction. Enter values from DZ(0) to DZ(IMAX-1). (Use uniform mesh size with second-order discretization methods.)</description>
      DZ(:DIM_K) = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the x-direction cyclic without pressure drop. no other boundary conditions for the x-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at x-boundary."/>
!  <valid value=".true." note="cyclic condition at x-boundary."/>
      CYCLIC_X = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the y-direction cyclic without pressure drop. no other boundary conditions for the y-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at y-boundary."/>
!  <valid value=".true." note="cyclic condition at x-boundary."/>
      CYCLIC_Y = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the z-direction cyclic without pressure drop. no other boundary conditions for the z-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at z-boundary."/>
!  <valid value=".true." note="cyclic condition at z-boundary."/>
      CYCLIC_Z = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the x-direction cyclic with pressure drop. if the keyword flux_g is given a value this becomes a cyclic boundary condition with specified mass flux. no other boundary conditions for the x-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at x-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at x-boundary."/>
      CYCLIC_X_PD = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the y-direction cyclic with pressure drop. if the keyword flux_g is given a value this becomes a cyclic boundary condition with specified mass flux. no other boundary conditions for the y-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at y-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at y-boundary."/>
      CYCLIC_Y_PD = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>flag for making the z-direction cyclic with pressure drop. if the keyword flux_g is given a value this becomes a cyclic boundary condition with specified mass flux. no other boundary conditions for the z-direction should be specified.</description>
!  <valid value=".false." note="no cyclic condition at z-boundary."/>
!  <valid value=".true." note="cyclic condition with pressure drop at z-boundary."/>
      CYCLIC_Z_PD = .FALSE.
!</keyword>


! Constants
      DO LC = 1, DIMENSION_C
         C(LC) = UNDEFINED
         C_NAME(LC) = '....................'
      ENDDO
!<keyword category="category name" required="false">
!  <description>gravitational acceleration. by default, the gravity force acts in the _ve y-direction. modify file b_force2.inc to change the body force term.</description>
      GRAVITY = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>coefficient of restitution for particle-particle collisions. (mfix 1.94 keyword _e_).</description>
      C_E = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>coefficient of friction between the particles of two solids phases.</description>
      C_F = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>angle of internal friction (in degrees). set this value to zero to turn off plastic regime stress calculations.</description>
      PHI = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specularity coefficient associated with particle-wall collisions.</description>
      PHIP = 0.6D0
!</keyword>

      k4phi = undefined

!<keyword category="category name" required="false">
!  <description>specify the value of specularity coefficient when the normalized slip velocity goes to zero when bc_jj_m is .true.. this variable is calculated internally in the code. do not modify unless an accurate number is known.</description>
!  <dependents>bc_jj_m</dependents>
      phip0 = undefined
!</keyword>

!<keyword category="category name" required="false">
!  <description>coefficient of restitution for particle-wall collisions.</description>
      E_W = 1.D0
!</keyword>

!<keyword category="category name" required="false">
!  <description>angle of internal friction (in degrees) at walls. set this value to non-zero (phi_w = 11.31 means tan_phi_w = mu = 0.2) when using jenkins or bc_jj_m boundary condition.</description>
      PHI_W = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>minimum solids fraction above which friction sets in.  (when friction = .true.)</description>
!  <dependents>friction</dependents>
      EPS_F_MIN = 0.5D0
!</keyword>

!<keyword category="category name" required="false">
!  <description>value of turbulent length initialized. this may be overwritten in specific regions with the keyword ic_l_scale.</description>
      L_SCALE0 = ZERO
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum value of the turbulent viscosity of the fluid.</description>
      MU_GMAX = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>excluded volume in boyle-massoudi stress.</description>
!  <valid value="0.0" note="b-m stress is turned off."/>
      V_EX = ZERO
!</keyword>

!<keyword category="category name" required="false">
!  <description>reference pressure.</description>
      P_REF = ZERO
!</keyword>

!<keyword category="category name" required="false">
!  <description>scale factor for pressure.</description>
      P_SCALE = ONE
!</keyword>

! gera: 08/15/03
!<keyword category="category name" required="false">
!  <description>used in calculating the initial slope of segregation: see gera et al. (2004) - recommended value 0.3. increasing this coefficient results in decrease in segregation of particles in binary mixtures.</description>
      SEGREGATION_SLOPE_COEFFICIENT=0.D0
!</keyword>

!<keyword category="category name" required="false">
!  <description>Maximum solids volume fraction at packing for polydisperse systems (more than one solids phase used). The value of EP_star may change during the computation if solids phases with different particle diameters are specified and Yu_Standish or Fedors_Landel correlations are used.</description>
!  <range min="0" max="1-ep_star" />
      EP_S_MAX(:DIM_M) = UNDEFINED
!</keyword>

! coefficient of restitution
      r_p(:DIM_M, :DIM_M) = UNDEFINED

! rong
!<keyword category="category name" required="false">
!  <description>success-factor for aggregation.</description>
      AGGREGATION_EFF=0.D0
!</keyword>

!<keyword category="category name" required="false">
!  <description>success-factor for breakage.</description>
      BREAKAGE_EFF=0.D0
!</keyword>

! numerics
!<keyword category="category name" required="false">
!  <description>factor to normalize the gas continuity equation residual.</description>
      NORM_G = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>factor to normalize the solids continuity equation residual.</description>
      NORM_S = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum residual at convergence (continuity+momentum).</description>
      TOL_RESID = 1.0D-3
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum residual at convergence (energy).</description>
      TOL_RESID_T = 1.0D-4
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum residual at convergence (species balance).</description>
      TOL_RESID_X = 1.0D-4
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum residual at convergence (scalar balances.)</description>
      TOL_RESID_Scalar = 1.0D-4
!</keyword>

      TOL_RESID_K_Epsilon = 1.0D-4

!<keyword category="category name" required="false">
!  <description>maximum residual at convergence (granular energy).</description>
      TOL_RESID_Th = 1.0D-4
!</keyword>

!<keyword category="category name" required="false">
!  <description>minimum residual for declaring divergence. when the fluid is incompressible, the velocity residuals take large values in the second iteration (e.g., 1e+8) and then drop down to a low value in the third iteration (e.g., 0.1). in such cases, it is desirable to increase this setting.</description>
      TOL_DIVERGE = 1.0D+4
!</keyword>

!<keyword category="category name" required="false">
!  <description>the code declares divergence if the velocity anywhere in the domain exceeds a maximum value.  this maximum value is automatically determined from the boundary values. the user may scale the maximum value by adjusting this scale factor.</description>
      MAX_INLET_VEL_FAC = ONE
!</keyword>

!<keyword category="category name" required="false">
!  <description>maximum number of iterations.</description>
      MAX_NIT = 500
!</keyword>

      LEQ_IT(1) = 20
      LEQ_IT(2) = 20
      LEQ_IT(3) = 5
      LEQ_IT(4) = 5
      LEQ_IT(5) = 5
      LEQ_IT(6) = 15
      LEQ_IT(7) = 15
      LEQ_IT(8) = 15
      LEQ_IT(9) = 15
      LEQ_METHOD(1) = 2
      LEQ_METHOD(2) = 2
      LEQ_METHOD(3) = 2
      LEQ_METHOD(4) = 2
      LEQ_METHOD(5) = 2
      LEQ_METHOD(6) = 2
      LEQ_METHOD(7) = 2
      LEQ_METHOD(8) = 2
      LEQ_METHOD(9) = 2
      LEQ_SWEEP(1) = 'RSRS'
      LEQ_SWEEP(2) = 'RSRS'
      LEQ_SWEEP(3) = 'RSRS'
      LEQ_SWEEP(4) = 'RSRS'
      LEQ_SWEEP(5) = 'RSRS'
      LEQ_SWEEP(6) = 'RSRS'
      LEQ_SWEEP(7) = 'RSRS'
      LEQ_SWEEP(8) = 'RSRS'
      LEQ_SWEEP(9) = 'RSRS'
      LEQ_TOL(1) = 1.0D-4
      LEQ_TOL(2) = 1.0D-4
      LEQ_TOL(3) = 1.0D-4
      LEQ_TOL(4) = 1.0D-4
      LEQ_TOL(5) = 1.0D-4
      LEQ_TOL(6) = 1.0D-4
      LEQ_TOL(7) = 1.0D-4
      LEQ_TOL(8) = 1.0D-4
      LEQ_TOL(9) = 1.0D-4
      LEQ_PC(1:9)  = 'LINE'

!<keyword category="category name" required="false">
!  <description>NA</description>
      DO_TRANSPOSE = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      icheck_bicgs = 1
!</keyword>

!<keyword category="category name" required="false">
!  <description>Print out additional statistics for parallel runs</description>
      solver_statistics = .FALSE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      opt_parallel = .FALSE.
!</keyword>

! AEOLUS: set default value for new variable to debug print whole
!         index layout
!<keyword category="category name" required="false">
!  <description>NA</description>
      DBGPRN_LAYOUT = .FALSE.
!</keyword>

! AEOLUS: set default value for enabling all processors write out
!         their *.LOG invidually
!<keyword category="category name" required="false">
!  <description>NA</description>
      ENABLE_DMP_LOG = .FALSE.
!</keyword>

      UR_FAC(1)  = 0.8D0             !pressure
      UR_FAC(2)  = 0.5D0             !rho, ep
      UR_FAC(3)  = 0.5D0             !U
      UR_FAC(4)  = 0.5D0             !V
      UR_FAC(5)  = 0.5D0             !W
      UR_FAC(6)  = 1.0D0             !T
      UR_FAC(7)  = 1.0D0             !X
      UR_FAC(8)  = 0.5D0             !Th
      UR_FAC(9)  = 0.8D0             !Scalar

!<keyword category="category name" required="false">
!  <description>the implicitness calculation of the gas-solids drag coefficient may be underrelaxed by changing ur_f_gs, which takes values between 0 to 1:</description>
!  <range min="0" max="1" />
      UR_F_gs    = 1.0D0             !drag coefficient update
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      UR_Kth_sml = 1.0D0            ! conductivity term in IA theory
!</keyword>

! INITIALIZE THE GAS PHASE SECTION
!<keyword category="category name" required="false">
!  <description>specified constant gas density. this value may be set to zero to make the drag zero and to simulate granular flow in a vacuum. for this case, users may turn off solving for gas momentum equations to accelerate convergence.</description>
      RO_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant gas viscosity.</description>
      MU_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant gas conductivity.</description>
      K_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant gas diffusivity.</description>
      DIF_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant gas specific heat.</description>
      C_PG0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>average molecular weight of gas.</description>
      MW_AVG = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Number of species in phase m. Note that the gas phase is indicated as m=0.</description>
      NMAX(0) = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>molecular weight of gas species n.</description>
      MW_G(:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant granular viscosity. if this value is specified, then the kinetic theory calculation is turned off and p_s = 0 and lambda_s = -2/3 mu_s0.</description>
      MU_S0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant solids conductivity.</description>
      K_S0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant solids diffusivity.</description>
      DIF_S0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>specified constant solids specific heat.</description>
      C_PS0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Initial particle diameters, same as the old D_P(m).</description>
      D_P0(:DIM_M) = UNDEFINED
!</keyword>

!--------------------------------------------------------------------------> JMusser.0 Start
      RO_S0 = UNDEFINED
      X_s0 = UNDEFINED
      RO_Xs0 = UNDEFINED
      INERT_SPECIES = UNDEFINED_I
!--------------------------------------------------------------------------> JMusser.0 End

!<keyword category="category name" required="false">
!  <description>Number of species in phase m. Note that the gas phase is indicated as m=0.</description>
      NMAX(1:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>indicates whether the solids phase forms a packed bed with a void fraction ep_star.</description>
      CLOSE_PACKED(:DIM_M) = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>Molecular weight of solids phase-m, species n.</description>
      MW_S = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>packed bed void fraction.</description>
      EP_STAR = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>Number of species comprising the gas phase.</description>
      NMAX_g = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>Name of gas phase species n as it appears in the materials database.</description>
      SPECIES_g(:) = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>User defined name for gas phase species n.</description>
      SPECIES_ALIAS_g(:) = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>Number of species comprising solids phase m.</description>
      NMAX_s(:) = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>Name of solids phase m, species n as it appears in the materials database.</description>
      SPECIES_s(:,:) = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>User defined name for solids phase m, species n</description>
      SPECIES_ALIAS_s(:,:) = UNDEFINED_C
!</keyword>

!<keyword category="category name" required="false">
!  <description>Use the automated reaction rate UDFS; usr_rates.f and usr_rates_des.f.</description>
      USE_RRATES = .FALSE.
!</keyword>

! NO_OF_RXNS is not a keyword. However, it is initialized here so that
! if there are no reactions, this value is assigned.
      NO_OF_RXNS = UNDEFINED_I

! INITIALIZE THE INITIAL CONDITIONS
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
!  <description>Type of initial condition. Mainly used in restart runs to overwrite values read from the .RES file by specifying it as _PATCH_. The user needs to be careful when using the _PATCH_ option, since the values from the .RES file are overwritten and no error checking is done for the patched values.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial void fraction in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas pressure in the IC region. If this quantity is not specified, MFIX will set up a hydrostatic pressure profile, which varies only in the y-direction.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_P_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids pressure in the IC region. Usually, this value is specified as zero.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_P_STAR(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Turbulence length scale in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_L_SCALE(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial gas phase temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Gas phase radiation coefficient in the IC region. Modify file radtn2.inc to change the source term.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_GAMA_RG(LC) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Gas phase radiation temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_T_RG(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial mass fraction of gas species n.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="species" min="1" max="DIM_N_g"/>
         IC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of gas velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_W_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial bulk density (rop_s = ro_s x ep_s) of solids phase-m in the IC region. Users need to specify this IC only for polydisperse flow (MMAX > 1). Users must make sure that summation of ( IC_ROP_s(ic,m) / RO_s(m) ) over all solids phases is equal to ( 1.0 – IC_EP_g(ic) ).</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial x-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial y-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial z-component of solids-phase velocity in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial solids phase-m granular temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Initial value of Scalar n.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_SCALAR(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

! sof: force users to set initial values for K and Epsilon.
         IC_K_Turb_G(LC) = UNDEFINED
         IC_E_Turb_G(LC) = UNDEFINED

!<keyword category="Initial Condition" required="false">
!  <description>Solids phase-m radiation coefficient in the IC region. Modify file radtn2.inc to change the source term.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_GAMA_RS(LC,:DIM_M) = ZERO
!</keyword>

!<keyword category="Initial Condition" required="false">
!  <description>Solids phase-m radiation temperature in the IC region.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
         IC_T_RS(LC,:DIM_M) = UNDEFINED
!</keyword>


! safe is set not inflate the inital lattice as the algorithm has not 
! been thoroughly tested 
          IC_DES_FIT_TO_REGION(LC) = .FALSE. 

!         IC_X_S(LC,1,1+:DIM_N_S+) = UNDEFINED
      ENDDO

!<keyword category="Initial Condition">
!  <description>Solids phase species mass fraction.</description>
!  <arg index="1" id="IC region" min="1" max="DIMENSION_IC"/>
!  <arg index="2" id="phase" min="1" max="DIM_M"/>
!  <arg index="3" id="species" min="1" max="DIM_N_s"/>
!  <range min="0.0" max="1.0" />
      IC_X_S = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>fluid pressure drop across xlength when a cyclic boundary condition with pressure drop is imposed in the x-direction.</description>
      DELP_X = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>fluid pressure drop across ylength when a cyclic boundary condition with pressure drop is imposed in the y-direction.</description>
      DELP_Y = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>fluid pressure drop across zlength when a cyclic boundary condition with pressure drop is imposed in the z-direction.</description>
      DELP_Z = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>if a value is specified (in units of g/cm2.s), the domain-averaged gas flux is held constant at that value in simulations over a periodic domain.  a pair of boundaries specified as periodic with fixed pressure drop is then treated as periodic with fixed mass flux.   even for this case a pressure drop must also be specified, which is used as the initial guess in the simulations.</description>
      Flux_g = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      U_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      V_G0 = UNDEFINED
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      W_G0 = UNDEFINED
!</keyword>

      U_S0(:DIM_M) = UNDEFINED
      V_S0(:DIM_M) = UNDEFINED
      W_S0(:DIM_M) = UNDEFINED

! Boundary Conditions
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
!  <description>Void fraction at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_EP_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas pressure at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_P_G(LC) = UNDEFINED
!</keyword>

         BC_ROP_G(LC) = UNDEFINED

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase temperature at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of gas species n at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The species diffusion boundary condition is implemented as dX/dn + hw (X - Xw) = c, where n is the normal pointing from the fluid into the wall. The coefficients hw, Xw, and c should be specified. Hw = 0 => specified species diffusion flux; hw = ∞ => specified species concentration at the boundary. To set hw = ∞ , leave it unspecified and give a value for Xw. Gas phase hw for mass transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Xw for mass transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_XW_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase C for mass transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_X_G(LC,:DIM_N_G) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>x-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_U_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_V_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z-component of gas velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_W_G(LC) = UNDEFINED
!</keyword>

         BC_VELMAG_G(LC) = UNDEFINED

!<keyword category="Boundary Condition" required="false">
!  <description>Type of boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <valid value='DUMMY' description='The specified boundary condition is ignored. This is useful for turning off some boundary conditions without having to delete them from the file.' />
!  <valid value='MASS_INFLOW' description='Mass inflow rates for gas and solids phases are specified at the boundary.' alias='MI' />
!  <valid value='MASS_OUTFLOW' description='The specified values of gas and solids mass outflow rates at the boundary are maintained, approximately. This condition should be used sparingly for minor outflows, when the bulk of the outflow is occurring through other constant pressure outflow boundaries.' alias='MO' />
!  <valid value='P_INFLOW' description='Inflow from a boundary at a specified constant pressure. To specify as the west, south, or bottom end of the computational region, add a layer of wall cells to the west, south, or bottom of the PI cells. Users need to specify all scalar quantities and velocity components. The specified values of fluid and solids velocities are only used initially as MFIX computes these values at this inlet boundary.' alias='PI' />
!  <valid value='P_OUTFLOW' description='Outflow to a boundary at a specified constant pressure. To specify as the west, south, or bottom end of the computational region, add a layer of wall cells to the west, south, or bottom of the PO cells.' alias='PO' />
!  <valid value='FREE_SLIP_WALL' description='Velocity gradients at the wall vanish. If BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is used for solids.  A FSW is equivalent to using a PSW with hw=0.' alias='FSW' />
!  <valid value='NO_SLIP_WALL' description='All components of the velocity vanish at the wall. If BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is used for solids.    A NSW is equivalent to using a PSW with vw=0 and hw undefined.' alias='NSW' />
!  <valid value='PAR_SLIP_WALL' description='Partial slip at the wall implemented as dv/dn + hw (v – vw) = 0, where n is the normal pointing from the fluid into the wall. The coefficients hw and vw should be specified. For free slip set hw = 0. For no slip leave hw undefined (hw=∞) and set vw = 0. To set hw = ∞, leave it unspecified. If BC_JJ_PS is equal to 1, the Johnson-Jackson boundary condition is used for solids.' alias='PSW' />
         BC_TYPE(LC) = UNDEFINED_C
!</keyword>

         BC_APPLY_TO_MPPIC(LC) = .true.

!<keyword category="Boundary Condition" required="false">
!  <description>Gas volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas mass flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval at the beginning when the normal velocity at the boundary is equal to BC_Jet_g0. When restarting, run this value and BC_Jet_g0 should be specified such that the transient jet continues correctly. MFIX does not store the jet conditions. For MASS_OUTFLOW boundary conditions, BC_DT_0 is the time period to average and print the outflow rates. The adjustment of velocities to get a specified mass or volumetric flow rate is based on the average outflow rate.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_DT_0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>The interval when normal velocity is equal to BC_Jet_gh.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_DT_H(LC) = UNDEFINED
!</keyword>

         BC_DT_L(LC) = UNDEFINED
!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the initial interval BC_DT_0.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_JET_G0(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Value of normal velocity during the interval BC_DT_h.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_JET_GH(LC) = UNDEFINED
!</keyword>

         BC_JET_GL(LC) = UNDEFINED

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_UW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_WW_G(LC) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase hw for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_T_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase Tw for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_TW_G(LC) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Gas phase C for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_T_G(LC) = UNDEFINED
!</keyword>


! Kapil & Anuj: 01/19/98
!<keyword category="Boundary Condition" required="false">
!  <description>Johnson and Jackson partial slip bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
!  <valid value='0' description='Do not use Johnson and Jackson partial slip bc. If granular energy transport equation is not solved: (GRANULAR_ENERGY=.FALSE.). />
!  <valid value='1' description='Use Johnson and Jackson partial slip bc. If granular energy transport equation is solved: (GRANULAR_ENERGY=.TRUE.). />
         BC_JJ_PS(LC) = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Bulk density of solids phase at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_ROP_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>x-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_U_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>y-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_V_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>z-component of solids-phase velocity at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_W_S(LC,:DIM_M) = UNDEFINED
!</keyword>

         BC_VELMAG_S(LC,:DIM_M) = UNDEFINED

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m temperature at the BC plane.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids volumetric flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VOLFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids mass flow rate through the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_MASSFLOW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Uw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_UW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Vw for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_VW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Ww for partial slip boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_WW_S(LC,:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Tw for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_TW_S(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase C for heat transfer.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_T_S(LC,:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>Hw for granular energy bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Tw for granular energy bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_THETAW_M(LC,:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>c for granular energy bc.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_THETA_M(LC,:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Boundary Condition" required="false">
!  <description>hw for scalar transfer at the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_HW_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Xw for scalar transfer at the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_ScalarW(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>C for scalar transfer at the boundary.</description>
!  <arg index="1" id="BC region" min="1" max="DIMENSION_BC"/>
         BC_C_Scalar(LC,:DIM_SCALAR) = UNDEFINED
!</keyword>

      ENDDO

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase-m granular temperature at the BC plane.</description>
      BC_THETA_M = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>NA</description>
      BC_Scalar = UNDEFINED
!</keyword>

! sof: force users to set inlet BC for K and Epsilon
!<keyword category="Boundary Condition" required="false">
!  <description>NA</description>
      BC_K_Turb_G = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>NA</description>
      BC_E_Turb_G = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Mass fraction of solids phase-m, species n at the BC plane.</description>
      BC_X_S = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase hw for mass transfer.</description>
      BC_HW_X_S = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase Xw for mass transfer.</description>
      BC_XW_S = UNDEFINED
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>Solids phase C for mass transfer.</description>
      BC_C_X_S = UNDEFINED
!</keyword>

! Point Source:
!<keyword category="Point Source" required="false">
!  <description>x coordinate of the west face or edge.</description>
      PS_X_W = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>x coordinate of the east face or edge.</description>
      PS_X_E = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y coordinate of the south face or edge.</description>
      PS_Y_S = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y coordinate of the north face or edge.</description>
      PS_Y_N = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z coordinate of the bottom face or edge.</description>
      PS_Z_B = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z coordinate of the top face or edge.</description>
      PS_Z_T = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>i index of the west-most cell.</description>
      PS_I_W = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>i index of the east-most cell.</description>
      PS_I_E = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>j index of the south-most cell.</description>
      PS_J_S = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>j index of the north-most cell.</description>
      PS_J_N = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>k index of the bottom-most cell.</description>
      PS_K_B = UNDEFINED_I
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>k index of the top-most cell.</description>
      PS_K_T = UNDEFINED_I
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>x-component of incoming gas velocity.</description>
      PS_U_G = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y-component of incoming gas velocity.</description>
      PS_V_G = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z-component of incoming gas velocity.</description>
      PS_W_G = UNDEFINED
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>Gas mass flow rate through the point source.</description>
      PS_MASSFLOW_G = UNDEFINED
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming gas.</description>
      PS_T_G = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Gas phase incoming species n mass fraction.</description>
      PS_X_G = UNDEFINED
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>x-component of incoming solids velocity.</description>
      PS_U_S = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>y-component of incoming solids velocity.</description>
      PS_V_S = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>z-component of incoming solids velocity.</description>
      PS_W_S = UNDEFINED
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>Solids mass flow rate through the point source.</description>
      PS_MASSFLOW_S = UNDEFINED
!</keyword>


!<keyword category="Point Source" required="false">
!  <description>Temperature of incoming solids.</description>
      PS_T_S = UNDEFINED
!</keyword>

!<keyword category="Point Source" required="false">
!  <description>Solids phase incoming species n mass fraction.</description>
      PS_X_S = UNDEFINED
!</keyword>



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
!  <description>permeability</description>
         IS_PC(LC,1) = LARGE_NUMBER
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Inertial resistance coefficient.</description>
         IS_PC(LC,2) = ZERO
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Type of internal surface</description>
!  <valid value="IMPERMEABLE" note="No gas or solids flow through the surface." alias="IP"/>
!  <valid value="SEMIPERMEABLE" note="Gas flows through the surface with an additional resistance. Solids velocity through the surface is set to zero or to a user-specified fixed value (i.e., solids momentum equation for this direction is not solved)." alias="SP"/>
         IS_TYPE(LC) = UNDEFINED_C
!</keyword>

!<keyword category="Internal Surface" required="false">
!  <description>Value of fixed solids velocity through semipermeable surfaces.</description>
         IS_VEL_S(LC,:DIM_M) = ZERO
!</keyword>
      ENDDO

      DO LC = 1, DIM_N_ALL
!<keyword category="category name" required="false">
!  <description>Names of gas and solids phase species as it appears in the materials database. The first NMAX(0) are the names of gas species. The next NMAX(1) are the names of solids phase-1 species, etc.</description>
         SPECIES_NAME(LC) = UNDEFINED_C
!</keyword>
      ENDDO

!<keyword category="category name" required="false">
!  <description>number of grid blocks in x-direction.</description>
      NODESI = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of grid blocks in y-direction.</description>
      NODESJ = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>number of grid blocks in z-direction.</description>
      NODESK = UNDEFINED_I
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      IS_SERIAL = .TRUE.
!</keyword>

!<keyword category="category name" required="false">
!  <description>NA</description>
      USE_DOLOOP = .FALSE.
!</keyword>

! Stiff Chemistry Solver:
      STIFF_CHEMISTRY = .FALSE.
      CALL_DI   = .FALSE.      ! Legacy Keyword
      CALL_GROW = .FALSE.      ! Legacy Keyword
      CALL_ISAT = .FALSE.      ! Legacy Keyword
      ISATdt    = UNDEFINED    ! Legacy Keyword

      bWrite_netCDF(:) = .false.


! Debug flags
      REPORT_NEG_DENSITY = .FALSE.

      CALL DES_INIT_NAMELIST

! Initialize QMOMK namelist
      CALL QMOMK_INIT_NAMELIST

      CALL USR_INIT_NAMELIST

! JFD: MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
      CALL CARTESIAN_GRID_INIT_NAMELIST


      RETURN
      END SUBROUTINE INIT_NAMELIST
