# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, unicode_literals, division

"""
Boundary Conditions Task Pane Window: This section allows a user to define the boundary
conditions for the described model. This section relies on regions named in the Regions section.
The top of the task pane is where users define/select BC regions


Icons to add/remove/duplicate boundary conditions are given at the top
Clicking the 'add' and 'duplicate' buttons triggers a popup window where the user must select a region to apply the boundary condition.
 Users cannot select inapplicable regions.
 BC regions must be planes or STLs (not volumes or points)
 No region can define more than one boundary condition.
Select boundary type
 Selection is required
 Available selections:
  Mass Inflow
  Plane regions set keyword BC_TYPE(#) to 'MI'
 STL regions set keyword BC_TYPE(#) to 'CG_MI'
 Pressure Outflow
 Plane regions set keyword BC_TYPE(#) to 'PO'
 STL regions set keyword BC_TYPE(#) to 'CG_PO'
 No Slip Wall
 Plane regions set keyword BC_TYPE(#) to 'NSW'
 STL regions set keyword BC_TYPE(#) to 'CG_NSW'
 Free Slip Wall
 Plane regions set keyword BC_TYPE(#) to 'FSW'
 STL regions set keyword BC_TYPE(#) to 'CG_FSW'
 Partial Slip Wall
 Plane regions set keyword BC_TYPE(#) to 'PSW'
 STL regions set keyword BC_TYPE(#) to 'CG_PSW'

Pressure Inflow
 Plane regions set keyword BC_TYPE(#) to 'PI'
 Not available for STL regions
Mass Outflow
 Plane regions set keyword BC_TYPE(#) to 'MO'
 STL regions set keyword BC_TYPE(#) to 'CG_MO'
 Specification always available
 DEFAULT - No slip wall
 Error check: mass fractions must sum to one

Tabs group boundary condition parameters for phases and additional equations. Tabs are
unavailable if no input is required from the user.
 Fluid tab - Unavailable if the fluid phase was disabled.
 Each solid phase will have its own tab. The tab name should be the name of the solid
 Group tab inputs by equation type (e.g., momentum, energy, species). Making the grouped inputs a 'collapsible list' may make navigation easier.
Subtask Pane Tab for Wall type (NSW, FSW, PSW, CG_NSW, CG_FSW, CG_PSW) Boundary
Condition Regions
Fluid (tab)
 Define transfer coefficient
 Specification only available with PSW
 Sets keyword BC_HW_G(#)
 DEFAULT value of 0.0
 Define Wall U-velocity
 Specification only available with PSW
 Sets keyword BC_UW_G(#)
 DEFAULT value of 0.0
 Define Wall V-velocity
 Specification only available with PSW
 Sets keyword BC_VW_G(#)
 DEFAULT value of 0.0
 Define Wall W-velocity
 Specification only available with PSW
 Sets keyword BC_WW_G(#)
 DEFAULT value of 0.0
 Select energy equation boundary type:
 Selection only available when solving energy equations
 Available selections:
 No-Flux (adiabatic) [DEFAULT]
 Sets keyword BC_HW_T_G(#) to 0.0
 Sets keyword BC_C_T_G(#) to 0.0
 Sets keyword BC_TW_G(#) to UNDEFINED
 Specified Temperature
 Sets keyword BC_HW_T_G(#) to UNDEFINED
 Sets keyword BC_C_T_G(#) to 0.0
 Requires BC_TW_G(#)
 Specified Flux
 Sets keyword BC_HW_T_G(#) to 0.0
 Requires BC_C_T_G(#)
 Sets keyword BC_TW_G(#) to UNDEFINED
 Convective Flux
 Requires BC_HW_T_G(#)
 Sets keyword BC_C_T_G(#) to 0.0
 Requires BC_TW_G(#)
Define wall temperature
 Specification only available with 'Specified Temperature' BC type
 Sets keyword BC_TW_G(#)
 DEFAULT value of 293.15
Define constant flux
 Specification only available with 'Specified Flux' BC type
 Sets keyword BC_C_T_G(#)
 DEFAULT value of 0.0
Define transfer coefficient
 Specification only available with 'Convective Flux' BC type
 Sets keyword BC_HW_T_G(#)
 DEFAULT value of 0.0
Define free stream temperature
 Specification only available with 'Convective Flux' BC type
 Sets keyword BC_TW_G(#)
 DEFAULT value of 0.0
Select species equation boundary type:
 Selection only available when solving species equations
 Available selections:
 No-Flux [DEFAULT]
 Sets keyword BC_HW_X_G(#) to 0.0
 Sets keyword BC_C_X_G(#) to 0.0
 Sets keyword BC_XW_G(#) to UNDEFINED
 Specified Mass Fraction
 Sets keyword BC_XW_T_G(#) to UNDEFINED
 Sets keyword BC_C_X_G(#) to 0.0
 Requires BC_XW_G(#)
 Specified Flux
 Sets keyword BC_HW_X_G(#) to 0.0
 Requires BC_C_X_G(#)
 Sets keyword BC_XW_G(#) to UNDEFINED
 Convective Flux
 Requires BC_HW_X_G(#)
 Sets keyword BC_C_X_G(#) to 0.0
 Requires BC_XW_G(#)
Define wall mass fraction
 Specification only available with 'Specified Mass Fraction' BC type
 Sets keyword BC_XW_G(#)
 DEFAULT value of 0.0
Define constant flux
 Specification only available with 'Specified Flux' BC type
 Sets keyword BC_C_X_G(#)
 DEFAULT value of 0.0
Define transfer coefficient
 Specification only available with 'Convective Flux' BC type
 Sets keyword BC_HW_X_G(#)
 DEFAULT value of 0.0
Define free stream mass fraction
 Specification only available with 'Convective Flux' BC type
 Sets keyword BC_XW_G(#)
 DEFAULT value of 0.0

Mockup of Task pane for specifying the Fluid properties for WALL boundary condition regions.
Solids-# (tab) - (Replace with phase name defined by the user)
 Enable Jackson-Johnson partial slip boundary
 Disabled (0.0) for CARTESIAN_GRID = .TRUE.
 Disabled (0.0) for KT_TYPE = 'ALGEBRAIC'
 Disabled (0.0) for KT_TYPE = 'GHD_2007'
 Sets keyword BC_JJ_PS(#)
 DEFAULT value of 1.0 when not disabled
 Select type of Jackson and Johnson BC:
 Selection only available BC_JJ_PS(#) = 1.0
 Available selections:
 Default Jackson-Johnson BC [DEFAULT]
 Sets keyword BC_JJ_M to .FALSE.
 Sets keyword JENKINS to .FALSE.
 Variable specularity coefficient
 Sets keyword BC_JJ_M to .TRUE.
 Sets keyword JENKINS to .FALSE.
Jenkins small frictional boundary
 Sets keyword BC_JJ_M to .FALSE.
 Sets keyword JENKINS to .TRUE.
Define restitution coefficient
 Specification only available with BC_JJ_PS(#) = 1.0
 Sets keyword E_W
 DEFAULT value of 1.0
 Required when available
Define specularity coefficient
 Specification only available with BC_JJ_PS(#)=1.0 and JENKINS=.FALSE.
 Sets keyword PHIP
 DEFAULT value of 0.6
 Required when available
Define specularity coefficient at zero slip
 Specification only available with BC_JJ_PS(#)=1.0 and BC_JJ_M=.TRUE.
 Sets keyword PHIP0
 DEFAULT -blank
 Optional when available
Define angle of internal friction
 Specification only available with BC_JJ_PS(#)=1.0 and (JENKINS=.TRUE.
FRICTION_MODEL=SRIVASTAVA)
 DEFAULT value of 11.31
 Required when available
Define transfer coefficient
 Specification only available with PSW
 Sets keyword BC_HW_S(#,#)
 DEFAULT value of 0.0
Define Wall U-velocity
 Specification only available with PSW or BC_JJ_PS(#) = 1.0
 Sets keyword BC_UW_S(#,#)
 DEFAULT value of 0.0
Define Wall V-velocity
 Specification only available with PSW or BC_JJ_PS(#) = 1.0
 Sets keyword BC_VW_S(#,#)
 DEFAULT value of 0.0
Define Wall W-velocity
 Specification only available with PSW or BC_JJ_PS(#) = 1.0
 Sets keyword BC_WW_S(#,#)
 DEFAULT value of 0.0
Select energy equation boundary type:
 Selection only available when solving energy equations
 Available selections:
 No-Flux (adiabatic) [DEFAULT]
  Sets keyword BC_HW_T_S(#,#) to 0.0
  Sets keyword BC_C_T_S(#,#) to 0.0
  Sets keyword BC_TW_S(#,#) to UNDEFINED
 Specified Temperature
  Sets keyword BC_HW_T_S(#,#) to UNDEFINED
  Sets keyword BC_C_T_S(#,#) to 0.0
  Requires BC_TW_S(#,#)
 Specified Flux
  Sets keyword BC_HW_T_S(#,#) to 0.0
  Requires BC_C_T_S(#)
  Sets keyword BC_TW_S(#,#) to UNDEFINED
 Convective Flux
  Requires BC_HW_T_S(#,#)
  Sets keyword BC_C_T_S(#,#) to 0.0
  Requires BC_TW_S(#,#)
Define wall temperature
  Specification only available with 'Specified Temperature' BC type
  Sets keyword BC_TW_S(#,#)
  DEFAULT value of 293.15
Define constant flux
  Specification only available with 'Specified Flux' BC type
  Sets keyword BC_C_T_S(#,#)
  DEFAULT value of 0.0
Define transfer coefficient
  Specification only available with 'Convective Flux' BC type
  Sets keyword BC_HW_T_S(#,#)
  DEFAULT value of 0.0
Define free stream temperature
  Specification only available with 'Convective Flux' BC type
  Sets keyword BC_TW_S(#,#)
  DEFAULT value of 0.0
Select granular energy equation boundary type:
  Selection only available with BC_JJ_PS(#)=0.0 and KT_TYPE /= 'ALGEBRAIC'
  Available selections:
  No-Flux [DEFAULT]
  Sets keyword BC_HW_THETA_M(#,#) to 0.0
  Sets keyword BC_C_THETA_M (#,#) to 0.0
  Sets keyword BC_THETAW_M(#,#) to UNDEFINED
  Specified Temperature
  Sets keyword BC_HW_THETA_M(#,#) to UNDEFINED
  Sets keyword BC_C_THETA_M(#,#) to 0.0
  Requires BC_THETAW_M(#,#)
  Specified Flux
  Sets keyword BC_HW_THETA_M(#,#) to 0.0
  Requires BC_C_THETA_M(#)
  Sets keyword BC_THETAW_M(#,#) to UNDEFINED
Define granular temperature
  Specification only available with 'Specified Temperature' BC type
  Sets keyword BC_THETAW_M(#,#)
  DEFAULT value of 0.0
Define constant flux
  Specification only available with 'Specified Flux' BC type
  Sets keyword BC_C_THETA_M(#,#)
  DEFAULT value of 0.0

When solving solids species equations:
  Set keyword BC_HW_X_S(#,#,#) to 0.0
  Set keyword BC_C_X_S(#,#,#) to 0.0
  Set keyword BC_XW_S(#,#,#) to UNDEFINED

Mockup of Task pane for specifying the Solid-# properties for WALL boundary condition regions.
Scalar (tab) - Tab only available if scalar equations are solved
  Select scalar boundary type:
  Available selections:
  No-Flux [DEFAULT]
  Sets keyword BC_HW_SCALAR(#,#) to 0.0
  Sets keyword BC_C_SCALAR(#,#) to 0.0
  Sets keyword BC_SCALARW(#,#) to UNDEFINED
  Specified Temperature
  Sets keyword BC_HW_T_S(#,#) to UNDEFINED
  Sets keyword BC_C_SCALAR (#,#) to 0.0
  Requires BC_SCALARW (#,#)
  Specified Flux
  Sets keyword BC_HW_T_S(#,#) to 0.0
  Requires BC_C_SCALAR (#)
  Sets keyword BC_SCALARW (#,#) to UNDEFINED
  Convective Flux
  Requires BC_HW_T_S(#,#)
  Sets keyword BC_C_SCALAR (#,#) to 0.0
  Requires BC_SCALARW (#,#)
  Define wall temperature
  Specification only available with 'Specified Temperature' BC type
  Sets keyword BC_SCALARW (#,#)
  DEFAULT value of 0.0
  Define constant flux
  Specification only available with 'Specified Flux' BC type
  Sets keyword BC_C_SCALAR (#,#)
  DEFAULT value of 0.0
  Define transfer coefficient
  Specification only available with 'Convective Flux' BC type
  Sets keyword BC_HW_SCALAR(#,#)
  DEFAULT value of 0.0
  Define free stream temperature
  Specification only available with 'Convective Flux' BC type
  Sets keyword BC_SCALARW (#,#)
  DEFAULT value of 0.0

Mockup of Task pane for specifying the Scalar properties for WALL boundary condition regions.
Note that this tab should only be available if scalar equations are being solved. Furthermore, the
number of scalars requiring input (here 2) comes from the number of scalar equations specified by
the user.
Subtask Pane Tab for INFLOW type (MI, PI, CG_MI) Boundary Condition Regions
Fluid (tab)
  Define volume fraction
  Specification always available
  Sets keyword BC_EP_G(#)
  DEFAULT value of 1.0 for MI and CG_MI; leave [UNDEFINED] for PI
  Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0
  Error Check: For PI - either all are defined and sum to 1, or all are undefined
  Define inflow properties: Mass inflow specification changes based on the BC_TYPE and
Region orientation (e.g., XZ-Plane)
  For BC_TYPE='MI' and XZ-Plane region
  Select mass inflow specification type:
  Available selections:
  Y-Axial Velocity (m/s) [DEFAULT]
  Sets keyword BC_V_G(#)
  DEFAULT value of 0.0
  Volumetric Flowrate (m3/s)
  Sets keyword BC_VOLFLOW_G(#)
  DEFAULT value of 0.0
  Mass Flowrate (kg/s)
  Sets keyword BC_MASSFLOW_G(#)
  DEFAULT value of 0.0
  Define Tangential Velocities:
  Define X-Axial Velocity
  Sets keyword BC_U_G(#)
  DEFAULT value of 0.0
  Define Z-Axial Velocity
  Sets keyword BC_W_G(#)
  DEFAULT value of 0.0
  For BC_TYPE='MI' and YZ-Plane region
  Select mass inflow specification type:
  Available selections:
  X-Axial Velocity (m/s) [DEFAULT]
  Sets keyword BC_U_G(#)
  DEFAULT value of 0.0
  Volumetric Flowrate (m3/s)
  Sets keyword BC_VOLFLOW_G(#)
  DEFAULT value of 0.0
  Mass Flowrate (kg/s)
  Sets keyword BC_MASSFLOW_G(#)
  DEFAULT value of 0.0
  Define Tangential Velocities:
  Define Y-Axial Velocity
  Sets keyword BC_V_G(#)
  DEFAULT value of 0.0
  Define Z-Axial Velocity
  Sets keyword BC_W_G(#)
  DEFAULT value of 0.0
  For BC_TYPE='MI' and XY-Plane region
  Select mass inflow specification type:
  Available selections:
  Z-Axial Velocity (m/s) [DEFAULT]
  Sets keyword BC_W_G(#)
  DEFAULT value of 0.0
  Volumetric Flowrate (m3/s)
  Sets keyword BC_VOLFLOW_G(#)
  DEFAULT value of 0.0
  Mass Flowrate (kg/s)
  Sets keyword BC_MASSFLOW_G(#)
  DEFAULT value of 0.0
  Define Tangential Velocities:
  Define X-Axial Velocity
  Sets keyword BC_U_G(#)
  DEFAULT value of 0.0
  Define Y-Axial Velocity
  Sets keyword BC_V_G(#)
  DEFAULT value of 0.0
  For BC_TYPE='CG_MI' or 'PI'
  Specify all velocity components:
  Define X-Axial Velocity
  Sets keyword BC_U_G(#)
  DEFAULT value of 0.0
  Define Y-Axial Velocity
  Sets keyword BC_V_G(#)
  DEFAULT value of 0.0
  Define Z-Axial Velocity
  Sets keyword BC_V_G(#)
  DEFAULT value of 0.0
Define temperature
  Specification always available
  Input required for any of the following
  Fluid density model: Ideal Gas Law
  Fluid viscosity model: Sutherland's Law
  Energy equations are solved
  Sets keyword BC_T_G(#)
  DEFAULT value of 293.15
Define pressure
  Specification always available
  Input required when combining ideal gas law and specified mass inflow rate
  Input required for BC_TYPE = PI
  Sets keyword BC_P_G(#)
  DEFAULT 101.325d3
Select species and set mass fractions (table format)
  Specification always available
  Input required for species equations
  Drop down menu of fluid species
  Sets keyword BC_X_G(#,#)
  DEFAULT - last defined species has mass fraction of 1.0
  Error check: mass fractions must sum to one
Turbulence: Define k-ε turbulent kinetic energy
  Specification only available with K-Epsilon turbulence model
  Sets keyword BC_K_TURB_G(#)
  DEFAULT value of 0.0
Turbulence: Define k-ε turbulent dissipation
  Specification only available with K-Epsilon turbulence model
  Sets keywords BC_E_TURB_G(#)
  DEFAULT value of 0.0

Mockup of Task pane for specifying the Fluid properties for MI, boundary condition regions.
Solid-# (tab) - Rename tab to user provided solids name.
Define volume fraction
  Specification always available
  Sets keyword BC_EP_S(#,#)
  DEFAULT value of 1.0 - (sum of previous tabs) for MI and CG_MI; leave
 UNDEFINED] for PI
  Error Check: For MI and CG_MI, BC_EP_G(#) + BC_EP_S(#,:) = 1.0
  Error Check: For PI - either all are defined and sum to 1, or all are undefined
Define inflow properties: Mass inflow specification changes based on the BC_TYPE and
Region orientation (e.g., XZ-Plane)
  For BC_TYPE='MI' and XZ-Plane region
  Select mass inflow specification type:
  Available selections:
  Y-Axial Velocity (m/s) [DEFAULT]
  Sets keyword BC_V_S(#,#)
  DEFAULT value of 0.0
  Volumetric Flowrate (m3/s)
  Sets keyword BC_VOLFLOW_S(#,#)
  DEFAULT value of 0.0
  Mass Flowrate (kg/s)
  Sets keyword BC_MASSFLOW_S(#,#)
  DEFAULT value of 0.0
  Define Tangential Velocities:
  Define X-Axial Velocity
  Sets keyword BC_U_S(#,#)
  DEFAULT value of 0.0
  Define Z-Axial Velocity
  Sets keyword BC_W_G(#,#)
  DEFAULT value of 0.0
  For BC_TYPE='MI' and YZ-Plane region
  Select mass inflow specification type:
  Available selections:
  X-Axial Velocity (m/s) [DEFAULT]
  Sets keyword BC_U_S(#,#)
  DEFAULT value of 0.0
  Volumetric Flowrate (m3/s)
  Sets keyword BC_VOLFLOW_S(#,#)
  DEFAULT value of 0.0
  Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_S#, (#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define Y-Axial Velocity
o Sets keyword BC_V_S(#,#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_W_S(#,#)
o DEFAULT value of 0.0
o For BC_TYPE='MI' and XY-Plane region
 Select mass inflow specification type:
 Available selections:
 Z-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_W_S(#,#)

Software Requirements Specification for MFIX-

o DEFAULT value of 0.0
Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_S(#,#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_S(#,#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define X-Axial Velocity
o Sets keyword BC_U_S(#,#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_S(#,#)
o DEFAULT value of 0.0
o For BC_TYPE='CG_MI' or 'PI'
 Specify all velocity components:
 Define X-Axial Velocity
o Sets keyword BC_U_S(#,#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_S(#,#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_V_S(#,#)
o DEFAULT value of 0.0
Define temperature
o Specification always available
o Input required when energy equations are solved
o Sets keyword BC_T_S(#,#)
o DEFAULT value of 293.15
Select species and set mass fractions (table format)
o Specification always available
o Input required for species equations
o Drop down menu of fluid species
o Sets keyword BC_X_S(#,#,#)
o DEFAULT - last defined species has mass fraction of 1.0
o Error check: mass fractions must sum to one






UI

Page 41

Software Requirements Specification for MFIX-

UI

Page 42

Mockup of Task pane for specifying the Solids properties for MI, CG_MI, and PI boundary
condition regions.
Scalar (tab) - Tab only available if scalar equations are solved
 Define initial scalar value
o Sets keyword BC_SCALAR(#,#)
o DEFAULT value of 0.0
Subtask Pane Tab for PRESSURE OUTFLOW type (PO) Boundary Condition Regions
Fluid (tab)
 Define pressure
o Specification always available
o Input required
o Sets keyword BC_P_G(#)
o DEFAULT 101.325d3
The remaining inputs are "optional." They do not have default values, because MFIX will calculate
appropriate values if they are unspecified and 'backflow' occurs at the outlet.
 Define volume fraction
o Specification always available
o Sets keyword BC_EP_G(#)
o No DEFAULT value
o Error Check: If any volume fraction for the BC region is specified, then all volume
factions for the BC region must be specified and must sum to one.
 Define temperature
o Specification always available
o NO DEFAULT value
o Sets keyword BC_T_G(#)
 Select species and set mass fractions (table format)
o Specification always available
o NO DEFAULT value
o Sets keyword BC_X_G(#,#)

Software Requirements Specification for MFIX-

UI

Page 43

o Error check: if specified, mass fractions must sum to one_
Solids-# (tab)
All inputs are optional. They do not have default values, because MFIX will calculate appropriate
values if they are unspecified and 'backflow' occurs at the outlet.
 Define volume fraction
o Specification always available
o Sets keyword BC_EP_S(#,#)
o No DEFAULT value
o Error Check: If any volume fraction for the BC region is specified, then all volume
factions for the BC region must be specified and must sum to one.
 Define temperature
o Specification always available
o NO DEFAULT value
o Sets keyword BC_T_S(#,#)
 Select species and set mass fractions (table format)
o Specification always available
o NO DEFAULT value
o Sets keyword BC_X_S(#,#,#)
o Error check: if specified, mass fractions must sum to one_
Scalar (tab) - Tab only available if scalar equations are solved
All inputs are optional. They do not have default values, because MFIX will calculate appropriate
values if they are unspecified and 'backflow' occurs at the outlet.


Define scalar value
o Sets keyword BC_SCALAR(#,#)
o NO DEFAULT value

Software Requirements Specification for MFIX-

UI

Page 44

Subtask Pane Tab for MASS OUTFLOW type (MO) Boundary Condition Regions
Fluid (tab)
 Define inflow properties: Mass inflow specification changes based on the BC_TYPE and
Region orientation (e.g., XZ-Plane)
o For BC_TYPE='MO' and XZ-Plane region
 Select mass inflow specification type:
 Available selections:
 Y-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_V_G(#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_G(#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_G(#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define X-Axial Velocity
o Sets keyword BC_U_G(#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_W_G(#)
o DEFAULT value of 0.0
o For BC_TYPE='MO' and YZ-Plane region
 Select mass inflow specification type:
 Available selections:
 X-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_U_G(#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_G(#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_G(#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define Y-Axial Velocity
o Sets keyword BC_V_G(#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_W_G(#)
o DEFAULT value of 0.0
o For BC_TYPE='MO' and XY-Plane region
 Select mass inflow specification type:
 Available selections:
 Z-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_W_G(#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_G(#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_G(#)

Software Requirements Specification for MFIX-

UI

Page 45

o DEFAULT value of 0.0
Define Tangential Velocities:
 Define X-Axial Velocity
o Sets keyword BC_U_G(#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_G(#)
o DEFAULT value of 0.0
o For BC_TYPE='CG_MO'
 Specify all velocity components:
 Define X-Axial Velocity
o Sets keyword BC_U_G(#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_G(#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_V_G(#)
o DEFAULT value of 0.0




Define duration to average outflow rate.
o Specification always available
o Input required
o Sets keyword BC_DT_0(#)
o DEFAULT value of 0.1
o Error Check: Value should be positive (non-negative)
o BC_DT_0 specification should persist across the gas and solids tabs. If the user sets
it in the gas phase tab, but then changes it under a solids tab, a warning messing
indicating that this value is 'constant' across all phases should be given.

The remaining inputs are only required when either the mass or the volumetric flowrates are
specified. They are not required if the velocities are given for the outlet.






Define volume fraction
o Specification only available with mass or volumetric flowrates.
o Input required
o Sets keyword BC_EP_G(#)
o DEFAULT value 1.0
o Error Check: If any volume fraction for the BC region is specified, then all volume
factions for the BC region must be specified and must sum to one.
Define temperature
o Specification only available with mass or volumetric flowrates and R_G0 is
UNDEFINED
o DEFAULT value 293.15
o Sets keyword BC_T_G(#)
Select species and set mass fractions (table format)
o Specification only available with mass or volumetric flowrates and R_G0 is
UNDEFINED
o DEFAULT value 1.0 of last defined species
o Sets keyword BC_X_G(#,#)
o Error check: if specified, mass fractions must sum to one_

Solids-# (tab)

Software Requirements Specification for MFIX-



UI

Page 46

Define inflow properties: Mass inflow specification changes based on the BC_TYPE and
Region orientation (e.g., XZ-Plane)
o For BC_TYPE='MO' and XZ-Plane region
 Select mass inflow specification type:
 Available selections:
 Y-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_V_S(#,#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_ S(#,#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_ S(#,#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define X-Axial Velocity
o Sets keyword BC_U_ S(#,#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_W_ S(#,#)
o DEFAULT value of 0.0
o For BC_TYPE='MO' and YZ-Plane region
 Select mass inflow specification type:
 Available selections:
 X-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_U_ S(#,#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_ S(#,#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_ S(#,#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define Y-Axial Velocity
o Sets keyword BC_V_ S(#,#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_W_ S(#,#)
o DEFAULT value of 0.0
o For BC_TYPE='MO' and XY-Plane region
 Select mass inflow specification type:
 Available selections:
 Z-Axial Velocity (m/s) [DEFAULT]
o Sets keyword BC_W_ S(#,#)
o DEFAULT value of 0.0
 Volumetric Flowrate (m3/s)
o Sets keyword BC_VOLFLOW_ S(#,#)
o DEFAULT value of 0.0
 Mass Flowrate (kg/s)
o Sets keyword BC_MASSFLOW_ S(#,#)
o DEFAULT value of 0.0
 Define Tangential Velocities:
 Define X-Axial Velocity

Software Requirements Specification for MFIX-

UI

Page 47

o Sets keyword BC_U_ S(#,#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_ S(#,#)
o DEFAULT value of 0.0
o For BC_TYPE='CG_MO'
 Specify all velocity components:
 Define X-Axial Velocity
o Sets keyword BC_U_ S(#,#)
o DEFAULT value of 0.0
 Define Y-Axial Velocity
o Sets keyword BC_V_ S(#,#)
o DEFAULT value of 0.0
 Define Z-Axial Velocity
o Sets keyword BC_V_ S(#,#)
o DEFAULT value of 0.0


Define duration to average outflow rate.
o Specification always available
o Input required
o Sets keyword BC_DT_0(#)
o DEFAULT value of 0.1
o Error Check: Value should be positive (non-negative)
BC_DT_0 specification should persist across the gas and solids tabs. If the user sets it in the gas
phase tab, but then changes it under a solids tab, a warning messing indicating that this value is
'constant' across all phases should be given.
"""

class BCS(object):
    def init_bcs(self):
        self.bcs_current_region = None

    def setup_bcs(self):
        pass
