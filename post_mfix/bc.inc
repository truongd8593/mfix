CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: bc.inc                                                 C
C  Purpose: Common block containing boundary conditions data           C
C                                                                      C
C  Author: M. Syamlal                                 Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
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
C                      x coordinate of the west face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_X_w (DIMENSION_BC) 
C
C                      x coordinate of the east face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_X_e (DIMENSION_BC)
C
C                      y coordinate of the south face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_Y_s (DIMENSION_BC)
C
C                      y coordinate of the north face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_Y_n (DIMENSION_BC)
C
C                      z coordinate of the bottom face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_Z_b (DIMENSION_BC)
C
C                      z coordinate of the top face of a region where
C                      boundary conditions are specified
      DOUBLE PRECISION BC_Z_t (DIMENSION_BC)
C
C                      i index of the west face of a region where
C                      boundary conditions are specified
      INTEGER          BC_I_w (DIMENSION_BC) 
C
C                      i index of the east face of a region where
C                      boundary conditions are specified
      INTEGER          BC_I_e (DIMENSION_BC)
C
C                      j index of the south face of a region where
C                      boundary conditions are specified
      INTEGER          BC_J_s (DIMENSION_BC)
C
C                      j index of the north face of a region where
C                      boundary conditions are specified
      INTEGER          BC_J_n (DIMENSION_BC)
C
C                      k index of the bottom face of a region where
C                      boundary conditions are specified
      INTEGER          BC_K_b (DIMENSION_BC)
C
C                      k index of the top face of a region where
C                      boundary conditions are specified
      INTEGER          BC_K_t (DIMENSION_BC)
C
C                      Void fraction in a specified boundary region
      DOUBLE PRECISION BC_EP_g (DIMENSION_BC)
C
C                      Gas pressure in a specified boundary region
      DOUBLE PRECISION BC_P_g (DIMENSION_BC)
C
C                      Microscopic density of gas in a specified
C                      boundary region
      DOUBLE PRECISION BC_RO_g (DIMENSION_BC)
C
C                      Macroscopic density of gas in a specified
C                      boundary region
      DOUBLE PRECISION BC_ROP_g (DIMENSION_BC)
C
C                      Macroscopic density of solids phases in a
C                      specified boundary region
      DOUBLE PRECISION BC_ROP_s (DIMENSION_BC, DIMENSION_M)
C
C                      Gas phase temperature in a specified boundary
C                      region
      DOUBLE PRECISION BC_T_g (DIMENSION_BC)
C
C                      Solids phase temperature in a specified 
C                      boundary region
      DOUBLE PRECISION BC_T_s (DIMENSION_BC, DIMENSION_M)
C
C                      Solids phase granular temperature in a specified 
C                      boundary region
      DOUBLE PRECISION BC_Theta_m (DIMENSION_BC, DIMENSION_M)
C
C                      x-component of gas velocity in a specified
C                      boundary region
      DOUBLE PRECISION BC_U_g (DIMENSION_BC)
C
C                      x-component of solids phase velocity in a
C                      specified boundary region
      DOUBLE PRECISION BC_U_s (DIMENSION_BC, DIMENSION_M)
C
C                      y-component of gas velocity in a specified
C                      boundary region
      DOUBLE PRECISION BC_V_g (DIMENSION_BC)
C
C                      y-component of solids phase velocity in a
C                      specified boundary region
      DOUBLE PRECISION BC_V_s (DIMENSION_BC, DIMENSION_M)
C
C                      z-component of gas velocity in a specified
C                      boundary region
      DOUBLE PRECISION BC_W_g (DIMENSION_BC)
C
C                      z-component of solids phase velocity in a
C                      specified boundary region
      DOUBLE PRECISION BC_W_s (DIMENSION_BC, DIMENSION_M)
C
C                      Type of boundary: MASS_INFLOW, MASS_OUTFLOW, 
C                      P_INFLOW, P_OUTFLOW, FREE_SLIP_WALL, NO_SLIP_WALL
      CHARACTER*16     BC_TYPE (DIMENSION_BC)
C
C                      Gas volumetric flow rate through the boundary
      DOUBLE PRECISION BC_VOLFLOW_g (DIMENSION_BC)
C
C                      Solids volumetric flow rate through the boundary
      DOUBLE PRECISION BC_VOLFLOW_s (DIMENSION_BC, DIMENSION_M)
C
C                      Gas mass flow rate through the boundary
      DOUBLE PRECISION BC_MASSFLOW_g (DIMENSION_BC)
C
C                      Solids mass flow rate through the boundary
      DOUBLE PRECISION BC_MASSFLOW_s (DIMENSION_BC, DIMENSION_M)
C
C                      Logical variable to determine whether a bc is defined
       LOGICAL         BC_DEFINED (DIMENSION_BC)
Cstart kapil&anuj 01/19/98
C
C                      Logical variable to determine whether the partial
C                      slip bc of Johnson and Jackson is used  
       INTEGER         BC_JJ_PS (DIMENSION_BC) 
Cend   kapil&anuj 01/19/98
C
C                      Character variable with values W, E, S, N, B, and T
C                      to determine the flow plane of a flow cell
       CHARACTER       BC_PLANE (DIMENSION_BC)
C
C                      The interval at the beginning when normal vel. is equal to
C                      BC_Jet_g0
      DOUBLE PRECISION BC_DT_0 (DIMENSION_BC)
C
C                      Stored value of normal velocity
      DOUBLE PRECISION BC_Jet_g (DIMENSION_BC)
C
C                      Value of normal vel. during the initial interval BC_DT_0
      DOUBLE PRECISION BC_Jet_g0 (DIMENSION_BC)
C
C                      The interval when normal vel. is equal to BC_Jet_gh
      DOUBLE PRECISION BC_DT_h (DIMENSION_BC)
C
C                      Value of normal vel. during the initial interval BC_DT_h
      DOUBLE PRECISION BC_Jet_gh (DIMENSION_BC)
C
C                      The interval when normal vel. is equal to BC_Jet_gl
      DOUBLE PRECISION BC_DT_l (DIMENSION_BC)
C
C                      Value of normal vel. during the initial interval BC_DT_l
      DOUBLE PRECISION BC_Jet_gl (DIMENSION_BC)
C
C                      Time to update a transient boundary condition
      DOUBLE PRECISION BC_TIME (DIMENSION_BC)
C
C                      Area of boundary surfaces
      DOUBLE PRECISION BC_AREA (DIMENSION_BC)
C
C                      Gas species mass fractions in a boundary region
      DOUBLE PRECISION BC_X_g (DIMENSION_BC, DIMENSION_N_g)
C
C                      Solids species mass fractions in a boundary region
      DOUBLE PRECISION BC_X_s (DIMENSION_BC, DIMENSION_M, DIMENSION_N_s)
C
C                      Accumulated or average mass outflow rate of gas
      DOUBLE PRECISION BC_MOUT_g(DIMENSION_BC)
C
C                      Accumulated or average mass outflow rate of solids
      DOUBLE PRECISION BC_MOUT_s(DIMENSION_BC, DIMENSION_M)
C
C                      Accumulated or average volumetric outflow rate of gas
      DOUBLE PRECISION BC_VOUT_g(DIMENSION_BC)
C
C                      Accumulated or average volumetric outflow rate of solids
      DOUBLE PRECISION BC_VOUT_s(DIMENSION_BC, DIMENSION_M)
C
C                      Number of outflow rate values accumulated
      INTEGER          BC_OUT_N (DIMENSION_BC)
C
C                      Pressure drop specified for cyclic b.c. in X
      DOUBLE PRECISION DELP_X
C
C                      Pressure drop specified for cyclic b.c. in Y
      DOUBLE PRECISION DELP_Y
C
C                      Pressure drop specified for cyclic b.c. in Z
      DOUBLE PRECISION DELP_Z
C
C                      Average gas velocity in X direction (for cyclic bc)
      DOUBLE PRECISION U_g0
C
C                      Average gas velocity in Y direction (for cyclic bc)
      DOUBLE PRECISION V_g0
C
C                      Average gas velocity in Z direction (for cyclic bc)
      DOUBLE PRECISION W_g0
C
C                      Average solids velocity in X direction (for cyclic bc)
      DOUBLE PRECISION U_s0 (DIMENSION_M)
C
C                      Average solids velocity in Y direction (for cyclic bc)
      DOUBLE PRECISION V_s0 (DIMENSION_M)
C
C                      Average solids velocity in Z direction (for cyclic bc)
      DOUBLE PRECISION W_s0 (DIMENSION_M)
C
C                      IJK location where P_g is fixed for cyclic b.c's
      INTEGER          IJK_P_g
C
C                      Coefficient in partial slip condition -- gas
      DOUBLE PRECISION BC_hw_g (DIMENSION_BC)
C
C                      Coefficient in partial slip condition -- solids
      DOUBLE PRECISION BC_hw_s (DIMENSION_BC, DIMENSION_M)
C
C                      Wall velocity for partial slip condition -- gas
      DOUBLE PRECISION BC_Uw_g (DIMENSION_BC)
C
C                      Wall velocity for partial slip condition -- gas
      DOUBLE PRECISION BC_Vw_g (DIMENSION_BC)
C
C                      Wall velocity for partial slip condition -- gas
      DOUBLE PRECISION BC_Ww_g (DIMENSION_BC)
C
C                      Wall velocity for partial slip condition -- solids
      DOUBLE PRECISION BC_Uw_s (DIMENSION_BC, DIMENSION_M)
C
C                      Wall velocity for partial slip condition -- solids
      DOUBLE PRECISION BC_Vw_s (DIMENSION_BC, DIMENSION_M)
C
C                      Wall velocity for partial slip condition -- solids
      DOUBLE PRECISION BC_Ww_s (DIMENSION_BC, DIMENSION_M)

C
C                      Coefficient in heat transfer boundary condition -- gas
      DOUBLE PRECISION BC_hw_T_g (DIMENSION_BC)
C
C                      Coefficient in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_hw_T_s (DIMENSION_BC, DIMENSION_M)
C
C                      Wall temperature in heat transfer boundary  condition -- gas
      DOUBLE PRECISION BC_Tw_g (DIMENSION_BC)
C
C                      Wall temperature in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_Tw_s (DIMENSION_BC, DIMENSION_M)
C
C                      Coefficient in heat transfer boundary condition -- gas
      DOUBLE PRECISION BC_C_T_g (DIMENSION_BC)
C
C                      Coefficient in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_C_T_s (DIMENSION_BC, DIMENSION_M)


C
C                      Coefficient in granular temp boundary condition -- solids
      DOUBLE PRECISION BC_hw_Theta_m (DIMENSION_BC, DIMENSION_M)
C
C                      Wall granular temperature -- solids
      DOUBLE PRECISION BC_Thetaw_m (DIMENSION_BC, DIMENSION_M)
C
C                      Coefficient in granualr tempearure boundary condition -- solids
      DOUBLE PRECISION BC_C_Theta_m (DIMENSION_BC, DIMENSION_M)

C
C                      Coefficient in heat transfer boundary condition -- gas
      DOUBLE PRECISION BC_hw_X_g (DIMENSION_BC, DIMENSION_N_g)
C
C                      Coefficient in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_hw_X_s (DIMENSION_BC, DIMENSION_M,
     &                            DIMENSION_N_s)
C
C                      Wall temperature in heat transfer boundary  condition -- gas
      DOUBLE PRECISION BC_Xw_g (DIMENSION_BC, DIMENSION_N_g)
C
C                      Wall temperature in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_Xw_s (DIMENSION_BC, DIMENSION_M,
     &                          DIMENSION_N_s)
C
C                      Coefficient in heat transfer boundary condition -- gas
      DOUBLE PRECISION BC_C_X_g (DIMENSION_BC, DIMENSION_N_g)
C
C                      Coefficient in heat transfer boundary condition -- solids
      DOUBLE PRECISION BC_C_X_s (DIMENSION_BC, DIMENSION_M, 
     &                           DIMENSION_N_s)
C
      COMMON / BC_DP /
     &       BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t,
     &       BC_EP_g, BC_P_g, BC_RO_g, BC_ROP_g, BC_ROP_s, 
     &       BC_T_g, BC_T_s,
     &       BC_U_g, BC_U_s, BC_V_g, BC_V_s, BC_W_g, BC_W_s,
     &       BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g, 
     &       BC_MASSFLOW_s, BC_Jet_g,
     &       BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_l, 
     &       BC_Jet_gl,
     &       BC_TIME, BC_AREA, BC_X_g, BC_X_s, 
     &       BC_MOUT_g, BC_MOUT_s, BC_VOUT_g, BC_VOUT_s, BC_OUT_N,
     &       DELP_X, DELP_Y, DELP_Z, 
     &       U_g0, U_s0, V_g0, V_s0, W_g0, W_s0,
     &       BC_hw_g, BC_hw_s,
     &       BC_Uw_g, BC_Uw_s, BC_Vw_g, BC_Vw_s, BC_Ww_g, BC_Ww_s,
     &       BC_hw_T_g, BC_hw_T_s, BC_Tw_g, BC_Tw_s, BC_C_T_g,
     &       BC_C_T_s,
     &       BC_hw_X_g, BC_hw_X_s, BC_Xw_g, BC_Xw_s, BC_C_X_g,
     &       BC_C_X_s,
     &       BC_hw_Theta_m, BC_Thetaw_m, BC_C_Theta_m,
     &       BC_Theta_m
      COMMON / BC_I /
     &       BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t,
     &       IJK_P_g, BC_JJ_PS
C
      COMMON / BC_L /
     &       BC_DEFINED  
C
      COMMON / BC_C /
     &       BC_TYPE, BC_PLANE
