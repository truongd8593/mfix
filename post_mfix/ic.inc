CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: ic.inc                                                 C
C  Purpose: Common block containing initial conditions data            C
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
C                      initial conditions are specified
      DOUBLE PRECISION IC_X_w (DIMENSION_IC) 
C
C                      x coordinate of the east face of a region where
C                      initial conditions are specified
      DOUBLE PRECISION IC_X_e (DIMENSION_IC)
C
C                      y coordinate of the south face of a region where
C                      initial conditions are specified
      DOUBLE PRECISION IC_Y_s (DIMENSION_IC)
C
C                      y coordinate of the north face of a region where
C                      initial conditions are specified
      DOUBLE PRECISION IC_Y_n (DIMENSION_IC)
C
C                      z coordinate of the bottom face of a region where
C                      initial conditions are specified
      DOUBLE PRECISION IC_Z_b (DIMENSION_IC)
C
C                      z coordinate of the top face of a region where
C                      initial conditions are specified
      DOUBLE PRECISION IC_Z_t (DIMENSION_IC)
C
C                      i index of the west face of a region where
C                      initial conditions are specified
      INTEGER          IC_I_w (DIMENSION_IC) 
C
C                      i index of the east face of a region where
C                      initial conditions are specified
      INTEGER          IC_I_e (DIMENSION_IC)
C
C                      j index of the south face of a region where
C                      initial conditions are specified
      INTEGER          IC_J_s (DIMENSION_IC)
C
C                      j index of the north face of a region where
C                      initial conditions are specified
      INTEGER          IC_J_n (DIMENSION_IC)
C
C                      k index of the bottom face of a region where
C                      initial conditions are specified
      INTEGER          IC_K_b (DIMENSION_IC)
C
C                      k index of the top face of a region where
C                      initial conditions are specified
      INTEGER          IC_K_t (DIMENSION_IC)
C
C                      Initial void fraction in a specified region
      DOUBLE PRECISION IC_EP_g (DIMENSION_IC)
C
C                      Initial gas pressure in a specified region
      DOUBLE PRECISION IC_P_g (DIMENSION_IC)
C
C                      Initial gas pressure in a specified region
      DOUBLE PRECISION IC_P_star (DIMENSION_IC)
C
C                      Initial turbulence length scale in a specified region
      DOUBLE PRECISION IC_L_scale (DIMENSION_IC)
C
C                      Initial macroscopic density of solids phases in a
C                      specified region
      DOUBLE PRECISION IC_ROP_s (DIMENSION_IC, DIMENSION_M)
C
C                      Initial gas phase temperature in a specified region
      DOUBLE PRECISION IC_T_g (DIMENSION_IC)
C
C                      Initial solids phase temperature in a specified 
C                      region
      DOUBLE PRECISION IC_T_s (DIMENSION_IC, DIMENSION_M)
C
C                      Initial granular temperature in a specified 
C                      region
      DOUBLE PRECISION IC_Theta_m (DIMENSION_IC, DIMENSION_M)
C
C                      Initial x-component of gas velocity in a specified
C                      region
      DOUBLE PRECISION IC_U_g (DIMENSION_IC)
C
C                      Initial x-component of solids phase velocity in a
C                      specified region
      DOUBLE PRECISION IC_U_s (DIMENSION_IC, DIMENSION_M)
C
C                      Initial y-component of gas velocity in a specified
C                      region
      DOUBLE PRECISION IC_V_g (DIMENSION_IC)
C
C                      Initial y-component of solids phase velocity in a
C                      specified region
      DOUBLE PRECISION IC_V_s (DIMENSION_IC, DIMENSION_M)
C
C                      Initial z-component of gas velocity in a specified
C                      region
      DOUBLE PRECISION IC_W_g (DIMENSION_IC)
C
C                      Type of initial condition: PATCH
      CHARACTER*16     IC_TYPE (DIMENSION_IC)
C
C                      Initial z-component of solids phase velocity in a
C                      specified region
      DOUBLE PRECISION IC_W_s (DIMENSION_IC, DIMENSION_M)
C
C                      Logical variable to determine whether an ic is defined
      LOGICAL          IC_DEFINED (DIMENSION_IC)
C
C                      Initial gas species mass fractions in a specified region
      DOUBLE PRECISION IC_X_g (DIMENSION_IC, DIMENSION_N_g)
C
C                      Initial solids species mass fractions in a 
C                      specified region
      DOUBLE PRECISION IC_X_s (DIMENSION_IC, DIMENSION_M, DIMENSION_N_s)
C
C                      Gas phase radiation coefficient 
      DOUBLE PRECISION IC_GAMA_Rg (DIMENSION_IC)
C
C                      Gas phase radiation temperature 
      DOUBLE PRECISION IC_T_Rg (DIMENSION_IC)
C
C                      Solids phase-1 radiation coefficient 
      DOUBLE PRECISION IC_GAMA_Rs (DIMENSION_IC, DIMENSION_M)
C
C                      Solids phase-1 radiation temperature 
      DOUBLE PRECISION IC_T_Rs (DIMENSION_IC, DIMENSION_M)
C
C
      COMMON / IC_DP /
     &       IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t,
     &       IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s,
     &       IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s,
     &       IC_X_g, IC_X_s, IC_P_star, IC_L_scale,
     &       IC_GAMA_Rg, IC_T_Rg, IC_GAMA_Rs, IC_T_Rs,
     &       IC_Theta_m
      COMMON / IC_I /
     &       IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t
      COMMON / IC_L /
     &       IC_DEFINED
C
      COMMON / IC_C /
     &       IC_TYPE
