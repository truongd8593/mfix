!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ic.inc                                                 C
!  Purpose: Common block containing initial conditions data            C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
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


      MODULE ic


      Use param
      Use param1


!
!                      x coordinate of the west face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_X_w (DIMENSION_IC)
!
!                      x coordinate of the east face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_X_e (DIMENSION_IC)
!
!                      y coordinate of the south face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_Y_s (DIMENSION_IC)
!
!                      y coordinate of the north face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_Y_n (DIMENSION_IC)
!
!                      z coordinate of the bottom face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_Z_b (DIMENSION_IC)
!
!                      z coordinate of the top face of a region where
!                      initial conditions are specified
      DOUBLE PRECISION IC_Z_t (DIMENSION_IC)
!
!                      i index of the west face of a region where
!                      initial conditions are specified
      INTEGER          IC_I_w (DIMENSION_IC)
!
!                      i index of the east face of a region where
!                      initial conditions are specified
      INTEGER          IC_I_e (DIMENSION_IC)
!
!                      j index of the south face of a region where
!                      initial conditions are specified
      INTEGER          IC_J_s (DIMENSION_IC)
!
!                      j index of the north face of a region where
!                      initial conditions are specified
      INTEGER          IC_J_n (DIMENSION_IC)
!
!                      k index of the bottom face of a region where
!                      initial conditions are specified
      INTEGER          IC_K_b (DIMENSION_IC)
!
!                      k index of the top face of a region where
!                      initial conditions are specified
      INTEGER          IC_K_t (DIMENSION_IC)
!
!                      Initial void fraction in a specified region
      DOUBLE PRECISION IC_EP_g (DIMENSION_IC)
!
!                      Initial gas pressure in a specified region
      DOUBLE PRECISION IC_P_g (DIMENSION_IC)
!
!                      Initial gas pressure in a specified region
      DOUBLE PRECISION IC_P_star (DIMENSION_IC)
!
!                      Initial turbulence length scale in a specified region
      DOUBLE PRECISION IC_L_scale (DIMENSION_IC)
!
!                      Initial macroscopic density of solids phases in a
!                      specified region
      DOUBLE PRECISION IC_ROP_s (DIMENSION_IC, DIM_M)
!
!                      Initial gas phase temperature in a specified region
      DOUBLE PRECISION IC_T_g (DIMENSION_IC)
!
!                      Initial solids phase temperature in a specified
!                      region
      DOUBLE PRECISION IC_T_s (DIMENSION_IC, DIM_M)
!
!                      Initial granular temperature in a specified
!                      region
      DOUBLE PRECISION IC_Theta_m (DIMENSION_IC, DIM_M)
!
!                      Initial x-component of gas velocity in a specified
!                      region
      DOUBLE PRECISION IC_U_g (DIMENSION_IC)
!
!                      Initial x-component of solids phase velocity in a
!                      specified region
      DOUBLE PRECISION IC_U_s (DIMENSION_IC, DIM_M)
!
!                      Initial y-component of gas velocity in a specified
!                      region
      DOUBLE PRECISION IC_V_g (DIMENSION_IC)
!
!                      Initial y-component of solids phase velocity in a
!                      specified region
      DOUBLE PRECISION IC_V_s (DIMENSION_IC, DIM_M)
!
!                      Initial z-component of gas velocity in a specified
!                      region
      DOUBLE PRECISION IC_W_g (DIMENSION_IC)
!
!                      Type of initial condition: PATCH
      CHARACTER*16     IC_TYPE (DIMENSION_IC)
!
!                      Initial z-component of solids phase velocity in a
!                      specified region
      DOUBLE PRECISION IC_W_s (DIMENSION_IC, DIM_M)
!
!                      Logical variable to determine whether an ic is defined
      LOGICAL          IC_DEFINED (DIMENSION_IC)
!
!                      Initial gas species mass fractions in a specified region
      DOUBLE PRECISION IC_X_g (DIMENSION_IC, DIM_N_g)
!
!                      Initial solids species mass fractions in a
!                      specified region
      DOUBLE PRECISION IC_X_s (DIMENSION_IC, DIM_M, DIM_N_s)
!
!                      Gas phase radiation coefficient
      DOUBLE PRECISION IC_GAMA_Rg (DIMENSION_IC)
!
!                      Gas phase radiation temperature
      DOUBLE PRECISION IC_T_Rg (DIMENSION_IC)
!
!                      Solids phase-1 radiation coefficient
      DOUBLE PRECISION IC_GAMA_Rs (DIMENSION_IC, DIM_M)
!
!                      Solids phase-1 radiation temperature
      DOUBLE PRECISION IC_T_Rs (DIMENSION_IC, DIM_M)
!
!


      END MODULE ic                                                                              
