CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: field_variables.inc                                    C
C  Purpose: Common block containing field variables data               C
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
C                      Void fraction
      DOUBLE PRECISION EP_g (DIMENSION_3)
C
C                      Previous-time-step value of Void fraction
      DOUBLE PRECISION EP_go (DIMENSION_3)
C
C                      Gas pressure
      DOUBLE PRECISION P_g (DIMENSION_3)
C
C                      Previous-time-step value of Gas pressure
      DOUBLE PRECISION P_go (DIMENSION_3)
C
C                      Gas density
      DOUBLE PRECISION RO_g (DIMENSION_3)
C
C                      Previous-time-step value of Gas density
      DOUBLE PRECISION RO_go (DIMENSION_3)
C
C                      Macroscopic gas density
      DOUBLE PRECISION ROP_g (DIMENSION_3)
C
C                      Previous-time-step value of macroscopic gas density
      DOUBLE PRECISION ROP_go (DIMENSION_3)
C
C                      Macroscopic density of solids phases
      DOUBLE PRECISION ROP_s (DIMENSION_3, DIMENSION_M)
C
C                      Previous-time-step value of macroscopic density of
C                      solids phases
      DOUBLE PRECISION ROP_so (DIMENSION_3, DIMENSION_M)
C
C                      Gas phase temperature
      DOUBLE PRECISION T_g (DIMENSION_3)
C
C                      Solids phase temperature
      DOUBLE PRECISION T_s (DIMENSION_3, DIMENSION_M)
C
C                      Previous-time-step value of Gas phase temperature
      DOUBLE PRECISION T_go (DIMENSION_3)
C
C                      Previous-time-step value of Solids phase temperature
      DOUBLE PRECISION T_so (DIMENSION_3, DIMENSION_M)
C
C                      Gas species mass fraction
      DOUBLE PRECISION X_g (DIMENSION_3, DIMENSION_N_g)
C
C                      Solids species mass fraction
      DOUBLE PRECISION X_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s)
C
C                      Previous-time-step value of Gas species mass fraction
      DOUBLE PRECISION X_go (DIMENSION_3, DIMENSION_N_g)
C
C                      Previous-time-step value of Solids species mass fraction
      DOUBLE PRECISION X_so (DIMENSION_3, DIMENSION_M, DIMENSION_N_s)
C
C                      x-component of gas velocity
      DOUBLE PRECISION U_g (DIMENSION_3)
C
C                      Previous-time-step value of x-component of gas velocity
      DOUBLE PRECISION U_go (DIMENSION_3)
C
C                      x-component of solids phase velocity
      DOUBLE PRECISION U_s (DIMENSION_3, DIMENSION_M)
C
C                      Previous-time-step value of x-component of solids phase
C                      velocity
      DOUBLE PRECISION U_so (DIMENSION_3, DIMENSION_M)
C
C                      y-component of gas velocity
      DOUBLE PRECISION V_g (DIMENSION_3)
C
C                      Previous time-step value of y-component of gas velocity
      DOUBLE PRECISION V_go (DIMENSION_3)
C
C                      y-component of solids phase velocity
      DOUBLE PRECISION V_s (DIMENSION_3, DIMENSION_M)
C
C                      Previous time-step value of y-component of solids phase
C                      velocity
      DOUBLE PRECISION V_so (DIMENSION_3, DIMENSION_M)
C
C                      z-component of gas velocity
      DOUBLE PRECISION W_g (DIMENSION_3)
C
C                      Previous time-step value of z-component of gas velocity
      DOUBLE PRECISION W_go (DIMENSION_3)
C
C                      z-component of solids phase velocity
      DOUBLE PRECISION W_s (DIMENSION_3, DIMENSION_M)
C
C                      Previous time-step value of z-component of solids phase
C                      velocity
      DOUBLE PRECISION W_so (DIMENSION_3, DIMENSION_M)
C
C                      Solids pressure as a result of granular motion
      DOUBLE PRECISION P_s (DIMENSION_3, DIMENSION_M)
C
C                      Solids pressure that maintains EP_g >= EP_star
      DOUBLE PRECISION P_star (DIMENSION_3)
C
C                      Previous-time-step value of Solids pressure that maintains EP_g >= EP_star
      DOUBLE PRECISION P_staro (DIMENSION_3)
C
C                      Granular temperature of mth phase
      DOUBLE PRECISION THETA_m (DIMENSION_3, DIMENSION_M)
C
C                      Previous-time-step value of Granular temperature of mth phase
      DOUBLE PRECISION THETA_mo (DIMENSION_3, DIMENSION_M)
C
      COMMON /FIELD_VARIABLES/
     &      EP_g, EP_go, P_g, P_go, RO_g, RO_go, ROP_g, ROP_go,
     &      ROP_s, ROP_so,
     &      T_g, T_s, T_go, T_so,
     &      U_g , U_s, V_g,  V_s, W_g, W_s,
     &      U_go, V_go, W_go,
     &      U_so, V_so, W_so,
     &      P_s, P_star, P_staro, THETA_m, THETA_mo,
     &      X_g, X_s, X_go, X_so
