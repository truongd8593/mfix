!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: field_variables.inc                                    C
!  Purpose: Common block containing field variables data               C
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
 
 
      MODULE fldvar
 
 
      Use param
      Use param1
 
 
!
!                      Void fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  EP_g 
!
!                      Previous-time-step value of Void fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  EP_go 
!
!                      Gas pressure
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_g 
!
!                      Previous-time-step value of Gas pressure
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_go 
!
!                      Gas density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  RO_g 
!
!                      Previous-time-step value of Gas density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  RO_go 
!
!                      Macroscopic gas density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_g 
!
!                      Previous-time-step value of macroscopic gas density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_go 
!
!                      Macroscopic density of solids phases
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_s 
!
!                      Previous-time-step value of macroscopic density of
!                      solids phases
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_so 
!
!                      Gas phase temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  T_g 
!
!                      Solids phase temperature
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  T_s 
!
!                      Previous-time-step value of Gas phase temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  T_go 
!
!                      Previous-time-step value of Solids phase temperature
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  T_so 
!
!                      Gas species mass fraction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  X_g 
!
!                      Solids species mass fraction
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  X_s 
!
!                      Previous-time-step value of Gas species mass fraction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  X_go 
!
!                      Previous-time-step value of Solids species mass fraction
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  X_so 
!
!                      x-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  U_g 
!
!                      Previous-time-step value of x-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  U_go 
!
!                      x-component of solids phase velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  U_s 
!
!                      Previous-time-step value of x-component of solids phase
!                      velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  U_so 
!
!                      y-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  V_g 
!
!                      Previous time-step value of y-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  V_go 
!
!                      y-component of solids phase velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  V_s 
!
!                      Previous time-step value of y-component of solids phase
!                      velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  V_so 
!
!                      z-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  W_g 
!
!                      Previous time-step value of z-component of gas velocity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  W_go 
!
!                      z-component of solids phase velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  W_s 
!
!                      Previous time-step value of z-component of solids phase
!                      velocity
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  W_so 
!
!                      Solids pressure as a result of granular motion
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  P_s 
!
!                      Solids pressure as a result of granular motion
!                      Collisional Contribution
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  P_s_c 
!
!                      Solids pressure that maintains EP_g >= EP_star
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_star 
!
!                      Previous-time-step value of Solids pressure that maintains EP_g >= EP_star
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_staro 
!
!                      Granular temperature of mth phase
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  THETA_m 
!
!                      Previous-time-step value of Granular temperature of mth phase
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  THETA_mo 
!
 
!HPF$ align EP_g(:) with TT(:)
!HPF$ align EP_go(:) with TT(:)
!HPF$ align P_g(:) with TT(:)
!HPF$ align P_go(:) with TT(:)
!HPF$ align RO_g(:) with TT(:)
!HPF$ align RO_go(:) with TT(:)
!HPF$ align ROP_g(:) with TT(:)
!HPF$ align ROP_go(:) with TT(:)
!HPF$ align ROP_s(:, *) with TT(:)
!HPF$ align ROP_so(:, *) with TT(:)
!HPF$ align T_g(:) with TT(:)
!HPF$ align T_s(:, *) with TT(:)
!HPF$ align T_go(:) with TT(:)
!HPF$ align T_so(:, *) with TT(:)
!HPF$ align X_g(:, *) with TT(:)
!HPF$ align X_s(:, *, *) with TT(:)
!HPF$ align X_go(:, *) with TT(:)
!HPF$ align X_so(:, *, *) with TT(:)
!HPF$ align U_g(:) with TT(:)
!HPF$ align U_go(:) with TT(:)
!HPF$ align U_s(:, *) with TT(:)
!HPF$ align U_so(:, *) with TT(:)
!HPF$ align V_g(:) with TT(:)
!HPF$ align V_go(:) with TT(:)
!HPF$ align V_s(:, *) with TT(:)
!HPF$ align V_so(:, *) with TT(:)
!HPF$ align W_g(:) with TT(:)
!HPF$ align W_go(:) with TT(:)
!HPF$ align W_s(:, *) with TT(:)
!HPF$ align W_so(:, *) with TT(:)
!HPF$ align P_s(:, *) with TT(:)
!HPF$ align P_s_c(:, *) with TT(:)
!HPF$ align P_star(:) with TT(:)
!HPF$ align P_staro(:) with TT(:)
!HPF$ align THETA_m(:, *) with TT(:)
!HPF$ align THETA_mo(:, *) with TT(:)

 
      END MODULE fldvar
