!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_IC                                                 !
!                                                                      !
!  Purpose: Common elements needed for DES initial conditions.         !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_IC

      USE param

! Logical variable to determine whether an initial condition is defined
      LOGICAL DES_IC_DEFINED (DIMENSION_IC)

! Set to true if any DES initial condition regions are specified. This
! logical is used to invoke the call to DES_SET_IC.
      LOGICAL DES_IC_EXIST

! Physical coordinates indicating initial condition region
      DOUBLE PRECISION DES_IC_X_w(DIMENSION_IC)
      DOUBLE PRECISION DES_IC_X_e(DIMENSION_IC)
      DOUBLE PRECISION DES_IC_Y_s(DIMENSION_IC)
      DOUBLE PRECISION DES_IC_Y_n(DIMENSION_IC)
      DOUBLE PRECISION DES_IC_Z_b(DIMENSION_IC)
      DOUBLE PRECISION DES_IC_Z_t(DIMENSION_IC)

! Initial solids phase temperature in a specified region
      DOUBLE PRECISION DES_IC_T_s (DIMENSION_IC, DIM_M)

! Initial solids species mass fractions in a specified region
      DOUBLE PRECISION DES_IC_X_s (DIMENSION_IC, DIM_M, DIM_N_s)

! Density of unreacted, shrinking core
      DOUBLE PRECISION DES_IC_CORE_Rho(DIMENSION_IC, DIM_M)

      END MODULE DES_IC

