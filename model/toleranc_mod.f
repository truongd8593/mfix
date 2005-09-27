!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: tolerance.inc                                          C
!  Purpose: Specify all tolerance parameters                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      MODULE toleranc


      Use param
      Use param1


!
!
!                      Minimum value of solids volume fraction tracked
      DOUBLE PRECISION, PARAMETER          ::  ZERO_EP_s = 1.0D-8
!
!                      Small value for species mass fraction for disregarding
!                      residual calculation
      DOUBLE PRECISION, PARAMETER          ::  ZERO_X_gs = 1.0D-7
 
!                      Dilute flow threshold.  When the volume fraction of a
!                      certain phase in a cell is smaller than this value the
!                      momentum equation for that phase is not solved in the cell.
      DOUBLE PRECISION, PARAMETER          ::  DIL_EP_s = 1.0D-4
!
!                      Tolerance used for comparing two numbers for equality
!                      in function compare(a, b)
      DOUBLE PRECISION, PARAMETER          ::  TOL_COM = 1.0D-4
!
!                      Upper bound for temperatures
      DOUBLE PRECISION, PARAMETER          ::  TMAX = 4000.D0
!
!                      Lower bound for temperatures
      DOUBLE PRECISION, PARAMETER          ::  TMIN = 250.D0
!
!                      Reciprocal of a maximum molecular weight
      DOUBLE PRECISION, PARAMETER          ::  oMW_MAX = (ONE/500.D0)
!
!                      Maximum value of velocities set to avoid divergence problems.
      DOUBLE PRECISION MAX_INLET_VEL
!
!                      User definable factor used to scale MAX_INLET_VEL.  Default value is 1.
      DOUBLE PRECISION MAX_INLET_VEL_FAC
!
!                      Maximum allowed velocity of gas or solids in case no inlet velocities
!                      (or zero velocities) are defined at inlet (see function check_vel_bound)
      DOUBLE PRECISION, PARAMETER          ::  MAX_ALLOWED_VEL = 500.0D+2
!
!  The following quantities can be specified through the input data file, with
!  namelist inputs of the same name.
!
!                      Tolerance in residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID
!
!                      Tolerance in energy eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_T
!
!                      Tolerance in species eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_X
!
!                      Tolerance in scalr eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_Scalar
!
!                      Tolerance in K & Epsilon eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_K_Epsilon
!
!                      Tolerance in Granular Temperature eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_Th
!
!                      Minimum residual for declaring divergence
      DOUBLE PRECISION TOL_DIVERGE
!
!                      Factor for normalizing the residual of gas cont. eq.
      DOUBLE PRECISION NORM_g
!
!                      Factor for normalizing the residual of solids cont. eq.
      DOUBLE PRECISION NORM_s
!
!


      END MODULE toleranc                                                                        
