CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: tolerance.inc                                          C
C  Purpose: Specify all tolerance parameters                           C
C                                                                      C
C  Author: M. Syamlal                                 Date: 24-JUL-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C
C                      Minimum value of solids volume fraction tracked
      DOUBLE PRECISION ZERO_EP_s
      PARAMETER (ZERO_EP_s = 1.0E-8)

C                      Dilute flow threshold.  When the volume fraction of a 
C                      certain phase in a cell is smaller than this value the
C                      momentum equation for that phase is not solved in the cell. 
      DOUBLE PRECISION DIL_EP_s
      PARAMETER (DIL_EP_s = 1.0E-4)
C
C                      Tolerance used for comparing two numbers for equality
C                      in function compare(a, b)
      DOUBLE PRECISION TOL_COM
      PARAMETER (TOL_COM = 1.0E-4)
C
C                      Upper bound for temperatures
      DOUBLE PRECISION TMAX
      PARAMETER (TMAX = 4000.)
C
C                      Lower bound for temperatures
      DOUBLE PRECISION TMIN
      PARAMETER (TMIN = 250.)
C
C                      Reciprocal of a maximum molecular weight
      DOUBLE PRECISION oMW_MAX
      PARAMETER (oMW_MAX = 1./500.)
C
C  The following quantities can be specified through the input data file, with
C  namelist inputs of the same name.
C
C                      Tolerance in residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID
C
C                      Tolerance in energy eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_T
C
C                      Tolerance in species eq residuals allowed for convergence
      DOUBLE PRECISION TOL_RESID_X
C
C                      Minimum residual for declaring divergence
      DOUBLE PRECISION TOL_DIVERGE
C
C                      Factor for normalizing the residual of gas cont. eq.
      DOUBLE PRECISION NORM_g
C
C                      Factor for normalizing the residual of solids cont. eq.
      DOUBLE PRECISION NORM_s
C
C
      COMMON / TOLERANCE_DP /
     &       TOL_RESID, TOL_RESID_T, TOL_RESID_X, TOL_DIVERGE,
     &       NORM_g, NORM_s
C
c      COMMON / TOLERANCE_I /
c     &      
C
c      COMMON / TOLERANCE_L /
c     &       
