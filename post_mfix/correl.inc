CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: correlations.inc                                       C
C  Purpose: Include file containing variables used for computing       C
C           averages, variances, and correlations.                     C
C                                                                      C
C  Author: M. Syamlal                                 Date: 01-AUG-92  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:                                               C
C  Variables modified:                                                 C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
C                      Array for accumulating the sum of EP_g
      DOUBLE PRECISION SUM_EP_g (DIMENSION_3)
C
C                      Average value of EP_g
      DOUBLE PRECISION AVG_EP_g (DIMENSION_3)
C
C                      Array for accumulating the sum of EP_g*EP_g
      DOUBLE PRECISION SUM_EPxEP_g (DIMENSION_3)
C
C                      Standard deviation of EP_g.
      DOUBLE PRECISION SDV_EP_g (DIMENSION_3)
C
C                      Array for accumulating the sum of V_g
      DOUBLE PRECISION SUM_V_g (DIMENSION_3)
C
C                      Average value of V_g
      DOUBLE PRECISION AVG_V_g (DIMENSION_3)
C
C                      Array for accumulating the sum of V_g*V_g
      DOUBLE PRECISION SUM_VxV_g (DIMENSION_3)
C
C                      Standard deviation of V_g.
      DOUBLE PRECISION SDV_V_g (DIMENSION_3)
C
C                      A flag to check whether the subroutine is called for 
C                      the first time
      INTEGER          STARTED
C
C                      Number of time step data summed
      INTEGER          NSUM
C
      COMMON /CORRELATIONS_DP/ 
     & SUM_EP_g, AVG_EP_g, SUM_EPxEP_g, SDV_EP_g,
     & SUM_V_g, AVG_V_g, SUM_VxV_g, SDV_V_g
C
      COMMON /CORRELATIONS_I/ STARTED, NSUM 
