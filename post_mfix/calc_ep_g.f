CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_EP_g(IJK)                                         C
C  Purpose: Calculate EP_g from known solids volume fractions          C
C                                                                      C
C  Author: M. Syamlal                                 Date: 29-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: MMAX, IJK                                     C
C  Variables modified: M, EP_g                                         C
C                                                                      C
C  Local variables: SUM_EPS                                            C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE CALC_EP_g(IJK)
C
      IMPLICIT NONE
C
C  Include param.inc file to specify parameter values
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C     Physical and Numerical Parameters Section
C
      INCLUDE 'physprop.inc'
C
C     Field Variables
C
      INCLUDE 'fldvar.inc'
      INCLUDE 'geometry.inc'
C
C     Indices
C
      INCLUDE 'indices.inc'
C
C  Local variables
C
C                      Sum of EP_s
      DOUBLE PRECISION SUM_EPS
      INTEGER          IJK, M
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      SUM_EPS = ZERO
      DO 50 M = 1, MMAX
        SUM_EPS = SUM_EPS + EP_s(IJK, M)
50    CONTINUE
      EP_g(IJK) = ONE - SUM_EPS
C
      RETURN
      END
