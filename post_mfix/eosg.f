CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: EOSG (MW, PG, TG)                                      C
C  Purpose: Equation of state for gas                                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 29-JAN-92  C
C  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: GAS_CONST                                     C
C  Variables modified: EOSG                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      DOUBLE PRECISION FUNCTION EOSG (MW, PG, TG)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'scales.inc'
C
C                      Molecular weight
      DOUBLE PRECISION MW
C
C                      Gas pressure
      DOUBLE PRECISION PG
C
C                      Gas temperature
      DOUBLE PRECISION TG
C
C                      dummy variable in sc_p_g2.inc
      DOUBLE PRECISION XXX
C
      INCLUDE 'sc_p_g1.inc'
      INCLUDE 'sc_p_g2.inc'
C
      EOSG = UNSCALE(PG) * MW / (GAS_CONST * TG)
C
      RETURN
      END

CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: dROodP_g(ROG, PG)                                      C
C  Purpose: derivative of gas density w.r.t pressure                   C
C                                                                      C
C  Author: M. Syamlal                                 Date: 14-AUG-96  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: GAS_CONST                                     C
C  Variables modified: EOSG                                            C
C                                                                      C
C  Local variables: None                                               C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      DOUBLE PRECISION FUNCTION dROodP_g (ROG, Pg)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'scales.inc'
C
C                      gas density
      DOUBLE PRECISION ROg
C
C                      Gas pressure
      DOUBLE PRECISION PG
C
      dROodP_g = ROg / (Pg + P_ref)
C
      RETURN
      END
