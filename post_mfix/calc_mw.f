CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_MW (X_g, DIM, L, NMAX, MW_g)                      C
C  Purpose: Calculate average molecular weight of gas                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 19-OCT-92  C
C  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced:None                                           C
C  Variables modified:None                                             C
C                                                                      C
C  Local variables: SUM, N                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      DOUBLE PRECISION FUNCTION CALC_MW(X_g, DIM, L, NMAX, MW_g)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'toleranc.inc'
C
C  passed variables
C
C                      Mass fraction array's Ist dimension
      INTEGER          DIM
C
C                      Mass fraction array
      DOUBLE PRECISION X_g(DIM, *)
C
C                      Moleculare weight array
      DOUBLE PRECISION MW_g(*)
C
C                      Mass fraction array Ist index
      INTEGER          L
C
C                      Max of X_g array 2nd index and MW_g array index
      INTEGER          NMAX
C
C  Local variable
C
C                      local sum
      DOUBLE PRECISION SUM
C
C                      local index
      INTEGER          N
C
      SUM = ZERO
      DO 10 N = 1, NMAX
        SUM = SUM + X_g(L, N) / MW_g(N)
10    CONTINUE
      CALC_MW = ONE/MAX(SUM,oMW_MAX)
C
      RETURN
      END
