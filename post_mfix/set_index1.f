CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,       C
C                            IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,       C
C                           IJKB, IJKT, IM, JM, KM)                    C
C  Purpose: Set the indices of the first neighbors of cell ijk         C
C           This version adds 'increments' stored in STORE_INCREMENTS  C
C           to IJK to find indices of the neighbors.                   C
C                                                                      C
C  Author: M. Syamlal, W. A. Rogers                   Date: 17-Dec-91  C
C  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
C           W. Rogers                                                  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Second index of store_increments changed to PARAMETER's.   C
C           Use IM_OF etc. instead of BOUND_FUNIJK for computing     C
C           IMJK etc.                                                  C
C  Author: M. Syamlal                                 Date: 18-FEB-92  C
C  Revision Number:2                                                   C
C  Purpose: change STORE_INCREMENTS to INCREMENT_FOR_xx.  Do only      C
C           calculation of the nearest neighbors.  Remove MIN and MAX  C
C           from IM, IP etc. calculations.                             C
C  Author: M. Syamlal                                 Date: 18-SEP-92  C
C  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: CELL_CLASS, INCREMENTS_FOR_xx, IJK, I, J, K   C
C                                                                      C
C  Variables modified: IJKN,IJKS,IJKE,IJKW,IJKT,IJKB, IJKM, IJMK, IMJK,C
C                      IPJK, IJPK, IJKP                                C
C                                                                      C
C  Local variables: ICLASS                                             C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SET_INDEX1(IJK, I, J, K, IMJK, IPJK, IJMK, IJPK,
     &                       IJKM, IJKP, IJKW, IJKE, IJKS, IJKN,
     &                       IJKB, IJKT, IM, JM, KM)
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
C
      INCLUDE 'geometry.inc'
      INCLUDE 'constant.inc'
C
C     Indices
C
      INCLUDE 'indices.inc'
C
C Local Variables
C
C                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP,
     &                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT,
     &                 IM, JM, KM
C
C                      Denotes cell class, column of STORE_INCREMENTS
      INTEGER          ICLASS
C
      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
C
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
C
C     Determine the effective indices of neighboring cells which are used to
C     define neighboring cell-centered quantities
C
      ICLASS = CELL_CLASS(IJK) !Determine the class of cell IJK
C
C     Determine the true indices of neighboring cells
C
C
      IJKW  = IJK + INCREMENT_FOR_w (ICLASS)
      IJKE  = IJK + INCREMENT_FOR_e (ICLASS)
      IJKS  = IJK + INCREMENT_FOR_s (ICLASS)
      IJKN  = IJK + INCREMENT_FOR_n (ICLASS)
      IJKB  = IJK + INCREMENT_FOR_b (ICLASS)
      IJKT  = IJK + INCREMENT_FOR_t (ICLASS)
      IMJK  = IJK + INCREMENT_FOR_im(ICLASS)
      IPJK  = IJK + INCREMENT_FOR_ip(ICLASS)
      IJMK  = IJK + INCREMENT_FOR_jm(ICLASS)
      IJPK  = IJK + INCREMENT_FOR_jp(ICLASS)
      IJKM  = IJK + INCREMENT_FOR_km(ICLASS)
      IJKP  = IJK + INCREMENT_FOR_kp(ICLASS)
C
      RETURN
      END
