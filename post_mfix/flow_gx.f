CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_gx(I, J, K, IJK)                                 C
C  Purpose: Calculate volumetric flow of gas in x direction            C
C                                                                      C
C  Author: M. Syamlal                                 Date: 22-NOV-93  C
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
      REAL FUNCTION VFLOW_gx(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IPJK
C
      INCLUDE 'function.inc'
C
      IF(U_g(IJK) .GT. ZERO) THEN
        VFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        VFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_gx(I, J, K, IJK)                                 C
C  Purpose: Calculate mass flow of gas in x direction                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 22-NOV-93  C
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
      REAL FUNCTION MFLOW_gx(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IPJK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(U_g(IJK) .GT. ZERO) THEN
        MFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        MFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)
     &             * CALC_RO_g(IPJK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_gx(IJK)                                           C
C  Purpose: Calculate mass flux of gas in x direction                  C
C                                                                      C
C  Author: M. Syamlal                                 Date: 22-NOV-93  C
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
      REAL FUNCTION FLUX_gx(IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, IPJK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(U_g(IJK) .GT. ZERO) THEN
        FLUX_gx = U_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IPJK = IP_OF(IJK)
        FLUX_gx = U_g(IJK) * EP_g(IPJK)
     &             * CALC_RO_g(IPJK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_gx(I, J, K, IJK, N)                              C
C  Purpose: Calculate gas species mass flow in x direction             C
C                                                                      C
C  Author: M. Syamlal                                 Date: 22-NOV-93  C
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
      REAL FUNCTION XFLOW_gx(I, J, K, IJK, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, N, IPJK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(U_g(IJK) .GT. ZERO) THEN
        XFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IPJK = IP_OF(IJK)
        XFLOW_gx = DY(J) * X_E(I) * DZ(K) * U_g(IJK) * EP_g(IPJK)
     &             * CALC_RO_g(IPJK) * X_g(IPJK, N)
      ENDIF
C
      RETURN
      END
