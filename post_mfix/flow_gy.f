CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_gy(I, J, K, IJK)                                 C
C  Purpose: Calculate volumetric flow of gas in y direction            C
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
      REAL FUNCTION VFLOW_gy(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IJPK
C
      INCLUDE 'function.inc'
C
      IF(V_g(IJK) .GT. ZERO) THEN
        VFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        VFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_gy(I, J, K, IJK)                                 C
C  Purpose: Calculate mass flow of gas in y direction                  C
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
      REAL FUNCTION MFLOW_gy(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IJPK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(V_g(IJK) .GT. ZERO) THEN
        MFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        MFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)
     &             * CALC_RO_g(IJPK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_gy(IJK)                                           C
C  Purpose: Calculate mass flux of gas in y direction                  C
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
      REAL FUNCTION FLUX_gy(IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, IJPK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(V_g(IJK) .GT. ZERO) THEN
        FLUX_gy = V_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IJPK = JP_OF(IJK)
        FLUX_gy = V_g(IJK) * EP_g(IJPK)
     &             * CALC_RO_g(IJPK)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_gy(I, J, K, IJK, N)                              C
C  Purpose: Calculate gas species mass flow in y direction             C
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
      REAL FUNCTION XFLOW_gy(I, J, K, IJK, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, N, IJPK
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(V_g(IJK) .GT. ZERO) THEN
        XFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IJPK = JP_OF(IJK)
        XFLOW_gy = DX(I) * X(I) * DZ(K) * V_g(IJK) * EP_g(IJPK)
     &             * CALC_RO_g(IJPK) * X_g(IJPK, N)
      ENDIF
C
      RETURN
      END
