CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_sy(I, J, K, IJK, M)                              C
C  Purpose: Calculate volumetric flow of solids in x direction         C
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
      REAL FUNCTION VFLOW_sy(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IJPK
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      IF(V_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IJPK = JP_OF(IJK)
        VFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * EP_s(IJPK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_sy(I, J, K, IJK, M)                              C
C  Purpose: Calculate mass flow of solids in y direction               C
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
      REAL FUNCTION MFLOW_sy(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IJPK
C
      INCLUDE 'function.inc'
C
      IF(V_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IJPK = JP_OF(IJK)
        MFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M) * ROP_s(IJPK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_sy(IJK, M)                                        C
C  Purpose: Calculate mass flux of solids in y direction               C
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
      REAL FUNCTION FLUX_sy(IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, M, IJPK
C
      INCLUDE 'function.inc'
C
      IF(V_s(IJK, M) .GT. ZERO) THEN
        FLUX_sy = V_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IJPK = JP_OF(IJK)
        FLUX_sy = V_s(IJK, M) * ROP_s(IJPK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_sy(I, J, K, IJK, M, N)                           C
C  Purpose: Calculate solids species mass flow in x direction          C
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
      REAL FUNCTION XFLOW_sy(I, J, K, IJK, M, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, N, IJPK
C
      INCLUDE 'function.inc'
C
      IF(V_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M)
     &             * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IJPK = JP_OF(IJK)
        XFLOW_sy = DX(I) * X(I) * DZ(K) * V_s(IJK, M)
     &             * ROP_s(IJPK, M) * X_s(IJPK, M, N)
      ENDIF
C
      RETURN
      END
