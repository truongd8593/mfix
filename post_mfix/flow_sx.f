CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_sx(I, J, K, IJK, M)                              C
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
      REAL FUNCTION VFLOW_sx(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IPJK
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      IF(U_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IPJK = IP_OF(IJK)
        VFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * EP_s(IPJK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_sx(I, J, K, IJK, M)                              C
C  Purpose: Calculate mass flow of solids in x direction               C
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
      REAL FUNCTION MFLOW_sx(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IPJK
C
      INCLUDE 'function.inc'
C
      IF(U_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IPJK = IP_OF(IJK)
        MFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M) * ROP_s(IPJK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_sx(IJK, M)                                        C
C  Purpose: Calculate mass flux of solids in x direction               C
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
      REAL FUNCTION FLUX_sx(IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, M, IPJK
C
      INCLUDE 'function.inc'
C
      IF(U_s(IJK, M) .GT. ZERO) THEN
        FLUX_sx = U_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IPJK = IP_OF(IJK)
        FLUX_sx = U_s(IJK, M) * ROP_s(IPJK, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_sx(I, J, K, IJK, M, N)                           C
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
      REAL FUNCTION XFLOW_sx(I, J, K, IJK, M, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, N, IPJK
C
      INCLUDE 'function.inc'
C
      IF(U_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M)
     &             * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IPJK = IP_OF(IJK)
        XFLOW_sx = DY(J) * X_E(I) * DZ(K) * U_s(IJK, M)
     &             * ROP_s(IPJK, M) * X_s(IPJK, M, N)
      ENDIF
C
      RETURN
      END
