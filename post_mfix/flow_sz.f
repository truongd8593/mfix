CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_sz(I, J, K, IJK, M)                              C
C  Purpose: Calculate volumetric flow of solids in z direction         C
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
      REAL FUNCTION VFLOW_sz(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IJKP
C
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
C
      IF(W_s(IJK, M) .GT. ZERO) THEN
        VFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * EP_s(IJK, M)
      ELSE
        IJKP = KP_OF(IJK)
        VFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * EP_s(IJKP, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_sz(I, J, K, IJK, M)                              C
C  Purpose: Calculate mass flow of solids in z direction               C
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
      REAL FUNCTION MFLOW_sz(I, J, K, IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, IJKP
C
      INCLUDE 'function.inc'
C
      IF(W_s(IJK, M) .GT. ZERO) THEN
        MFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IJKP = KP_OF(IJK)
        MFLOW_sz = DX(I) * DY(J) * W_s(IJK, M) * ROP_s(IJKP, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_sz(IJK, M)                                        C
C  Purpose: Calculate mass flux of solids in z direction               C
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
      REAL FUNCTION FLUX_sz(IJK, M)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, M, IJKP
C
      INCLUDE 'function.inc'
C
      IF(W_s(IJK, M) .GT. ZERO) THEN
        FLUX_sz = W_s(IJK, M) * ROP_s(IJK, M) 
      ELSE
        IJKP = KP_OF(IJK)
        FLUX_sz = W_s(IJK, M) * ROP_s(IJKP, M)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_sz(I, J, K, IJK, M, N)                           C
C  Purpose: Calculate solids species mass flow in z direction          C
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
      REAL FUNCTION XFLOW_sz(I, J, K, IJK, M, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, M, N, IJKP
C
      INCLUDE 'function.inc'
C
      IF(W_s(IJK, M) .GT. ZERO) THEN
        XFLOW_sz = DX(I) * DY(J) * W_s(IJK, M)
     &             * ROP_s(IJK, M) * X_s(IJK, M, N)
      ELSE
        IJKP = KP_OF(IJK)
        XFLOW_sz = DX(I) * DY(J) * W_s(IJK, M)
     &             * ROP_s(IJKP, M) * X_s(IJKP, M, N)
      ENDIF
C
      RETURN
      END
