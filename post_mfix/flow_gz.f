CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: VFLOW_gz(I, J, K, IJK)                                 C
C  Purpose: Calculate volumetric flow of gas in z direction            C
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
      REAL FUNCTION VFLOW_gz(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IJKP
C
      INCLUDE 'function.inc'
C
      IF(W_g(IJK) .GT. ZERO) THEN
        VFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        VFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: MFLOW_gz(I, J, K, IJK)                                 C
C  Purpose: Calculate mass flow of gas in z direction                  C
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
      REAL FUNCTION MFLOW_gz(I, J, K, IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, IJKP
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(W_g(IJK) .GT. ZERO) THEN
        MFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        MFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)
     &             * CALC_RO_g(IJKP)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: FLUX_gz(IJK)                                           C
C  Purpose: Calculate mass flux of gas in z direction                  C
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
      REAL FUNCTION FLUX_gz(IJK)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER IJK, IJKP
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(W_g(IJK) .GT. ZERO) THEN
        FLUX_gz = W_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK)
      ELSE
        IJKP = KP_OF(IJK)
        FLUX_gz = W_g(IJK) * EP_g(IJKP)
     &             * CALC_RO_g(IJKP)
      ENDIF
C
      RETURN
      END
CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: XFLOW_gz(I, J, K, IJK, N)                              C
C  Purpose: Calculate gas species mass flow in z direction             C
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
      REAL FUNCTION XFLOW_gz(I, J, K, IJK, N)
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
C
      INTEGER I, J, K, IJK, N, IJKP
C
C  Function subroutines
C
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C
      IF(W_g(IJK) .GT. ZERO) THEN
        XFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJK) 
     &             * CALC_RO_g(IJK) * X_g(IJK, N)
      ELSE
        IJKP = KP_OF(IJK)
        XFLOW_gz = DX(I) * DY(J) * W_g(IJK) * EP_g(IJKP)
     &             * CALC_RO_g(IJKP) * X_g(IJKP, N)
      ENDIF
C
      RETURN
      END
