CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SET_GEOMETRY                                           C
C  Purpose: Calculate X, X_E,  oX, oX_E                                C
C                                                                      C
C  Author: M. Syamlal                                 Date: 21-JAN-92  C
C  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
C           W. Rogers                                                  C
C                                                                      C
C  Revision Number: 1                                                  C
C  Purpose: Fix bugs                                                   C
C  Author: M. Syamlal                                 Date: 10-FEB-92  C
C  Revision Number: 2                                                  C
C  Purpose: Include logic for Variable Grid Spacing capability         C
C  Author: W. Rogers                                  Date: 06-APR-92  C
C  Revision Number: 3                                                  C
C  Purpose: Add oX, and oX_E calculations                              C
C  Author: M. Syamlal                                 Date: 8-MAY-92   C
C  Revision Number: 4                                                  C
C  Purpose: Add FX, FX_bar, FX_E, FX_E_bar, FY_N, FY_N_bar, FZ_T, and  C
C           FZ_T_bar calculations                                      C
C  Author: M. Syamlal                                 Date: 31-AUG-92  C
C  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: COORDINATES, IMAX2, DX, JMAX2, DY, KMAX2,     C
C                        DZ,                                           C
C                                                                      C
C  Variables modified: X, X_E, I,                                      C
C                      J, K,  oX,                                      C
C                      oX_E, FX, FX_bar, FX_E, FX_E_bar, FY_N,         C
C                      FY_N_bar, FZ_T, FZ_T_bar                        C
C                                                                      C
C  Local variables: DX_E, DY_N, DZ_T                                   C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SET_GEOMETRY
C
      IMPLICIT NONE
C
C  Include param.inc file to specify parameter values
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
C
C     Run control data
C
      INCLUDE 'run.inc'
C
C     Geometry and Discretization Section
C
      INCLUDE 'geometry.inc'
C
C
C  Function calls
C
      LOGICAL COMPARE
C
C  Local variables
C
C                      Indices
      INTEGER          I, J, K
C
C     X-direction dimension of U-momentum cell
      DOUBLE PRECISION DX_E
C
C     Y-direction dimension of V-momentum cell
      DOUBLE PRECISION DY_N
C
C     Z-direction dimension of W-momentum cell
      DOUBLE PRECISION DZ_T
C
      IF(CYCLIC_X_PD)CYCLIC_X = .TRUE.
      IF(CYCLIC_Y_PD)CYCLIC_Y = .TRUE.
      IF(CYCLIC_Z_PD)CYCLIC_Z = .TRUE.
      CYCLIC = CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z

      IF(CYLINDRICAL
     &  .AND.  COMPARE(ZLENGTH, 8.D0*ATAN(ONE)) .AND. DO_K ) THEN
        CYCLIC_Z = .TRUE.
      ENDIF
C
      IF(CYCLIC_X) THEN
        DX(1)     = DX(IMAX1)
        DX(IMAX2) = DX(IMIN1)
      ENDIF
      IF(CYCLIC_Y) THEN
        DY(1)     = DY(JMAX1)
        DY(JMAX2) = DY(JMIN1)
      ENDIF
      IF(CYCLIC_Z) THEN
        DZ(1)     = DZ(KMAX1)
        DZ(KMAX2) = DZ(KMIN1)
      ENDIF
C
      IF( COORDINATES .EQ. 'CARTESIAN')THEN
        DO 100 I = 1, IMAX2
           X(I)      = ONE
           X_E(I)     = ONE
           oX(I)     = ONE
           oX_E(I)    = ONE
           oDX(I)   = ONE / DX(I)
100     CONTINUE
      ELSEIF( CYLINDRICAL)THEN
        IF(XMIN .EQ. ZERO)THEN
          oDX(1)   = ONE / DX(1)
          oX(1)      = UNDEFINED
          oX_E(1)     = UNDEFINED
          IF(DO_I) THEN
            X(1)  = - HALF * DX(1)
            X_E(1) = 0.0
          ELSE
            X(1)  = HALF * DX(1)
            X_E(1) = DX(1)
          ENDIF
        ELSE
          IF(DO_I) THEN
            X_E(1) = XMIN
            X(1)   = XMIN - HALF * DX(1)
          ELSE
            X_E(1) = XMIN +  DX(1)
            X(1)   = XMIN + HALF * DX(1)
          ENDIF
          oX(1)      = ONE / X(1)
          oX_E(1)     = ONE / X_E(1)
          oDX(1)   = ONE / DX(1)
        ENDIF
C                                                #1 add the DO_I IF block
        IF(DO_I) THEN
          DO 200 I = IMIN1, IMAX2
            X(I)       = X(I-1) + (DX(I-1) + DX(I) ) / 2.
            X_E(I)      = X_E(I-1) + DX(I)
            oX(I)      = ONE / X(I)
            oX_E(I)     = ONE / X_E(I)
            oDX(I)   = ONE / DX(I)
200       CONTINUE
        ENDIF
      ENDIF
C
      DO 300 J = 1, JMAX2
        oDY(J) = ONE / DY(J)
300   CONTINUE
C
      DO 400 K = 1, KMAX2
        oDZ(K)   = ONE / DZ(K)
400   CONTINUE
C
C================================================#2 
C       Look at first U-, V-, and W-momentum cells
        DX_E = HALF * (DX(1) + DX(IMIN1))
        DY_N = HALF * (DY(1) + DY(JMIN1))
        DZ_T = HALF * (DZ(1) + DZ(KMIN1))
C
C
        oDX_E(1)    = ONE / DX_E
        oDY_N(1)    = ONE / DY_N
        oDZ_T(1)    = ONE / DZ_T
        FX(1)       = HALF
        FX_bar(1)   = HALF
        FX_E(1)     = HALF
        FX_E_bar(1) = HALF
        FY_N(1)     = HALF
        FY_N_bar(1) = HALF
        FZ_T(1)     = HALF
        FZ_T_bar(1) = HALF
C       ..........................................
C       Look at 2 through IMAX1 U-momentum cells
        IF(DO_I) THEN
          DO 210 I = IMIN1, IMAX1
            DX_E = HALF * ( DX(I+1) + DX(I) )
            oDX_E(I)    = ONE / DX_E
c            FX(I)       = DX_E * X_E(I) /
c     &             (DX_E * X_E(I) + HALF * (DX(I) + DX(I-1)) * X_E(I-1))
            FX(I)       = HALF
            FX_bar(I)   = ONE - FX(I)
c            FX_E(I)     = X(I+1) * DX(I+1) /
c     &                  ( X(I+1) * DX(I+1) + X(I) * DX(I) )
            FX_E(I)     = DX(I+1) / ( DX(I+1) + DX(I) )
            FX_E_bar(I) = ONE - FX_E(I)
210       CONTINUE
        ENDIF
C
C       Look at 2 through JMAX1 V-momentum cells
        IF(DO_J) THEN
          DO 211 J = JMIN1, JMAX1
            DY_N = HALF * ( DY(J+1) + DY(J) )
            oDY_N(J)   = ONE / DY_N
            FY_N(J)  = DY(J+1) / ( DY(J+1) + DY(J) )
            FY_N_bar(J) = ONE - FY_N(J)
211       CONTINUE
        ENDIF
C
C       Look at 2 through KMAX1 W-momentum cells
        IF(DO_K) THEN
          DO 212 K = KMIN1, KMAX1
            DZ_T = HALF * ( DZ(K+1) + DZ(K) )
            oDZ_T(K)   = ONE / DZ_T
            FZ_T(K) = DZ(K+1) / ( DZ(K+1) + DZ(K) )
            FZ_T_bar(K) = ONE - FZ_T(K)
212       CONTINUE
        ENDIF
C       ..........................................
C       Look at last U-, V-, and W-momentum cells
        DX_E = DX(IMAX2)
        DY_N = DY(JMAX2)
        DZ_T = DZ(KMAX2)
        oDX_E(IMAX2)    = ONE / DX_E
        oDY_N(JMAX2)    = ONE / DY_N
        oDZ_T(KMAX2)    = ONE / DZ_T
        FX(IMAX2)       = HALF
        FX_bar(IMAX2)   = HALF
        FX_E(IMAX2)     = HALF
        FX_E_bar(IMAX2) = HALF
        FY_N(JMAX2)     = HALF
        FY_N_bar(JMAX2) = HALF
        FZ_T(KMAX2)     = HALF
        FZ_T_bar(KMAX2) = HALF
        IF(CYCLIC_X) THEN
          FX_E(1)     = FX_E(IMAX1)
          FX_E_bar(1) = FX_E_bar(IMAX1)
        ENDIF
        IF(CYCLIC_Y) THEN
          FY_N(1)     = FY_N(JMAX1)
          FY_N_bar(1) = FY_N_bar(JMAX1)
        ENDIF
        IF(CYCLIC_Z) THEN
          FZ_T(1)     = FZ_T(KMAX1)
          FZ_T_bar(1) = FZ_T_bar(KMAX1)
        ENDIF
C=====================================================================
C
      RETURN
      END
