CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: INTERP_RES                                             C
C  Purpose: Interpolate from old RES file to create a new RES file     C
C                                                                      C
C  Author: M. Syamlal                                 Date: 03-DEC-93  C
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
      SUBROUTINE INTERP_RES
C
      IMPLICIT NONE
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'run.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'xforms.inc'
C
C  Function subroutines
C
      LOGICAL OPEN_FILEP
      INTEGER GET_INDEX
C
C  Local variables
C
      DOUBLE PRECISION
     &  EP_g_OLD(DIMENSION_3), P_g_OLD(DIMENSION_3),
     &  P_star_OLD(DIMENSION_3), RO_g_OLD(DIMENSION_3),
     &  ROP_g_OLD(DIMENSION_3), T_g_OLD(DIMENSION_3),
     &  T_s_OLD(DIMENSION_3, DIMENSION_M), 
     &  THETA_M_OLD(DIMENSION_3, DIMENSION_M),
     &  X_g_OLD(DIMENSION_3, DIMENSION_N_g), U_g_OLD(DIMENSION_3),
     &  V_g_OLD(DIMENSION_3), W_g_OLD(DIMENSION_3),
     &  ROP_s_OLD(DIMENSION_3, DIMENSION_M),
     &  U_s_OLD(DIMENSION_3, DIMENSION_M),
     &  V_s_OLD(DIMENSION_3, DIMENSION_M),
     &  W_s_OLD(DIMENSION_3, DIMENSION_M),
     &  X_s_OLD(DIMENSION_3, DIMENSION_M, DIMENSION_N_s),
     &  TIME_OLD
      REAL
     &  XDIST_SC_OLD(DIMENSION_I), XDIST_VEC_OLD(DIMENSION_I), 
     &  YDIST_SC_OLD(DIMENSION_J), YDIST_VEC_OLD(DIMENSION_J),
     &  ZDIST_SC_OLD(DIMENSION_K), ZDIST_VEC_OLD(DIMENSION_K)
      INTEGER
     &  NMAX_OLD(0:DIMENSION_M), IMAX2_OLD, JMAX2_OLD, KMAX2_OLD,
     &  IJMAX2_OLD, MMAX_OLD, FLAG_OLD(DIMENSION_3)
      INTEGER
     &  I_OLD, J_OLD, K_OLD, IV_OLD, JV_OLD, KV_OLD, IJK_OLD,
     &  IM_OLD, JM_OLD, KM_OLD, IP_OLD, JP_OLD, KP_OLD, 
     &  IVJK_OLD, IJVK_OLD, IJKV_OLD, I1SAVE , J1SAVE , K1SAVE,
     &  NSTEP_OLD, L
      LOGICAL EXT_I, EXT_J, EXT_K, DONE, SHIFT

      INTEGER I, J, K, IJK, M, N
C
      INCLUDE 'function.inc'
C
      WRITE(*,*)' Processing data. Please wait. '
      EXT_I = .FALSE.
      EXT_J = .FALSE.
      EXT_K = .FALSE.
      I1SAVE = I1
      J1SAVE = J1
      K1SAVE = K1
C
C  Read old RES file
C
      CALL READ_RES1
C
C  Save old values
C
      DO 100 K = 1, KMAX2
      DO 100 J = 1, JMAX2
      DO 100 I = 1, IMAX2
        IJK = FUNIJK(I, J, K)
        FLAG_OLD(IJK)    = FLAG(IJK)
        EP_g_OLD(IJK)    = EP_g(IJK)
        P_g_OLD(IJK)     = P_g(IJK)
        P_star_OLD(IJK)  = P_star(IJK)
        RO_g_OLD(IJK)    = RO_g(IJK)
        ROP_g_OLD(IJK)   = ROP_g(IJK)
        T_g_OLD(IJK)     = T_g(IJK)
        DO 80 N = 1, NMAX(0)
          X_g_OLD(IJK, N) = X_g(IJK, N)
80      CONTINUE
        U_g_OLD(IJK)     = U_g(IJK)
        V_g_OLD(IJK)     = V_g(IJK)
        W_g_OLD(IJK)     = W_g(IJK)
        DO 95 M = 1, MMAX
          ROP_s_OLD(IJK, M) = ROP_s(IJK, M)
          T_s_OLD(IJK, M)    = T_s(IJK, M)
          U_s_OLD(IJK, M)   = U_s(IJK, M)
          V_s_OLD(IJK, M)   = V_s(IJK, M)
          W_s_OLD(IJK, M)   = W_s(IJK, M)
	  Theta_m_OLD(IJK, M) = Theta_m(IJK, M)
          DO 85 N = 1, NMAX(M)
            X_s_OLD(IJK, M, N) = X_s(IJK, M, N)
85        CONTINUE
95      CONTINUE
100   CONTINUE
      DO 105 I = 1, IMAX2
        XDIST_SC_OLD(I)  = XDIST_SC(I)
        XDIST_VEC_OLD(I) = XDIST_VEC(I)
105   CONTINUE
      DO 110 J = 1, JMAX2
        YDIST_SC_OLD(J)  = YDIST_SC(J)
        YDIST_VEC_OLD(J) = YDIST_VEC(J)
110   CONTINUE
      DO 115 K = 1, KMAX2
        ZDIST_SC_OLD(K)  = ZDIST_SC(K)
        ZDIST_VEC_OLD(K) = ZDIST_VEC(K)
115   CONTINUE
      IMAX2_OLD  = IMAX2
      JMAX2_OLD  = JMAX2
      KMAX2_OLD  = KMAX2
      IJMAX2_OLD = IJMAX2
      MMAX_OLD   = MMAX
      DO 120 M = 0, MMAX
        NMAX_OLD(M) = NMAX(M)
120   CONTINUE
      TIME_OLD = TIME
      NSTEP_OLD = NSTEP
C
C  Read the new data file
C
      CALL INIT_NAMELIST
      CALL READ_NAMELIST(1)
C
C  Do initial calculations
C
      SHIFT = .TRUE.
      CALL CHECK_DATA_03(SHIFT)
      CALL SET_GEOMETRY
      CALL CALC_DISTANCE (XMIN,DX,IMAX2,XDIST_SC,XDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DY,JMAX2,YDIST_SC,YDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DZ,KMAX2,ZDIST_SC,ZDIST_VEC)
C
C  Open new RES files
C
      CLOSE(UNIT_RES)
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/A,$)') ' Enter a new RUN_NAME > '
         READ(*,'(A)')RUN_NAME
      ELSE
         RUN_NAME = TEMP_FILE(1:60)
         DO I = 1,60
            IF (RUN_NAME(I:I).EQ.CHAR(0)) RUN_NAME(I:I) = ' '
         END DO
      END IF
      CALL MAKE_UPPER_CASE(RUN_NAME, 60)
200   DONE = OPEN_FILEP(RUN_NAME, 'NEW', 0)
      IF(.NOT.DONE) THEN
        WRITE(*,'(/A,$)')
     &    ' Unable to open RES file. Enter a new RUN_NAME > '
        READ(*,'(A)')RUN_NAME
        GOTO 200
      ENDIF
C
C  Interpolate values for the new grid
C
      DO 500 K = KMIN1, KMAX1
      DO 500 J = JMIN1, JMAX1
      DO 500 I = IMIN1, IMAX1
        IJK = FUNIJK(I, J, K)
C
C  compute I, J, and K for the old coordinate system
C
        I_OLD = GET_INDEX
     &         (XDIST_SC(I), XDIST_SC_OLD, IMAX2_OLD, EXT_I, I1,'X')
        J_OLD = GET_INDEX
     &         (YDIST_SC(J), YDIST_SC_OLD, JMAX2_OLD, EXT_J, J1,'Y')
        K_OLD = GET_INDEX
     &         (ZDIST_SC(K), ZDIST_SC_OLD, KMAX2_OLD, EXT_K, K1,'Z')
        IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
C
C  If the old IJK location is a wall cell, search the near by cells
C  for a non-wall cell.  Although this interpolation will not be
C  accurate, it is essential for restarting a run, since non-zero values
C  are required for quantities such as void fraction, temperature etc.
C
        IF(FLAG_OLD(IJK_OLD) .GE. 100) THEN
          DO 380 L = 1, 1000
            IM_OLD = MAX((I_OLD - L), 1)
            IJK_OLD  = IM_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
            IP_OLD = MIN((I_OLD + L), IMAX2_OLD)
            IJK_OLD  = IP_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
            JM_OLD = MAX((J_OLD - L), 1)
            IJK_OLD  = I_OLD + (JM_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
            JP_OLD = MIN((J_OLD + L), JMAX2_OLD)
            IJK_OLD  = I_OLD + (JP_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
            KM_OLD = MAX((K_OLD - L), 1)
            IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (KM_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
            KP_OLD = MIN((K_OLD + L), KMAX2_OLD)
            IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (KP_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)GOTO 390
380       CONTINUE
        ENDIF
390     IV_OLD= GET_INDEX
     &         (XDIST_VEC(I), XDIST_VEC_OLD, IMAX2_OLD, EXT_I, I1,'X_E')
        JV_OLD= GET_INDEX
     &         (YDIST_VEC(J), YDIST_VEC_OLD, JMAX2_OLD, EXT_J, J1,'Y_N')
        KV_OLD= GET_INDEX
     &         (ZDIST_VEC(K), ZDIST_VEC_OLD, KMAX2_OLD, EXT_K, K1,'Z_T')
        IVJK_OLD = IV_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
        IJVK_OLD = I_OLD + (JV_OLD - 1) * IMAX2_OLD
     &            + (K_OLD - 1) * IJMAX2_OLD
        IJKV_OLD = I_OLD + (J_OLD - 1) * IMAX2_OLD
     &            + (KV_OLD - 1) * IJMAX2_OLD
C
C  Set the values for the new arrays
C
        EP_g(IJK)    = EP_g_OLD(IJK_OLD)
        P_g(IJK)     = P_g_OLD(IJK_OLD)
        P_star(IJK)  = P_star_OLD(IJK_OLD)
        RO_g(IJK)    = RO_g_OLD(IJK_OLD)
        ROP_g(IJK)   = ROP_g_OLD(IJK_OLD)
        T_g(IJK)     = T_g_OLD(IJK_OLD)
        DO 400 N = 1, NMAX(0)
          IF(N .LE. NMAX_OLD(0))THEN
            X_g(IJK, N) = X_g_OLD(IJK_OLD, N)
          ELSE
            X_g(IJK, N) = ZERO
          ENDIF
400     CONTINUE
        IF(U_g_OLD(IVJK_OLD) .NE. UNDEFINED) THEN
          U_g(IJK)     = U_g_OLD(IVJK_OLD)
        ELSE
          U_g(IJK)     = ZERO
        ENDIF
        IF(V_g_OLD(IJVK_OLD) .NE. UNDEFINED) THEN
          V_g(IJK)     = V_g_OLD(IJVK_OLD)
        ELSE
          V_g(IJK)     = ZERO
        ENDIF
        IF(W_g_OLD(IJKV_OLD) .NE. UNDEFINED) THEN
          W_g(IJK)     = W_g_OLD(IJKV_OLD)
        ELSE
          W_g(IJK)     = ZERO
        ENDIF
        DO 450 M = 1, MMAX
          IF(M .LE. MMAX_OLD) THEN
            ROP_s(IJK, M) = ROP_s_OLD(IJK_OLD, M)
            T_s(IJK, M)    = T_s_OLD(IJK_OLD, M)
            IF(U_s_OLD(IVJK_OLD, M) .NE. UNDEFINED) THEN
              U_s(IJK, M)   = U_s_OLD(IVJK_OLD, M)
            ELSE
              U_s(IJK, M)   = ZERO
            ENDIF
            IF(V_s_OLD(IJVK_OLD, M) .NE. UNDEFINED) THEN
              V_s(IJK, M)   = V_s_OLD(IJVK_OLD, M)
            ELSE
              V_s(IJK, M)   = ZERO
            ENDIF
            IF(W_s_OLD(IJKV_OLD, M) .NE. UNDEFINED) THEN
              W_s(IJK, M)   = W_s_OLD(IJKV_OLD, M)
            ELSE
              W_s(IJK, M)   = ZERO
            ENDIF
            IF(Theta_m_OLD(IJKV_OLD, M) .NE. UNDEFINED) THEN
              Theta_m(IJK, M)   = Theta_m_OLD(IJKV_OLD, M)
            ELSE
              Theta_m(IJK, M)   = ZERO
            ENDIF
            DO 420 N = 1, NMAX(M)
              IF(N .LE. NMAX_OLD(M)) THEN
                X_s(IJK, M, N) = X_s_OLD(IJK_OLD, M, N)
              ELSE
                X_s(IJK, M, N) = ZERO
              ENDIF
420         CONTINUE
          ELSE
            ROP_s(IJK, M)  = ZERO
            U_s(IJK, M)    = ZERO
            V_s(IJK, M)    = ZERO
            W_s(IJK, M)    = ZERO
            Theta_m(IJK, M)   = ZERO
            X_s(IJK, M, 1) = ONE
            DO 440 N = 2, NMAX(M)
              X_s(IJK, M, N) = ZERO
440         CONTINUE
          ENDIF
450     CONTINUE
500   CONTINUE
C
C  Write the new RES file
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/A,$)')
     *         ' Do you need time to be reset to 0.0 (Y/N) ?'
         READ(*,'(A)')RUN_NAME
      ELSE
         RUN_NAME = 'N'
         IF (RESET_TIME) RUN_NAME = 'Y'
      END IF
      IF(RUN_NAME(1:1) .EQ. 'Y' .OR. RUN_NAME(1:1) .EQ. 'y') THEN 
        TIME = ZERO
        NSTEP = 0
      ELSE
        TIME  = TIME_OLD
        NSTEP = NSTEP_OLD
       ENDIF
      CALL WRITE_RES0
      CALL WRITE_RES1
C
      WRITE(*,*)' New RES file written.  Start the new run by setting'
      WRITE(*,*)' RUN_TYPE = RESTART_2'
C
C SET IN CASE CODE IS CHANGED TO RETURN INSTEAD OF STOP
C
      I1 = I1SAVE
      J1 = J1SAVE
      K1 = K1SAVE
C
      STOP
      END
