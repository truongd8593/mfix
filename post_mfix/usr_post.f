CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: USR_POST                                               C
C                                                                      C
C  Purpose: Do user defined calculations from TIME_START to TIME_LAST  C
C                                                                      C
C  Author: M. Syamlal                               Date: 26-OCT-93    C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: V_g, EP_g, IMIN1, IMAX1, JMIN1, JMAX1         C
C                        KMIN1, KMAX1, IJKMAX2                         C
C  Variables modified: PLOT_TYPE, VAR_INDEX, LOC_X, LOC_Y, LOC_Z       C
C                      I, J, K, IJK                                    C
C                                                                      C
C  Local variables: FILE_NAME, NX, NY, NZ    L, NT                     C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR_POST
C
      IMPLICIT NONE
      INTEGER MAX_COUNT
      PARAMETER (MAX_COUNT=1000)
C
      REAL a1, a2, a3, Re_c, EP_c
c      PARAMETER (a1 = 250.)
      PARAMETER (a1 = 1500.)
      PARAMETER (a2 = 0.005)
      PARAMETER (a3 = 90.0)
      PARAMETER (Re_c = 5.)
      PARAMETER (EP_c = 0.92)
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'run.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'physprop.inc'
C
      REAL              TIME_START , TIME_REAL(N_SPX), TIME_FOUND
      REAL              TIME_LAST, TIME_NOW, SUM, FC_DIST, Re, VREL, FC
      CHARACTER         FILE_NAME*60
      INTEGER           NX , NY , NZ, NSTEP_1
      INTEGER           REC_POINTER(N_SPX) , L , NT
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
C
      INTEGER           FC_COUNT(MAX_COUNT)
      INTEGER           IJK
C
      INCLUDE 'function.inc'
C
      WRITE (*,*) ' Enter start-time and end-time'
      READ  (*,*) TIME_START, TIME_LAST
C
      MU_g0 = 1.8e-4
      RO_g0 = 1.8e-3
      DO 10 L = 1, MAX_COUNT
        FC_COUNT(L) = 0
10    CONTINUE
C
      CALL GET_FILE_NAME(FILE_NAME)
      OPEN (UNIT=40,FILE=FILE_NAME,STATUS='UNKNOWN')
C
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
C
C   Enable the required files
C
      READ_SPX(1) = .TRUE.    ! EP_g
c      READ_SPX(2) = .TRUE.    ! P_g, P_star
      READ_SPX(3) = .TRUE.    ! U_g, V_g, W_g
      READ_SPX(4) = .TRUE.    ! U_s, V_s, W_s
c      READ_SPX(5) = .TRUE.    ! ROP_s
c      READ_SPX(6) = .TRUE.    ! T_g, T_s1, T_s2
c      READ_SPX(7) = .TRUE.    ! X_g, X_s
      CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
C
      NT = 0
C
C  Time loop -- Start
C
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,
     &                    TIME_NOW, TIME_REAL, NSTEP_1)
C
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_LAST) GOTO 500
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      NT = NT + 1
C
C  DO required computations
C
      DO 200 IJK = 1,IJKMAX2
        IF(FLUID_AT(IJK)) THEN
          VREL = SQRT(  (U_g(IJK) - U_s(IJK,1))**2
     &                + (V_g(IJK) - V_s(IJK,1))**2
     &                + (W_g(IJK) - W_s(IJK,1))**2 )
          Re =  D_p(1) * VREL * RO_g0 / MU_g0
          FC = (ONE + a1 * exp(-a2*(Re - Re_c)**2 
     &                          - a3*(EP_g(IJK)-ep_c)**2)
     &               * Re * (1. - EP_g(IJK))               )

          L = MIN(MAX_COUNT, NINT(FC))
          FC_COUNT(L) = FC_COUNT(L) + 1
        ENDIF
200   CONTINUE
C
      GOTO 100
C
C  Time loop -- End
C
500   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
C
C  Do final data processing printing
C
      SUM = ZERO
      DO 600 L = 1, MAX_COUNT
        SUM = SUM + FC_COUNT(L)
600   CONTINUE
C
      DO 650 L = 1, MAX_COUNT
        FC_DIST = REAL(FC_COUNT(L))/SUM
        WRITE(40,*)L, FC_DIST
650   CONTINUE
C
      CLOSE (UNIT=40)
C
      RETURN
      END
