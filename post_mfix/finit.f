CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: POST_MFIX                                              C
C  Purpose: main routine for postprocessing MFIX results               C
C                                                                      C
C  Author: P. Nicoletti                               Date: 14-FEB-92  C
C  Reviewer:                                                           C
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
      SUBROUTINE F_INIT
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'xforms.inc'
C
      CHARACTER*30  FILE_NAME
      INTEGER L
      INTEGER REC_POINTER(N_SPX)
      LOGICAL READ_SPX(N_SPX)
      LOGICAL AT_EOF(N_SPX)
      REAL    TIME_REAL(N_SPX), TIME_LAST
      LOGICAL OPEN_FILEP
C
      INTEGER M_PASS,N_PASS
C
C set up machine constants
C
      CALL MACHINE_CONS
C
C get the RUN_NAME from the user
C
      IF (DO_XFORMS) THEN
         L = INDEX(DIR_NAME,'.RES')
         IF (L.EQ.0) GOTO 10
         RUN_NAME = DIR_NAME(1:L-1)
         GOTO 20
      END IF
C
10    WRITE (*,'(A,$)') ' Enter the RUN_NAME to post_process > '
      READ  (*,'(A)') RUN_NAME
      CALL MAKE_UPPER_CASE(RUN_NAME,60)
 20   IF(.NOT. OPEN_FILEP(RUN_NAME,'RESTART_1',N_SPX))GOTO 10
      CALL INIT_NAMELIST
      CALL READ_RES0
      CALL READ_RES1
C
C  Initializations
C
      IF(COORDINATES .EQ. 'CYLINDRICAL') THEN
        CYLINDRICAL = .TRUE.
      ELSE
        CYLINDRICAL = .FALSE.
      ENDIF
C
      IF(IMAX2 .GE. 3) THEN
        DO_I = .TRUE.
        NO_I = .FALSE.
      ELSE
        DO_I = .FALSE.
        NO_I = .TRUE.
      ENDIF
C
      IF(JMAX2 .GE. 3) THEN
        DO_J = .TRUE.
        NO_J = .FALSE.
      ELSE
        DO_J = .FALSE.
        NO_J = .TRUE.
      ENDIF
C
      IF(KMAX2 .GE. 3) THEN
        DO_K = .TRUE.
        NO_K = .FALSE.
      ELSE
        DO_K = .FALSE.
        NO_K = .TRUE.
      ENDIF
C
      CALL SET_MAX2
      CALL SET_GEOMETRY
      CALL SET_CONSTANTS
      CALL SET_INCREMENTS
      TIME_LAST = TIME
C
C READ INITIAL RECORDS OF THE .SPx files ...
C
      DO L = 1,N_SPX
         READ_SPX(L)    = .TRUE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      CALL READ_SPX0(READ_SPX)
      DO L = 1,N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
C
C CALCULATE DISTANCES AND VOLUMES
C
      CALL CALC_DISTANCE (XMIN,DX,IMAX2,XDIST_SC,XDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DY,JMAX2,YDIST_SC,YDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DZ,KMAX2,ZDIST_SC,ZDIST_VEC)
      CALL CALC_VOL
C
      IF (DO_XFORMS) RETURN
50    CALL HEADER_MAIN
      IF (SELECTION.EQ.0) THEN
        STOP
      ELSEIF (SELECTION.EQ.1) THEN
        CALL EXAMINE_DATA
      ELSEIF (SELECTION.EQ.2) THEN
        CALL RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
      ELSEIF (SELECTION.EQ.3) THEN
        CALL INTERP_RES
      ELSEIF (SELECTION.EQ.4) THEN
        CALL CALC_QUANTITIES
      ELSEIF (SELECTION.EQ.5) THEN
        CALL PRINT_OUT
      ELSEIF (SELECTION.EQ.6) THEN
        CALL USR_POST
      ELSEIF (SELECTION.EQ.7) THEN
        CALL SELECT_SPX_REC
      ENDIF
      GOTO 50
      END
C
      SUBROUTINE CLOSE_OLD_RUN()
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'funits.inc'
C
      INTEGER L
C
      CLOSE (UNIT=UNIT_RES)
      DO L = 1,N_SPX
         CLOSE (UNIT=UNIT_SPX+L)
      END DO
C
      RETURN
      END
