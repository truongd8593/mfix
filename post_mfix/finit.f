!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: POST_MFIX                                              C
!  Purpose: main routine for postprocessing MFIX results               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 14-FEB-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE F_INIT
!
      Use param
      Use param1
      Use run
      Use post3d
      Use geometry
      Use indices
      Use fldvar
      Use physprop
      Use constant
      Use funits
      Use parallel_mpi
      Use gridmap
!
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      CHARACTER*30  FILE_NAME
      INTEGER L
      INTEGER REC_POINTER(N_SPX)
      LOGICAL READ_SPX(N_SPX)
      LOGICAL AT_EOF(N_SPX)
      REAL    TIME_REAL(N_SPX), TIME_LAST
      LOGICAL OPEN_FILEP
!
      INTEGER M_PASS,N_PASS

      nodesi = 1
      nodesj = 1
      nodesk = 1
      call parallel_init()
!// Partition the domain and set indices
!
!
! set up machine constants
!
      CALL MACHINE_CONS
!
! get the RUN_NAME from the user
!
      IF (DO_XFORMS) THEN
         L = INDEX(DIR_NAME,'.RES')
         IF (L.EQ.0) GOTO 10
         RUN_NAME = DIR_NAME(1:L-1)
         GOTO 20
      END IF
!
10    WRITE (*,'(A,$)') ' Enter the RUN_NAME to post_process > '
      READ  (*,'(A)') RUN_NAME
      CALL MAKE_UPPER_CASE(RUN_NAME,60)
 20   IF(.NOT. OPEN_FILEP(RUN_NAME,'RESTART_1',N_SPX))GOTO 10
      CALL INIT_NAMELIST
      nodesi = 1
      nodesj = 1
      nodesk = 1


      imin2 = 1
      jmin2 = 1
      kmin2 = 1
      CALL READ_RES0
      if(no_k .and. dz(1) .eq. UNDEFINED)DZ(1)=ONE
      call SET_MAX2
      call GRIDMAP_INIT
      CALL SET_INCREMENTS

      CALL READ_RES1
!
!  Initializations
!
      IF(COORDINATES .EQ. 'CYLINDRICAL') THEN
        CYLINDRICAL = .TRUE.
      ELSE
        CYLINDRICAL = .FALSE.
      ENDIF
!
      IF(IMAX2 .GE. 3) THEN
        DO_I = .TRUE.
        NO_I = .FALSE.
      ELSE
        DO_I = .FALSE.
        NO_I = .TRUE.
      ENDIF
!
      IF(JMAX2 .GE. 3) THEN
        DO_J = .TRUE.
        NO_J = .FALSE.
      ELSE
        DO_J = .FALSE.
        NO_J = .TRUE.
      ENDIF
!
      IF(KMAX2 .GE. 3) THEN
        DO_K = .TRUE.
        NO_K = .FALSE.
      ELSE
        DO_K = .FALSE.
        NO_K = .TRUE.
      ENDIF
!
      CALL SET_MAX2
      CALL SET_GEOMETRY
      CALL SET_CONSTANTS
      CALL SET_INCREMENTS
      TIME_LAST = TIME
!
! READ INITIAL RECORDS OF THE .SPx files ...
!
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
!
! CALCULATE DISTANCES AND VOLUMES
!
      CALL CALC_DISTANCE (XMIN,DX,IMAX2,XDIST_SC,XDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DY,JMAX2,YDIST_SC,YDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DZ,KMAX2,ZDIST_SC,ZDIST_VEC)
      CALL CALC_VOL
!
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
!
      SUBROUTINE CLOSE_OLD_RUN()
!
      Use param
      Use param1
      Use funits
!
      IMPLICIT NONE
!
      INTEGER L
!
      CLOSE (UNIT=UNIT_RES)
      DO L = 1,N_SPX
         CLOSE (UNIT=UNIT_SPX+L)
      END DO
!
      RETURN
      END
