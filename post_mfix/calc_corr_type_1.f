CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: CALC_CORR_TYPE_1                                       C
C  Purpose: Driver routine for correlation calculations : type # 1     C
C           avg(ep_g),sdv(ep_g),                                       C
C           avg(v_g) , sdv(v_g)                                        C
C           Uses N_SPX                                                 C
C                                                                      C
C  Author: P. Nicoletti                               Date: 02-AUG-92  C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1,NSUM C
C  Variables modified: STARTED, PLOT_TYPE, VAR_INDEX, LOC_X, LOC_Y     C
C                      LOC_Z, PLOT_TYPE2, I, J, K, IJK                 C
C                                                                      C
C  Local variables: TIME_1, TIME_2, FILE1_INDEX,FILE2_INDEX            C
C                   NX, NY, NZ, L, NT, FINISH                          C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CORR_TYPE_1
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'run.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'correl.inc'
      INCLUDE 'xforms.inc'
C
      REAL              TIME_REAL(N_SPX) , TIME_1 , TIME_2, TIME_FOUND
      REAL              TIME_NOW
      INTEGER           NSTEP_1 
      INTEGER           NX , NY , NZ
      INTEGER           REC_POINTER(N_SPX) , L , NT
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX) , FINISH
      LOGICAL           INTER
      INTEGER           I, J, K, IJK
C
      INCLUDE 'function.inc'
C
      INTER = .FALSE.
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,*) ' ENTER TIME_START,TIME_END'
         READ  (*,*) TIME_1 , TIME_2
      ELSE
         TIME_1 = TIME_START
         TIME_2 = TIME_END
      END IF
C
      DO L = 1,N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(3) = .TRUE.    !  V_g
      READ_SPX(1) = .TRUE.    !  EP_g
      CALL SEEK_TIME(READ_SPX, TIME_1, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
      NT = 0
      FINISH = .FALSE.
      STARTED = 0
C
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,
     &                    TIME_NOW, TIME_REAL, NSTEP_1)
      IF (DO_XFORMS) THEN
         CALL CHECK_INTER(INTER)
         IF (INTER) RETURN
      END IF
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_2) GOTO 200
      IF (TIME_NOW .LT. TIME_1) GOTO 100
      NT = NT + 1
      CALL CALC_CORR_01 (FINISH,INTER)
      GOTO 100
C
200   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
C
      FINISH = .TRUE.
      CALL CALC_CORR_01(FINISH,INTER)
      IF (INTER) RETURN
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)')' Average EP_g'
         CALL GET_FILE_NAME(TEMP_FILE4)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE4,STATUS='UNKNOWN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE (40,*)' Average EP_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')
     &   'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')
     &  'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) AVG_EP_g(IJK) 
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)')' Standard Deviation of EP_g'
         CALL GET_FILE_NAME(TEMP_FILE2)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE2,STATUS='UNKNOWN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Standard Deviation of EP_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')
     &   'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')
     &  'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) SDV_EP_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)') ' Average V_g'
         CALL GET_FILE_NAME(TEMP_FILE3)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE3,STATUS='UNKNOWN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Average V_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')
     &   'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')
     &  'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) AVG_V_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)') ' Standard deviation of V_g'
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Standard Deviation of V_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')
     &   'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')
     &  'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) SDV_V_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
C
      RETURN
      END
