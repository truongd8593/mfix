CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SELECT_SPX_REC                                         C
C  Purpose: Select records from SPx files interactively and write to a C
C           new SPx file                                               C
C           NOTE : USES N_SPX                                          C
C                                                                      C
C  Author: M. Syamlal                                 Date: 6-APR-94   C
C  Reviewer:                                                           C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IJKMAX2                                       C
C  Variables modified: IJK                                             C
C                                                                      C
C  Local variables:                                                    C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SELECT_SPX_REC
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'run.inc'
      INCLUDE 'machine.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'xforms.inc'
C
      REAL     TIME_FOR_RES, TIME_FOUND
      LOGICAL  AT_EOF(N_SPX), READ_SPX(N_SPX),SELECT
      INTEGER  REC_POINTER(N_SPX), REC_POINTER_t(N_SPX)
      INTEGER  NSTEP_1 , ERROR_CODE
      CHARACTER*1  IANS 
      CHARACTER*13 LINE
      REAL     TIME_REAL(N_SPX)
      LOGICAL  ERROR
C
      INTEGER L, L_SPX , LL
C
      ERROR  = .FALSE.
      SELECT = .TRUE.
      L_SPX  = SPX_NUM
C
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A,$)') 'Enter the number of Spx file > '
         READ  (*,*) L_SPX
10       WRITE (*,'(A,$)') 'Enter the name of new file > '
         READ  (*,'(A60)') TEMP_FILE
         OPEN (UNIT=UNIT_SPX+L_SPX,FILE=TEMP_FILE,STATUS='NEW',
     &      RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=10)
      END IF
C
c      IF (DO_XFORMS.AND.GET_TIMES) THEN
c         OPEN (UNIT=UNIT_SPX+L_SPX,FILE=TEMP_FILE,STATUS='UNKNOWN',
c     &      RECL=OPEN_N1,ACCESS='DIRECT',FORM='UNFORMATTED',ERR=10)
c      END IF
C
      CALL WRITE_SPX0(L_SPX)
C
      DO 20 L = 1, N_SPX
        READ_SPX(L) = .FALSE.
        REC_POINTER(L) = 4
        AT_EOF(L) = .FALSE.
20    CONTINUE
      READ_SPX(L_SPX) = .TRUE.
C
      L = 0
100   CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL,NSTEP_1)
      IF (.NOT.AT_EOF(L_SPX)) THEN
         L = L + 1
         IF (.NOT.DO_XFORMS) THEN
            WRITE(*,'(A,G12.5,A,$)')
     &           'Write time ', TIME_REAL(L_SPX), '(Y/N) [Y]'
            READ(*,'(A1)')IANS
            IF (IANS .NE. 'n' .AND. IANS .NE. 'N') THEN
               TIME = DBLE(TIME_REAL(L_SPX))
               CALL WRITE_SPX1(L_SPX)
            END IF
         ENDIF
         IF (DO_XFORMS) THEN
            IF (GET_TIMES) THEN
               WRITE (LINE,'(G12.5)') TIME_REAL(L_SPX)
               LINE(13:13) = CHAR(0)
               IF (EVERY_OTHER_TIME) THEN
                  SELECT = .TRUE.
                  IF (MOD(L,2).EQ.0) SELECT = .FALSE.
               END IF
               CALL ADD_TO_SPX_BR(SELECT,LINE)
            ELSE
               CALL SPX_TIME_SELECTED(L,LL)
               IF (LL.NE.0) THEN
                  TIME = DBLE(TIME_REAL(L_SPX))
                  CALL WRITE_SPX1(L_SPX)
                  CALL SPX_DESELECT_TIME(L)
               END IF
            END IF
         END IF
         GOTO 100
      ENDIF
      IF (.NOT.GET_TIMES) CLOSE(UNIT_SPX+L_SPX)
      WRITE(*,*)
      WRITE(*,'(A,I1,A)')' A New SP',L_SPX, ' file written '
      RETURN
      END
