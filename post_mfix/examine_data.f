CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: EXAMINE_DATA                                           C
C  Purpose: Examine/print selected data                                C
C                                                                      C
C  Author: M. Syamlal                                 Date: 03-NOV-93  C
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
      SUBROUTINE EXAMINE_DATA
C
      IMPLICIT NONE
      INTEGER  N_VAR
      PARAMETER (N_VAR=48)
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'run.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'post3d.inc'
      INCLUDE 'xforms.inc'
C
      CHARACTER*80 LINE
      CHARACTER*120 STRING, SUBSTR
      CHARACTER*8  VAR, VAR_DAT(N_VAR)
      CHARACTER*120    FILE_NAME
      INTEGER      L, L3, L4, LMAX, IANS, NSTEP_1
      REAL         DX_E, DY_N, DZ_T
      REAL         DIST(DIMENSION_3), VALUE(DIMENSION_3)
      REAL         TIME_IN_RES
      INTEGER      DISPLAY, DIRECTION, NT
      INTEGER      I2d, J2d, K2d
      LOGICAL      FILE_EXIST,END_AVERAGE
      LOGICAL      SUM
      LOGICAL      STRCMP,INTER
      REAL         XTMP, YTMP, ZTMP, VALUE_TMP
      INTEGER      IJK1
      INTEGER      M_LOCAL, mIJK, lIJK, IER
      INTEGER      I, J, K, IJK, M, N   
      REAL         DELm, DELl, FAC1, FAC2 
C
      REAL              TIME_REAL(N_SPX), TIME_FOUND, TIME_NOW, TIME_OLD
      INTEGER           REC_POINTER(N_SPX)
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
C
C  Function subroutines
C
      REAL XFLOW_gx, XFLOW_gy, XFLOW_gz
      REAL VFLOW_gx, VFLOW_gy, VFLOW_gz, MFLOW_gx, MFLOW_gy, MFLOW_gz
      REAL XFLOW_sx, XFLOW_sy, XFLOW_sz
      REAL VFLOW_sx, VFLOW_sy, VFLOW_sz, MFLOW_sx, MFLOW_sy, MFLOW_sz
      REAL FLUX_gx, FLUX_gy, FLUX_gz
      REAL FLUX_sx, FLUX_sy, FLUX_sz
      REAL CALC_RO_g
C
      INCLUDE 'function.inc'
C                   1       2      3         4      5      6
      DATA VAR_DAT/'EP_g', 'P_g', 'P_star', 'U_g', 'V_g', 'W_g',
      
C                   7      8      9      10       11     12 
     &             'U_s', 'V_s', 'W_s', 'ROP_s', 'T_g', 'T_s',
     
C                   13      14     15     16          17
     &             'T_s2', 'X_g', 'X_s', 'XFLOW_gx', 'XFLOW_gy', 
     
C                   18          19          20          21     
     &             'XFLOW_gz', 'XFLOW_sx', 'XFLOW_sy', 'XFLOW_sz',
     
C                   22          23          24          25     
     &             'MFLOW_gx', 'MFLOW_gy', 'MFLOW_gz', 'MFLOW_sx', 
     
C                   26          27         28           29     
     &             'MFLOW_sy', 'MFLOW_sz', 'VFLOW_gx', 'VFLOW_gy', 
     
C                   30          31          32          33
     &             'VFLOW_gz', 'VFLOW_sx', 'VFLOW_sy', 'VFLOW_sz',
     
C                   34        35        36         37     
     &             'MASS_g', 'MASS_s', 'FLUX_gx', 'FLUX_gy', 
     
C                   38         39         40         41
     &             'FLUX_gz', 'FLUX_sx', 'FLUX_sy', 'FLUX_sz' ,
     
C                   42      43     44    45     46     47
     &             'KE_g', 'KE_s','P_s','PE_g','PE_s','BERN_s', 
     
C                   48
     &   	   'Theta_m' /


      DZ_T(K) = HALF * (DZ(K) + DZ(Kp1(K)))
      DY_N(J) = HALF * (DY(J) + DY(Jp1(J)))
      DX_E(I) = HALF * (DX(I) + DX(Ip1(I)))
C
      CALL READ_RES1
      TIME_IN_RES = TIME
C
      SUM  = .FALSE.
      LMAX = LEN(STRING)
      INTER = .FALSE.
C
      IF (.NOT.DO_XFORMS) THEN
         TIME_START = 0.
         TIME_END   = 0.
         TIME_AVERAGE = .FALSE.
         I1         = 1
         I2         = 1
         I_AVERAGE  = .FALSE.
         J1         = 1
         J2         = 1
         J_AVERAGE  = .FALSE.
         K1         = 1
         K2         = 1
         K_AVERAGE  = .FALSE.
         M = 1
         N = 1
         VAR     = 'EP_g'
      ELSE
         VAR     = VAR_DAT(VAR_NO)
      END IF
C
C
      FILE_NAME  = '*'
      IF (.NOT.DO_XFORMS) THEN
         RES_P_g    = .FALSE.
         RES_T_g    = .FALSE.
         RES_X_g    = .FALSE.
         MINMAX     = -1
      END IF
      M          = M_USE
      N          = N_USE
      DO L = 1, N_SPX
         REC_POINTER(L) = 4
      END DO
C
C  
      IF (DO_XFORMS) GOTO 5500
C
C
      WRITE(*,*)
      WRITE(*,*)
     &' Interactive data retrieval program. Type ? any time for help,'
      WRITE(*,*)
     &' or press RETURN to select default values shown in parenthesis.'
C
C  Read time
C
10    WRITE(*,*)
      WRITE(*,'(A,F7.3,A,F7.3,A,$)')
     & ' Time: (',TIME_START,',',TIME_END,') > '
      READ(*,'(1A60)',ERR=10) STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(10)
        GOTO 10
      ELSEIF(STRING(1:1) .EQ. 'e' .OR. STRING(1:1) .EQ. 'E' .OR.
     &       STRING(1:1) .EQ. 'q' .OR. STRING(1:1) .EQ. 'Q') THEN
        IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
        RETURN
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)TIME_START
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=10)TIME_END
      IF(TIME_START .LT. ZERO .OR. TIME_END .LT. ZERO) THEN
        IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
        RETURN
      ENDIF
      IF(TIME_START .GE. TIME_IN_RES) THEN
        TIME_START = TIME_IN_RES
        TIME_END   = TIME_IN_RES
      ENDIF
      IF(TIME_END .LT. TIME_START)GOTO 10
      IF(TIME_START .NE. TIME_END) THEN
        IF(TIME_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
11      WRITE(*, '(A,1A1,A,$)')' Time average ? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=11) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(11)
          GOTO 11
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          TIME_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          TIME_AVERAGE = .FALSE.
        ENDIF
      ENDIF
C
C  Read variable name
C
20    CONTINUE
      IF(MINMAX .EQ. 1) THEN
        SUBSTR(1:1) = '1'
        SUBSTR(2:9) = VAR
      ELSEIF(MINMAX .EQ. 0) THEN
        SUBSTR(1:1) = '0'
        SUBSTR(2:9) = VAR
      ELSE
        SUBSTR(1:8) = VAR
        SUBSTR(9:9) = ' '
      ENDIF
      WRITE(*,'(A,1A9,A,$)')
     & ' Variable: (', SUBSTR(1:9), ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(20)
        GOTO 20
      ENDIF
      L3 = 1
C
      IF(STRING(1:1) .EQ. '1') THEN
        MINMAX = 1
        L3 =2
      ELSEIF(STRING(1:1) .EQ. '0')THEN
        MINMAX = 0
        L3 = 2
      ELSEIF(STRING(1:1) .NE. ' ')THEN
        MINMAX = -1
      ENDIF
C
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,'(1A8)',ERR=20) VAR
C
C     Identify variable number
C
      DO 22 L = 1, N_VAR
        IF (STRCMP(VAR,VAR_DAT(L)))THEN
          VAR_NO = L
          GOTO 23
        ENDIF
22    CONTINUE
      WRITE(*,'(A,1A8,A)')' Variable ', VAR, ' not found'
      GOTO 20
23    CONTINUE
C
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE. 10) .OR.
     &   (VAR_NO .EQ. 15                     ) .OR.
     &   (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.
     &   (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR.
     &   (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.
     &   (VAR_NO .EQ. 35                     ) .OR.
     &   (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.
     &   (VAR_NO .EQ. 43                     ) .OR.
     &   (VAR_NO .EQ. 44                     ) .OR.
     &   (VAR_NO .EQ. 46                     ) .OR.
     &   (VAR_NO .EQ. 47                     ) .OR.
     &   (VAR_NO .EQ. 48                     ) 
     &                                             )THEN
        IF(MMAX .GT. 1) THEN
24        WRITE(*,'(A,I2,A,$)') ' Solids phase: (', M, ') > '
          READ(*,'(1A60)',ERR=24) STRING
          IF(STRING(1:1) .EQ. '?') THEN
            CALL HELP(24)
            GOTO 24
          ENDIF
          L3 = 1
          CALL GET_SUBSTR(STRING, L3, SUBSTR)
          IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=24)M
          IF(M .GT. MMAX) THEN
            WRITE(*,*)' Value should not exceed ', MMAX
            M = MMAX
            GOTO 24
          ENDIF
        ELSE
          M = 1
        ENDIF
      ENDIF
C
      IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21)THEN
25      WRITE(*,'(A,I2,A,$)') ' Species: (', N, ') > '
        READ(*,'(1A60)',ERR=25) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(25)
          GOTO 25
        ENDIF
        L3 = 1
        CALL GET_SUBSTR(STRING, L3, SUBSTR)
        IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=25)N
        IF((VAR_NO .EQ. 14 .OR. VAR_NO .EQ. 16 .OR. VAR_NO .EQ. 17 .OR.
     &     VAR_NO .EQ. 18 ) .AND. N .GT. NMAX(0)) THEN
          WRITE(*,*)' Value should not exceed ', NMAX(0)
          N = NMAX(0)
          GOTO 25
        ELSEIF((VAR_NO .EQ. 15 .OR. VAR_NO .EQ. 19 .OR. VAR_NO .EQ. 20
     &          .OR. VAR_NO .EQ. 21) .AND. N .GT. NMAX(M)) THEN
          WRITE(*,*)' Value should not exceed ', NMAX(M)
          N = NMAX(M)
          GOTO 25
        ENDIF
      ENDIF
C
C
 5500 CONTINUE
C
C
      IF(VAR_NO .EQ.  4 .OR. VAR_NO .EQ.  7 .OR. VAR_NO .EQ. 16 .OR.
     &   VAR_NO .EQ. 19 .OR. VAR_NO .EQ. 22 .OR. VAR_NO .EQ. 25 .OR.
     &   VAR_NO .EQ. 28 .OR. VAR_NO .EQ. 31 .OR. VAR_NO .EQ. 36 .OR.
     &   VAR_NO .EQ. 39)THEN
        DIRECTION = 1
      ELSEIF
     &  (VAR_NO .EQ.  5 .OR. VAR_NO .EQ.  8 .OR. VAR_NO .EQ. 17 .OR.
     &   VAR_NO .EQ. 20 .OR. VAR_NO .EQ. 23 .OR. VAR_NO .EQ. 26 .OR.
     &   VAR_NO .EQ. 29 .OR. VAR_NO .EQ. 32 .OR. VAR_NO .EQ. 37 .OR.
     &   VAR_NO .EQ. 40)THEN
        DIRECTION = 2
      ELSEIF
     &  (VAR_NO .EQ.  6 .OR. VAR_NO .EQ.  9 .OR. VAR_NO .EQ. 18 .OR.
     &   VAR_NO .EQ. 21 .OR. VAR_NO .EQ. 24 .OR. VAR_NO .EQ. 27 .OR.
     &   VAR_NO .EQ. 30 .OR. VAR_NO .EQ. 33 .OR. VAR_NO .EQ. 38 .OR.
     &   VAR_NO .EQ. 41)THEN
        DIRECTION = 3
      ELSE
c        DIRECTION = 0
      ENDIF
C
      IF(VAR_NO .GE. 16 .AND. VAR_NO .LE. 41) THEN
        SUM = .TRUE.
      ELSE
        SUM = .FALSE.
      ENDIF
C
C  Enable the required SPX file
C
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         AT_EOF(L)      = .FALSE.
      END DO
C
C
      IF(VAR_NO .EQ. 1 .OR. (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.
     &   (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR. 
     &   (VAR_NO .GE. 28 .AND. VAR_NO .LE. 30) .OR.
     &    VAR_NO .EQ. 34                       .OR.
     &   (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.
     &   (VAR_NO .GE. 42 .AND. VAR_NO .LE. 47)    )            THEN
        READ_SPX(1) = .TRUE.    ! EP_g
      ENDIF
      IF(VAR_NO .EQ. 2 .OR. VAR_NO .EQ. 3  .OR.
     &   VAR_NO .EQ. 44) THEN
        READ_SPX(2) = .TRUE.    ! P_g, P_star
      ENDIF
      IF((VAR_NO .GE.  4 .AND. VAR_NO .LE.  6) .OR.
     &   (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.
     &   (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR. 
     &   (VAR_NO .GE. 28 .AND. VAR_NO .LE. 30) .OR.
     &   (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.
     &   (VAR_NO .EQ. 42                     ) 
     &                                           ) THEN  
        READ_SPX(3) = .TRUE.    ! U_g, V_g, W_g
      ENDIF
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE.  9) .OR.
     &   (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.
     &   (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. 
     &   (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.
     &   (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.
     &   (VAR_NO .EQ. 43                     ) .OR.
     &   (VAR_NO .EQ. 44                     ) .OR.
     &   (VAR_NO .EQ. 47                     ) 
     &                                             ) THEN
        READ_SPX(4) = .TRUE.    ! U_s, V_s, W_s
      ENDIF
      IF (VAR_NO .EQ. 44) THEN
        READ_SPX(5) = .TRUE.
      END IF
      IF (VAR_NO .EQ. 10 .OR.
     &   (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.
     &   (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. 
     &   (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.
     &    VAR_NO .EQ. 35                       .OR.
     &   (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41)    ) THEN
        IF(MMAX .EQ. 1) THEN
          READ_SPX(1) = .TRUE.    ! EP_g
        ELSE
          READ_SPX(5) = .TRUE.    ! ROP_s
        ENDIF
      ENDIF
      IF(VAR_NO .GE. 11 .AND. VAR_NO .LE. 13) THEN
        READ_SPX(6) = .TRUE.    ! T_g, T_s, T_s2
      ENDIF
      IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21) THEN
        READ_SPX(7) = .TRUE.    ! X_g, X_s
      ENDIF
      
      IF(VAR_NO .EQ. 48 ) THEN
        READ_SPX(8) = .TRUE.    ! Theta_m
      ENDIF
C
C  Open P_g, T_g, and X_g files, if gas density needs to be determined
C
      IF (DO_XFORMS) GOTO 1125
      IF( RO_g0 .EQ. UNDEFINED .AND. 
     &   ( (VAR_NO .GE. 16 .AND. VAR_NO .LE. 18) .OR.
     &     (VAR_NO .GE. 22 .AND. VAR_NO .LE. 24) .OR.
     &     (VAR_NO .EQ. 34                     ) .OR.
     &     (VAR_NO .GE. 36 .AND. VAR_NO .LE. 38) .OR.
     &     (VAR_NO .EQ. 42                     ) .OR.
     &     (VAR_NO .EQ. 45                     ) 
     &                                              ) ) THEN
          IF (.NOT.DO_XFORMS) THEN
             WRITE(*,*)
     &     ' To calculate gas density P_g, T_g, and X_g are needed'
C
           IF(RES_P_g)THEN
             SUBSTR(1:1) = 'Y'
           ELSE
             SUBSTR(1:1) = 'N'
           ENDIF
26         WRITE(*, '(A,1A1,A,$)')
     &      ' P_g from RES file? (',SUBSTR(1:1),') > '
           READ(*,'(1A60)',ERR=26) STRING
           IF(STRING(1:1) .EQ. '?') THEN
             CALL HELP(26)
             GOTO 26
           ENDIF
           IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
             RES_P_g = .TRUE.
           ELSEIF(STRING(1:1) .NE. ' ')THEN
             RES_P_g = .FALSE.
           ENDIF
C
           IF(RES_T_g)THEN
             SUBSTR(1:1) = 'Y'
           ELSE
             SUBSTR(1:1) = 'N'
           ENDIF
27         WRITE(*, '(A,1A1,A,$)')
     &       ' T_g from RES file? (',SUBSTR(1:1),') > '
           READ(*,'(1A60)',ERR=27) STRING
           IF(STRING(1:1) .EQ. '?') THEN
             CALL HELP(27)
             GOTO 27
           ENDIF
           IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
             RES_T_g = .TRUE.
           ELSEIF(STRING(1:1) .NE. ' ')THEN
             RES_T_g = .FALSE.
           ENDIF
        ELSE
           CONTINUE
        END IF
C
        RES_X_G = .FALSE.
        IF(MW_avg .EQ. UNDEFINED) THEN
          IF (.NOT.DO_XFORMS) THEN
             IF(RES_X_g)THEN
               SUBSTR(1:1) = 'Y'
             ELSE
               SUBSTR(1:1) = 'N'
             ENDIF
28           WRITE(*, '(A,1A1,A,$)')
     &         ' X_g from RES file? (',SUBSTR(1:1),') > '
             READ(*,'(1A60)',ERR=28) STRING
             IF(STRING(1:1) .EQ. '?') THEN
               CALL HELP(28)
               GOTO 28
             ENDIF
             IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
               RES_X_g = .TRUE.
             ELSEIF(STRING(1:1) .NE. ' ')THEN
               RES_X_g = .FALSE.
             ENDIF
          ELSE
             CONTINUE
          END IF
        ENDIF
C
C
        IF(.NOT.RES_P_g) READ_SPX(2) = .TRUE.    ! P_g, P_star
        IF(.NOT.RES_X_g) READ_SPX(7) = .TRUE.    ! X_g, X_s
        IF(.NOT.RES_T_g) READ_SPX(6) = .TRUE.    ! T_g, T_s, T_s2
        IF(RES_P_g .OR. RES_T_g .OR. RES_X_g) CALL READ_RES1
C
C
      ENDIF
C
1125  IF (DO_XFORMS) GOTO 5501
C
C
C  Read I range
C
30    WRITE(*,'(A,I3,A,I3,A,$)')
     & ' I range: (', I1, ',', I2, ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(30)
        GOTO 30
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) I1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=30) I2
C
C     Check bounds
C
      IF(I2 .LT. I1) I2 = I1
      IF(I1 .LT. 1 .OR. I1 .GT. IMAX2 .OR.
     &   I2 .LT. 1 .OR. I2 .GT. IMAX2     ) THEN
        WRITE(*,'(A,I3)')
     &    ' I1 and I2 should be in the range 1 to ', IMAX2
        GOTO 30
      ENDIF
C
      IF(I1 .NE. I2) THEN
        IF(I_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
31      WRITE(*, '(A,1A1,A,$)')
     &  ' Average or sum over I? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=31) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(31)
          GOTO 31
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          I_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          I_AVERAGE = .FALSE.
        ENDIF
      ENDIF
C
C  Read J range
C
40    WRITE(*,'(A,I3,A,I3,A,$)')
     & ' J range: (', J1,',',J2, ') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(40)
        GOTO 40
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=40) J1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=40) J2
C
C     Check bounds
C
      IF(J2 .LT. J1) J2 = J1
      IF(J1 .LT. 1 .OR. J1 .GT. JMAX2 .OR.
     &   J2 .LT. 1 .OR. J2 .GT. JMAX2     ) THEN
        WRITE(*,'(A,I3)')
     &    ' J1 and J2 should be in the range 1 to ', JMAX2
        GOTO 40
      ENDIF
C
      IF(J1 .NE. J2) THEN
        IF(J_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
41      WRITE(*, '(A,1A1,A,$)')
     &  ' Average or sum over J? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=41) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(41)
          GOTO 41
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          J_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          J_AVERAGE = .FALSE.
        ENDIF
      ENDIF
C
C  Read K range
C
50    WRITE(*,'(A,I3,A,I3,A,$)')
     & ' K range: (', K1,',',K2,') >'
      READ(*,'(1A60)') STRING
      IF(STRING(1:1) .EQ. '?') THEN
        CALL HELP(50)
        GOTO 50
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=50) K1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF(SUBSTR(1:1) .NE. ' ')READ(SUBSTR,*,ERR=50) K2
C
C     Check bounds
C
      IF(K2 .LT. K1) K2 = K1
      IF(K1 .LT. 1 .OR. K1 .GT. KMAX2 .OR.
     &   K2 .LT. 1 .OR. K2 .GT. KMAX2     ) THEN
        WRITE(*,'(A,I3)')
     &    ' K1 and K2 should be in the range 1 to ', KMAX2
        GOTO 50
      ENDIF
C
      IF(K1 .NE. K2) THEN
        IF(K_AVERAGE)THEN
          SUBSTR(1:1) = 'Y'
        ELSE
          SUBSTR(1:1) = 'N'
        ENDIF
51      WRITE(*, '(A,1A1,A,$)')
     &  ' Average or sum over K? (',SUBSTR(1:1),') > '
        READ(*,'(1A60)',ERR=51) STRING
        IF(STRING(1:1) .EQ. '?') THEN
          CALL HELP(51)
          GOTO 51
        ENDIF
        IF(STRING(1:1) .EQ. 'Y' .OR. STRING(1:1) .EQ. 'y')THEN
          K_AVERAGE = .TRUE.
        ELSEIF(STRING(1:1) .NE. ' ')THEN
          K_AVERAGE = .FALSE.
        ENDIF
      ENDIF
C
C
 5501 CONTINUE
C
C
C  Read file name
C
c      IF (DO_XFORMS) THEN
c         L3 = INDEX(TEMP_FILE,'*')
c         IF (L3.EQ.0) THEN
c            IF (APPEND_MODE) THEN
c               OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN',
c     &                       ACCESS='APPEND')
c            ELSE
c               OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN')
c            END IF
c            FILE_NAME(1:1) = 'A'
c         ELSE
c            FILE_NAME(1:1) = '*'
c         END IF
c         GOTO 5502
c      END IF
C
C
70    WRITE(*,'(A,1A30,A,$)') ' File: (', FILE_NAME,') >'
      READ(*,'(1A60)') STRING
      IF (STRING(1:1) .EQ. '?') THEN
         CALL HELP(70)
         GOTO 70
      ELSE IF (STRING(1:1) .EQ. '!') THEN
         GOTO 10
      ENDIF
      L3 = 1
      CALL GET_SUBSTR(STRING, L3, SUBSTR)
      IF (SUBSTR(1:1) .NE. ' ')THEN
         IF (.NOT.STRCMP(STRING(1:30),FILE_NAME)) THEN
            IF (FILE_NAME(1:1) .NE. '*') CLOSE(40)
            CALL STREQS(FILE_NAME,STRING(1:30))
            IF (FILE_NAME(1:1) .NE. '*') THEN
               INQUIRE (FILE=FILE_NAME,EXIST=FILE_EXIST)
               IF (FILE_EXIST .AND. .NOT.DO_XFORMS) THEN
                  WRITE(*,'(A,$)')' File exists.  Over write? (1=Yes) >'
                  READ(*,*)IANS
                  IF(IANS .NE. 1)GOTO 70
               ENDIF
               OPEN (UNIT=40,FILE=FILE_NAME,STATUS='UNKNOWN')
            ENDIF
         ENDIF
      ENDIF
C
C
 5502 CONTINUE
C
      IF (TIME_START .LT. TIME_IN_RES) THEN
         CALL SEEK_TIME(READ_SPX,TIME_START,REC_POINTER,TIME_FOUND)
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) THEN
               RETURN
            END IF
         END IF
        IF(TIME_FOUND .LT. ZERO) THEN
          WRITE(*,*)' Could not find record for TIME_START'
          GOTO 10
        ENDIF
      ENDIF
C
C  write initial data
C
      CALL WRITE_LINE(FILE_NAME,' ',1)
      CALL WRITE_LINE(FILE_NAME,' ',1)
      DISPLAY = 15
      IF(TIME_START .EQ. TIME_END) THEN
        DISPLAY = DISPLAY - 8
      ENDIF
      IF(I1 .EQ. I2) THEN
        IF(DIRECTION .EQ. 1)THEN
          XTMP = XDIST_VEC(I1)
        ELSE
          XTMP = XDIST_SC(I1)
        ENDIF
        DISPLAY = DISPLAY - 4
        WRITE(LINE,'(A,G12.5)')' X = ', XTMP
        CALL WRITE_LINE(FILE_NAME,LINE,17)
      ELSEIF(I_AVERAGE) THEN
        IF(DIRECTION .EQ. 1)THEN
          XTMP = XDIST_VEC(I1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for X = ', XTMP, ' to ',XDIST_VEC(I2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for X = ', XTMP, ' to ',XDIST_VEC(I2)
          ENDIF
        ELSE
          XTMP = XDIST_SC(I1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for X = ', XTMP, ' to ',XDIST_SC(I2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for X = ', XTMP, ' to ',XDIST_SC(I2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 4
        CALL WRITE_LINE(FILE_NAME,LINE,51)
      ENDIF
C
      IF(J1 .EQ. J2) THEN
        IF(DIRECTION .EQ. 2)THEN
          YTMP = YDIST_VEC(J1)
        ELSE
          YTMP = YDIST_SC(J1)
        ENDIF
        DISPLAY = DISPLAY - 2
        WRITE(LINE,'(A,G12.5)')' Y = ', YTMP
        CALL WRITE_LINE(FILE_NAME,LINE,17)
      ELSEIF(J_AVERAGE) THEN
        IF(DIRECTION .EQ. 2)THEN
          YTMP = YDIST_VEC(J1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for Y = ', YTMP, ' to ',YDIST_VEC(J2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for Y = ', YTMP, ' to ',YDIST_VEC(J2)
          ENDIF
        ELSE
          YTMP = YDIST_SC(J1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for Y = ', YTMP, ' to ',YDIST_SC(J2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for Y = ', YTMP, ' to ',YDIST_SC(J2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 2
        CALL WRITE_LINE(FILE_NAME,LINE,51)
      ENDIF
C
      IF(K1 .EQ. K2) THEN
        IF(DIRECTION .EQ. 3)THEN
          ZTMP = ZDIST_VEC(K1)
        ELSE
          ZTMP = ZDIST_SC(K1)
        ENDIF
        DISPLAY = DISPLAY - 1
        WRITE(LINE,'(A,G12.5)')' Z = ', ZTMP
        CALL WRITE_LINE(FILE_NAME,LINE,17)
      ELSEIF(K_AVERAGE) THEN
        IF(DIRECTION .EQ. 3)THEN
          ZTMP = ZDIST_VEC(K1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for Z = ', ZTMP, ' to ',ZDIST_VEC(K2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for Z = ', ZTMP, ' to ',ZDIST_VEC(K2)
          ENDIF
        ELSE
          ZTMP = ZDIST_SC(K1)
          IF(SUM) THEN
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Sum of values for Z = ', ZTMP, ' to ',ZDIST_SC(K2)
          ELSE
            WRITE(LINE,'(A,G12.5, A, G12.5)')
     &      ' Average value for Z = ', ZTMP, ' to ',ZDIST_SC(K2)
          ENDIF
        ENDIF
        DISPLAY = DISPLAY - 1
        CALL WRITE_LINE(FILE_NAME,LINE,51)
      ENDIF
C
      IF(DISPLAY .EQ. 8 .AND. .NOT.TIME_AVERAGE)THEN
        WRITE(LINE,'(5X,A,10X,1A8)')'Time', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 27)
      ENDIF
      IF((VAR_NO .GE.  7 .AND. VAR_NO .LE. 10) .OR.
     &   (VAR_NO .EQ. 15                     ) .OR.
     &   (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.
     &   (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR.
     &   (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.
     &   (VAR_NO .EQ. 35                     ) .OR.
     &   (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.
     &   (VAR_NO .EQ. 43                     ) .OR.
     &   (VAR_NO .EQ. 44                     ) .OR.
     &   (VAR_NO .EQ. 46                     ) .OR.
     &   (VAR_NO .EQ. 47                     ) .OR.
     &   (VAR_NO .EQ. 48                     ) 
     &                                             )THEN
        WRITE(LINE,'(A,I2)') ' Solids phase = ', M
        CALL WRITE_LINE(FILE_NAME, LINE, 18)
      ENDIF
      IF(VAR_NO .GE. 14 .AND. VAR_NO .LE. 21) THEN
        WRITE(LINE,'(A,I2)') ' Species = ', N
        CALL WRITE_LINE(FILE_NAME, LINE, 13)
      ENDIF
C
C  Read data
C
      END_AVERAGE = .FALSE.
      TIME_OLD = -1.
      IF(TIME_AVERAGE) THEN
        NT = 0
        DO 90 K = K1, K2
        DO 90 J = J1, J2
        DO 90 I = I1, I2
          IJK = FUNIJK(I,J,K)
          VALUE(IJK) = ZERO
90      CONTINUE
      ENDIF
C
      IF(MINMAX .EQ. 0) THEN
        WRITE(LINE,'(6X,A, 10X, A,13X,A,13X,A,5X,A,1A8)')
     &      'Time', 'X', 'Y', 'Z', 'Minimum ', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 70)
      ELSEIF(MINMAX .EQ. 1) THEN
        WRITE(LINE,'(6X,A, 10X, A,13X,A,13X,A,5X,A,1A8)')
     &      'Time', 'X', 'Y', 'Z', 'Maximum ', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 70)
      ENDIF
C
100   CONTINUE
      IF (TIME_START .LT. TIME_IN_RES) THEN
         CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,
     &                    TIME_NOW, TIME_REAL, NSTEP_1)
         IF (DO_XFORMS) THEN
            CALL CHECK_INTER(INTER)
            IF (INTER) THEN
               RETURN
            END IF
         END IF
      ELSE
         CALL READ_RES1
         TIME_NOW = TIME_IN_RES
      ENDIF
C
      IF (TIME_NOW .LT. ZERO) THEN
        IF(.NOT.TIME_AVERAGE) THEN
           IF (.NOT.DO_XFORMS) THEN
              GOTO 10
           ELSE
              IF(FILE_NAME(1:1) .NE. '*') CLOSE(40)
              RETURN
           END IF
        END IF
        TIME_NOW = TIME_OLD
        END_AVERAGE = .TRUE.
        GOTO 106
      ENDIF
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      TIME_OLD = TIME_NOW
C
      IF(MMAX .EQ. 1) THEN
        IF (VAR_NO .EQ. 10 .OR.
     &   (VAR_NO .GE. 19 .AND. VAR_NO .LE. 21) .OR.
     &   (VAR_NO .GE. 25 .AND. VAR_NO .LE. 27) .OR. 
     &   (VAR_NO .GE. 31 .AND. VAR_NO .LE. 33) .OR.
     &   (VAR_NO .EQ. 35                     ) .OR.
     &   (VAR_NO .GE. 39 .AND. VAR_NO .LE. 41) .OR.    
     &   (VAR_NO .EQ. 43                     ) .OR.
     &   (VAR_NO .EQ. 44                     ) .OR.
     &   (VAR_NO .EQ. 46                     ) .OR.
     &   (VAR_NO .EQ. 47                     ) 
     &                                            ) THEN
          DO 102 K = K1, K2
          DO 102 J = J1, J2
          DO 102 I = I1, I2
            IJK = FUNIJK(I, J, K)
            ROP_s(IJK, 1) = (ONE - EP_g(IJK)) * RO_s(1)
102       CONTINUE
        ENDIF
      ENDIF
C
C     FIND THETA IF CALCULATING SOLIDS PRESSURE, P_s
C
      IF(VAR_NO .EQ. 44 .OR. VAR_NO .EQ. 47) THEN
        M_LOCAL=M
        DO M = 1, MMAX
          CALL CALC_MU_s(M, IER)
        ENDDO
        M=M_LOCAL
      ENDIF
C
      NT = NT + 1
      DO 105 K = K1, K2
      DO 105 J = J1, J2
      DO 105 I = I1, I2
        IJK = FUNIJK(I, J, K)
        IF(VAR_NO .EQ. 1) THEN
          VALUE_TMP = EP_g(IJK)
        ELSEIF(VAR_NO .EQ. 2)THEN
          VALUE_TMP = P_g(IJK)
        ELSEIF(VAR_NO .EQ. 3)THEN
          VALUE_TMP = P_star(IJK)
        ELSEIF(VAR_NO .EQ. 4)THEN
          VALUE_TMP = U_g(IJK)
        ELSEIF(VAR_NO .EQ. 5)THEN
          VALUE_TMP = V_g(IJK)
        ELSEIF(VAR_NO .EQ. 6)THEN
          VALUE_TMP = W_g(IJK)
        ELSEIF(VAR_NO .EQ. 7)THEN
          VALUE_TMP = U_s(IJK, M)
        ELSEIF(VAR_NO .EQ. 8)THEN
          VALUE_TMP = V_s(IJK, M)
        ELSEIF(VAR_NO .EQ. 9)THEN
          VALUE_TMP = W_s(IJK, M)
        ELSEIF(VAR_NO .EQ. 10)THEN
          VALUE_TMP = ROP_s(IJK, M)
        ELSEIF(VAR_NO .EQ. 11)THEN
          VALUE_TMP = T_g(IJK)
        ELSEIF(VAR_NO .EQ. 12)THEN
          VALUE_TMP = T_s(IJK, M)
        ELSEIF(VAR_NO .EQ. 13)THEN
          VALUE_TMP = T_s(IJK, 2)
        ELSEIF(VAR_NO .EQ. 14)THEN
          VALUE_TMP = X_g(IJK, N)
        ELSEIF(VAR_NO .EQ. 15)THEN
          VALUE_TMP = X_s(IJK, M, N)
        ELSEIF(VAR_NO .EQ. 16)THEN
          VALUE_TMP = XFLOW_gx(I, J, K, IJK, N)
        ELSEIF(VAR_NO .EQ. 17)THEN
          VALUE_TMP = XFLOW_gy(I, J, K, IJK, N)
        ELSEIF(VAR_NO .EQ. 18)THEN
          VALUE_TMP = XFLOW_gz(I, J, K, IJK, N)
        ELSEIF(VAR_NO .EQ. 19)THEN
          VALUE_TMP = XFLOW_sx(I, J, K, IJK, M, N)
        ELSEIF(VAR_NO .EQ. 20)THEN
          VALUE_TMP = XFLOW_sy(I, J, K, IJK, M, N)
        ELSEIF(VAR_NO .EQ. 21)THEN
          VALUE_TMP = XFLOW_sz(I, J, K, IJK, M, N)
        ELSEIF(VAR_NO .EQ. 22)THEN
          VALUE_TMP = MFLOW_gx(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 23)THEN
          VALUE_TMP = MFLOW_gy(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 24)THEN
          VALUE_TMP = MFLOW_gz(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 25)THEN
          VALUE_TMP = MFLOW_sx(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 26)THEN
          VALUE_TMP = MFLOW_sy(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 27)THEN
          VALUE_TMP = MFLOW_sz(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 28)THEN
          VALUE_TMP = VFLOW_gx(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 29)THEN
          VALUE_TMP = VFLOW_gy(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 30)THEN
          VALUE_TMP = VFLOW_gz(I, J, K, IJK)
        ELSEIF(VAR_NO .EQ. 31)THEN
          VALUE_TMP = VFLOW_sx(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 32)THEN
          VALUE_TMP = VFLOW_sy(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 33)THEN
          VALUE_TMP = VFLOW_sz(I, J, K, IJK, M)
        ELSEIF(VAR_NO .EQ. 34)THEN
          VALUE_TMP = EP_g(IJK) * CALC_RO_g(IJK) * VOL_SC(IJK)
        ELSEIF(VAR_NO .EQ. 35)THEN
          VALUE_TMP = ROP_s(IJK, M) * VOL_SC(IJK)
        ELSEIF(VAR_NO .EQ. 36)THEN
          VALUE_TMP = FLUX_gx(IJK)
        ELSEIF(VAR_NO .EQ. 37)THEN
          VALUE_TMP = FLUX_gy(IJK)
        ELSEIF(VAR_NO .EQ. 38)THEN
          VALUE_TMP = FLUX_gz(IJK)
        ELSEIF(VAR_NO .EQ. 39)THEN
          VALUE_TMP = FLUX_sx(IJK, M)
        ELSEIF(VAR_NO .EQ. 40)THEN
          VALUE_TMP = FLUX_sy(IJK, M)
        ELSEIF(VAR_NO .EQ. 41)THEN
          VALUE_TMP = FLUX_sz(IJK, M)
        ELSEIF(VAR_NO .EQ. 42)THEN
          VALUE_TMP = 0.5*CALC_RO_g(IJK)*EP_g(IJK)
     &              *(U_g(IJK)**2 + V_g(IJK)**2 + W_g(IJK)**2)
        ELSEIF(VAR_NO .EQ. 43)THEN
          VALUE_TMP = 0.5*ROP_s(IJK,M)
     &        *(U_s(IJK,M)**2 + V_s(IJK,M)**2 + W_s(IJK,M)**2)
        ELSEIF(VAR_NO .EQ. 44)THEN
          VALUE_TMP = P_s(IJK,M)
            IF(EP_g(IJK) .LT. EP_star) THEN
              VALUE_TMP = P_star(IJK)
            ENDIF
        ELSEIF(VAR_NO .EQ. 45)THEN
          VALUE_TMP = GRAVITY*YDIST_SC(J)*CALC_RO_g(IJK)*EP_g(IJK)
        ELSEIF(VAR_NO .EQ. 46)THEN
          VALUE_TMP = GRAVITY*YDIST_SC(J)*ROP_s(IJK,M)
        ELSEIF(VAR_NO .EQ. 47)THEN
          mIJK= FUNIJK(I,J+1,K)
          lIJK= FUNIJK(I,J-1,K)
          DELm=YDIST_SC(J+1)-YDIST_SC(J)
          DELl=YDIST_SC(J)-YDIST_SC(J-1)
          FAC1= (ROP_s(mIJK,M)-ROP_s(IJK,M))/DELm
     &         +(ROP_s(IJK,M)-ROP_s(lIJK,M))/DELl
          FAC2=P_s(IJK,M)/(ROP_s(IJK,M)**2)
C          IF(EP_g(IJK) .LT. EP_star) THEN
C             FAC2=P_star(IJK)/(ROP_s(IJK,M)**2)
C          ENDIF
          VALUE_TMP = DY(J)*FAC1*FAC2/2
        ELSEIF(VAR_NO .EQ. 48)THEN
          VALUE_TMP = Theta_m(IJK, M)
        ENDIF
	
	
        IF(TIME_AVERAGE)THEN
          VALUE(IJK) = VALUE(IJK) + VALUE_TMP
        ELSE
          VALUE(IJK) = VALUE_TMP
        ENDIF
105   CONTINUE
C
106   IF(TIME_AVERAGE)THEN
        IF(TIME_NOW .GE. TIME_END .OR. END_AVERAGE)THEN
          IF(NT .EQ. 0) THEN
            WRITE(*,*)' Could not do time averaging'
            GOTO 10
          ENDIF
          DO 110 K = K1, K2
          DO 110 J = J1, J2
          DO 110 I = I1, I2
            IJK = FUNIJK(I,J,K)
            VALUE(IJK) = VALUE(IJK) / REAL(NT)
110       CONTINUE
          WRITE(LINE,'(1X,A,1A8,A,G12.5,A,G12.5)')
     &     'Time average of ',VAR, ' from Time = ',TIME_START,
     &     ' to ', TIME_NOW
          CALL WRITE_LINE(FILE_NAME,LINE,66)
        ELSE
          GOTO 100
        ENDIF
      ENDIF
C
C  DO I, J, or K averaging
C
      IF(K_AVERAGE)THEN
        K = K1
        DO 115 J = J1, J2
        DO 115 I = I1, I2
          IJK = FUNIJK(I, J, K)
          IF(WALL_AT(IJK))THEN
            VALUE(IJK) = ZERO
            DIST(IJK)  = ZERO
          ELSEIF(.NOT. SUM) THEN
            IF(DIRECTION .EQ. 3)THEN
              VALUE(IJK) = VALUE(IJK) * DZ_T(K)
              DIST(IJK)  = DZ_T(K)
            ELSE
              VALUE(IJK) = VALUE(IJK) * DZ(K)
              DIST(IJK)  = DZ(K)
            ENDIF
          ENDIF
115     CONTINUE
C
        DO 116 K = K1+1, K2
        DO 116 J = J1, J2
        DO 116 I = I1, I2
          IJK = FUNIJK(I, J, K)
          IJK1 = FUNIJK(I, J, K1)
          IF(.NOT.WALL_AT(IJK)) THEN
            IF(SUM) THEN
              VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
            ELSE
              IF(DIRECTION .EQ. 3)THEN
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DZ_T(K)
                DIST(IJK1)  = DIST(IJK1)  + DZ_T(K)
              ELSE
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DZ(K)
                DIST(IJK1)  = DIST(IJK1)  + DZ(K)
              ENDIF
            ENDIF
          ENDIF
116     CONTINUE
C
        IF(.NOT. SUM) THEN
          K = K1
          DO 117 J = J1, J2
          DO 117 I = I1, I2
            IJK = FUNIJK(I, J, K)
            IF(DIST(IJK) .NE. ZERO) THEN
              VALUE(IJK) = VALUE(IJK) / DIST(IJK)
            ELSEIF(VALUE(IJK) .NE. ZERO) THEN
              WRITE(*,*)' Error in K-averaging'
            ENDIF
117       CONTINUE
        ENDIF
        K2d = K1
      ELSE
        K2d = K2
      ENDIF
C
      IF(J_AVERAGE)THEN
        J = J1
        DO 118 K = K1, K2d
        DO 118 I = I1, I2
          IJK = FUNIJK(I, J, K)
          IF(WALL_AT(IJK) .AND. .NOT.K_AVERAGE)THEN
            VALUE(IJK) = ZERO
            DIST(IJK)  = ZERO
          ELSEIF(.NOT. SUM) THEN
            IF(DIRECTION .EQ. 2)THEN
              VALUE(IJK) = VALUE(IJK) * DY_N(J)
              DIST(IJK)  = DY_N(J)
            ELSE
              VALUE(IJK) = VALUE(IJK) * DY(J)
              DIST(IJK)  = DY(J)
            ENDIF
          ENDIF
118     CONTINUE
C
        DO 119 K = K1, K2d
        DO 119 J = J1+1, J2
        DO 119 I = I1, I2
          IJK = FUNIJK(I, J, K)
          IJK1 = FUNIJK(I, J1, K)
          IF(.NOT.WALL_AT(IJK) .OR. K_AVERAGE) THEN
            IF(SUM) THEN
              VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
            ELSE
              IF(DIRECTION .EQ. 2)THEN
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DY_N(J)
                DIST(IJK1)  = DIST(IJK1)  + DY_N(J)
              ELSE
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DY(J)
                DIST(IJK1)  = DIST(IJK1)  + DY(J)
              ENDIF
            ENDIF
          ENDIF
119     CONTINUE
C
        IF(.NOT.SUM) THEN
          J = J1
          DO 120 K = K1, K2d
          DO 120 I = I1, I2
            IJK = FUNIJK(I, J, K)
            IF(DIST(IJK) .NE. ZERO) THEN
              VALUE(IJK) = VALUE(IJK) / DIST(IJK)
            ELSEIF(VALUE(IJK) .NE. ZERO) THEN
              WRITE(*,*)' Error in J-averaging'
            ENDIF
120       CONTINUE
        ENDIF
        J2d = J1
      ELSE
        J2d = J2
      ENDIF
C
      IF(I_AVERAGE)THEN
        I = I1
        DO 121 K = K1, K2d
        DO 121 J = J1, J2d
          IJK = FUNIJK(I, J, K)
          IF(WALL_AT(IJK) .AND. .NOT.J_AVERAGE .AND. .NOT.K_AVERAGE)THEN
            VALUE(IJK) = ZERO
            DIST(IJK)  = ZERO
          ELSEIF(.NOT. SUM) THEN
            IF(DIRECTION .EQ. 1)THEN
              VALUE(IJK) = VALUE(IJK) * DX_E(I) * X_E(I)
              DIST(IJK)  = DX_E(I) * X_E(I)
            ELSE
              VALUE(IJK) = VALUE(IJK) * DX(I) * X(I)
              DIST(IJK)  = DX(I) * X(I)
            ENDIF
          ENDIF
121     CONTINUE
C
        DO 122 K = K1, K2d
        DO 122 J = J1, J2d
        DO 122 I = I1+1, I2
          IJK = FUNIJK(I, J, K)
          IJK1 = FUNIJK(I1, J, K)
          IF(.NOT.WALL_AT(IJK) .OR. J_AVERAGE .OR. K_AVERAGE) THEN
            IF(SUM) THEN
              VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK)
            ELSE
              IF(DIRECTION .EQ. 1)THEN
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DX_E(I)* X_E(I)
                DIST(IJK1)  = DIST(IJK1)  + DX_E(I) * X_E(I)
              ELSE
                VALUE(IJK1) = VALUE(IJK1) + VALUE(IJK) * DX(I) * X(I)
                DIST(IJK1)  = DIST(IJK1)  + DX(I) * X(I)
              ENDIF
            ENDIF
          ENDIF
122     CONTINUE
C
        IF(.NOT.SUM) THEN
          I = I1
          DO 123 K = K1, K2d
          DO 123 J = J1, J2d
            IJK = FUNIJK(I, J, K)
            IF(DIST(IJK) .NE. ZERO) THEN
              VALUE(IJK) = VALUE(IJK) / DIST(IJK)
            ELSEIF(VALUE(IJK) .NE. ZERO) THEN
              WRITE(*,*)' Error in I-averaging'
            ENDIF
123       CONTINUE
        ENDIF
        I2d = I1
      ELSE
        I2d = I2
      ENDIF
C
C  Display data or write data to file
C
      IF(MINMAX .GE. 0) THEN
        IF(MINMAX .EQ. 0)THEN
          VALUE_TMP = 1E32
          DO 128 K = K1, K2d
          DO 128 J = J1, J2d
          DO 128 I = I1, I2d
            IJK = FUNIJK(I,J,K)
            IF(WALL_AT(IJK) .AND. 
     &        .NOT. (I_AVERAGE .OR. J_AVERAGE .OR. K_AVERAGE) )GOTO 128
            IF(VALUE(IJK) .GE. VALUE_TMP)GOTO 128
            XTMP = XDIST_SC(I)
            YTMP = YDIST_SC(J)
            ZTMP = ZDIST_SC(K)
            IF(DIRECTION .EQ. 1) THEN
              XTMP = XDIST_VEC(I)
            ELSEIF(DIRECTION .EQ. 2) THEN
              YTMP = YDIST_VEC(J)
            ELSEIF(DIRECTION .EQ. 3) THEN
              ZTMP = ZDIST_VEC(K)
            ENDIF
            VALUE_TMP = VALUE(IJK)
128       CONTINUE
        ELSEIF(MINMAX .EQ. 1)THEN
          VALUE_TMP = -1E32
          DO 129 K = K1, K2d
          DO 129 J = J1, J2d
          DO 129 I = I1, I2d
            IJK = FUNIJK(I,J,K)
            IF(WALL_AT(IJK) .AND. 
     &        .NOT. (I_AVERAGE .OR. J_AVERAGE .OR. K_AVERAGE) )GOTO 129
            IF(VALUE(IJK) .LE. VALUE_TMP)GOTO 129
            XTMP = XDIST_SC(I)
            YTMP = YDIST_SC(J)
            ZTMP = ZDIST_SC(K)
            IF(DIRECTION .EQ. 1) THEN
              XTMP = XDIST_VEC(I)
            ELSEIF(DIRECTION .EQ. 2) THEN
              YTMP = YDIST_VEC(J)
            ELSEIF(DIRECTION .EQ. 3) THEN
              ZTMP = ZDIST_VEC(K)
            ENDIF
            VALUE_TMP = VALUE(IJK)
129       CONTINUE
        ENDIF
        WRITE(LINE,'(1X,5(G12.5,2X))')TIME_NOW, XTMP, YTMP, ZTMP,
     &                                VALUE_TMP
        CALL WRITE_LINE(FILE_NAME, LINE, 71)
      ELSEIF(DISPLAY .EQ. 8) THEN
        IJK = FUNIJK(I1, J1, K1)
        IF(TIME_AVERAGE) THEN
          WRITE(LINE,'(1X,G12.5)')VALUE(IJK)
          CALL WRITE_LINE(FILE_NAME,LINE,13)
        ELSE
          WRITE(LINE,'(1X,G12.5,2X,G12.5)')TIME_NOW, VALUE(IJK)
          CALL WRITE_LINE(FILE_NAME,LINE,27)
        ENDIF
      ELSEIF(DISPLAY .EQ. 4 .OR. DISPLAY .EQ. 12) THEN
        WRITE(LINE,'(A,G12.5)')' Time = ',TIME_NOW
        CALL WRITE_LINE(FILE_NAME,LINE,20)
        WRITE(LINE,'(6X,A,13X,1A8)')'X', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 28)
        DO 130 I = I1, I2d
          IJK = FUNIJK(I, J1, K1)
          IF(DIRECTION .EQ. 1)THEN
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')XDIST_VEC(I), VALUE(IJK)
          ELSE
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')XDIST_SC(I), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, 27)
130     CONTINUE
      ELSEIF(DISPLAY .EQ. 2 .OR. DISPLAY .EQ. 10) THEN
        WRITE(LINE,'(A,G12.5)')' Time = ',TIME_NOW
        CALL WRITE_LINE(FILE_NAME,LINE,20)
        WRITE(LINE,'(6X,A,13X,1A8)')'Y', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 28)
        DO 140 J = J1, J2d
          IJK = FUNIJK(I1, J, K1)
          IF(DIRECTION .EQ. 2)THEN
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')YDIST_VEC(J), VALUE(IJK)
          ELSE
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')YDIST_SC(J), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, 27)
140     CONTINUE
      ELSEIF(DISPLAY .EQ. 1 .OR. DISPLAY .EQ. 9) THEN
        WRITE(LINE,'(A,G12.5)')' Time = ',TIME_NOW
        CALL WRITE_LINE(FILE_NAME,LINE,20)
        WRITE(LINE,'(6X,A,13X,1A8)')'Z', VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 28)
        DO 150 K = K1, K2d
          IJK = FUNIJK(I1, J1, K)
          IF(DIRECTION .EQ. 3)THEN
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')ZDIST_VEC(K), VALUE(IJK)
          ELSE
            WRITE(LINE,'(1X,G12.5,2X,G12.5)')ZDIST_SC(K), VALUE(IJK)
          ENDIF
          CALL WRITE_LINE(FILE_NAME, LINE, 27)
150     CONTINUE
      ELSEIF(DISPLAY .EQ. 0) THEN
        WRITE(LINE,'(A,G12.5)')' Time = ',TIME_NOW
        CALL WRITE_LINE(FILE_NAME,LINE,20)
        IJK = FUNIJK(I1, J1, K1)
        WRITE(LINE,'(1X,1A8,A,G12.5)')VAR,' = ',VALUE(IJK)
        CALL WRITE_LINE(FILE_NAME,LINE,24)
      ELSE
        WRITE(LINE,'(A,G12.5)')' Time = ',TIME_NOW
        CALL WRITE_LINE(FILE_NAME,LINE,20)
        WRITE(LINE,'(6X,A,13X,A,13X,A,13X,1A8)')'X','Y','Z',VAR
        CALL WRITE_LINE(FILE_NAME, LINE, 56)
        DO 250 K = K1, K2d
        DO 250 J = J1, J2d
        DO 250 I = I1, I2d
          IJK = FUNIJK(I,J,K)
          XTMP = XDIST_SC(I)
          YTMP = YDIST_SC(J)
          ZTMP = ZDIST_SC(K)
          IF(DIRECTION .EQ. 1) THEN
            XTMP = XDIST_VEC(I)
          ELSEIF(DIRECTION .EQ. 2) THEN
            YTMP = YDIST_VEC(J)
          ELSEIF(DIRECTION .EQ. 3) THEN
            ZTMP = ZDIST_VEC(K)
          ENDIF
          WRITE(LINE,'(1X,4(G12.5,2X))')XTMP, YTMP, ZTMP, VALUE(IJK)
          CALL WRITE_LINE(FILE_NAME, LINE, 57)
250     CONTINUE
      ENDIF
C
      IF (TIME_AVERAGE .AND. END_AVERAGE) THEN
        IF (DO_XFORMS) THEN
           IF (FILE_NAME(1:1) .NE. '*') CLOSE(40)
           RETURN
        ELSE
           GOTO 10
        END IF
      END IF
C
      IF(TIME_NOW .GE. TIME_END) THEN
        IF (DO_XFORMS) THEN
           IF (FILE_NAME(1:1) .NE. '*') CLOSE(UNIT=40)
           RETURN
        ELSE
           GOTO 10
        END IF
      END IF
C
      GOTO 100
      END
c
      SUBROUTINE WRITE_LINE(FILE_NAME,LINE,NCHARS)
      IMPLICIT NONE
      CHARACTER*(*) FILE_NAME, LINE
      CHARACTER*81  LINE2
      INTEGER       NCHARS
C
      INCLUDE 'xforms.inc'
C
      IF (DO_XFORMS) THEN
         LINE2 = LINE
         LINE2(NCHARS+1:NCHARS+1) = CHAR(0)
         CALL ADD_TO_RESULTS_BROWSER(LINE2(1:NCHARS+1))
      END IF
      IF(FILE_NAME(1:1) .EQ. '*')THEN
        IF (.NOT.DO_XFORMS) WRITE(*,*) LINE(1:NCHARS)
      ELSE
        WRITE(40,*) LINE
      ENDIF
      RETURN
      END
C
      SUBROUTINE HELP(N)
      IMPLICIT NONE
      INTEGER N
C
      WRITE(*,*)
C
      IF(N .EQ. 10)THEN
        WRITE(*,*)
     &  ' Enter start time and end time for data retrieval.  If the'
        WRITE(*,*)
     &  ' time entered is negative, control returns to the main menu.'
      ELSEIF(N .EQ. 11)THEN
        WRITE(*,*)
     &  ' Enter Y or y for time averaging.'
      ELSEIF(N .EQ. 20)THEN
        WRITE(*,'(42(A,/))')
     &  ' Valid variable names:',
     &  '   EP_g      - Void fraction',
     &  '   P_g       - Gas pressure, dyne/cm^2',
     &  '   P_star    - Solids pressure (frictional regime), dyne/cm^2',
     &  '   Gas velocity, cm/s',
     &  '   U_g       - X component',
     &  '   V_g       - Y component',
     &  '   W_g       - Z component',
     &  '   Solids velocity, cm/s',
     &  '   U_s       - X component',
     &  '   V_s       - Y component',
     &  '   W_s       - Z component',
     &  '   ROP_s     - Solids density x volume fraction, g/cm^3',
     &  '   T_g       - Gas temperature, K',
     &  '   T_s       - Solids temperature, K',
     &  '   Theta_m   - Granular temperature, cm^2/s^2',
     &  '   X_g       - Gas species mass fraction',
     &  '   X_s       - Solids species mass fraction',
     &  '   Gas species mass flow rates, g/s',
     &  '   XFLOW_gx  - in X direction',
     &  '   XFLOW_gy  - in Y direction',
     &  '   XFLOW_gz  - in Z direction',
     &  '   Solids species mass flow rates, g/s',
     &  '   XFLOW_sx  - in X direction',
     &  '   XFLOW_sy  - in Y direction',
     &  '   XFLOW_sz  - in Z direction',
     &  '   Gas total mass flow rates, g/s',
     &  '   MFLOW_gx  - in X direction',
     &  '   MFLOW_gy  - in Y direction',
     &  '   MFLOW_gz  - in Z direction',
     &  '   Solids total mass flow rates, g/s',
     &  '   MFLOW_sx  - in X direction',
     &  '   MFLOW_sy  - in Y direction',
     &  '   MFLOW_sz  - in Z direction',
     &  '   Gas total volumetric flow rates, cm^3/s',
     &  '   VFLOW_gx  - in X direction',
     &  '   VFLOW_gy  - in Y direction',
     &  '   VFLOW_gz  - in Z direction',
     &  '   Solids total volumetric flow rates, cm^3/s',
     &  '   VFLOW_sx  - in X direction',
     &  '   VFLOW_sy  - in Y direction',
     &  '   VFLOW_sz  - in Z direction',
     &  '   Gas total mass flux, g/(cm^2.s)', 
     &  '   FLUX_gx   - in X direction',
     &  '   FLUX_gy   - in Y direction',
     &  '   FLUX_gz   - in Z direction',
     &  '   Solids total mass flux, g/(cm^2.s)', 
     &  '   FLUX_sx   - in X direction',
     &  '   FLUX_sy   - in Y direction',
     &  '   FLUX_sz   - in Z direction',
     &  '   MASS_g    - Total mass of gas in specified volume, g',
     &  '   MASS_s    - Total mass of solids in specified volume, g',
     &  ' To determine the minimum or maximum preceed the variable',
     &  ' name with a 0 or 1.'
      ELSEIF(N .EQ. 24)THEN
        WRITE(*,*)' Enter the solids phase number'
      ELSEIF(N .EQ. 25)THEN
        WRITE(*,*)' Enter the species number'
      ELSEIF(N .EQ. 26)THEN
        WRITE(*,*)' Enter Y or y to use P_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of P_g'
      ELSEIF(N .EQ. 27)THEN
        WRITE(*,*)' Enter Y or y to use T_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of T_g'
      ELSEIF(N .EQ. 28)THEN
        WRITE(*,*)' Enter Y or y to use X_g values from the RES file'
        WRITE(*,*)' The values may differ from current values of X_g'
      ELSEIF(N .EQ. 30)THEN
        WRITE(*,*)' Enter the starting and ending values of I'
      ELSEIF(N .EQ. 31)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 40)THEN
        WRITE(*,*)' Enter the starting and ending values of J'
      ELSEIF(N .EQ. 41)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 50)THEN
        WRITE(*,*)' Enter the starting and ending values of K'
      ELSEIF(N .EQ. 51)THEN
        WRITE(*,*)' Intensive variables, such as EP_g, are averaged'
        WRITE(*,*)' Extensive variables, such as VFLOW_gy, are summed'
      ELSEIF(N .EQ. 70)THEN
        WRITE(*,*)' Enter the file name.  Enter * for displaying'
        WRITE(*,*)' the values at the terminal.  Enter ! for going'
        WRITE(*,*)' to the begining of this menu.'
      ENDIF
C
      WRITE(*,*)
C
      RETURN
      END
