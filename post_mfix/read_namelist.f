CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: READ_NAMELIST(POST)                                    C
C  Purpose: Read in the NAMELIST variables                             C
C                                                                      C
C  Author: P. Nicoletti                               Date: 25-NOV-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: None                                          C
C  Variables modified: RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME,   C
C                      TSTOP, DT, RES_DT, SPX_DT , OUT_DT, USR_DT,     C
C                      USR_TYPE, USR_VAR, USR_FORMAT, USR_EXT, USR_X_w C
C                      USR_X_e , USR_Y_s , USR_Y_n , USR_Z_b , USR_Z_t,C
C                      NLOG, COORDINATES, IMAX, DX, XLENGTH, JMAX, DY, C
C                      YLENGTH, KMAX, DZ, ZLENGTH, MMAX, D_p, RO_s,    C
C                      EP_star, MU_g0, MW_AVG, IC_X_w, IC_X_e, IC_Y_s, C
C                      IC_Y_n, IC_Z_b, IC_Z_t, IC_I_w, IC_I_e, IC_J_s, C
C                      IC_J_n, IC_K_b, IC_K_t, IC_EP_g, IC_P_g,        C
C                      IC_ROP_s, IC_T_g, IC_T_s,  IC_U_g,     C
C                      IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s, BC_X_w, C
C                      BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t, BC_I_w, C
C                      BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t, BC_EP_g,C
C                      BC_P_g, BC_RO_g, BC_ROP_g, BC_ROP_s, BC_T_g,    C
C                      BC_T_s,  BC_U_g, BC_U_s,BC_V_g, BC_V_s,C
C                      BC_W_g, BC_W_s, BC_TYPE, BC_VOLFLOW_g,          C
C                      BC_VOLFLOW_s, BC_MASSFLOW_g, BC_MASSFLOW_s,     C
C                      BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_l,C
C                      BC_Jet_gl, RO_g0                                C
C                                                                      C
C  Local variables: LINE_STRING, MAXCOL, LINE_TOO_BIG, COMMENT_INDEX   C
C                   SEEK_COMMENT                                       C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE READ_NAMELIST(POST)
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE 'run.inc'
      INCLUDE 'output.inc'
      INCLUDE 'physprop.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'ic.inc'
      INCLUDE 'is.inc'
      INCLUDE 'bc.inc'
      INCLUDE 'fldvar.inc'
      INCLUDE 'constant.inc'
      INCLUDE 'indices.inc'
      INCLUDE 'toleranc.inc'
      INCLUDE 'funits.inc'
      INCLUDE 'scales.inc'
      INCLUDE 'ur_facs.inc'
      INCLUDE 'leqsol.inc'
      INCLUDE 'usrnlst.inc'
      INCLUDE 'residual.inc'
      INCLUDE 'rxns.inc'
C
C               LINE_STRING(1:MAXCOL) has valid input data
      INTEGER   MAXCOL
      PARAMETER (MAXCOL = 80)
C
C               holds one line in the input file
      CHARACTER LINE_STRING*132
C
C               integer function ... if > 0 ... data past column
C               MAXCOL in input file
      INTEGER   LINE_TOO_BIG 
C
C               Length of noncomment string
      INTEGER   LINE_LEN
C
C               integer function which returns COMMENT_INDEX
      INTEGER   SEEK_COMMENT 
C
C               This routine is called from: 0 -  mfix; 1 - post_mfix
      INTEGER   POST 
C
C                      Coefficient of restitution (old symbol)
      DOUBLE PRECISION e

C
C                      Indicate whether currently reading rxns or rate
      LOGICAL          RXN_FLAG
C
C                      Indicate whether to do a namelist read on the line
      LOGICAL          READ_FLAG

C
C                      blank line function
      LOGICAL          BLANK_LINE
      INTEGER   L

      INCLUDE 'namelist.inc'


      e = UNDEFINED

      RXN_FLAG     = .FALSE.

      READ_FLAG = .TRUE.

      NO_OF_RXNS = 0

C
C OPEN MFIX ASCII INPUT FILE, AND LOOP THRU IT
C
      OPEN (UNIT=UNIT_DAT, FILE='mfix.dat', STATUS='OLD', ERR=910)

100   CONTINUE
      READ (UNIT_DAT, 1100, END=500) LINE_STRING

      LINE_LEN = SEEK_COMMENT(LINE_STRING, LEN(LINE_STRING)) - 1 
      CALL REMOVE_COMMENT(LINE_STRING, LINE_LEN+1, LEN(LINE_STRING))

      IF(LINE_LEN .LE. 0)GOTO 100           !comment line
      IF(BLANK_LINE(LINE_STRING))GOTO 100   !blank line

      IF (LINE_TOO_BIG(LINE_STRING, LINE_LEN, MAXCOL) .GT. 0) then
        WRITE(*, 1300) LINE_STRING
        CALL EXIT
      END IF


      CALL MAKE_UPPER_CASE(LINE_STRING, LINE_LEN)

C
C  Complete arithmetic operations and expand line
C

      CALL PARSE_LINE(LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG)

C
C Write the current line to a scratch file
C and read the scratch file in NAMELIST format
C
      IF(READ_FLAG) THEN
        OPEN (UNIT=UNIT_TMP,STATUS='SCRATCH',ERR=900)
        WRITE (UNIT_TMP,1000)
        WRITE(UNIT_TMP,1100) LINE_STRING(1:LINE_LEN)
        WRITE (UNIT_TMP,1200)
        REWIND (UNIT=UNIT_TMP)
C       READ (UNIT_TMP,NML=INPUT_DATA,ERR=930,END=930)  ! Use this for FPS
        READ (UNIT_TMP,INPUT_DATA,ERR=400,END=930)
        GOTO 410
400     IF(POST .EQ. 1) GOTO 410  !If called from POST_MFIX ignore the error
        REWIND (UNIT=UNIT_TMP)
        WRITE (UNIT_TMP,1010)
        WRITE(UNIT_TMP,1100) LINE_STRING(1:LINE_LEN)
        WRITE (UNIT_TMP,1200)
        REWIND (UNIT=UNIT_TMP)
        READ (UNIT_TMP,USR_INPUT_DATA,ERR=930,END=930)
410     CLOSE (UNIT=UNIT_TMP)
      ENDIF

      GOTO 100
C
500   CLOSE (UNIT=UNIT_DAT)
C
      IF(e .NE. UNDEFINED)C_e = e
C
      DO 510 L = 1, NO_OF_RXNS
c
c  Do all the rxns have a rate defined?
c
        IF(.NOT.GOT_RATE(L))THEN

          WRITE(*,1610)RXN_NAME(L)
          STOP
   
        ENDIF
c
c  Do all the rxns have a stoichiometry defined?
c
        IF(.NOT.GOT_RXN(L))THEN

          WRITE(*,1620)RXN_NAME(L)
          STOP
   
        ENDIF
510   CONTINUE

      RETURN
C
C HERE IF AN ERROR OCCURED OPENNING/READING THE FILES
C
900   WRITE(*,1400)
      CALL EXIT
910   WRITE(*,1500)
      CALL EXIT
930   WRITE(*,1600)LINE_STRING(1:LINE_LEN)
      CALL EXIT
C
1000  FORMAT(1X,'$INPUT_DATA')
1010  FORMAT(1X,'$USR_INPUT_DATA')
1100  FORMAT(A)
1200  FORMAT(1X,'$END')
1300  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'The following line in mfix.dat has data past column 80',
     &/1X, A80 ,/1X, 70('*')/)
1400  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'Unable to open scratch file' ,/1X, 70('*')/)
1500  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'Unable to open mfix.dat file' ,/1X, 70('*')/)
1600  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'Error while reading the following line in mfix.dat',/1X, A80,
     &/1X,'Possible causes are 1. incorrect format, 2. unknown name,',
     &/1X,'3. the item is dimensioned too small (see PARAM.INC file).',
     &/1X, 70('*')/)
1610  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'No rxn rate defined for rxn: ', A, /1X, 70('*')/)
1620  FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',
     &'No stoichiometry defined for rxn: ', A, /1X, 70('*')/)
C
      END


      LOGICAL FUNCTION blank_line(LINE)

      IMPLICIT none

      character*(*) LINE
      integer l

      blank_line = .false.
      do l = 1, len(line)
        if(line(l:l) .ne. ' ' .and. line(l:l) .ne. '	' )return
      end do

      blank_line = .true.
      return

      end
