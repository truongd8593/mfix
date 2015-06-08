!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_NAMELIST(POST)

      USE bc
      USE cdist
      USE compar
      USE constant
      USE cutcell
      USE dashboard
      USE des_bc
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE error_manager
      USE fldvar
      USE funits
      USE geometry
      USE ic
      USE indices
      USE is
      USE leqsol
      USE mfix_pic
      USE output
      USE parallel
      USE param
      USE param1
      USE particle_filter
      USE physprop
      USE pic_bc
      USE polygon
      USE ps
      USE qmom_kinetic_equation
      USE quadric
      USE residual
      USE run
      USE rxns
      USE scalars
      USE scales
      USE stiff_chem
      USE toleranc
      USE ur_facs
      USE usr
      USE utilities, ONLY: blank_line, line_too_big, seek_comment
      USE vtk
      Use stl

      IMPLICIT NONE

! Dummy Arguments:
!------------------------------------------------------------------------//
!  This routine is called from: 0 -  mfix; 1 - post_mfix
      INTEGER :: POST

! Local Variables:
!------------------------------------------------------------------------//
! LINE_STRING(1:MAXCOL) has valid input data
      INTEGER, PARAMETER :: MAXCOL = 80
! Holds one line in the input file
      CHARACTER(LEN=512) :: LINE_STRING
! Length of noncomment string
      INTEGER :: LINE_LEN
! Line number
      INTEGER :: LINE_NO
! Coefficient of restitution (old symbol)
      DOUBLE PRECISION e
! Indicates whether currently reading rxns or rate
      LOGICAL :: RXN_FLAG
! Indicate whether to do a namelist read on the line
      LOGICAL :: READ_FLAG
! Logical to check if file exits.
      LOGICAL :: lEXISTS
! Error flag
      LOGICAL :: ERROR

      CHARACTER(len=256) :: STRING
      INTEGER :: IOS, II

! External namelist files:
!---------------------------------------------------------------------//
      INCLUDE 'usrnlst.inc'
      INCLUDE 'namelist.inc'
      INCLUDE 'desnamelist.inc'
      INCLUDE 'cartesian_grid_namelist.inc'
      INCLUDE 'qmomknamelist.inc'

      E = UNDEFINED
      RXN_FLAG = .FALSE.
      READ_FLAG = .TRUE.
      NO_OF_RXNS = 0
      LINE_NO = 0

! Open the mfix.dat file. Report errors if the file is not located or
! there is difficulties opening it.
      inquire(file='mfix.dat',exist=lEXISTS)
      IF(.NOT.lEXISTS) THEN
         IF(myPE == PE_IO) WRITE(*,1000)
         CALL MFIX_EXIT(myPE)

 1000 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1000: ',    &
         'The input data file, mfix.dat, is missing. Aborting.',/1x,   &
         70('*'),2/)

      ELSE
         OPEN(UNIT=UNIT_DAT, FILE='mfix.dat', STATUS='OLD', IOSTAT=IOS)
         IF(IOS /= 0) THEN
            IF(myPE == PE_IO) WRITE (*,1100)
            CALL MFIX_EXIT(myPE)
         ENDIF

 1001 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1001: ',    &
         'Unable to open the mfix.dat file. Aborting.',/1x,70('*'),2/)
      ENDIF


! Loop through the mfix.dat file and process the input data.
      READ_LP: DO
         READ (UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP
         IF(IOS > 0) THEN
         ENDIF

         LINE_NO = LINE_NO + 1

         LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(LINE_STRING)) - 1
         CALL REMOVE_COMMENT(LINE_STRING, LINE_LEN+1, LEN(LINE_STRING))

         IF(LINE_LEN <= 0) CYCLE READ_LP           ! comment line
         IF(BLANK_LINE(LINE_STRING)) CYCLE READ_LP ! blank line

         IF(LINE_TOO_BIG(LINE_STRING,LINE_LEN,MAXCOL) > 0) THEN
            WRITE (*, 1100) trim(iVAL(LINE_NO)), trim(ival(MAXCOL)), &
               LINE_STRING(1:MAXCOL)
            CALL MFIX_EXIT(myPE)
         ENDIF

 1100 FORMAT(//1X,70('*')/1x,'From: READ_NAMELIST',/1x,'Error 1100: ', &
         'Line ',A,' in mfix.dat has is too long. Input lines should', &
         /1x,'not pass column ',A,'.',2/3x,A,2/1x,'Please correct ',   &
         'the mfix.dat file.',/1X,70('*'),2/)

! All subsequent lines are thermochemical data
         IF(LINE_STRING(1:11) == 'THERMO DATA') EXIT READ_LP

         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
! At this point, the keyword was not identified therefore it is
! either deprecated or unknown.
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF

      ENDDO READ_LP

      DO II=1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(ii,LINE_STRING)
         LINE_LEN = len(line_string)
         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF
      ENDDO

      CLOSE(UNIT=UNIT_DAT)
      IF (E /= UNDEFINED) C_E = E

      RETURN

CONTAINS


! returns true if there is an error
  SUBROUTINE SET_KEYWORD(ERROR)

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) ::ERROR

    ERROR = .FALSE.

! Make upper case all except species names
    if(index(LINE_STRING,'SPECIES_NAME') == 0 .AND. &
         index(LINE_STRING,'species_name') == 0 .AND. &
         index(LINE_STRING,'Species_Name') == 0 .AND. &
         index(LINE_STRING,'SPECIES_g') == 0 .AND.    &
         index(LINE_STRING,'Species_g') == 0 .AND.    &
         index(LINE_STRING,'species_g') == 0 .AND.    &
         index(LINE_STRING,'SPECIES_s') == 0 .AND.    &
         index(LINE_STRING,'Species_s') == 0 .AND.    &
         index(LINE_STRING,'species_s') == 0)         &
         CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)

    CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
    CALL REMOVE_PAR_BLANKS(LINE_STRING)

! Complete arithmetic operations and expand line
    CALL PARSE_LINE (LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
    IF(.NOT.READ_FLAG) RETURN

! Standard model input parameters.
    STRING=''; STRING = '&INPUT_DATA '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
    READ(STRING, NML=INPUT_DATA, IOSTAT=IOS)
    IF(IOS == 0)  RETURN

! Stop processing keyword inputs if runing POST_MFIX
    IF(POST == 1) RETURN

! Discrete Element model input parameters.
    STRING=''; STRING = '&DES_INPUT_DATA '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
    READ(STRING, NML=DES_INPUT_DATA, IOSTAT=IOS)
    IF(IOS == 0)  RETURN

! User defined input parameters.
    STRING=''; STRING = '&USR_INPUT_DATA '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
    READ(STRING, NML=USR_INPUT_DATA, IOSTAT=IOS)
    IF(IOS == 0)  RETURN

! Cartesian grid cut-cell input parameters.
    STRING=''; STRING = '&CARTESIAN_GRID_INPUT_DATA '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
    READ(STRING, NML=CARTESIAN_GRID_INPUT_DATA, IOSTAT=IOS)
    IF(IOS == 0)  RETURN

! QMOMK input parameters.
    STRING=''; STRING = '&QMOMK_INPUT_DATA '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
    READ(STRING, NML=QMOMK_INPUT_DATA, IOSTAT=IOS)
    IF(IOS == 0)  RETURN

    ERROR = .TRUE.

    RETURN

  END SUBROUTINE SET_KEYWORD

END SUBROUTINE READ_NAMELIST
