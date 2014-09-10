!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_NAMELIST(POST)

      USE param
      USE param1
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE is
      USE bc
      USE ps
      USE fldvar
      USE constant
      USE indices
      USE toleranc
      USE funits
      USE scales
      USE ur_facs
      USE leqsol
      USE residual
      USE rxns
      USE scalars
      USE compar
      USE parallel
      USE discretelement
      USE mfix_pic
      USE usr
      USE des_bc
      USE pic_bc
      USE des_thermo
      USE des_rxns
      USE stiff_chem
      USE cdist
      USE quadric
      USE cutcell
      USE vtk
      USE polygon
      USE dashboard
      Use stl
      USE qmom_kinetic_equation

      use error_manager

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
      CHARACTER :: LINE_STRING*132
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

      CHARACTER(len=256) :: STRING
      INTEGER :: IOS

! External Functions
!---------------------------------------------------------------------//
! Returns integer if data past column MAXCOL.
      INTEGER, EXTERNAL :: LINE_TOO_BIG
! Integer function which returns COMMENT_INDEX
      INTEGER, EXTERNAL :: SEEK_COMMENT
! Blank line function
      LOGICAL, EXTERNAL :: BLANK_LINE

! External namelist files:
!---------------------------------------------------------------------//
      INCLUDE 'usrnlst.inc'
      INCLUDE 'namelist.inc'
      INCLUDE 'des/desnamelist.inc'
      INCLUDE 'cartesian_grid/cartesian_grid_namelist.inc'
      INCLUDE 'qmomk/qmomknamelist.inc'



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
         'the mfix.dat file.',/1X,70('*')2/)

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

! All subsequent lines are thermochemical data
         IF(LINE_STRING(1:11) == 'THERMO DATA') EXIT READ_LP

         CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
         CALL REMOVE_PAR_BLANKS(LINE_STRING)

! Complete arithmetic operations and expand line
         CALL PARSE_LINE (LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
         IF(.NOT.READ_FLAG) CYCLE READ_LP

! Standard model input parameters.
         STRING=''; STRING = '&INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=INPUT_DATA, IOSTAT=IOS)
         IF(IOS ==0) CYCLE READ_LP

! Stop processing keyword inputs if runing POST_MFIX
         IF(POST == 1) CYCLE READ_LP

! Discrete Element model input parameters.
         STRING=''; STRING = '&DES_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=DES_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0) CYCLE READ_LP

! User defined input parameters.
         STRING=''; STRING = '&USR_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=USR_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0) CYCLE READ_LP

! Cartesian grid cut-cell input parameters.
         STRING=''; STRING = '&CARTESIAN_GRID_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=CARTESIAN_GRID_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0) CYCLE READ_LP

! QMOMK input parameters.
         STRING=''; STRING = '&QMOMK_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=QMOMK_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0) CYCLE READ_LP

! At this point, the keyword was not identified therefore it is 
! either deprecated or unknown.
        CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))

      ENDDO READ_LP

      CLOSE(UNIT=UNIT_DAT)
      IF (E /= UNDEFINED) C_E = E

      RETURN

      END SUBROUTINE READ_NAMELIST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Funcation: BLANK_LINE                                                !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Return .TRUE. if a line contains no input or only spaces.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION BLANK_LINE (line)

      IMPLICIT NONE

      CHARACTER :: LINE*(*)

      INTEGER :: L

      BLANK_LINE = .FALSE.
      DO L=1, len(line)
         IF(line(L:L)/=' ' .and. line(L:L)/='    ')RETURN
      ENDDO

      BLANK_LINE = .TRUE.
      RETURN
      END FUNCTION BLANK_LINE
