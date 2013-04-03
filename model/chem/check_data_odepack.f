!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_CHEM                                        C
!  Purpose: check the chemical rxns namelist variables for             C
!           CALL_DI or CALL_ISAT                                     C
!                                                                      C
!  Author: NAN XIE                              Date: 02-Aug-04        C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Variables referenced:  CALL_DI , CALL_ISAT , ISATdt               C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_ODEPACK

! Global Variables:
!---------------------------------------------------------------------//
      use compar,     only : myPe
      use funits,     only : UNIT_LOG
      use param1,     only : ZERO, UNDEFINED
      use physprop,   only : RO_G0
      use run,        only : ENERGY_EQ, SPECIES_EQ
      use rxns,       only : SUM_R_g, SUM_R_s
      use stiff_chem, only : CALL_ISAT, CALL_DI, STIFF_CHEMISTRY, ISATDT

      implicit none

! Error - ISAT is no longer available.
      IF(CALL_ISAT) THEN
!         IF(DMP_LOG) THEN
            WRITE(*,1001); WRITE(*,1000)
            WRITE(UNIT_LOG,1001); WRITE(UNIT_LOG,1000)
!         ENDIF
         CALL MFIX_EXIT(myPE)
      END IF

      IF(CALL_DI .AND. .NOT.STIFF_CHEMISTRY) THEN
         WRITE(*,1003); WRITE(*,1000)
         WRITE(UNIT_LOG,1003); WRITE(UNIT_LOG,1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

! If the stiff solver is not being used, there is no need for the
! following checks.
      IF(.NOT.STIFF_CHEMISTRY) RETURN

! Message - ISATDT is no longer used.
      IF(ISATdt /= UNDEFINED) WRITE(UNIT_LOG,1002)

! Energy equations must be solved.
      IF(.NOT.ENERGY_EQ) THEN
!         IF(DMP_LOG) THEN
            WRITE(*,1004)'ENERGY_EQ = .FALSE.'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1004)'ENERGY_EQ = .FALSE.'
            WRITE(UNIT_LOG,1000)
!         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF


      IF(RO_G0 /= UNDEFINED) THEN
!         IF(DMP_LOG) THEN
            WRITE(*,1003)'RO_G0 /= UNDEFINED'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1003)'RO_G0 /= UNDEFINED'
            WRITE(UNIT_LOG,1000)
!         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(.NOT.SPECIES_EQ(0)) THEN
!         IF(DMP_LOG) THEN
            WRITE(*,1003)'SPECIES_EQ(0) = .FALSE.'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1003)'SPECIES_EQ(0) = .FALSE.'
            WRITE(UNIT_LOG,1000)
!         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF





! Initialize ODEPACK operating parameters.
      IF(STIFF_CHEMISTRY) THEN
         if(allocated(SUM_R_g)) SUM_R_g = ZERO
         if(allocated(SUM_R_s)) SUM_R_s = ZERO
         CALL ODEPACK_INIT
      ENDIF

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1001 FORMAT(//1X,70('*')/' From: CHECK_DATA_ODEPACK',/' Error 1001:', &
         ' ISAT functionality has been disabled. Chemcial reactions',/ &
         ' may be included using the urs_rates.f UDF and if desired',  &
         ' the stiff',/' chemistry solve may be invoked setting the',  &
         ' keyword STIFF_CHEMISTRY',/' in the mfix.dat file.')

 1002 FORMAT(//1X,70('*')/' From: CHECK_DATA_ODEPACK',/                &
         ' Message 1002: ISATDT is a legacy variable that is no',      &
         ' longer necessary.',/1X,70('*')//)

 1003 FORMAT(//1X,70('*')/' From: CHECK_DATA_ODEPACK',/                &
         ' Error 1003: CALL_DI is a legacy variable. This keyword was',&
         ' replaced',/' by the keyword STIFF_CHEMISTRY.')

 1004 FORMAT(//1X,70('*')/' From: CHECK_DATA_ODEPACK',/' Error 1004:', &
         ' Invalid parameters for stiff chemistry solver!',//          &
         ' The following criteria must be satisfied:',/                &
         '   > Solving the energy equations.',/                        &
         '   > Compressible gas phase.',/                              &
         '   > Solving gas phase species equations.',//                &
         ' >>> Invalid Parameter: ',A)

      RETURN
      END SUBROUTINE CHECK_DATA_ODEPACK



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM_ODEPACK_INIT                                     C
!     Purpose: controlling values for ODEAPCK(reference to ODEPACK manual)C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ODEPACK_INIT

      use physprop,   only : MMAX, NMAX
      use run,        only : SOLID_RO_V

      use stiff_chem

      implicit none

      INTEGER :: LRN, LRS
      INTEGER :: N_VAR

! Dimension of ODEs solved in ISAT
      ODE_DIMN = 2 + NMAX(0) ! Gas Phase density, temperature, species

      N_VAR = 2 ! Solids Phase volume fraction, temperature
      IF(CALL_GROW)  N_VAR = N_VAR + 1  ! diameter
      IF(SOLID_RO_V) N_VAR = N_VAR + 1  ! density

      ODE_DIMN = ODE_DIMN + MMAX*N_VAR + sum(NMAX(1:MMAX))

! Indicates type of Error control.
      ODE_ITOL = 2 ! :: EWT(i) = RTOL * ABS(Y(i)) * ATOL(i)

! Relative error tolerance paramter.
      ODE_RTOL(1) = 1.0D-3

! Absolue error tolerance parameter.
      IF(.NOT.(allocated(ODE_ATOL))) allocate(ODE_ATOL(ODE_DIMN))
      ODE_ATOL(:) = 1.0D-5  ! All Equations

! Declared length of RWORK.
      LRN = 20 + 16*ODE_DIMN
      LRS = 22 + 9*ODE_DIMN + (ODE_DIMN**2)
      ODE_LRW = max(LRN, LRS)

! Declared length of IWORK.
      ODE_LIW = 20 + ODE_DIMN

! Jacobian type indicator.
      ODE_JT = 2 ! Internally generated.

      return
      END SUBROUTINE ODEPACK_INIT
