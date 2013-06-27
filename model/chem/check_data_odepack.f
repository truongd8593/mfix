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
! File unit of log file.
      use funits, only: UNIT_LOG
! Dimension of IJK arrays.
      use param, only: DIMENSION_3
! Double precision zero.
      use param1, only: ZERO
! Double precision value for undefined variables.
      use param1, only: UNDEFINED
! Constant gas phase density
      use physprop, only : RO_G0
! Runtime logical for solving energy equations.
      use run, only: ENERGY_EQ
! Run time logical for solving species equations.
      use run, only: SPECIES_EQ
! Net rate of gas phase production/consumption
      use rxns, only: SUM_R_g
! Net rate of solids phase production/consumption
      use rxns, only: SUM_R_s
! Run time logical for using stiff chemistry solver
      use stiff_chem, only: STIFF_CHEMISTRY
! Run time logicals for identifying cells owned by myPE
      use stiff_chem, only: notOwner

! Legacy Global Variables:
!---------------------------------------------------------------------//
! Run time logical for using ISAT (legacy)
      use stiff_chem, only: CALL_ISAT
! Run time logical for using direct integration (legacy)
      use stiff_chem, only: CALL_DI
! Time step for using ISAT (legacy)
      use stiff_chem, only: ISATDT

! Full access to the following modules:
!---------------------------------------------------------------------//
      use compar
      use geometry
      use indices

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK

      include 'function.inc'

! Error - ISAT is no longer available.
      IF(CALL_ISAT) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,1001); WRITE(*,1000)
            WRITE(UNIT_LOG,1001); WRITE(UNIT_LOG,1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      END IF

      IF(CALL_DI .AND. .NOT.STIFF_CHEMISTRY) THEN
         IF(myPE == PE_IO) then
            WRITE(*,1003); WRITE(*,1000)
            WRITE(UNIT_LOG,1003); WRITE(UNIT_LOG,1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

! If the stiff solver is not being used, there is no need for the
! following checks.
      IF(.NOT.STIFF_CHEMISTRY) RETURN

! Message - ISATDT is no longer used.
      IF(ISATdt /= UNDEFINED) WRITE(UNIT_LOG,1002)

! Energy equations must be solved.
      IF(.NOT.ENERGY_EQ) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,1004)'ENERGY_EQ = .FALSE.'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1004)'ENERGY_EQ = .FALSE.'
            WRITE(UNIT_LOG,1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF


      IF(RO_G0 /= UNDEFINED) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,1003)'RO_G0 /= UNDEFINED'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1003)'RO_G0 /= UNDEFINED'
            WRITE(UNIT_LOG,1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(.NOT.SPECIES_EQ(0)) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,1003)'SPECIES_EQ(0) = .FALSE.'
            WRITE(*,1000)
            WRITE(UNIT_LOG,1003)'SPECIES_EQ(0) = .FALSE.'
            WRITE(UNIT_LOG,1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

! The stiff chemistry solver only needs to loop over physical cells
! owned by a process (e.g., not ghost cells). To avoid having a 
! triple do loop, this array is populated to identify the cells that
! are not owned.
      ALLOCATE( notOwner(DIMENSION_3) ); notOwner = .TRUE.
      do k=kstart, kend
      do j=jstart, jend
      do i=istart, iend
         ijk = funijk(i,j,k)
         notOwner(IJK) = .FALSE.
      enddo
      enddo
      enddo



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

      use physprop,       only : MMAX, NMAX
      use run,            only : SOLID_RO_V, SPECIES_EQ
      use stiff_chem_dbg, only : INIT_ODE_STATS0

      use stiff_chem

      implicit none

      INTEGER :: LRN, LRS
      INTEGER :: M, N_VAR


      NEQ_DIMN = 2 + MMAX

! Dimension of ODEs solved if only the gas phase is present:
! Gas density, temperature, and species
      ODE_DIMN_g = 2 + NMAX(0)
! Solids temperature excluding for all phases.
      ODE_DIMN_g = ODE_DIMN_g + MMAX


! Calculate the total number of ODEs that are solve.
      ODE_DIMN_all = ODE_DIMN_g
      DO M=1, MMAX
! Solids bulk density and species.
         IF(SPECIES_EQ(M)) ODE_DIMN_all = ODE_DIMN_all + (1 + NMAX(M))
      ENDDO

! Indicates type of Error control.
      ODE_ITOL = 2 ! :: EWT(i) = RTOL * ABS(Y(i)) * ATOL(i)

! Relative error tolerance paramter.
      ODE_RTOL(1) = 1.0D-5

! Absolue error tolerance parameter.
      IF(.NOT.(allocated(ODE_ATOL))) allocate(ODE_ATOL(ODE_DIMN_all))
      ODE_ATOL(:) = 1.0D-6  ! All Equations

! Declared length of RWORK.
      LRN = 20 + 16*ODE_DIMN_all
      LRS = 22 + 9*ODE_DIMN_all + (ODE_DIMN_all**2)
      ODE_LRW = max(LRN, LRS)

! Declared length of IWORK.
      ODE_LIW = 20 + ODE_DIMN_all

! Jacobian type indicator.
      ODE_JT = 2 ! Internally generated.

!
      CALL INIT_ODE_STATS0(ODE_DIMN_all)

      return
      END SUBROUTINE ODEPACK_INIT
