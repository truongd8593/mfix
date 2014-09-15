!----------------------------------------------------------------------!
!                                                                      !
!         --->>> THIS IS A DUMMY ROUTINE FOR POST_MFIX <<<---          !
!                                                                      !
!----------------------------------------------------------------------!
      MODULE STIFF_CHEM

! Runtime Flags:
!---------------------------------------------------------------------//
! Flag to invoke stiff chemistry solver.
      LOGICAL :: STIFF_CHEMISTRY
! Flag to invoke the variable solids diameter model.
      LOGICAL :: CALL_GROW
! Flag indicating if cell IJK is own by myPE.
      LOGICAL, dimension(:), allocatable :: notOwner

! ODEPACK Controlling parameters:
!---------------------------------------------------------------------//
! Dimension of ODEs solved in ISAT or DI
      INTEGER :: ODE_DIMN
! Indicates type of Error control.
      INTEGER :: ODE_ITOL
! Relative error tolerance paramter.
      DOUBLE PRECISION, DIMENSION(1) :: ODE_RTOL
! Absolue error tolerance parameter. (Dimension (ODE_DIMN))
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ODE_ATOL
! Declared length of RWORK.
      INTEGER :: ODE_LRW
! Declared length of IWORK.
      INTEGER :: ODE_LIW
! Jacobian type indicator.
      INTEGER :: ODE_JT
! The maximum number of steps ODEPACK may use to integrate.
      INTEGER :: STIFF_CHEM_MAX_STEPS


      END MODULE STIFF_CHEM
