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


! Legacy Variables:
!---------------------------------------------------------------------//
! Former keyword for invoking stiff solver.
      LOGICAL :: CALL_DI
! Keyword for using ISAT tables with stiff solver. (disabled)
      LOGICAL :: CALL_ISAT
! Time step for isat calculation. (disabled)
      DOUBLE PRECISION :: ISATdt


      END MODULE STIFF_CHEM
