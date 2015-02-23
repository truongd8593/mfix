      MODULE leqsol

      use param, only: DIM_EQS

! Maximum number of outer iterations
      INTEGER :: MAX_NIT

! Automatic adjustment of leq parameters possible (set in iterate after
! the completion of first iteration).
      LOGICAL :: LEQ_ADJUST

! Maximum number of linear equation solver iterations
      INTEGER :: LEQ_IT(DIM_EQS)

! Linear equation solver method
      INTEGER :: LEQ_METHOD(DIM_EQS)

! Total Iterations
      INTEGER :: ITER_TOT(DIM_EQS) = 0

! Linear equation solver sweep direction
      CHARACTER(LEN=4) :: LEQ_SWEEP(DIM_EQS)

! Linear equation solver tolerance
      DOUBLE PRECISION :: LEQ_TOL(DIM_EQS)

! Preconditioner option
      CHARACTER(LEN=4) :: LEQ_PC(DIM_EQS)

! Option to minimize dot products
      LOGICAL :: MINIMIZE_DOTPRODUCTS

! Option to transpose A_m
      LOGICAL :: DO_TRANSPOSE

! Frequency of convergence check in BiCGStab
      INTEGER :: ICHECK_BICGS

! Optimize for massively parallel machine
      LOGICAL :: OPT_PARALLEL

! Linear and non-linear solver statistics
      LOGICAL :: SOLVER_STATISTICS

      END MODULE leqsol
