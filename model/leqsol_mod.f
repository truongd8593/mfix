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

      CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE REPORT_SOLVER_STATS(TNIT, STEPS)

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: TNIT, STEPS

      INTEGER :: LC

      WRITE(ERR_MSG,1100) iVal(TNIT), iVal(TNIT/STEPS)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1100 FORMAT(/2x,'Total number of non-linear iterations: ', A,/2x,&
         'Average number per time-step: ',A)

      WRITE(ERR_MSG,1200)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1200 FORMAT(2x,'|',10('-'),'|',13('-'),'|',14('-'),'|',/&
         2x,'| Equation |  Number of  |  Avg Solves  |',/&
         2x,'|  Number  |   Solves    |   for NIT    |',/&
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')

      DO LC = 1, DIM_EQS
         WRITE(ERR_MSG,1201) LC, ITER_TOT(LC), ITER_TOT(LC)/TNIT
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

 1201 FORMAT(2x,'|',3x,I3,4x,'|',2x,I9,2x,'|',2x,I10,2x,'|',/ &
         2x,'|',10('-'),'|',13('-'),'|',14('-'),'|')


      RETURN
      END SUBROUTINE REPORT_SOLVER_STATS


      END MODULE leqsol
