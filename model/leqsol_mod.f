      MODULE leqsol


      Use param
      Use param1


!     leqsol.inc
!     1   p
!     2   rho, ep
!     3   u
!     4   v
!     5   w
!     6   T
!     7   X
!     8   Th
!     9   S
!
!
!                      Maximum number of outer iterations
      INTEGER          MAX_NIT
!
!                      Automatic adjustment of leq parameters possible (set in
!                      iterate after the completion of first iteration)
      LOGICAL          LEQ_ADJUST
!
!                      Maximum number of linear equation solver iterations
      INTEGER          LEQ_IT(9)
!
!                      linear equation solver method
      INTEGER          LEQ_METHOD(9)
!
!                      linear equation solver sweep direction
      CHARACTER*4 ::   LEQ_SWEEP(9)
!
!                      linear equation solver tolerance
      DOUBLE PRECISION LEQ_TOL(9)
!
!		       Preconditioner option
      CHARACTER*4 ::   LEQ_PC
!
!      COMMON / ITERS_DP /
!     &
!
!


      END MODULE leqsol                                                                          
