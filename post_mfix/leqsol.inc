C     leqsol.inc
C     1   p
C     2   rho, ep
C     3   u
C     4   v
C     5   w
C     6   T
C     7   X
C     8   Th
C
C
C                      Maximum number of outer iterations
      INTEGER          MAX_NIT
C
C                      Automatic adjustment of leq parameters possible (set in
C                      iterate after the completion of first iteration)
      LOGICAL          LEQ_ADJUST
CC
C                      Maximum number of linear equation solver iterations
      INTEGER          LEQ_IT(8)
C
C                      linear equation solver method
      INTEGER          LEQ_METHOD(8)
C
C
c      COMMON / ITERS_DP /
c     &      
C
C
      COMMON / ITERS_L /
     &      LEQ_ADJUST
C
      COMMON / ITERS_I /
     &       MAX_NIT, LEQ_IT, LEQ_METHOD
C
