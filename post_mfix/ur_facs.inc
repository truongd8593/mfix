C     ur_facs.inc
C     1   p
C     2   rho, ep
C     3   u
C     4   v
C     5   w
C     6   T
C     7   X
C     8   Th
C
C                      Under relaxation factors
      DOUBLE PRECISION UR_FAC(8)
C
C
      INTEGER           STEPS_LAST
C
      COMMON /UR_FACS_DP/
     &                  UR_FAC
C
      COMMON /UR_FACS_I/
     &                  STEPS_LAST 
