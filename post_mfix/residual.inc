C
C     residual.inc
C
      INTEGER          MAX_RESID_INDEX
      PARAMETER        (MAX_RESID_INDEX = 8)    !for printing. do not change
C
      INTEGER          NRESID
      INTEGER          NPREFIX
      INTEGER          RESID_p, RESID_ro, RESID_u, RESID_v, RESID_w,
     &                 RESID_t, RESID_x, RESID_th
      PARAMETER        (RESID_p  = 1)     !pressure
      PARAMETER        (RESID_ro = 2)     !density, volume fraction
      PARAMETER        (RESID_u  = 3)     !u-velocity
      PARAMETER        (RESID_v  = 4)     !v-velocity
      PARAMETER        (RESID_w  = 5)     !w-velocity
      PARAMETER        (RESID_t  = 6)     !temperature
      PARAMETER        (RESID_th = 7)     !granular temperature
      PARAMETER        (NRESID   = 7 +MAX(DIMENSION_N_g, DIMENSION_N_s))
      PARAMETER        (RESID_x  = 8)     !mass fraction
      PARAMETER        (NPREFIX  = 8)
C
C                      prefix of Residuals string
      CHARACTER        RESID_PREFIX(NPREFIX)

C
C                      Average residual
      DOUBLE PRECISION RESID(NRESID, 0:DIMENSION_M)
C
C                      Maximum residual
      DOUBLE PRECISION MAX_RESID(NRESID, 0:DIMENSION_M) 
C
C                      sum of residuals every 5 iterations
      DOUBLE PRECISION SUM5_RESID 
C
C                      IJK location of maximum residual
      INTEGER          IJK_RESID(NRESID, 0:DIMENSION_M) 

C
C                      Residuals to be printed out
      CHARACTER*4      RESID_STRING(MAX_RESID_INDEX)
C
C                      Indices of residuals to be printed out
      INTEGER          RESID_INDEX(MAX_RESID_INDEX, 2) 
C
      COMMON /RESIDUAL_DP/
     &                  RESID, MAX_RESID, SUM5_RESID
C
      COMMON /RESIDUAL_I/
     &                   IJK_RESID, RESID_INDEX
C
      COMMON /RESIDUAL_C/
     &                   RESID_STRING, RESID_PREFIX
