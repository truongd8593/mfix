      MODULE residual


      Use param
      Use param1


!
!     residual.inc
!
      INTEGER          MAX_RESID_INDEX
      PARAMETER        (MAX_RESID_INDEX = 8)    !for printing. do not change
!
      INTEGER          NRESID
      INTEGER          NPREFIX
      INTEGER          RESID_p, RESID_ro, RESID_u, RESID_v, RESID_w,&
                       RESID_t, RESID_x, RESID_th
      PARAMETER        (RESID_p  = 1)     !pressure
      PARAMETER        (RESID_ro = 2)     !density, volume fraction
      PARAMETER        (RESID_u  = 3)     !u-velocity
      PARAMETER        (RESID_v  = 4)     !v-velocity
      PARAMETER        (RESID_w  = 5)     !w-velocity
      PARAMETER        (RESID_t  = 6)     !temperature
      PARAMETER        (RESID_th = 7)     !granular temperature
      PARAMETER        (NRESID   = 7 + DIM_N)
      PARAMETER        (RESID_x  = 8)     !mass fraction
      PARAMETER        (NPREFIX  = 8)
!
!                      prefix of Residuals string
      CHARACTER        RESID_PREFIX(NPREFIX)
 
!
!                      Average residual
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RESID
!
!                      Maximum residual
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  MAX_RESID
!
!                      sum of residuals every 5 iterations
      DOUBLE PRECISION SUM5_RESID
!
!                      IJK location of maximum residual
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IJK_RESID
 
!
!                      Residuals to be printed out
      CHARACTER*4      RESID_STRING(MAX_RESID_INDEX)
!
!                      Indices of residuals to be printed out
      INTEGER          RESID_INDEX(MAX_RESID_INDEX, 2)
!


      END MODULE residual                                                                        
