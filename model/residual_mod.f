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
                       RESID_t, RESID_x, RESID_th, RESID_sc,RESID_ke 
      PARAMETER        (RESID_p  = 1)     !pressure
      PARAMETER        (RESID_ro = 2)     !density, volume fraction
      PARAMETER        (RESID_u  = 3)     !u-velocity
      PARAMETER        (RESID_v  = 4)     !v-velocity
      PARAMETER        (RESID_w  = 5)     !w-velocity
      PARAMETER        (RESID_t  = 6)     !temperature
      PARAMETER        (RESID_th = 7)     !granular temperature
      PARAMETER        (RESID_sc = 8)     !user-defined scalar
      PARAMETER        (NRESID   = 8 + DIM_N)
      PARAMETER        (RESID_ke = 9)     !k-epsilon equations
      PARAMETER        (RESID_x  = 10)     !mass fraction (keep this the last)
      PARAMETER        (NPREFIX  = 10)
!
!                      prefix of Residuals string
      CHARACTER, PARAMETER, DIMENSION(NPREFIX) :: RESID_PREFIX = &
        (/ 'P', 'R', 'U', 'V', 'W', 'T', 'G', 'S', 'K', 'X' /)
 
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

!                        fluid and solids accumulation, for checking the over-all fluid mass balance
        DOUBLE PRECISION accum_resid_g, accum_resid_s(DIM_M)

      END MODULE residual                                                                        
