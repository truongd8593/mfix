      MODULE tau_g

      Use param
      Use param1

!
!                      cross terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  TAU_W_g

!!!HPF$ align TAU_U_g(:) with TT(:)
!!!HPF$ align TAU_V_g(:) with TT(:)
!!!HPF$ align TAU_W_g(:) with TT(:)


      END MODULE tau_g
