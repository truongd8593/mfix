      MODULE tau_s
 
 
      Use param
      Use param1
 
 
!
!                      cross terms
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  TAU_U_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  TAU_V_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  TAU_W_s
 
 
!!!HPF$ align TAU_U_s(:, *) with TT(:)
!!!HPF$ align TAU_V_s(:, *) with TT(:)
!!!HPF$ align TAU_W_s(:, *) with TT(:)

 
      END MODULE tau_s
