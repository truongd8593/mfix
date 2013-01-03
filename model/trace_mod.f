      MODULE trace
 
 
      Use param
      Use param1
 
 
!                      Trace of D_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_C 
!
!                      Trace of the square of D_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s2 
 
!	               Trace of D_s at previous timestep
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_Co 
!QX
!	               Trace of D_s at previous time
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  trD_s_Co2
!
!     JEG Added 
!     University of Colorado, Hrenya Research Group
!                      Trace of the dot of D_s (M solid phase) and
!                      D_sl (L solid phase)
      DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: trD_s2_ip
!     END JEG
 
 
 
!!!HPF$ align trD_s_C(:, *) with TT(:)
!!!HPF$ align trD_s2(:, *) with TT(:)
!!!HPF$ align trD_s_Co(:, *) with TT(:)

      END MODULE trace
