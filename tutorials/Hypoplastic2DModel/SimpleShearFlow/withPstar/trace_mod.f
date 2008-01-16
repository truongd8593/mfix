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
!
!     JEG Added 
!     University of Colorado, Hrenya Research Group
!                      Trace of the dot of D_s (M solid phase) and
!                      D_sl (L solid phase)
      DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: trD_s2_ip
!     END JEG
 
!  Store components of solids strain rate tensor, note that D_sij = D_sji
!  SOF 
!	               D_s(1,1)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s11 
 
!	               D_s(2,2)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s22 
 
!	               D_s(3,3)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s33 
 
!	               D_s(1,2)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s12 
 
!	               D_s(1,3)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s13 
 
!	               D_s(2,3)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  D_s23 
 
!  computes vorticity components of phase m 1/2 (delVs - delVsTranspose)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  vort_s12 
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  vort_s13  
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  vort_s23 
!  End SOF 
 
 
 
!!!HPF$ align trD_s_C(:, *) with TT(:)
!!!HPF$ align trD_s2(:, *) with TT(:)
!!!HPF$ align trD_s_Co(:, *) with TT(:)

      END MODULE trace
