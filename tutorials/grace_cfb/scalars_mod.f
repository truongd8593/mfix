
      MODULE scalars


      Use param
      Use param1


!
!                      Number of scalar equations solved
      INTEGER          NScalar
!
!                      Index of phase associated with scalar n
      INTEGER, DIMENSION(1:DIM_Scalar) :: Phase4Scalar
!
!
!                      Source term for User-defined Scalars is linearized as
!                      S = Scalar_c - Scalar_p * Scalar
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_c 
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_p
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_c_O
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Scalar_p_O
!
!                      Diffusion coefficient for User-defined Scalars
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Dif_Scalar 

      END MODULE scalars                                                                             
