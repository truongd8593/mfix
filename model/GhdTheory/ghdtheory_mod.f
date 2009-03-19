  MODULE ghdtheory
 
 
      Use param
      Use param1
 
!
!     Zeroth order dissipation term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zeta0
!
!     cooling rate transport coefficient (1st order)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZetaU
!
!     Thermal diffusivity DiT
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DiT
!
!     Mass mobility coefficient
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: DijF
!
!     Thermal mobility
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Lij
!
!     Ordinary diffusion
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Dij
!
!     Dufour coefficient
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: DijQ
!
!     Species mass flux in X-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiX
!
!     Species mass flux in Y-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiY
!
!     Species mass flux in Z-direction
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: JoiZ


      END MODULE ghdtheory
