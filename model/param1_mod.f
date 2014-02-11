      MODULE param1


      Use param


!
! File: param1.inc
! 1.
!
!                      Number of single precision .SPx files
      INTEGER, PARAMETER          :: N_SPX = 11
!
!                      Maximum number of cell classes, increased to 1 Million due to re-indexing
      INTEGER, PARAMETER          :: MAX_CLASS = 1000000
!
!                      Maximum number of corner cells
      INTEGER, PARAMETER          :: MAX_NCORN = 4000
!
!
! 2.  Parameter calculation section.  This contains parameters calculated
!     from those defined in Section 1.
!
!                       Maximum dimension for the upper triangle of
!                       an MxM matrix
      INTEGER           DIMENSION_LM
!
!                      Total number of species
      INTEGER          DIMENSION_N_all
!
! 3.  Parameter definition section.  This contains parameters defined
!       for programming convenience.
!
!
!                      A large number
      DOUBLE PRECISION, PARAMETER          :: LARGE_NUMBER = 1.0D32
!
!                      A real number to indicate that a certain variable is
!                      not initialized
      DOUBLE PRECISION, PARAMETER          :: UNDEFINED = 9.87654321D31
!
!                      An integer to indicate that a certain variable is
!                      not initialized
      INTEGER, PARAMETER          :: UNDEFINED_I = 987654321
!
!                      A character to indicate that a certain variable is
!                      not initialized
      CHARACTER, PARAMETER          :: UNDEFINED_C = ' '
!
!                      A small number
      DOUBLE PRECISION, PARAMETER          :: SMALL_NUMBER = 1.0D-15
!
!                      Zero, nil, zilch, naught, nought, nix, ...
      DOUBLE PRECISION, PARAMETER          :: ZERO = 0.0d0
!
!                      One
      DOUBLE PRECISION, PARAMETER          :: ONE = 1.0d0
!
!                      1/2
      DOUBLE PRECISION, PARAMETER          :: HALF = 0.5d0


      END MODULE param1                                                                          
