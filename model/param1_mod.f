      MODULE param1


      Use param


!
! File: param1.inc
! 1.
!
!                      Number of single precision .SPx files
      INTEGER          N_SPX
      PARAMETER        (N_SPX=9)
!
!                      Maximum number of cell classes
      INTEGER          MAX_CLASS
      PARAMETER (MAX_CLASS = 1000)
!
!                      Maximum number of corner cells
      INTEGER          MAX_NCORN
      PARAMETER (MAX_NCORN = 1000)
!
!
! 2.  Parameter calculation section.  This contains parameters calculated
!     from those defined in Section 1.
!
!                       Maximum dimension for the upper triangle of
!                       an MxM matrix
      INTEGER           DIMENSION_LM
!      PARAMETER&
!        (DIMENSION_LM = (DIMENSION_M * (DIMENSION_M-1) / 2)+1 )
!
!                      Total number of species
      INTEGER          DIMENSION_N_all
!      PARAMETER        (DIMENSION_N_all = DIMENSION_N_g&
!                        + DIMENSION_M * DIMENSION_N_s )
!
! 3.  Parameter definition section.  This contains parameters defined
!       for programming convenience.
!
!
!                      A large number
      DOUBLE PRECISION LARGE_NUMBER
      PARAMETER (LARGE_NUMBER = 1.0D32)
!
!                      A real number to indicate that a certain variable is
!                      not initialized
      DOUBLE PRECISION UNDEFINED
      PARAMETER (UNDEFINED = 9.87654321D31)
!
!                      An integer to indicate that a certain variable is
!                      not initialized
      INTEGER          UNDEFINED_I
      PARAMETER (UNDEFINED_I = 987654321 )
!
!                      A character to indicate that a certain variable is
!                      not initialized
      CHARACTER        UNDEFINED_C
      PARAMETER (UNDEFINED_C = ' ')
!
!                      A small number
      DOUBLE PRECISION SMALL_NUMBER
      PARAMETER (SMALL_NUMBER = 1.0D-15)
!
!                      Zero, nil, zilch, naught, nought, nix, ...
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0d0)
!
!                      One
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0d0)
!
!                      1/2
      DOUBLE PRECISION HALF
      PARAMETER (HALF = 0.5d0)


      END MODULE param1                                                                          
