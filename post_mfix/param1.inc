C
C File: param1.inc
C 1.
C
C                      Number of single precision .SPx files
      INTEGER          N_SPX
      PARAMETER        (N_SPX=8)
C
C                      Maximum number of cell classes
      INTEGER          MAX_CLASS
      PARAMETER (MAX_CLASS = 1000)
C
C
C 2.  Parameter calculation section.  This contains parameters calculated
C     from those defined in Section 1.
C
C                       Maximum dimension for the upper triangle of
C                       an MxM matrix
      INTEGER           DIMENSION_LM
      PARAMETER 
     &  (DIMENSION_LM = (DIMENSION_M * (DIMENSION_M-1) / 2)+1 ) 
C
C                      Total number of species
      INTEGER          DIMENSION_N_all
      PARAMETER        (DIMENSION_N_all = DIMENSION_N_g 
     &                  + DIMENSION_M * DIMENSION_N_s )
C
C 3.  Parameter definition section.  This contains parameters defined
C       for programming convenience.
C
C
C                      A large number
      DOUBLE PRECISION LARGE_NUMBER
      PARAMETER (LARGE_NUMBER = 1.0D32)
C
C                      A real number to indicate that a certain variable is
C                      not initialized
      DOUBLE PRECISION UNDEFINED
      PARAMETER (UNDEFINED = 9.87654321D31)
C
C                      An integer to indicate that a certain variable is
C                      not initialized
      INTEGER          UNDEFINED_I
      PARAMETER (UNDEFINED_I = 987654321 )
C
C                      A character to indicate that a certain variable is
C                      not initialized
      CHARACTER        UNDEFINED_C
      PARAMETER (UNDEFINED_C = ' ')
C
C                      A small number
      DOUBLE PRECISION SMALL_NUMBER
      PARAMETER (SMALL_NUMBER = 1.0D-15)
C
C                      Zero, nil, zilch, naught, nought, nix, ...
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0)
C
C                      One
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.0)
C
C                      1/2
      DOUBLE PRECISION HALF
      PARAMETER (HALF = 0.5)
