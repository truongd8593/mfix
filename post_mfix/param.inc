C
C======Begin file param.inc ====================================================
C
C 1.  Problem size dependent parameters:
C
C                      Maximum of the number of cells in the x direction
      INTEGER          DIMENSION_I
      PARAMETER (DIMENSION_I = 502)
C
C                      Maximum of the number of cells in the y direction
      INTEGER          DIMENSION_J
      PARAMETER (DIMENSION_J = 502)
C
C                      Maximum of the number of cells in the z direction
      INTEGER          DIMENSION_K
      PARAMETER (DIMENSION_K = 502)
C
C                      Maximum number of computational cells
      INTEGER          DIMENSION_3
      PARAMETER (DIMENSION_3 =75000)
C
C                      Maximum number of solids phases
      INTEGER          DIMENSION_M
      PARAMETER (DIMENSION_M = 3)
C
C                      Maximum number of gas species
      INTEGER          DIMENSION_N_g
      PARAMETER (DIMENSION_N_g = 8)
C
C                      Maximum number of solids species. For multiple solids
C                      phase use the maximum for all phases
      INTEGER          DIMENSION_N_s
      PARAMETER (DIMENSION_N_s = 8)
C
C                      Maximum number of reactions
      INTEGER          DIMENSION_RXN
      PARAMETER        (DIMENSION_RXN = 10)
C
C===============================================================================
C
C 2.  Parameters for controlling data input:  This section seldom needs
C     to be changed.
C
C                      Number of user defined constants
      INTEGER          DIMENSION_C
      PARAMETER (DIMENSION_C = 500)
C
C                      Maximum number of items for specifying initial
C                      conditions
      INTEGER          DIMENSION_IC
      PARAMETER (DIMENSION_IC = 500)
C
C                      Maximum number of items for specifying boundary
C                      conditions
      INTEGER          DIMENSION_BC
      PARAMETER (DIMENSION_BC = 500)
C
C                      Maximum number of items for specifying internal
C                      surfaces.  Because FLAG_E, FLAG_N, and FLAG_T also
C                      points to the IS number, their manipulation should be
C                      modified if this dimension is made greater than 999.
      INTEGER          DIMENSION_IS
      PARAMETER (DIMENSION_IS = 500)
C
C                      Maximum number of user-defined output files
      INTEGER          DIMENSION_USR
      PARAMETER (DIMENSION_USR = 5)
C======End file param.inc ======================================================
C
