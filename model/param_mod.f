      MODULE param


!
!======Begin file param.inc ============================================
!
! 1.  Problem size dependent parameters:
!
!
!                      Maximum of the number of cells in the x direction
      INTEGER          DIMENSION_I
!      PARAMETER (DIMENSION_I = 9)
!
!                      Maximum of the number of cells in the y direction
      INTEGER          DIMENSION_J
!      PARAMETER (DIMENSION_J = 102)
!
!                      Maximum of the number of cells in the z direction
      INTEGER          DIMENSION_K
!      PARAMETER (DIMENSION_K = 1)
!
!                      Maximum number of computational cells
      INTEGER          DIMENSION_3
!      PARAMETER (DIMENSION_3 =918)
!//375 1025 added two new global and local subdomain size variables
      INTEGER          DIMENSION_3L, DIMENSION_3G
!
!                      Maximum number of solids phases
      INTEGER          DIMENSION_M
!      PARAMETER (DIMENSION_M = 1)
!
!                      Maximum number of gas species
      INTEGER          DIMENSION_N_g
!      PARAMETER (DIMENSION_N_g = 1)
!
!                      Maximum number of solids species. For multiple solids
!                      phase use the maximum for all phases
      INTEGER          DIMENSION_N_s
!      PARAMETER (DIMENSION_N_s = 1)
!
!===============================================================================
!
! 2.  Parameters for controlling data input:  This section seldom needs
!     to be changed.
!
!                      Maximum number of reactions
      INTEGER          DIMENSION_RXN
      PARAMETER        (DIMENSION_RXN = 100)
!
!                      Number of user defined constants
      INTEGER          DIMENSION_C
      PARAMETER (DIMENSION_C = 500)
!
!                      Maximum number of items for specifying initial
!                      conditions
      INTEGER          DIMENSION_IC
      PARAMETER (DIMENSION_IC = 500)
!
!                      Maximum number of items for specifying boundary
!                      conditions
      INTEGER          DIMENSION_BC
      PARAMETER (DIMENSION_BC = 500)
!
!                      Maximum number of items for specifying internal
!                      surfaces.  Because FLAG_E, FLAG_N, and FLAG_T also
!                      points to the IS number, their manipulation should be
!                      modified if this dimension is made greater than 999.
      INTEGER          DIMENSION_IS
      PARAMETER (DIMENSION_IS = 500)
!
!                      Maximum number of solids phases that can be read
      INTEGER          DIM_M
      PARAMETER        (DIM_M = 10)
!
!                      Maximum number of gas species  that can be read
      INTEGER          DIM_N_g
      PARAMETER        (DIM_N_g = 100)
!
!                      Maximum number of solids species. For multiple solids
!                      phase use the maximum for all phases that can be read
      INTEGER          DIM_N_s
      PARAMETER       (DIM_N_s = 100)
!
!                      maximum of DIM_N_g and DIM_N_s
      INTEGER          DIM_N
      PARAMETER       (DIM_N = MAX(DIM_N_g, DIM_N_s))
!
!                      Total number of species allowed
      INTEGER          DIM_N_ALL
      PARAMETER       (DIM_N_ALL = 2*DIM_N)
!
!                      Maximum of the number of cells in the x direction  that can be read
      INTEGER          DIM_I
      PARAMETER (DIM_I = 5000)
!
!                      Maximum of the number of cells in the y direction  that can be read
      INTEGER          DIM_J
      PARAMETER (DIM_J = 5000)
!
!                      Maximum of the number of cells in the z direction  that can be read
      INTEGER          DIM_K
      PARAMETER (DIM_K = 5000)
!!
!                      Maximum number of user-defined output files
      INTEGER          DIMENSION_USR
      PARAMETER (DIMENSION_USR = 5)
!======End file param.inc ==============================================
!


!HPF$ TEMPLATE TT(1000)
!HPF$ DISTRIBUTE TT(BLOCK)
      END MODULE param                                                                           
