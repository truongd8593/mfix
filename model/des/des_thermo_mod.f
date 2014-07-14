!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO                                             !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM heat transfer.                !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_THERMO

      USE param

! Run time logicals
!--------------------------
! Slove DES Energy Equations
!      LOGICAL DES_ENERGY_EQ

! Specifies that the convection model is to be solved.
! Requires DES_ENERGY_EQ = .TRUE.
! Default [.TRUE.]
!      LOGICAL DES_CONV_EQ
! Specifies that the particle conduction models are to be solved.
! Requires DES_ENERGY_EQ = .TRUE.
! Default [.TRUE.]
!      LOGICAL DES_COND_EQ
! These logicals can be set in the mfix.dat file to omit specific
! particle-particle heat conduction models from being solved.
! If the DEM conduction model is not solved (DES_COND_EQ), the values
! are automatically set to .FALSE. in check_des_thermo.f.
! Default [.TRUE.]
!      LOGICAL DES_COND_EQ_PFP ! particle-fluid-particle
! Default [.TRUE.]
!      LOGICAL DES_COND_EQ_PP  ! particle-particle
! Specifies that the particle-environment radiation model is to be
! solved. Requires DES_ENERGY_EQ = .TRUE.
! Default [.TRUE.]
!      LOGICAL DES_RADI_EQ

! Heat transfer correlation specified in mfix.dat
! Default [RANZ_1952]
      CHARACTER*24 :: DES_CONV_CORR

      INTEGER :: DES_CONV_CORR_ENUM
      INTEGER, PARAMETER :: RANZ_1952 = 0

! Particle properties
!--------------------------
! Particle temperature at current time step (S_TIME)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_T_s_OLD !(PARTICLES) 
! Particle temperature at previous time step (S_TIME - DT_SOLID)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_T_s_NEW !(PARTICLES) 

! DES Specified constant solids thermal conductivity by solids phase
!      DOUBLE PRECISION DES_K_s0(DIM_M)
! DES solids thermal conductivity by particle ! (not currently used)
!      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_K_s

! DES Specified constant solids specific heat by solids phase
!      DOUBLE PRECISION DES_C_ps0(DIM_M)
! DES specific heat of particles by particle
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_C_ps

! Emissivity of particles
      DOUBLE PRECISION DES_Em(DIM_M)
! Stefan-Boltzmann Constant
      DOUBLE PRECISION SB_CONST
! Bulk solids temperature for radiative heat transfer
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: avgDES_T_s


! Rates of heat transfer
!--------------------------
! Internal heat generation resulting from a chemical reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Qint !(PARTICLES)
! Heat transfer source to the particle.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q_Source
! Previous time step's rate of heat transfer. Used for Adams-Bashforth
! time integration scheme
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q_Source0

! Thermodynamic Neighborhood
!--------------------------
! Fluid Lens Proportion Constant used to calculate the radius of the
! fluid lens that surrounds the particle for particle-fluid-particle
! conduction.  Default [ 0.2 ]
      DOUBLE PRECISION FLPC
! Mininum separation distance between the surface of two contacting
! particles. This value is used to remove the singluarity that the
! particle-fluid-particle conduciton model develops at the contact
! interface. [4.0x10^(-10) meters]
      DOUBLE PRECISION DES_MIN_COND_DIST


      END MODULE DES_THERMO
