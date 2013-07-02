!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  MODULE: COEFF                                                       !
!  Purpose: Contains logic flags that tells the code whether to        !
!           perform the indicated type of calculation when the         !
!           value is true                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE coeff

! Flags used by PHYSICAL_PROP :: (0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: DENSITY(:)  ! Density
      LOGICAL, ALLOCATABLE :: SP_HEAT(:)  ! Specific heat 
      LOGICAL, ALLOCATABLE :: PSIZE(:)    ! Particle diameter


! Flags used by TRANSPORT_PROP :: (0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: VISC(:)      ! Viscosity
      LOGICAL, ALLOCATABLE :: COND(:)      ! Conductivity
      LOGICAL, ALLOCATABLE :: DIFF(:)      ! Diffusivity
      LOGICAL, ALLOCATABLE :: GRAN_DISS(:) ! Granular energy dissipation


! Flags used by EXCHANGE :: (0:DIMENSION_M)x(0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: DRAGCOEF(:,:) ! Drag coefficient
      LOGICAL, ALLOCATABLE :: HEAT_TR(:,:)  ! Heat transfer coeff


      contains

!**********************************************************************!
!  SUBROUTINE: INIT_COEFF                                              !
!                                                                      !
!  Purpose: Initialize logical flags.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_COEFF(IER)

! Global Variables:
!-----------------------------------------------------------------------
! Number of solids phases.
      use param, only: DIMENSION_M
! Kinetic theory model.
      use run, only: KT_TYPE
! Run-time flag for invoking DQMOM
      use run, only: CALL_DQMOM
! Real number of solids phases (GHD theory)
      use physprop, only: SMAX
! Run-time flag for to solve energy equations.
      use run, only: ENERGY_EQ
! Run-time flag for to solve species equations.
      use run, only: SPECIES_EQ
! Flag to recalculate gas viscosity.
      use visc_g, only: RECALC_VISC_G
! Run-time flag for invoking discrete element model
      use discretelement, only: DISCRETE_ELEMENT
! Run-time flag for invoking TFM/DEM hybrid model
      use discretelement, only: DES_CONTINUUM_HYBRID
! Run-time flag invoking QMOM theory
      use qmom_kinetic_equation, only: QMOMK
! Specified constant gas phase density (incompressible)
      use physprop, only: RO_G0
! Specified constant specific heat.
      use physprop, only: C_PG0, C_PS0
! Specified constant thermal conductivity.
      use physprop, only: K_G0,  K_S0
! Specified number of solids phases.
      use physprop, only: MMAX
! Number to determined unspecified double percision entries.
      use param1, only: UNDEFINED

      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Error flag.
      INTEGER, intent(inout) :: IER

! Invoke debug routine:
!-----------------------------------------------------------------------
      LOGICAL, parameter :: dbg_coeffs = .FALSE.


! Allocate and initialize:
!```````````````````````````````````````````````````````````````````````
      allocate( DENSITY(0:DIMENSION_M)); DENSITY = .FALSE.
      allocate( SP_HEAT(0:DIMENSION_M)); SP_HEAT = .FALSE.
      allocate( PSIZE(0:DIMENSION_M)); PSIZE   = .FALSE.

      allocate( VISC(0:DIMENSION_M)); VISC = .FALSE.
      allocate( COND(0:DIMENSION_M)); COND = .FALSE.
      allocate( DIFF(0:DIMENSION_M)); DIFF = .FALSE.
      allocate( GRAN_DISS(0:DIMENSION_M)); GRAN_DISS = .FALSE.

      allocate( DRAGCOEF(0:DIMENSION_M,0:DIMENSION_M))
      DRAGCOEF = .FALSE.

      allocate( HEAT_TR(0:DIMENSION_M,0:DIMENSION_M))
      HEAT_TR = .FALSE.


! Coefficients for gas phase parameters.
!```````````````````````````````````````````````````````````````````````
! Compressible flow.
      if(RO_G0 == UNDEFINED) DENSITY(0) = .TRUE.
! Viscosity is recalculated iteration-to-iteration if:
! 1) the energy equations are solved
! 2) a turbulace length scale is defined (L_SCALE0 /= ZERO)
! 3) K-Epsilon model is used.
      VISC(0) = RECALC_VISC_G
! Specific heat and thermal conductivity.
      if(ENERGY_EQ) then
         if(C_PG0 == UNDEFINED) SP_HEAT(0) = .TRUE. 
         if(K_G0  == UNDEFINED) COND(0) = .TRUE. 
      endif
! Species diffusivity.
      if(SPECIES_EQ(0)) DIFF(0) = .TRUE.


! Interphase transfer terms.
!```````````````````````````````````````````````````````````````````````
      if(.NOT.QMOMK) DRAGCOEF(0:MMAX,0:MMAX) = .TRUE. 

 
! Coefficients for solids phase parameters.
!```````````````````````````````````````````````````````````````````````
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
! Interphase heat transfer coefficient (GAMA)
         if(ENERGY_EQ) HEAT_TR(0:MMAX,0:MMAX) = .TRUE.

! Solids viscosity. CALC_MU_s must be invoked every iteration, even if
! MU_s0 /= UNDEFINED, so that initialization of global variables occurs.
         VISC(1:MMAX) = .TRUE. 

! Specific heat and thermal conductivity.
         if(ENERGY_EQ) THEN
            if(C_PS0 == UNDEFINED) SP_HEAT(1:MMAX) = .TRUE. 
            if(K_S0  == UNDEFINED) COND(1:MMAX) = .TRUE. 
         endif

! Species diffusivity. There is no reason to invoke this routine as the
! diffusion coefficient for solids is always zero.
         DIFF(1:MMAX) = .FALSE.

! Particle-Particle Energy Dissipation
         IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP' .OR. &
             TRIM(KT_TYPE) .EQ. 'GD_99') THEN
            GRAN_DISS(:MMAX) = .TRUE.
         ENDIF

! Particle diameter.
         if(Call_DQMOM) PSIZE(1:SMAX)=.TRUE.

      ENDIF   ! end if (.not.discrete_element .or des_continuum_hybrid)

      if(dbg_coeffs) CALL DEBUG_COEFF

      END SUBROUTINE INIT_COEFF

!**********************************************************************!
!  SUBROUTINE: DEBUG_COEFF                                             !
!                                                                      !
!  Purpose: Dump the coefficient arrays for debugging.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEBUG_COEFF

      use compar
      use physprop, only: MMAX

      implicit none

      INTEGER :: M, MM

      if(myPE /= PE_IO) return

      write(*,"(/3x,'From DEBUG_COEFF:')")

      write(*,"(/3x,'Gas phase coefficients:')")
      write(*,"( 5x,'Density (RO_g):',1x,1L)") DENSITY(0)
      write(*,"( 5x,'Specific heat (C_pg):',1x,1L)") SP_HEAT(0)
      write(*,"( 5x,'Viscosity: (MU_g)',1x,1L)") VISC(0)
      write(*,"( 5x,'Thermal conductivity (K_g):',1x,1L)") COND(0)
      write(*,"( 5x,'Species diffusivity: (DIF_G)',1x,1L)") DIFF(0)


      DO M=1, MMAX
         write(*,"(/3x,'Solids ',I1,' phase coefficients:')") M
         write(*,"( 5x,'Density: (RO_s)',1x,1L)") DENSITY(M)
         write(*,"( 5x,'Specific heat (C_ps):',1x,1L)") SP_HEAT(M)
         write(*,"( 5x,'Viscosity (MU_s):',1x,1L)") VISC(M)
         write(*,"( 5x,'Thermal conductivity (K_s):',1x,1L)") COND(M)
         write(*,"( 5x,'Species diffusivity (DIF_s):',1x,1L)") DIFF(M)
         write(*,"( 5x,'Gran. Dissipation (D_p):',1x,1L)") GRAN_DISS(M)
         write(*,"( 5x,'Diameter (D_p):',1x,1L)") PSIZE(M)
      ENDDO


      write(*,"(/3x,'Interphase drag:')")
      write(*,"( 5x,'ref',$)")
      DO M=0, MMAX
         write(*,"(2x,I3,$)")M
      ENDDO
      write(*,"('')")

      DO M=0, MMAX
         write(*,"( 5x,I3,$)") M
         DO MM=0, MMAX
            write(*,"(2x,L3,$)")DRAGCOEF(M, MM)
         ENDDO
         write(*,"('')")
      ENDDO

      write(*,"(/3x,'Interphase heat transfer:')")
      write(*,"( 5x,'ref',$)")
      DO M=0, MMAX
         write(*,"(2x,I3,$)")M
      ENDDO
      write(*,"('')")
      DO M=0, MMAX
         write(*,"( 5x,I3,$)") M
         DO MM=0, MMAX
            write(*,"(2x,L3,$)")HEAT_TR(M, MM)
         ENDDO
         write(*,"('')")
      ENDDO

      write(*,"(/3x,'DEBUG_COEFF - Exit',3/)")

      END SUBROUTINE DEBUG_COEFF

      END MODULE coeff
