!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS_DEM                            !
!                                                                      !
!  Purpose: check the initial conditions input section for DEM model   !
!     
!     - calculate the number of particles needed to initalize the      !
!        DEM model                                                     !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!  Comments: Most of the code in this routine is a consolidation of    !
!            codes from several existing routines                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_INITIAL_CONDITIONS_DEM
      

! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only : gener_part_config
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s

      USE discretelement, only: DES_CONTINUUM_HYBRID

! Computed volume of IC region for seeding
      USE discretelement, only: VOL_IC_REGION
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles to read from input file.
      USE discretelement, only: PARTICLES


! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! IC Region gas volume fraction.
      USE ic, only: IC_THETA_M 
      
      USE ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE
! Maximum number of IC regions and solids phases
      USE param, only: DIMENSION_IC, DIM_M 
! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: xlength, ylength, zlength, dx, dy, dz 
! Constant: 3.14159...
      USE constant, only: PI

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

      INTEGER :: ICV, ICV2, M, IDIM,  I,J,K
      INTEGER :: COUNT_IC, COUNT_IC_WITH_SOLS
      INTEGER :: FIRST_DEF_IC 

      DOUBLE PRECISION :: IC_ORIG(3), IC_END(3), IC2_ORIG(3) , IC2_END(3)
      DOUBLE PRECISION :: IC_MIN, IC_MAX, IC2_MIN, IC2_MAX , TOL_IC_REG
      
      LOGICAL :: SEP_AXIS, first_ic_ok

      DOUBLE PRECISION :: EP_SM
! Temp variable for storing number of particles.
      DOUBLE PRECISION :: PARTS_TEMP 
!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL, EXTERNAL :: COMPARE 


      IF (.NOT.GENER_PART_CONFIG) RETURN 

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_INITIAL_CONDITIONS_DEM")

      IF (DES_CONTINUUM_HYBRID) THEN
         WRITE(ERR_MSG, 999) 
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF     

 999  format('Error # 999: Gener_part_config set to', & 
      ' true for DES_continuum hybrid', /, &
      ' This is not allowed, specify the initial particle', &
      ' configuration explicitly', /, &
      ' See MFIX readme', /,  &
      ' Please correct the data file.')

      
! Determine the domain volume which is used to calculate the total
! number of particles and the number of particles in each phase. 
! Values of DZ(1) or zlength are guaranteed at this point due to
! check_data_03. If the user left both undefined and NO_K = .T., then
! they are set to ONE. If dz(1) is undefined but zlength is defined,
! then dz(1) is set to zlength (and vice versa).  If both are defined
! they must be equal.
      IF(DIMN.EQ.2) THEN 
         IF (DES_MMAX.EQ.1) THEN
! Warn the user if the domain depth is not equal to the particle
! diameter as it may cause problems for coupled simulations.
! The user should also be aware of this when interpreting
! volume/void fraction calculations (including bulk density).
            IF(.NOT.COMPARE(ZLENGTH,DES_D_P0(1))) THEN
               WRITE(ERR_MSG, 1000)
               CALL FLUSH_ERR_MSG
            ENDIF
         ELSE
! Let the user know basis of depth dimension for calculating number of
! particles. this will also be important when considering volume/void
! fraction calculations.
            WRITE(ERR_MSG, 1001)
            CALL FLUSH_ERR_MSG
         ENDIF
      ENDIF


 1000 FORMAT(' Message: ',&
      'WARNING: zlength or dz(1) is used to calculate the ',&
      'number of particles in the 2D simulation when ',&
      'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
      'does not equal D_P0(1).') 

 1001 FORMAT(' Message: ',&
      'WARNING: zlength or dz(1) is used to calculate the ',&
      'number of particles in the 2D simulation when ',&
      'GENER_PART_CONFIG is T and DIMN = 2.') 



! Now compute the total number of particles that will be initialized 
! in each IC region. PARTICLES calculated here is used to allocate
! DEM arrays. 

      IF(.NOT.allocated(vol_ic_region)) & 
      ALLOCATE(VOL_IC_REGION(DIMENSION_IC))

      if(.not.allocated(part_mphase_byic)) & 
      ALLOCATE(PART_MPHASE_BYIC(DIMENSION_IC, DIM_M))    

      PARTICLES = 0
      COUNT_IC = 0 

      DO ICV = 1, DIMENSION_IC 
         VOL_IC_REGION(ICV) = zero 
         PART_MPHASE_BYIC(ICV,:) = 0
         
         IF(.not.ic_defined(icv)) cycle         

         DO K = IC_K_B(ICV), IC_K_T(ICV) 
            DO J = IC_J_S(ICV), IC_J_N(ICV) 
               DO I = IC_I_W(ICV), IC_I_E(ICV) 
                  VOL_IC_REGION(ICV) = VOL_IC_REGION(ICV)+DX(I)*DY(J)*DZ(K)
               ENDDO 
            ENDDO
         ENDDO

         IF (IC_EP_G(ICV).lt.ONE) THEN 
            COUNT_IC = COUNT_IC + 1 
            DO M = 1, DES_MMAX
               EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
               PARTS_TEMP = FLOOR((6.D0*EP_SM*VOL_IC_REGION(ICV))/&
               (PI*(DES_D_P0(M)**3)))
               PART_MPHASE_BYIC(ICV, M) = PARTS_TEMP
               PARTICLES = PARTICLES + PARTS_TEMP
            ENDDO
         ENDIF
      ENDDO
    

! Local reporting for the user 
      WRITE(ERR_MSG, 2015)  COUNT_IC, DES_MMAX
      CALL FLUSH_ERR_MSG(FOOTER = .false.)

      DO ICV = 1, DIMENSION_IC 
         IF(.not.ic_defined(icv)) cycle 

         IF (IC_EP_G(ICV).lt.ONE) THEN 
            
            WRITE(ERR_MSG,2016) ICV, VOL_IC_REGION(ICV)
            CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)

            DO M = 1, DES_MMAX
               EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
               WRITE(ERR_MSG,2017) M,  &
               DES_D_P0(M), EP_SM, PART_MPHASE_BYIC(ICV, M) 
               CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)
            ENDDO

         ENDIF

      ENDDO

      WRITE(ERR_MSG,2018)
      CALL FLUSH_ERR_MSG(HEADER = .false.)


 2015 FORMAT('Gener_Part_Config specified as true', /,5x, &
      'Total IC region count with non zero solids = ', I5,2X, /5x, & 
      'Total discrete solid phases = ', I5,2X, /5x, & 
      'Particle configuration will be automatically generated: ')

 2016 FORMAT('IC number and Volume        = ', I4, 2x, (ES15.7,2X) )
      
 2017 FORMAT( &
      'PHASE Index, M              =  ', I5,2X, /5x, & 
      'D_P0(M)                     =  ', (ES15.7,2X), /5x, &
      'Solids vol fraction for M   =  ', (G15.8,2X), /5x, & 
      '# of particles for phase M  = ',  (I10,2X))
      
 2018 FORMAT( 'End of Message' ) 


      CALL FINL_ERR_MSG
      
      END SUBROUTINE CHECK_INITIAL_CONDITIONS_DEM
      
