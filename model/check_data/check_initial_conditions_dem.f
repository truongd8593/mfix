!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS_DEM                            !
!                                                                      !
!  Purpose: check the initial conditions input section for DEM model   !
!     - ensure the first IC is defined over the entire domain with     ! 
!        ep_g = 1 when more than one IC has solids                     !
!     - ensure the ICs are non-overlapping                             !
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
      
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
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

! First check if multiple IC regions are defined for non-zero solids volume 
! fraction, then the first IC is specified over the whole domain with IC_EP_g = 1 

      !total count of defined ICs
      COUNT_IC           = 0 
      !total count of defined IC's with solids 
      COUNT_IC_WITH_SOLS = 0
      FIRST_DEF_IC = UNDEFINED_I
      DO ICV = 1, DIMENSION_IC 
         
         IF (IC_DEFINED(ICV)) THEN 
            COUNT_IC = COUNT_IC + 1 
            FIRST_DEF_IC = MIN(FIRST_DEF_IC, ICV) 

            IF(IC_EP_G(ICV).LT.ONE) COUNT_IC_WITH_SOLS & 
            = COUNT_IC_WITH_SOLS  + 1
            
         ENDIF                  ! if(ic_defined(icv))
      end DO
      
      IF(COUNT_IC_WITH_SOLS.GE.1.AND.COUNT_IC.GE.COUNT_IC_WITH_SOLS+1) THEN
! If the number of IC's with solids is greater then one 
! and if such number of IC's are greater than total number of ICs then 
! make sure the first IC spans the entire domain with voidage of 1
! This is done so that users can specify arbitrary non-overlapping IC's with solids 
! without having to break up and specify IC for rest of the regions. 
         ICV = FIRST_DEF_IC 
         FIRST_IC_OK = .FALSE. 
         IF(IC_EP_G(ICV).EQ.ONE & 
              .AND.IC_X_W(ICV).LE.ZERO.AND.IC_X_E(ICV).GE.XLENGTH &
              .AND.IC_Y_S(ICV).LE.ZERO.AND.IC_Y_N(ICV).GE.YLENGTH) FIRST_IC_OK = .TRUE.
         !check for z also 
         IF (DIMN.EQ.3.AND.IC_Z_B(ICV).LE.ZERO.AND.IC_Z_T(ICV).GE.ZLENGTH) FIRST_IC_OK = .TRUE.
         IF(.NOT.FIRST_IC_OK) THEN 
            
            WRITE(ERR_MSG, 1003) ICV, IC_EP_G(ICV), & 
            IC_X_W(ICV), IC_X_E(ICV), &
            IC_Y_S(ICV), IC_Y_N(ICV), & 
            IC_Z_B(ICV), IC_Z_T(ICV)
            
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            
         ENDIF
      ENDIF

 1003 format('Error number 1003:',/,     &
      'For initial particle seeding with more than once IC having solids', /, &
      'first defined IC should span the entire domain and', /, & 
      'and gas voidage should be set to one.', /, & 
      'Not ensuring this will exit the code' / & 
      'IC number and IC_EP_g', I10, 2x, g17.8, /, & 
      'IC spans X direction                :', 2(2x, g17.8) , /, & 
      'IC spans Y direction                :', 2(2x, g17.8) , /, & 
      'IC spans Z direction (ignore for 2D):', 2(2x, g17.8) , /, & 
      'Please check mfix.dat file. Exitting')
      
! Now check if the ICs are non-overlapping. 

      TOL_IC_REG  = 1E-04
      ICVLOOP : DO ICV = 1, DIMENSION_IC 
         
         IF (.not.IC_DEFINED(ICV).or.IC_EP_G(ICV).eq.1.d0) cycle ICVLOOP
         IC_ORIG(1) = IC_X_W(ICV)
         IC_ORIG(2) = IC_Y_S(ICV)
         IC_ORIG(3) = IC_Z_B(ICV) 
         IC_END(1)  = IC_X_E(ICV) 
         IC_END(2)  = IC_Y_N(ICV) 
         IC_END(3)  = IC_Z_T(ICV) 
         ICVTWOLOOP : DO ICV2 = ICV+1, DIMENSION_IC 
            
            IF (.not.IC_DEFINED(ICV2).or.IC_EP_G(ICV2).eq.1.d0) cycle ICVTWOLOOP
            
            IC2_ORIG(1) = IC_X_W(ICV2)
            IC2_ORIG(2) = IC_Y_S(ICV2)
            IC2_ORIG(3) = IC_Z_B(ICV2) 
            IC2_END(1)  = IC_X_E(ICV2) 
            IC2_END(2)  = IC_Y_N(ICV2) 
            IC2_END(3)  = IC_Z_T(ICV2) 

            sep_axis  = .false. 
            DO idim = 1, dimn 

               ic_min = IC_ORIG(idim) 
               ic_max = IC_END(idim)
               ic2_min = IC2_ORIG(idim) 
               ic2_max = ic2_END(idim) 
               !Check for separating axis. If the separating axis exists, then 
               !the IC regions can't overlap 
               !generally equality implies lack of sep_axis, and thus, overlapping
               !However, doing so will flag all IC's as overlapping since 
               !IC's have to share common edges. So here the equality is considered 
               !as existence of a separating axis, and hence, no overlap 
               !equality is also considered as separating axis which is
               if ((ic_min .ge. ic2_max)  .or. (ic_max .le. ic2_min) ) then 
                  sep_axis = .true. 
                  exit 
               endif
            end DO

            if(.not.sep_axis) then
               !implies the IC regions could not find a separating axis and are thereofre 
               !overlapping 
               
               write(err_msg, 1004) ICV, ICV2
               CALL FLUSH_ERR_MSG(footer = .false.)
               
               DO IDIM = 1, DIMN 
                     
                  write(err_msg, 1005) 'IC1', IDIM, IC_ORIG(IDIM), IC_END(IDIM) 
                  CALL FLUSH_ERR_MSG(header = .false., footer = .false.)
                     
                  write(err_msg, 1005) 'IC2', IDIM, IC2_ORIG(IDIM), IC2_END(IDIM) 
                  CALL FLUSH_ERR_MSG(header = .false., footer = .false.)

               ENDDO
               write(err_msg, 1006) 

               CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
               
            endif
         end DO ICVTWOLOOP
      end DO ICVLOOP
      

 1004 FORMAT('Error # 1004 for DEM Solids IC:',/5x, & 
      'Overlapping IC regions with non zero', /, & 
      'solids volume fraction  not allowed if gener_part_config is true.', /, &
      'Overlapping ICs are', 2(2x, i4))

 1005 FORMAT('Spans of ', A, ' in dir ', I2, /5x, 2(2x, g17.8))
      

 1006 Format('Please correct the data file. Exiting.')


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
      
