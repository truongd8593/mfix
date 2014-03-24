!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS_MPPIC                          !
!                                                                      !
!  Purpose: check the initial conditions input section for MPPIC model !
!     - ensure the first IC is defined over the entire domain with     ! 
!        ep_g = 1 when more than one IC has solids                     !
!     - ensure the ICs are non-overlapping                             !
!     - calculate the number of particles needed to initalize the      !
!        MPPIC model                                                   !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!  Comments: Most of the code in this routine is a consolidation of    !
!            codes from several existing routines                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_INITIAL_CONDITIONS_MPPIC
      

! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only : gener_part_config
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Number of particles seeded per phase 
      USE discretelement, only: PART_MPHASE
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles to read from input file.
      USE discretelement, only: PARTICLES

! Implied total number of physical particles 
      USE mfix_pic, only: rnp_pic 
! Total number of computational particles/parcels 
      USE mfix_pic, only: cnp_pic 

      USE mfix_pic, only: cnp_array
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! IC Region solid volume fraction.
      USE ic, only: IC_EP_S
      
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
      USE ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE

! MPPIC specific IC region specification. 
      USE ic, only: IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT
! Maximum number of IC regions and solids phases
      USE param, only: DIMENSION_IC, DIM_M 

! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: xlength, ylength, zlength, dx, dy, dz 

! Constant: 3.14159...
      USE constant, only: PI

      USE mpi_utility
      
! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

      DOUBLE PRECISION :: EP_SM
! Temp variable for storing number of particles.
      DOUBLE PRECISION :: PARTS_TEMP 
! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M)
! Volume of the cell 
      DOUBLE PRECISION :: VOLIJK
      INTEGER :: I, J, K, IJK, ICV, M, COUNT_IC
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'function.inc'

      IF (.NOT.GENER_PART_CONFIG) RETURN 

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_INITIAL_CONDITIONS_MPPIC")

! First check if either a constant npc or constant statwt 
! is specified for each IC 
      DO ICV = 1, DIMENSION_IC 
         
         IF(.not.ic_defined(icv)) cycle         

         IF (IC_EP_G(ICV).lt.ONE) THEN 
            DO M = 1, DES_MMAX           
               CONST_NPC    = (IC_PIC_CONST_NPC   (ICV, M) .ne. 0)
               CONST_STATWT = (IC_PIC_CONST_STATWT(ICV, M) .ne. ZERO  )
               IF(CONST_NPC.and.CONST_STATWT.and.ic_ep_s(icv,m).gt.zero) then
                  WRITE(ERR_MSG, 1100) ICV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               
               IF(.not.CONST_NPC.and.(.not.CONST_STATWT).and. & 
               ic_ep_s(icv,m).gt.zero) then
                  WRITE(ERR_MSG, 1101) ICV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

 1100 FORMAT('Error 1100: In MPPIC model for IC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'Non zero Values specified for both ', &
      'IC_PIC_CONST_NPC and IC_PIC_CONST_STATWT.', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')


 1101 FORMAT('Error 1101: In MPPIC model for IC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'A non-zero value not specified for ', &
      'IC_PIC_CONST_NPC or IC_PIC_CONST_STATWT. ', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')

      if(.not.allocated(part_mphase_byic)) & 
      ALLOCATE(PART_MPHASE_BYIC(DIMENSION_IC, DES_MMAX))    

! cnp_array(ijk, 0) will contain the cumulative number of real 
! particles later in the handling of inflow BC for MPPIC. See 
! the mppic_mi_bc in mppic_wallbc_mod.f 
      ALLOCATE(CNP_ARRAY(DIMENSION_3, DES_MMAX))

      ALLOCATE(RNP_PIC(DES_MMAX))
      ALLOCATE(CNP_PIC(DES_MMAX))

      CNP_ARRAY(:, :) = 0
      PART_MPHASE_BYIC(:,:) = 0 

      RNP_PIC(:) = ZERO
      CNP_PIC(:) = ZERO 

      COUNT_IC = 0
      DO ICV = 1, DIMENSION_IC 
         IF(.not.ic_defined(icv)) cycle         
         
         IF (IC_EP_G(ICV).lt.ONE) THEN 
            COUNT_IC = COUNT_IC + 1 

            DO K = IC_K_B(ICV), IC_K_T(ICV) 
               DO J = IC_J_S(ICV), IC_J_N(ICV) 
                  DO I = IC_I_W(ICV), IC_I_E(ICV) 

                     IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle 
                     IJK = FUNIJK(I,J,K)

                    ! IF(.not.FLUID_AT(IJK)) cycle 

                     VOLIJK = Dx(I)*Dy(J)*Dz(K)
                     DO M = 1, DES_MMAX
                        
                        REAL_PARTS(M) = 0.
                        COMP_PARTS(M) = 0

                        EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)

                        REAL_PARTS(M) = 6.d0*EP_SM*VOLIJK/ &
                        (PI*(Des_D_P0(M)**3.d0))

                        CONST_NPC    = (IC_PIC_CONST_NPC   (ICV, M) .ne. 0) &
                        .AND. (EP_SM.gt.Zero)

                        CONST_STATWT = (IC_PIC_CONST_STATWT(ICV, M) .ne. ZERO  ) &
                        .AND. (EP_SM.gt.Zero)
                        
                        IF(CONST_NPC) THEN 
                           COMP_PARTS(M) = IC_PIC_CONST_NPC(ICV,M)
                        ELSEIF(CONST_STATWT) THEN
                           COMP_PARTS(M) = MAX(1, & 
                           INT(REAL_PARTS(M)/REAL(IC_PIC_CONST_STATWT(ICV,M))))
                        ENDIF
                        
                        RNP_PIC(M) = RNP_PIC(M) + REAL_PARTS(M)
                        CNP_PIC(M) = CNP_PIC(M) + COMP_PARTS(M)
                        
                        PART_MPHASE_BYIC(ICV, M) =  PART_MPHASE_BYIC(ICV, M) & 
                        + COMP_PARTS(M)
                                               
                        CNP_ARRAY(IJK,M) = COMP_PARTS(M)
                     ENDDO      ! M = 1, DES_MMAX
                  ENDDO
               ENDDO
            ENDDO
         ENDIF                  !ic defined
      ENDDO                     !ICV = 1, DIMENSION_IC


      CALL global_all_sum(RNP_PIC)
      CALL global_all_sum(CNP_PIC)
      CALL global_all_sum(PART_MPHASE_BYIC)

      PART_MPHASE(1:DES_MMAX) = CNP_PIC(1:DES_MMAX)
      PARTICLES = SUM(PART_MPHASE(1:DES_MMAX))
      
      
! Local reporting for the user 
      WRITE(ERR_MSG, 2015)  COUNT_IC, DES_MMAX
      CALL FLUSH_ERR_MSG(FOOTER = .false.)

      DO ICV = 1, DIMENSION_IC 
         IF(.not.ic_defined(icv)) cycle 

         IF (IC_EP_G(ICV).lt.ONE) THEN 
            
            WRITE(ERR_MSG,2016) ICV
            CALL FLUSH_ERR_MSG(HEADER = .false., FOOTER = .false.)

            DO M = 1, DES_MMAX
               
               CONST_NPC    = (IC_PIC_CONST_NPC   (ICV, M) &
               .ne. 0)
               CONST_STATWT = (IC_PIC_CONST_STATWT(ICV, M) & 
               .ne. ZERO  )

               EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
               WRITE(ERR_MSG,2017) M, DES_D_P0(M), EP_SM, & 
               CONST_NPC, CONST_STATWT, PART_MPHASE_BYIC(ICV, M) 
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

 2016 FORMAT('IC number:         = ', I4 )
      
 2017 FORMAT( &
      'PHASE Index, M              =  ', I5,2X, /5x, & 
      'D_P0(M)                     =  ', (ES15.7,2X), /5x, &
      'Solids vol fraction for M   =  ', (G15.8,2X), /5x, & 
      'Const Npc?', 2x, L4, 'Const Statwt?', 2x, L4, /5x, & 
      '# of computational particles or parcels for phase M  = ',  (I10,2X))
      
 2018 FORMAT( 'End of Message' ) 



      CALL FINL_ERR_MSG
      
      END SUBROUTINE CHECK_INITIAL_CONDITIONS_MPPIC
      
