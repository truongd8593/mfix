!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_COMMON_DISCRETE                                !
!                                                                      !
!  Purpose: check the initial conditions input section common to both  !
!           DEM and MPPIC models 
!     - ensure the first IC is defined over the entire domain with     ! 
!        ep_g = 1 when more than one IC has solids                     !
!     - ensure the ICs are non-overlapping                             !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!  Comments: Most of the code in this routine is a consolidation of    !
!            codes from several existing routines                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_COMMON_DISCRETE 
      

! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only : gener_part_config
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! IC Region solid volume fraction.
      USE ic, only: IC_EP_S
! IC Region gas volume fraction.
      USE ic, only: IC_THETA_M 
      
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t

      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE

! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: xlength, ylength, zlength


! Maximum number of IC regions 
      USE param, only: DIMENSION_IC
      
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

      IF (.NOT.GENER_PART_CONFIG) RETURN 

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_COMMON_DISCRETE")
      
! First check if multiple IC regions are defined for non-zero solids volume 
! fraction, then check if the first IC is specified over the whole domain with IC_EP_g = 1 

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
             
         IF (FIRST_IC_OK.AND.DIMN.EQ.3.AND.IC_Z_B(ICV).LE.ZERO.AND.IC_Z_T(ICV).GE.ZLENGTH) FIRST_IC_OK = .TRUE.
         
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



! Check if IC_theta_M is specified for solids phases wherever IC_EP_g lt 1
      DO ICV = 1, DIMENSION_IC 
         
         IF (IC_DEFINED(ICV).and.IC_EP_G(ICV).LT.ONE) THEN 
            DO M = 1, DES_MMAX
               IF(IC_THETA_M(ICV,M)==UNDEFINED) THEN
                  IF(IC_EP_S(ICV,M).gt.Zero) then 
                     WRITE(ERR_MSG, 2019) ICV,M
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ELSE
                     IC_Theta_M(ICV,M) = ZERO 
                  ENDIF
               ENDIF 
            ENDDO
         ENDIF
      ENDDO

 2019 FORMAT('Error 2019: For IC #', 2x, i4, /, & 
      'Initial value of granular temperature (IC_Theta_M) not ',&
      'defined for solid phase (M) : ', 2x, i3, /, & 
      'Exitting')

      CALL FINL_ERR_MSG
      
      END SUBROUTINE CHECK_IC_COMMON_DISCRETE 
