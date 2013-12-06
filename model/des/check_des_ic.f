!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_IC                                           !
!                                                                      !
!  Purpose: Check the data provided for des initial conditions.        !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Feb-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_IC

      USE des_ic
      USE discretelement 
      USE param
      USE param1
      USE des_thermo
      USE des_rxns
      USE compar
      USE constant
      USE funits  
      USE geometry
      USE indices
      USE physprop
      USE run

      IMPLICIT NONE


      INTEGER ICV

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 

!-----------------------------------------------           

      IF(DMP_LOG) WRITE(*,'(1X,A)') &
         '---------- START CHECK_DES_IC ---------->'

! Initialize the initial condition logical.
      DES_IC_EXIST = .FALSE.

! Check for des inlet/outlet information:
      CHECK_IC: DO ICV = 1, DIMENSION_IC 
         DES_IC_DEFINED(ICV) = .FALSE. 
         IF (DES_IC_X_w(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_X_e(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Y_s(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Y_n(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Z_b(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 
         IF (DES_IC_Z_t(ICV) /= UNDEFINED) DES_IC_DEFINED(ICV) = .TRUE. 

! If an initial condition is specified, verify necessary data
! Check that all dimensions have been defined
         IF (DES_IC_DEFINED(ICV)) THEN
            IF(DES_IC_X_w(ICV) == UNDEFINED .OR. &
               DES_IC_X_e(ICV) == UNDEFINED .OR. &
               DES_IC_Y_s(ICV) == UNDEFINED .OR. &
               DES_IC_Y_n(ICV) == UNDEFINED) THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1000) ICV
                  WRITE (*, 1000) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF 
            IF(DIMN == 3)THEN
               IF(DES_IC_Z_b(ICV) == UNDEFINED .OR. &
                 DES_IC_Z_t(ICV) == UNDEFINED) THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1000) ICV
                  WRITE (*, 1000) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
              ENDIF
            ENDIF
! Check that the region does not overlap
            IF(DES_IC_X_w(ICV) .LT. ZERO .OR. &
               DES_IC_Y_s(ICV) .LT. ZERO .OR. &
               DES_IC_X_e(ICV) .GT. XLENGTH .OR. &
               DES_IC_Y_n(ICV) .GT. YLENGTH .OR. &
               DES_IC_X_w(ICV) .GT. DES_IC_X_e(ICV) .OR. &
               DES_IC_Y_s(ICV) .GT. DES_IC_Y_n(ICV))THEN
               IF(DMP_LOG) THEN
                  WRITE (UNIT_LOG, 1001) ICV
                  WRITE (*, 1001) ICV
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF

            IF(DIMN == 3)THEN
               IF(DES_IC_Z_b(ICV) .LT. 0 .OR. &
                  DES_IC_Z_t(ICV) .GT. ZLENGTH .OR. &
                  DES_IC_Z_b(ICV) .GT. DES_IC_Z_t(ICV))THEN
                  IF(DMP_LOG) THEN
                     WRITE (UNIT_LOG, 1001) ICV
                     WRITE (*, 1001)ICV
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF

! This is used to invoke the call to DES_SET_IC which is only done
! if there are valid initial conditions specified.
            DES_IC_EXIST = .TRUE.

         ENDIF   ! end if DES_IC_DEFINED(ICV)

      ENDDO CHECK_IC


! If the DES energy equations or DES species equation are being solve
! then initial conditions must be specified to provide the intial 
! sate of the particles.
      IF(.NOT.DES_IC_EXIST) THEN
         IF(DES_ENERGY_EQ) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1002)'temperature','energy equation'
               WRITE(UNIT_LOG,1002)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ELSEIF(ANY_DES_SPECIES_EQ) THEN
            IF(DMP_LOG) THEN
               WRITE(*,1002)'species mass fractions','species equation'
               WRITE(UNIT_LOG,1002)
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      WRITE(*,'(1X,A)') '<---------- END CHECK_DES_IC ----------'

 1000 FORMAT(/1X,70('*')/, ' From: CHECK_DES_IC',/, ' Message: ',&
         'Insufficient DEM initial condition infomation',/10X,&
         'Check initial condition number: ',I3,/1X,70('*')/)

 1001 FORMAT(/1X,70('*')/, ' From: CHECK_DES_IC',/, ' Messsage: ',&
         'Improper DEM initial condition information',/10X,&
         'Check initial condition number: ',I3,/1X,70('*')/)

 1002 FORMAT(/1X,70('*')/,' From: CHECK_DES_IC',/'Message: Initial',   &
         ' condition data specifying the paritcle',/A,' is requied',   &
         ' when solving the ',A,'.',/' Check mfix.dat.',/1X,70('*')/)


      RETURN
      END SUBROUTINE CHECK_DES_IC
      
      SUBROUTINE CHECK_DES_IC_EPS

      USE param 
      USE param1 
      USE geometry
      USE ic
      USE fldvar
      USE physprop
      USE run
      USE indices
      USE funits 
      USE scalars
      USE compar      
      USE mpi_utility      
      USE sendrecv    
      USE discretelement      
      
      IMPLICIT NONE


      INTEGER ICV, M 

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 

      INTEGER :: I,J,K, count_ic, count_ic_with_sols
      integer :: first_def_ic 
      DOUBLE PRECISION SUM, SUM_EP
! Solids phase density in IC region.
      DOUBLE PRECISION :: IC_ROs

      logical :: first_ic_ok
!-----------------------------------------------           

      IF(DMP_LOG) WRITE(UNIT_LOG,'(1X,A,/, 1X,A)') &
      '---------- GENER PART CONFIG TRUE----------', &
      '---------- So Calling CHECK_DES_IC_EPS ------->'
!      IF(mype.eq.pe_IO) WRITE(*,'(1X,A,/, 1X,A)') &
!      '---------- GENER PART CONFIG TRUE----------', &
!      '---------- So Calling CHECK_DES_IC_EPS ------->'

      IF(.NOT.allocated(vol_ic_region)) & 
      ALLOCATE(VOL_IC_REGION(DIMENSION_IC))

      if(.not.allocated(part_mphase_byic)) & 
      ALLOCATE(PART_MPHASE_BYIC(DIMENSION_IC, DIM_M))
      
      IF (DES_CONTINUUM_HYBRID) THEN
         IF(DMP_LOG) write(UNIT_LOG, 1001)
         IF(mype.eq.pe_io) write(*, 1001)
         call  mfix_exit(mype)
      ENDIF     
      !total count of defined ICs
      count_ic           = 0 
      !total count of defined IC's with solids 
      count_ic_with_sols = 0
      first_def_ic = undefined_i
      DO ICV = 1, DIMENSION_IC 

         IF (IC_DEFINED(ICV)) THEN 
            VOL_IC_REGION(ICV) = zero 
            count_ic = count_ic + 1 
            first_def_ic = min(first_def_ic, icv) 
!----------------------------------------------------------------->>>
! For restart runs IC is defined only if IC_TYPE='PATCH'
            IF (RUN_TYPE/='NEW' .AND. IC_TYPE(ICV)/='PATCH') THEN 
               IC_DEFINED(ICV) = .FALSE. 
               CYCLE  
            ENDIF 


            SUM_EP = IC_EP_G(ICV) 
            IF(IC_EP_G(ICV).LT.ONE) count_ic_with_sols & 
                 = count_ic_with_sols  + 1

            DO M = 1, SMAX
! Bulk density must be explicitly defined if there are more than one
! solids phase or using hybrid model.
               IF(IC_ROP_S(ICV,M) == UNDEFINED) THEN
                  IF(IC_EP_G(ICV) == ONE) THEN
                     IC_ROP_S(ICV,M) = ZERO
                  ELSEIF(SMAX > 1) THEN
                     IF(DMP_LOG) WRITE(UNIT_LOG, 1100)           &
                     'IC_ROP_s', ICV, M 
                     CALL MFIX_EXIT(myPE) 
                  ENDIF
               ENDIF
! Calculate the solid density.

               IC_ROs = DES_RO_s(M)

! Sanity check on solids phase density.
               IF(IC_ROs <= ZERO .OR. IC_ROs==UNDEFINED) THEN
                  IF(DMP_LOG)THEN
                     WRITE(*,1401) M, ICV
                     WRITE(UNIT_LOG,1401) M, ICV
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF

! IC_ROP_S may still be undefined if there is only one solids phase
! and the hybrid model is not in use. Back out the bulk density with
! the assigned solids density.
               IF(IC_ROP_S(ICV,M) == UNDEFINED)                  &
               IC_ROP_S(ICV,M) = (ONE - IC_EP_G(ICV))*IC_ROs

! at this point ic_rop_s must be defined
! sum of void fraction and solids volume fractions                  
               SUM_EP = SUM_EP + IC_ROP_S(ICV,M)/IC_ROs 


! check that solids phase m velocity components are initialized
               IF (IC_U_S(ICV,M) == UNDEFINED) THEN 
                  IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_I) THEN 
                     IC_U_S(ICV,M) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                          'IC_U_s', ICV, M 
                     call mfix_exit(myPE) 
                  ENDIF
               ENDIF

               IF (IC_V_S(ICV,M) == UNDEFINED) THEN 
                  IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_J) THEN 
                     IC_V_S(ICV,M) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                          'IC_V_s', ICV, M 
                     call mfix_exit(myPE) 
                  ENDIF
               ENDIF

               IF (IC_W_S(ICV,M) == UNDEFINED) THEN 
                  IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_K) THEN 
                     IC_W_S(ICV,M) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                          'IC_W_s', ICV, M 
                     call mfix_exit(myPE) 
                  ENDIF
               ENDIF
            
               IF( IC_THETA_M(ICV,M) == UNDEFINED) THEN 
                  IF (IC_ROP_S(ICV,M) == ZERO) THEN 
                     IC_THETA_M(ICV,M) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                          'IC_Theta_m', ICV, M 
                     call mfix_exit(myPE) 
                  ENDIF
               ENDIF
               
            ENDDO               ! end loop over (m=1,smax)

            IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
               IF(DMP_LOG)WRITE (UNIT_LOG, 1125) ICV 
               IF(mype.eq.pe_IO)WRITE (*, 1125) ICV 
               call mfix_exit(myPE)
            ENDIF 

            DO K = IC_K_B(ICV), IC_K_T(ICV) 
               DO J = IC_J_S(ICV), IC_J_N(ICV) 
                  DO I = IC_I_W(ICV), IC_I_E(ICV) 
                     VOL_IC_REGION(ICV) = VOL_IC_REGION(ICV)+DX(I)*DY(J)*DZ(K)
                  ENDDO 
               ENDDO
            ENDDO

            
! SOLIDS PHASE Quantities
! --------------------------------------------<<<
         ENDIF                  ! if(ic_defined(icv))
      end DO
      
      if(count_ic_with_sols.ge.1.and.count_ic.ge.count_ic_with_sols+1) then
! If the number of IC's with solids is greater then one 
! and if such number of IC's are greater than total number of ICs then 
! make sure the first IC spans the entire domain with voidage of 1
! This is done so that users can specify arbitrary non-overlapping IC's with solids 
! without having to break up and specify IC for rest of the regions. 
         ICV = first_def_ic 
         first_ic_ok = .false. 
         if(IC_EP_G(ICV).eq.ONE & 
              .and.IC_X_W(ICV).le.zero.AND.IC_X_E(ICV).ge.Xlength &
              .AND.IC_Y_S(ICV).le.zero.and.IC_Y_N(ICV).ge.Ylength) first_ic_ok = .true.
         !check for z also 
         if (dimn.eq.3.AND.IC_Z_B(ICV).le.zero.and.IC_Z_T(ICV).ge.Zlength) first_ic_ok = .true.
         IF(.not.first_ic_ok) then 
               
            if(mype.eq.pe_io) write(*,1402) ICV, ic_ep_g(icv), & 
                 IC_X_W(ICV), IC_X_E(ICV), &
                 IC_Y_S(icv), ic_y_n(icv), & 
                 ic_z_b(icv), ic_z_t(icv)

            if(DMP_LOG) write(UNIT_LOG,1402) ICV, ic_ep_g(icv), & 
                 IC_X_W(ICV), IC_X_E(ICV), &
                 IC_Y_S(icv), ic_y_n(icv), & 
                 ic_z_b(icv), ic_z_t(icv)
            call mfix_exit(mype)
         endif
      endif
   
1001  format(/1X,70('*')/' From: CHECK_DES_IC_EPS',/, &
           ' Error 1001:',       &
           ' Gener_part_config set to true for DES_continuum hybrid', /, &
           ' This is not allowed and specify the particle initial', &
           ' configuration explicitly', /, &
           ' See MFIX readme', /,  &
           ' Please correct the data file. Exiting.',/, 1X,70('*')/)
1100  FORMAT(/1X,70('*')//' From: CHECK_DES_IC_EPS',/ &
           ' Message: ',A,'(',I2,',',I1,&
           ') not specified',/1X,70('*')/) 
      
      
1125  FORMAT(/1X,70('*')//' From: CHECK_DES_IC_EPS',/, & 
           ' Message: IC number:',I2,&
           ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
      
1401  FORMAT(//1X,70('*')/' From: CHECK_DES_IC_EPS',/,' Error 1401:' , &
           ' Solids phase ',I2,' failed sanity check in IC region ',I3,  &
           '. ',/' Please check mfix.dat file.',/1X,70('*')//)
      
1402  format(//1X,70('*')/' From: CHECK_DES_IC_EPS',/,' Error 1402:'      &
           'For initial particle seeding with more than once IC having solids', /, &
           'first defined IC should span the entire domain and', /, & 
           'and gas voidage should be set to one.', /, & 
           'Not respecting this will exit the code' / & 
           'IC number and IC_EP_g', I10, 2x, g17.8, /, & 
           'IC spans X direction                :', 2(2x, g17.8) , /, & 
           'IC spans Y direction                :', 2(2x, g17.8) , /, & 
           'IC spans Z direction (ignore for 2D):', 2(2x, g17.8) , /, & 
           'Please check mfix.dat file. Exitting',/1X,70('*')//)
      
      IF(DMP_LOG) WRITE(UNIT_LOG,'(1X,A)') '<---------- END CHECK_DES_IC_EPS ----------'
      !IF(mype.eq.pe_IO) WRITE(*,'(1X,A)') '<---------- END CHECK_DES_IC_EPS ----------'
      end SUBROUTINE CHECK_DES_IC_EPS
 
      subroutine des_check_ic_overlap

      USE param 
      USE param1 
      USE geometry
      USE ic
      USE fldvar
      USE physprop
      USE run
      USE indices
      USE funits 
      USE scalars
      USE compar      
      USE mpi_utility      
      USE sendrecv    
      USE discretelement      
      
      IMPLICIT NONE

      INTEGER ICV, ICV2, M, idim
      double precision :: ic_orig(3), ic_end(3), ic2_orig(3) , ic2_end(3)
      double precision :: ic_min, ic_max, ic2_min, ic2_max , tol_ic_reg
      
      logical :: sep_axis

      tol_ic_reg  = 1e-04
      ICVLOOP : DO ICV = 1, DIMENSION_IC 
         
         IF (.not.IC_DEFINED(ICV).or.IC_EP_G(ICV).eq.1.d0) cycle ICVLOOP
         IC_ORIG(1) = XE(IC_I_W(ICV)) - DX(IC_I_W(ICV))
         IC_ORIG(2) = YN(IC_J_S(ICV)) - DY(IC_J_S(ICV))
         IC_ORIG(3) = ZT(IC_K_B(ICV)) - DZ(IC_K_B(ICV))
         IC_END(1)  = XE(IC_I_E(ICV))
         IC_END(2)  = YN(IC_J_N(ICV))
         IC_END(3)  = ZT(IC_K_T(ICV))
         ICVTWOLOOP : DO ICV2 = ICV+1, DIMENSION_IC 
            
            IF (.not.IC_DEFINED(ICV2).or.IC_EP_G(ICV2).eq.1.d0) cycle ICVTWOLOOP
            
            
            IC2_ORIG(1) = XE(IC_I_W(ICV2)) - DX(IC_I_W(ICV2))
            IC2_ORIG(2) = YN(IC_J_S(ICV2)) - DY(IC_J_S(ICV2))
            IC2_ORIG(3) = ZT(IC_K_B(ICV2)) - DZ(IC_K_B(ICV2))
            IC2_END(1)  = XE(IC_I_E(ICV2))
            IC2_END(2)  = YN(IC_J_N(ICV2))
            IC2_END(3)  = ZT(IC_K_T(ICV2))

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
               if(dmp_log) then 
                  write(unit_log, 1001) ICV, ICV2
                  DO IDIM = 1, DIMN 
                     write(unit_log, 1002) 'IC1', IDIM, IC_ORIG(IDIM), IC_END(IDIM) 
                     write(unit_log, 1002) 'IC2', IDIM, IC2_ORIG(IDIM), IC2_END(IDIM) 
                  ENDDO
                  write(unit_log, 1003) 
               endif
               if(mype.eq.pe_io) then 
                  write(*, 1001) ICV, ICV2
                  DO IDIM = 1, DIMN 
                     write(*, 1002) 'IC1', IDIM, IC_ORIG(IDIM), IC_END(IDIM) 
                     write(*, 1002) 'IC2', IDIM, IC2_ORIG(IDIM), IC2_END(IDIM) 
                  ENDDO
                  write(*, 1003) 
               end if
               
               call mfix_exit(mype)
            endif 
         end DO ICVTWOLOOP
      end DO ICVLOOP
      

1001  FORMAT(/1X,70('*'),/5x, 'From: DES_CHECK_IC_OVERLAP: Error 1001',/5x, & 
           'Overlapping IC regions with non zero', /5x, & 
           'solids volume fraction  not allowed if gener_part_config is true.', /5x, &
           'Overlapping ICs are', 2(2x, i4))

1002  FORMAT(5x, &
           'Spans of ', A, ' in dir ', I2, /5x, 2(2x, g17.8))
      

1003  Format(/5x, 'Please correct the data file. Exiting.',/, 1X,70('*')/)

    end subroutine des_check_ic_overlap
    
      
