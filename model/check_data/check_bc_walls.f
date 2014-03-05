!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_P_INFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!     Unlike the MI boundary, for the PI boundary the velocities at    !
!     the inflow face are calculated by solving the momentum eqns      !
!     and are not fixed. In this way, the PI is similar to the PO      !
!     except that the flow is into the domain and hence all other      !
!     scalars (e.g., mass fractions, void fraction, temperature,       !
!     etc.,) at the inflow cells need to be specified. To satisfy      !
!     the error routines at the start of the simulation, both the      !
!     tangential and normal components at the inflow also need to      !
!     be specified. The velocities values essentially serve as IC.     !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS(M_TOT, SKIP, BCV)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE bc
      USE indices
      USE funits 
      USE scalars
      USE compar
      USE sendrecv
      USE discretelement
      USE mfix_pic
      USE cutcell

      use error_manager

      IMPLICIT NONE


      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT
      LOGICAL, INTENT(in) :: SKIP(DIM_M)


      INTEGER :: I, J ,K, IJK
      INTEGER :: M, N
! Solids phase density in BC region.
      DOUBLE PRECISION :: BC_ROs(MMAX)
! Index of inert species
      INTEGER :: INERT

      DOUBLE PRECISION SUM, SUM_EP

      DOUBLE PRECISION, EXTERNAL :: EOSS
      LOGICAL , EXTERNAL :: COMPARE 


      CALL INIT_ERR_MSG("CHECK_BC_WALLS")

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


! Following section to check quantities for wall type boundaries
! ---------------------------------------------------------------->>>      

! Default specification of Johnson-Jackson bc
            IF (GRANULAR_ENERGY) THEN 
               IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 1 
            ELSE 
               IF (BC_JJ_PS(BCV) == UNDEFINED_I) BC_JJ_PS(BCV) = 0 
            ENDIF 

! Set Jenkins default specification, modify BC_JJ accordingly
            IF (GRANULAR_ENERGY .AND. JENKINS) BC_JJ_PS(BCV) = 1

            IF (BC_TYPE(BCV)=='FREE_SLIP_WALL' .OR. &
                BC_TYPE(BCV)=='NO_SLIP_WALL' .OR. &
                BC_TYPE(BCV)=='PAR_SLIP_WALL') THEN 

! The wall velocities are not needed for no-slip or free-slip                
               IF (BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN 
                  IF (BC_UW_G(BCV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Uw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_VW_G(BCV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Vw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (.NOT.NO_K) THEN 
                     IF (BC_WW_G(BCV) == UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Ww_g', BCV
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE 
                     BC_WW_G(BCV) = ZERO 
                  ENDIF 
               ENDIF   ! end if (bc_type(bcv)=='par_slip_wall')

               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                   .OR.MPPIC) THEN
                  IF (BC_TYPE(BCV)=='PAR_SLIP_WALL' .OR. BC_JJ_PS(BCV)==1) THEN
                     DO M = 1, MMAX
                        IF (BC_UW_S(BCV,M) == UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Uw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_VW_S(BCV,M) == UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Vw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (.NOT.NO_K) THEN 
                           IF (BC_WW_S(BCV,M) == UNDEFINED) THEN 
                             IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Ww_s', BCV, M 
                             call mfix_exit(myPE)  
                           ENDIF 
                        ELSE 
                           BC_WW_S(BCV,M) = ZERO 
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (bc_type(bcv)=='par_slip_wall' or bc_jj_ps(bcv)==1)
               ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid .or.
                       !          mppic)

               IF (ENERGY_EQ) THEN 
                  IF (BC_HW_T_G(BCV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1003) 'BC_hw_T_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_HW_T_G(BCV)/=ZERO .AND. &
                      BC_TW_G(BCV)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_Tw_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 
                  IF (BC_HW_T_G(BCV)/=UNDEFINED .AND. &
                      BC_C_T_G(BCV)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_C_T_g', BCV 
                     call mfix_exit(myPE)  
                  ENDIF 

                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID &
                      .OR. MPPIC) THEN
                     DO M = 1, SMAX 
                        IF (BC_HW_T_S(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1103) 'BC_hw_T_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_T_S(BCV,M)/=ZERO .AND.&
                            BC_TW_S(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Tw_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_T_S(BCV,M)/=UNDEFINED .AND. &
                            BC_C_T_S(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_C_T_s', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid)
               ENDIF   ! end if (energy_eq)

               IF (.NOT. DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID  &
                   .OR. MPPIC) THEN
                  IF (GRANULAR_ENERGY .AND. BC_JJ_PS(BCV)==0) THEN 
                     DO M = 1, MMAX 
                        IF (BC_HW_THETA_M(BCV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1103) &
                              'BC_hw_Theta_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_THETA_M(BCV,M)/=ZERO .AND. &
                            BC_THETAW_M(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_Thetaw_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                        IF (BC_HW_THETA_M(BCV,M)/=UNDEFINED .AND. &
                            BC_C_THETA_M(BCV,M)==UNDEFINED) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'BC_C_Theta_m', BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDDO 
                  ENDIF   ! end if (granular_energy .and. bc_jj_ps(bcv)==0)
               ENDIF    ! end if (.not.discrete_element .or. des_continuum_hybrid 
                        !         .or. mppic)

               IF (SPECIES_EQ(0)) THEN 
                  DO N = 1, NMAX(0) 
                     IF (BC_HW_X_G(BCV,N) < ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_X_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                     IF (BC_HW_X_G(BCV,N)/=ZERO .AND. &
                         BC_XW_G(BCV,N)==UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Xw_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                     IF (BC_HW_X_G(BCV,N)/=UNDEFINED .AND. &
                         BC_C_X_G(BCV,N)==UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_X_g', BCV, N 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ENDDO 
               ENDIF   ! end if (species_eq(0))

               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID .OR.&
                   MPPIC) THEN
                  DO M = 1, SMAX 
                     IF (SPECIES_EQ(M)) THEN 
                        DO N = 1, NMAX(M) 
                           IF (BC_HW_X_S(BCV,M,N) < ZERO) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1105) &
                                 'BC_hw_X_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_HW_X_S(BCV,M,N)/=ZERO .AND. &
                               BC_XW_S(BCV,M,N)==UNDEFINED) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1104) &
                                 'BC_Xw_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                           IF (BC_HW_X_S(BCV,M,N)/=UNDEFINED .AND. &
                               BC_C_X_S(BCV,M,N)==UNDEFINED) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1104) &
                                 'BC_C_X_s', BCV, M, N 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ENDDO 
                     ENDIF    ! end if (species_eq(m))
                  ENDDO    ! end loop over (m=1,smax)
               ENDIF   ! end if (.not.discrete_element .or. des_continuum_hybrid)
                       !         .or. mppic)

               DO N = 1, NScalar
                  IF (BC_HW_Scalar(BCV,N) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1005) 'BC_hw_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
                  IF (BC_HW_Scalar(BCV,N)/=ZERO .AND. &
                      BC_Scalarw(BCV,N)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_ScalarW', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
                  IF (BC_HW_Scalar(BCV,N)/=UNDEFINED .AND. &
                      BC_C_Scalar(BCV,N)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_C_Scalar', BCV, N 
                     CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO    ! end loop over (n=1,nscalar)

            ENDIF  ! end if (bc_type(bcv)=='free_slip_wall' .or. 
                   !        bc_type(bcv)=='no_slip_wall' .or.
                   !        bc_type(bcv)=='par_slip_wall')



      CALL FINL_ERR_MSG


      RETURN  





!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      CALL FINL_ERR_MSG

      RETURN

! 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
!         'correct the mfix.dat file.')


 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ') 
 1002 FORMAT(5X,A16) 
 1003 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1004 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ') not specified',/1X,70('*')/) 
 1005 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_P_g( ',I2,&
         ') = ',G12.5,/&
         ' Pressure should be greater than zero for compressible flow',/1X,70(&
         '*')/) 
 1011 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified.',/1X,'These serve as initial values in ',&
         'the boundary region. Set to 0 by default',/1X,70('*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/) 
 1060 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1065 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1101 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: For ', &
         'BC_TYPE= ', A, ' BC_EP_G(',I2,') should not be defined',/1X,&
         'or the sum of the volume fractions may not equal one',&
         /1X,70('*')/)
 1102 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: For ',& 
         'BC_TYPE= ', A, ' Since BC_EP_G('I2,') is defined,',/1X,&
         'BC_ROP_S should also be defined to ensure that the ',&
         'sum of their',/1X, 'volume fractions will equal one',&
         /1X,70('*')/)

 1103 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/) 
 1104 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/) 
 1105 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1111 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') not specified.',/1X,'These serve as initial values in ',&
         'the boundary region. Set to 0 by default',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 

 !1125 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
 !        ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 

 1130 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - void fraction is unphysical (>1 or <0)',/1X,70('*')/) 
 1150 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: BC number:',I2,&
         ' - ',A,I1,' should be ',A,' zero.',/1X,70('*')/) 
 1160 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: Boundary condition no', &
         I2,' is a second outflow condition.',/1X,&
         '  Only one outflow is allowed.  Consider using P_OUTFLOW.',/1X, 70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K') 
 1410 FORMAT(I5,3X,I5,3X,I5) 
 1420 FORMAT(/1X,70('*')/) 

 1500 FORMAT(//1X,70('*')/' From: CHECK_DATA_07',/,' Error 1500:'      &
         ' Solids phase ',I2,' failed sanity check in IC region ',I3,  &
         '. ',/' Please check mfix.dat file.',/1X,70('*')//)



      END SUBROUTINE CHECK_BC_WALLS
