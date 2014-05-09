!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS                                           !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Driver routine to call checks for WALL BCs.                 !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS(M_TOT, SKIP, BCV)


! Global Variables:
!---------------------------------------------------------------------//
! Flag: Identifies solids model (TFM,DEM,PIC)
      use run, only: SOLIDS_MODEL
! Flag: Solve granular energy PDE
      use run, only: GRANULAR_ENERGY
! Flag: Use Jenkins small small frication BC
      use run, only: JENKINS
! Flag: Use revised phihp for JJ BC.
      use bc, only: BC_JJ_PS
! User-input: solids kinetic-theory model.
      use run, only: KT_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases
      use param, only: DIM_M
! Parameter constant.
      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
! Total number of solids phases.
      INTEGER, INTENT(in) :: M_TOT
! Flag. Solids not present at this BC (used for flow BCs).
      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! Local Variables:
!---------------------------------------------------------------------//
! Loop/counter variable.
      INTEGER :: M
! Local total number of solids phases
      INTEGER :: MTOT_L
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS")

! Input checks for gas phase.
      CALL CHECK_BC_WALLS_GAS(BCV)

      MTOT_L = merge( M_TOT+1, M_TOT, KT_TYPE(1:3) == 'GHD')

! Input checks for solid phases.
      DO M=1, MTOT_L
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM'); CALL CHECK_BC_WALLS_TFM(BCV, M)
         CASE ('DEM'); CALL CHECK_BC_WALLS_DISCRETE(BCV, M)
         CASE ('PIC'); CALL CHECK_BC_WALLS_DISCRETE(BCV, M)
         END SELECT
      ENDDO

! Input checks for user-defined scalar equations.
      CALL CHECK_BC_WALLS_SCALAR_EQ(BCV)

! Set the default specification of Johnson-Jackson BC
      IF(BC_JJ_PS(BCV) == UNDEFINED_I)                                 &
         BC_JJ_PS(BCV) = merge(1,0,GRANULAR_ENERGY)

! Set Jenkins default specification, modify BC_JJ accordingly
      IF(GRANULAR_ENERGY .AND. JENKINS) BC_JJ_PS(BCV) = 1

      CALL FINL_ERR_MSG

      RETURN  
      END SUBROUTINE CHECK_BC_WALLS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_GAS                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for gas phase WALL BC parameters.          !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_GAS(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc, only: BC_TYPE
! User-Input: gas velocity at wall BCs.
      use bc, only: BC_UW_G, BC_VW_G, BC_WW_G
! User-Input: gas energy eq BCs.
      use bc, only: BC_HW_T_G, BC_TW_G, BC_C_T_G
! User-Input: gas species eq BCs.
      use bc, only: BC_HW_X_G, BC_XW_G, BC_C_X_G
! Total number of speices in each phase.
      use physprop, only: NMAX
! Flag: Solve energy equations.
      use run, only: ENERGY_EQ
! Flag: Solve species equations.
      use run, only: SPECIES_EQ
! Flag: Solve K-th direction (3D)
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
      INTEGER, INTENT(in) :: BCV

      INTEGER :: N
!......................................................................!


! Initialize the error manger.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_GAS")

! The wall velocities are not needed for no-slip or free-slip                
      IF(BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN 
         IF(BC_UW_G(BCV) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_VW_G(BCV) == UNDEFINED) THEN 
            WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_WW_G(BCV) == UNDEFINED) THEN
            IF(DO_K)THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_g',BCV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSE 
               BC_WW_G(BCV) = ZERO 
            ENDIF
         ENDIF
      ENDIF

! Check energy equation input.
      IF(ENERGY_EQ) THEN
         IF(BC_HW_T_G(BCV) < ZERO) THEN 
            WRITE(ERR_MSG,1001) trim(iVar('BC_HW_T_g',BCV)),           &
               trim(iVal(BC_HW_T_G(BCV)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_T_G(BCV)/=ZERO .AND.                                 &
            BC_TW_G(BCV)==UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_Tw_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_T_G(BCV)/=UNDEFINED .AND.                            &
            BC_C_T_G(BCV)==UNDEFINED) THEN 
            WRITE(ERR_MSG,1000) trim(iVar('BC_C_T_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF


! Check species equation input.
      IF(SPECIES_EQ(0)) THEN
         DO N=1, NMAX(0)
            IF(BC_HW_X_G(BCV,N) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_X_g',BCV,N)),      &
                  trim(iVal(BC_HW_X_G(BCV,N)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_X_G(BCV,N)/=ZERO .AND.                            &
               BC_XW_G(BCV,N)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Xw_g',BCV,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_X_G(BCV,N)/=UNDEFINED .AND.                       &
               BC_C_X_G(BCV,N)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_X_g',BCV,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDDO 
      ENDIF


! Clear the error manager.
      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_TFM                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for TFM solids WALL BC parameters.         !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_TFM(BCV,M)



! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc, only: BC_TYPE
! User-Input: solids velocity at wall BCs.
      use bc, only: BC_UW_s, BC_VW_s, BC_WW_s
! User-Input: solids energy eq BCs.
      use bc, only: BC_HW_T_s, BC_TW_s, BC_C_T_s
! User-Input: solids species eq BCs.
      use bc, only: BC_HW_X_s, BC_XW_s, BC_C_X_s
! User-Input: granular energy eq BCs.
      use bc, only: BC_HW_THETA_M, BC_ThetaW_M, BC_C_Theta_M
! Total number of solids phases
      use physprop, only: MMAX
! Total number of speices in each phase.
      use physprop, only: NMAX
! Flag: Solve energy equations.
      use run, only: ENERGY_EQ
! Flag: Solve species equations.
      use run, only: SPECIES_EQ
! Flag: Solve Granular energy PDE
      use run, only: GRANULAR_ENERGY
! Flag: Use revised phihp for JJ BC.
      use bc, only: BC_JJ_PS
! Flag: Solve K-th direction (3D)
      use geometry, only: DO_K
! User-input: solids kinetic-theory model.
      use run, only: KT_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants.
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC being checked.
      INTEGER, INTENT(in) :: BCV
! Index of solids phase.
      INTEGER, INTENT(in) :: M

! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable counter.
      INTEGER :: N
! Flag to check momentum eq input.
      LOGICAL :: CHECK_MOMENTUM
! Flag to check scalar eq input.
      LOGICAL :: CHECK_SCALARS
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_TFM")

! Toggle the momentum and scalar input variable checks.
      SELECT CASE(trim(adjustl(KT_TYPE)))
      CASE ('GHD')
         CHECK_MOMENTUM = (M == MMAX)
         CHECK_SCALARS  = (M /= MMAX)
      CASE DEFAULT
         CHECK_MOMENTUM = .TRUE.
         CHECK_SCALARS  = .TRUE.
      END SELECT

! The wall velocities are not needed for no-slip or free-slip
      IF(CHECK_MOMENTUM) THEN
         IF(BC_TYPE(BCV) == 'PAR_SLIP_WALL') THEN 
            IF(BC_UW_S(BCV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Uw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(BC_VW_S(BCV,M) == UNDEFINED) THEN 
               WRITE(ERR_MSG,1000) trim(iVar('BC_Vw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(BC_WW_S(BCV,M) == UNDEFINED) THEN
               IF(DO_K)THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Ww_s',BCV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ELSE 
                  BC_WW_S(BCV,M) = ZERO 
               ENDIF
            ENDIF
         ENDIF

         IF(GRANULAR_ENERGY .AND. BC_JJ_PS(BCV)==0) THEN
            IF(BC_HW_THETA_M(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_Theta_M',BCV,M)),  &
                  trim(iVal(BC_HW_Theta_M(BCV,M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
            IF(BC_HW_THETA_M(BCV,M)/=ZERO .AND.                        &
               BC_THETAW_M(BCV,M)==UNDEFINED) THEN 
               WRITE(ERR_MSG,1000) trim(iVar('BC_ThetaW_M',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
            IF(BC_HW_THETA_M(BCV,M)/=UNDEFINED .AND.                   &
               BC_C_THETA_M(BCV,M)==UNDEFINED) THEN 
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_THETA_M',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ELSE
      ENDIF

      IF(CHECK_SCALARS)THEN
         IF(ENERGY_EQ) THEN
            IF(BC_HW_T_S(BCV,M) < ZERO) THEN 
               WRITE(ERR_MSG,1001) trim(iVar('BC_HW_T_s',BCV,M)),      &
                  trim(iVal(BC_HW_T_S(BCV,M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_T_S(BCV,M)/=ZERO .AND.                            &
               BC_TW_S(BCV,M)==UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('BC_Tw_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(BC_HW_T_S(BCV,M)/=UNDEFINED .AND.                       &
               BC_C_T_S(BCV,M)==UNDEFINED) THEN 
               WRITE(ERR_MSG,1000) trim(iVar('BC_C_T_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(SPECIES_EQ(M)) THEN
            DO N=1, NMAX(M)
               IF(BC_HW_X_S(BCV,M,N) < ZERO) THEN
                  WRITE(ERR_MSG,1001) trim(iVar('BC_HW_X_s',BCV,M,N)), &
                     trim(iVal(BC_HW_X_S(BCV,M,N)))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               IF(BC_HW_X_S(BCV,M,N)/=ZERO .AND.                       &
                  BC_XW_S(BCV,M,N)==UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_Xw_s',BCV,M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               IF(BC_HW_X_S(BCV,M,N)/=UNDEFINED .AND.                  &
                  BC_C_X_S(BCV,M,N)==UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('BC_C_X_s',BCV,M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF 
            ENDDO 
         ENDIF ! Species Equation
      ELSE
      ENDIF ! Check Scalars

      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_TFM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_DISCRETE                                  !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for DEM/PIC solids WALL BC parameters.     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_DISCRETE(BCV,M)


! Global Variables:
!---------------------------------------------------------------------//
! User-input: type of BC
      use bc, only: BC_TYPE
! User-Input: solids velocity at wall BCs.
      use bc, only: BC_UW_s, BC_VW_s, BC_WW_s
! User-Input: solids energy eq BCs.
      use bc, only: BC_HW_T_s, BC_TW_s, BC_C_T_s
! User-Input: solids species eq BCs.
      use bc, only: BC_HW_X_s, BC_XW_s, BC_C_X_s
! Total number of solids phases
      use physprop, only: MMAX
! Total number of speices in each phase.
      use physprop, only: NMAX
! Flag: Solve energy equations.
      use run, only: ENERGY_EQ
! Flag: Solve species equations.
      use run, only: SPECIES_EQ

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of possible species.
      use param, only: DIM_N_S
! Parameter constants.
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC getting checked.
      INTEGER, INTENT(in) :: BCV
! Index of solid phase getting checked.
      INTEGER, INTENT(in) :: M

! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable counter.
      INTEGER :: N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_DISCRETE")

! DEM and PIC are restricted to adibatic walls.
      IF(BC_HW_T_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_HW_T_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_TW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_Tw_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_C_T_S(BCV,M) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_C_T_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: ',A,' should not specified for DEM/PIC',/    &
         'simulations as they are currently limited to adiabatic BCs.',&
         /'Please correct the mfix.dat file.')


! The following checks verify that TFM solids parameters are not
! specified for discrete solids.


! The wall velocities are not needed DEM/PIC solids
      IF(BC_UW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Uw_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_VW_S(BCV,M) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Vw_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(BC_WW_S(BCV,M) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Ww_s',BCV,M))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! DEM cannot have a species flux at the walls.
      DO N=1, DIM_N_s
         IF(BC_HW_X_S(BCV,M,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_HW_X_s',BCV,M,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_XW_S(BCV,M,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_Xw_s',BCV,M,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(BC_C_X_S(BCV,M,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1101) BCV, trim(iVar('BC_C_X_s',BCV,M,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDDO 

 1101 FORMAT('Error 1101: Illegal input for boundary condition: 'I3,/  &
         A,' should not be specified for DEM/PIC simulations.',/       &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_WALLS_DISCRETE




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_WALLS_SCALAR_EQ                                 !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Check user-input for generic scalar eq WALL BC parameters.  !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_WALLS_SCALAR_EQ(BCV)


! Global Variables:
!---------------------------------------------------------------------//
! User-input: number of generic scalar equations to solve.
      use scalars, only: NSCALAR
! User-Input: generic scalar eq at wall BCs.
      use bc, only: BC_HW_SCALAR, BC_SCALARw, BC_C_SCALAR

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments.
!---------------------------------------------------------------------//
! Index of BC getting checked.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Loop/counter variable.
      INTEGER :: N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_WALLS_SCALAR_EQ")

      DO N=1, NSCALAR
         IF(BC_HW_Scalar(BCV,N) < ZERO) THEN 
            WRITE(ERR_MSG,1001) trim(iVar('BC_HW_SCALAR',BCV,N)),      &
               trim(iVal(BC_HW_Scalar(BCV,N)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
         IF(BC_HW_SCALAR(BCV,N) /= ZERO .AND.                          &
            BC_SCALARw(BCV,N) == UNDEFINED) THEN 
            WRITE(ERR_MSG,1000) trim(iVar('BC_SCALARw',BCV,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_HW_SCALAR(BCV,N) /= UNDEFINED .AND.                     &
            BC_C_Scalar(BCV,N) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('BC_C_SCALAR',BCV,N))
         ENDIF 
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unphysical input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_WALLS_SCALAR_EQ
