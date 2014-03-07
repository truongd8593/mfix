!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's (GET_FLOW_BC)          !
!     - compute area of boundary surfaces (GET_BC_AREA)                !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     ! 
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BOUNDARY_CONDITIONS

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use physprop, only: SMAX
! Total number of discrete solids.
      use discretelement, only: DES_MMAX
! Type of run: NEW/RESTART
      use run, only: RUN_TYPE
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc, only: BC_TYPE
! User specifed BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED
! Maximum number of BCs
      use param, only: DIMENSION_BC
! Maximum number of disperse phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for BCs
      INTEGER :: BCV
! Total number of solids phases (continuum + discrete)
      INTEGER :: MMAX_TOT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BOUNDARY_CONDITIONS")

! Determine which BCs are DEFINED
      CALL CHECK_BC_GEOMETRY

! Find and validate I/J/K locations for wall BCs.
      CALL GET_WALLS_BC 

! Find and validate I/J/K locations of flow BCs
      CALL GET_FLOW_BC 

! Compute area of boundary surfaces
      CALL GET_BC_AREA 

! Total number of solids.
      MMAX_TOT = SMAX + DES_MMAX

! Loop over each defined BC and check the user data.
      DO BCV = 1, DIMENSION_BC 

         IF (BC_DEFINED(BCV)) THEN

! Determine which solids phases are present.
            SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
               .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

            IF(MMAX_TOT == 1 .AND. BC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

            SELECT CASE (TRIM(BC_TYPE(BCV)))

            CASE ('MASS_INFLOW')
               CALL CHECK_BC_MASS_INFLOW(MMAX_TOT, SKIP, BCV)

            CASE ('P_INFLOW')
               CALL CHECK_BC_P_INFLOW(MMAX_TOT, SKIP, BCV)

            CASE ('OUTFLOW')
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, SKIP, BCV)

            CASE ('MASS_OUTFLOW')
               CALL CHECK_BC_MASS_OUTFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, SKIP, BCV)

            CASE ('P_OUTFLOW')
               CALL CHECK_BC_P_OUTFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, SKIP, BCV)

            CASE ('FREE_SLIP_WALL')
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE ('NO_SLIP_WALL')
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE ('PAR_SLIP_WALL')
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            END SELECT

! Check whether BC values are specified for undefined BC locations
         ELSEIF(BC_TYPE(BCV) /= 'DUMMY' .AND.                          &
            BC_TYPE(BCV)(1:2) /= 'CG') THEN

            CALL CHECK_BC_RANGE(BCV)

         ENDIF
      ENDDO

! Verify that ICBC flags are set for all fluid cells.
      IF (RUN_TYPE(1:3) == 'NEW') CALL CHECK_ICBC_FLAG

! Cleanup and exit.
      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_BOUNDARY_CONDITIONS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Flag: BC contains geometric data and/or specified type
      use bc, only: BC_DEFINED
! User specified BC 
      use bc, only: BC_TYPE
! User specifed: BC geometry
      use bc, only: BC_X_e, BC_X_w, BC_I_e, BC_I_w
      use bc, only: BC_Y_n, BC_Y_s, BC_J_n, BC_J_s
      use bc, only: BC_Z_t, BC_Z_b, BC_K_t, BC_K_b
! User specified: System geometry
      use geometry, only: NO_I, XLENGTH
      use geometry, only: NO_J, YLENGTH
      use geometry, only: NO_K, ZLENGTH

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC
! Parameter constants
      use param1, only: ZERO, UNDEFINED, UNDEFINED_I, UNDEFINED_C

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: BCV, I
! Error flag
      LOGICAL :: RECOGNIZED_BC_TYPE
! Total number of valid BC types
      INTEGER, PARAMETER :: DIM_BCTYPE = 21 
! Valid boundary condition types
      CHARACTER*16, DIMENSION(1:DIM_BCTYPE) ::VALID_BC_TYPE = (/&
           'MASS_INFLOW     ', 'MI              ',&
           'MASS_OUTFLOW    ', 'MO              ',&
           'P_INFLOW        ', 'PI              ',&
           'P_OUTFLOW       ', 'PO              ',&
           'FREE_SLIP_WALL  ', 'FSW             ',&
           'NO_SLIP_WALL    ', 'NSW             ',&
           'PAR_SLIP_WALL   ', 'PSW             ',&
           'OUTFLOW         ', 'OF              ',&
           'CG_NSW          ', 'CG_FSW          ',&
           'CG_PSW          ', 'CG_MI           ',&
           'CG_PO           '/)
!......................................................................!

      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY")


      L50: DO BCV = 1, DIMENSION_BC 

         BC_DEFINED(BCV) = .FALSE.
         IF(BC_X_W(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_X_E(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_Y_S(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_Y_N(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_Z_B(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_Z_T(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_I_W(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_I_E(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_J_S(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_J_N(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_K_B(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_K_T(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE. 
         IF(BC_TYPE(BCV) == 'CG_NSW')   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_TYPE(BCV) == 'CG_FSW')   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_TYPE(BCV) == 'CG_PSW')   BC_DEFINED(BCV) = .TRUE. 
         IF(BC_TYPE(BCV) == 'CG_MI')    BC_DEFINED(BCV) = .TRUE. 
         IF(BC_TYPE(BCV) == 'CG_PO')    BC_DEFINED(BCV) = .TRUE. 

         IF (BC_TYPE(BCV) == 'DUMMY') BC_DEFINED(BCV) = .FALSE. 

         IF(BC_TYPE(BCV)/=UNDEFINED_C .AND. BC_TYPE(BCV)/='DUMMY')THEN

            RECOGNIZED_BC_TYPE = .FALSE.
            DO I = 1, DIM_BCTYPE
                VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
                IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN
                   RECOGNIZED_BC_TYPE = .TRUE.
                   EXIT
                ENDIF 
            ENDDO 
            
            IF(.NOT.RECOGNIZED_BC_TYPE) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)), &
                  trim(BC_TYPE(BCV)), VALID_BC_TYPE
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(.NOT.BC_DEFINED(BCV)) CYCLE
         IF(BC_TYPE(BCV)(1:2) == 'CG') CYCLE

         IF(BC_X_W(BCV)==UNDEFINED .AND. BC_I_W(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN 
               BC_X_W(BCV) = ZERO
            ELSE 
               WRITE(ERR_MSG,1101) BCV, 'BC_X_w and BC_I_w'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_X_E(BCV)==UNDEFINED .AND. BC_I_E(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN 
               BC_X_E(BCV) = XLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1101) BCV, 'BC_X_e and BC_I_e'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Y_S(BCV)==UNDEFINED .AND. BC_J_S(BCV)==UNDEFINED_I) THEN
            IF(NO_J) THEN 
               BC_Y_S(BCV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1101) BCV, 'BC_Y_s and BC_J_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Y_N(BCV)==UNDEFINED .AND. BC_J_N(BCV)==UNDEFINED_I) THEN 
            IF(NO_J) THEN 
               BC_Y_N(BCV) = YLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1101) BCV, 'BC_Y_n and BC_J_n'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Z_B(BCV)==UNDEFINED .AND. BC_K_B(BCV)==UNDEFINED_I) THEN 
            IF(NO_K) THEN 
               BC_Z_B(BCV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1101) BCV, 'BC_Z_b and BC_K_b'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Z_T(BCV)==UNDEFINED .AND. BC_K_T(BCV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               BC_Z_T(BCV) = ZLENGTH
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_Z_t and BC_K_t'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

 1101 FORMAT('Error 1101: Boundary condition ',I3,' is ill-defined.',/ &
         ' > ',A,' are not specified.',/'Please correct the mfix.dat ',&
         'file.')

! Swap BC aliases for the "full name" complement.
         DO I = 1, DIM_BCTYPE 
            VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
            IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN 
               IF(MOD(I,2) == 0) BC_TYPE(BCV) = VALID_BC_TYPE(I-1)
               CYCLE  L50
            ENDIF 
         ENDDO 

         WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)),               &
            trim(BC_TYPE(BCV)), VALID_BC_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ENDDO L50   ! end loop over (bcv=1,dimension_bc)

      CALL FINL_ERR_MSG

      RETURN


 1100 FORMAT('Error 1101: Illegal entry: ',A,' = ',A,/'Valid entries ',&
         'include:'10(/5X,A,2x,A),/5X,A)

      END SUBROUTINE CHECK_BC_GEOMETRY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_RANGE(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase BC varaibles
      use bc, only: BC_EP_g, BC_T_g, BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g
! Sslids phase BC variables.
      USE bc, only: BC_EP_s, BC_ROP_s, BC_T_s, BC_X_s
      use bc, only: BC_U_s, BC_V_s, BC_W_s
! Scalar equation BC variables.
      USE bc, only: BC_SCALAR


! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecifed values.
      use param1, only: UNDEFINED
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIMENSION_N_G, DIMENSION_N_S
! Maximum number of scalar equations.
      use param, only: DIM_SCALAR


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop varaibles.
      INTEGER :: M, N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_RANGE")


! Check gas phase variables.
      IF(BC_U_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_U_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF(BC_V_G(BCV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_V_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
      ENDIF 
      IF (BC_W_G(BCV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_W_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
      ENDIF 
      IF (BC_EP_G(BCV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_EP_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
      ENDIF 
      IF (BC_P_G(BCV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_P_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
      ENDIF 
      IF (BC_T_G(BCV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG,1100) trim(iVar('BC_T_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
      ENDIF 

      DO N = 1, DIMENSION_N_G 
         IF(BC_X_G(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',BCV,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 
      ENDDO 

! Check solids phase variables.
      DO M = 1, DIM_M
         IF(BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG,1100) trim(iVar('BC_ROP_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 
         IF(BC_EP_S(BCV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG,1100) trim(iVar('BC_EP_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 
         IF(BC_U_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_U_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 
         IF(BC_V_S(BCV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG,1100) trim(iVar('BC_V_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 

         IF(BC_W_S(BCV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG,1100) trim(iVar('BC_W_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 
         IF(BC_T_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_T_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
         ENDIF 

         DO N = 1, DIMENSION_N_S 
            IF(BC_X_S(BCV,M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1100) trim(iVar('BC_X_s',BCV,M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
            ENDIF 
         ENDDO 

      ENDDO 

! Check scalar equation variables.
      DO N = 1, DIM_SCALAR
         IF(BC_Scalar(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_Scalar',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDDO 

  
      CALL FINL_ERR_MSG


      RETURN  

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

      END SUBROUTINE CHECK_BC_RANGE




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_ICBC_FLAG                                         !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!  Note that the error message may be incomplete 
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ICBC_FLAG

      use mpi_utility
      use sendrecv

      use error_manager

      IMPLICIT NONE

      LOGICAL :: ERROR = .FALSE.

      INTEGER :: I, J ,K, IJK, IER

      INCLUDE 'function.inc'


      CALL INIT_ERR_MSG("CHECK_ICBC_FLAG")


! First check for any errors.
      DO K = kStart2, kEnd2 
      DO J = jStart2, jEnd2 
      DO I = iStart2, iEnd2
         IF(ICBC_FLAG(FUNIJK(I,J,K)) == '   ') ERROR = .TRUE. 
      ENDDO 
      ENDDO 
      ENDDO 

! Sync up the error flag across all processes.
      CALL GLOBAL_ALL_OR(ERROR)

! If an error is detected, have each rank open a log file and write
! it's own message. Otherwise, we need to send all the data back to
! PE_IO and that's too much work!
      IF(ERROR) THEN

         CALL OPEN_PE_LOG(IER)

         WRITE(ERR_MSG, 1100) trim(iVal(myPE))
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

         DO K = kStart2, kEnd2
         DO J = jStart2, jEnd2 
         DO I = iStart2, iEnd2

            IF(ICBC_FLAG(FUNIJK(I,J,K)) == '   ') THEN
               WRITE(ERR_MSG,1101) I, J, K
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

         ENDDO 
         ENDDO 
         ENDDO 

         WRITE(ERR_MSG, 1102)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

      ELSE
! If no erros, sync up the ghost cell layers.
         CALL SEND_RECV(ICBC_FLAG,2)
      ENDIF

! Clean up and return.
      CALL FINL_ERR_MSG

      RETURN  

 1100 FORMAT('Error 1100 (PE ',A,') : No initial or boundary ',        &
         'condtions specified in','the following cells:',/             &
         '    I       J       K')

 1101 FORMAT(I5,3X,I5,3X,I5) 

 1102 FORMAT('Please correct the mfix.dat file.') 

      END SUBROUTINE CHECK_ICBC_FLAG
