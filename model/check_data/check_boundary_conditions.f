!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DATA_07                                           !
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

      use physprop, only: SMAX

      use run, only: RUN_TYPE
      use discretelement, only: DES_MMAX
      use bc

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER :: BCV

! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
      INTEGER :: MMAX_TOT



      CALL INIT_ERR_MSG("CHECK_BOUNDARY_CONDITIONS")

      CALL CHECK_BC_GEOMETRY

      CALL GET_WALLS_BC 

! Find and validate i, j, k locations of flow BC's
      CALL GET_FLOW_BC 

! Compute area of boundary surfaces
      CALL GET_BC_AREA 

! Total number of solids.
      MMAX_TOT = SMAX + DES_MMAX

      DO BCV = 1, DIMENSION_BC 

! Skip BCs that are not defined.
         IF (BC_DEFINED(BCV)) THEN

! Determine which solids phases are present.
            SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
               .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

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


      CALL CHECK_ICBC_FLAG

      RETURN  


      END SUBROUTINE CHECK_BOUNDARY_CONDITIONS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_GEOMETRY


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

!-----------------------------------------------
! Local Parameters
!-----------------------------------------------
! JFD: for cartesian grid implementation
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! error flag
      LOGICAL :: ERROR
! loop/variable indices
      INTEGER :: BCV
      INTEGER :: I, J ,K, IJK
      INTEGER :: M, N

!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 


! loop/variable indices
      INTEGER, PARAMETER :: DIM_BCTYPE = 21 
! error flag
      LOGICAL :: RECOGNIZED_BC_TYPE


! valid boundary condition types
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
           'CG_PO           '                     &
            /)


      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY")

! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
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

         IF(BC_TYPE(BCV) /= UNDEFINED_C.AND.BC_TYPE(BCV) /= 'DUMMY') THEN

            RECOGNIZED_BC_TYPE = .FALSE.
            DO I = 1, DIM_BCTYPE 
                VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
                IF (VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN 
                   RECOGNIZED_BC_TYPE = .TRUE.
                   EXIT
                ENDIF 
            ENDDO 
            
            IF(.NOT.RECOGNIZED_BC_TYPE) THEN
               WRITE(UNIT_LOG, 1001) BCV, BC_TYPE(BCV) 
               WRITE(UNIT_LOG, 1002) VALID_BC_TYPE 
               call mfix_exit(myPE)
            ENDIF
         ENDIF

! Validate the BC postion for all non cut-cell boundaries:
         IF(.NOT.BC_DEFINED(BCV)) CYCLE
         IF(BC_TYPE(BCV)(1:2) == 'CG') CYCLE

         IF(BC_X_W(BCV)==UNDEFINED .AND. BC_I_W(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN 
               BC_X_W(BCV) = ZERO
            ELSE 
               WRITE(UNIT_LOG,1100) BCV, 'BC_X_w and BC_I_w'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_X_E(BCV)==UNDEFINED .AND. BC_I_E(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN 
               BC_X_E(BCV) = XLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1100) BCV, 'BC_X_e and BC_I_e'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Y_S(BCV)==UNDEFINED .AND. BC_J_S(BCV)==UNDEFINED_I) THEN
            IF(NO_J) THEN 
               BC_Y_S(BCV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1100) BCV, 'BC_Y_s and BC_J_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Y_N(BCV)==UNDEFINED .AND. BC_J_N(BCV)==UNDEFINED_I) THEN 
            IF(NO_J) THEN 
               BC_Y_N(BCV) = YLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1100) BCV, 'BC_Y_n and BC_J_n'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(BC_Z_B(BCV)==UNDEFINED .AND. BC_K_B(BCV)==UNDEFINED_I) THEN 
            IF(NO_K) THEN 
               BC_Z_B(BCV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1100) BCV, 'BC_Z_b and BC_K_b'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Z_T(BCV)==UNDEFINED .AND. BC_K_T(BCV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               BC_Z_T(BCV) = ZLENGTH
            ELSE
               WRITE(ERR_MSG, 1100) BCV, 'BC_Z_t and BC_K_t'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

 1100 FORMAT('Error 1100: Boundary condition ',I3,' is ill-defined.',/ &
         ' > ',A,' are not specified.',/'Please correct the mfix.dat ',&
         'file.')

         DO I = 1, DIM_BCTYPE 
            VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
            IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN 
               IF(MOD(I,2) == 0) BC_TYPE(BCV) = VALID_BC_TYPE(I-1)
               CYCLE  L50
            ENDIF 
         ENDDO 

         IF(DMP_LOG)WRITE (UNIT_LOG, 1001) BCV, BC_TYPE(BCV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1002) VALID_BC_TYPE 
         call mfix_exit(myPE)  


      ENDDO L50   ! end loop over (bcv=1,dimension_bc)

      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_07',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ') 
 1002 FORMAT(5X,A16) 

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

      USE param 
      USE param1 
      USE bc
      USE scalars

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER, INTENT(in) :: BCV

      INTEGER :: M, N


      CALL INIT_ERR_MSG("CHECK_BC_RANGE")


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

      DO N = 1, NScalar
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
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_ICBC_FLAG

      use run, only: RUN_TYPE

      use mpi_utility
      use sendrecv

      use error_manager

      IMPLICIT NONE

      LOGICAL :: ERROR

      INTEGER :: I, J ,K, IJK

      INCLUDE 'function.inc'


      CALL INIT_ERR_MSG("CHECK_ICBC_FLAG")


      ERROR = .FALSE. 

      IF (RUN_TYPE(1:3) == 'NEW') THEN

         DO K = kstart2, kend2 
         DO J = jstart2, Jend2 
         DO I = istart2, Iend2
            IJK = FUNIJK(I,J,K) 
            IF(ICBC_FLAG(IJK) == '   ') THEN
               IF(.NOT.ERROR) THEN
                  WRITE(ERR_MSG, 1100)
                  CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               ENDIF
               WRITE(ERR_MSG,1101) I, J, K
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               ERROR = .TRUE. 
            ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 

! Sync up the error across all processes.
         CALL GLOBAL_ALL_OR(ERROR)

! Finalize any error messages and exit.
         IF(ERROR) THEN
            WRITE(ERR_MSG,"('Aborting.')")
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)
         ENDIF 

! If no erros, sync up the ghost cell layers.
         CALL SEND_RECV(ICBC_FLAG,2)

      ENDIF

      CALL FINL_ERR_MSG

      RETURN  

 1100 FORMAT('Error 1100: No Initial or boundary conditions ',          &
         'specified:',/'    I       J       K')

 1101 FORMAT(I5,3X,I5,3X,I5) 

      END SUBROUTINE CHECK_ICBC_FLAG
