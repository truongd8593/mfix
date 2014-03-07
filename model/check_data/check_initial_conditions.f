!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_INITIAL_CONDITIONS

      use ic, only: IC_DEFINED
      use ic, only: IC_TYPE

      use ic, only: IC_I_W, IC_I_E
      use ic, only: IC_J_S, IC_J_N
      use ic, only: IC_K_B, IC_K_T

      use sendrecv
      use mpi_utility      
      use error_manager

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: ICV
      INTEGER :: I, J, K, IJK

      INCLUDE 'function.inc'

      CALL INIT_ERR_MSG("CHECK_INITIAL_CONDITIONS")

      CALL INIT_ICBC_FLAGS

      CALL CHECK_IC_GEOMETRY

      IC_LP: DO ICV=1, DIMENSION_IC

         IF(IC_DEFINED(ICV)) THEN

! Skip checks for PATCH restarts.
            IF (IC_TYPE(ICV) == 'PATCH') CYCLE IC_LP

            CALL CHECK_IC_GAS_PHASE(ICV)

            CALL CHECK_IC_SOLIDS_PHASES(ICV)

!  Set ICBC flag
            DO K = IC_K_B(ICV), IC_K_T(ICV)
            DO J = IC_J_S(ICV), IC_J_N(ICV) 
            DO I = IC_I_W(ICV), IC_I_E(ICV) 
               IF(.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF(DEAD_CELL_AT(I,J,K)) CYCLE
               IJK = FUNIJK(I,J,K)
               WRITE(ICBC_FLAG(IJK)(1:3),"('.',I2.2)") MOD(ICV,100)
            ENDDO 
            ENDDO 
            ENDDO 

         ELSE

! Check whether physical quantities are specified for undefined
! initial conditions and if so flag error
!            CALL CHECK_IC_OVERFLOW(ICV)

         ENDIF

      ENDDO IC_LP

! Update the ICBC flag on ghost cells.
      CALL SEND_RECV(ICBC_FLAG, 2)

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_INITIAL_CONDITIONS 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: INIT_ICBC_FLAGS                                          !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE INIT_ICBC_FLAGS

      use run, only: RUN_TYPE

      use mpi_utility

      implicit none

      INTEGER :: I, J, K, IJK

      include 'function.inc'

! Initialize the icbc_flag array.
      DO K = kStart3, kEnd3       
      DO J = jStart3, jEnd3 
      DO I = iStart3, iEnd3

         IJK = FUNIJK(I,J,K)

! Initialize the ICBC Flag
         ICBC_FLAG(IJK) = merge('   ', '.--', RUN_TYPE == 'NEW')


! If at domain boundaries then set default values (wall or, if
! specified, cyclic)
         IF (DO_K) THEN 
            IF(K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)THEN
               IF (CYCLIC_Z_PD) THEN 
                  ICBC_FLAG(IJK) = 'C--' 
               ELSEIF (CYCLIC_Z) THEN 
                  ICBC_FLAG(IJK) = 'c--' 
               ELSE 
                  ICBC_FLAG(IJK) = 'W--' 
               ENDIF 
            ENDIF 
         ENDIF 

         IF(DO_J)THEN 
            IF(J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMAX3)THEN 
               IF (CYCLIC_Y_PD) THEN 
                  ICBC_FLAG(IJK) = 'C--' 
               ELSEIF (CYCLIC_Y) THEN 
                  ICBC_FLAG(IJK) = 'c--' 
               ELSE 
                 ICBC_FLAG(IJK) = 'W--' 
               ENDIF 
            ENDIF
         ENDIF 

         IF(DO_I)THEN 
            IF(I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3)THEN 
               IF (CYCLIC_X_PD) THEN 
                  ICBC_FLAG(IJK) = 'C--' 
               ELSEIF (CYCLIC_X) THEN 
                  ICBC_FLAG(IJK) = 'c--' 
               ELSE 
                  ICBC_FLAG(IJK) = 'W--' 
               ENDIF 
            ENDIF 
            IF (I==1 .AND. CYLINDRICAL .AND. XMIN==ZERO) &
               ICBC_FLAG(IJK) = 'S--'
         ENDIF 

! corner cells are wall cells
         IF ((I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3) .AND. &
             (J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMIN3) .AND. &
             (K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)) THEN 
            IF (ICBC_FLAG(IJK) /= 'S--') ICBC_FLAG(IJK) = 'W--' 
         ENDIF 
       
      ENDDO ! end do loop (i=istart3, iend3)
      ENDDO ! end do loop (j=jstart3, jend3)
      ENDDO ! end do loop (k=kstart3, kend3)

      RETURN

      END SUBROUTINE INIT_ICBC_FLAGS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GEOMETRY

      use run, only: RUN_TYPE
      use param, only: DIMENSION_IC

      use ic

      use mpi_utility
      use error_manager

      implicit none


      INTEGER :: ICV
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t


      CALL INIT_ERR_MSG("CHECK_IC_GEOMETRY")

! Check geometry of any specified IC region      
      DO ICV = 1, DIMENSION_IC 

         IC_DEFINED(ICV) = .FALSE. 

         IF (IC_X_W(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_X_E(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_S(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_N(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_B(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_T(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 

         IF (IC_I_W(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_I_E(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_S(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_N(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_B(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_T(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 


! For restart runs IC is defined only if IC_TYPE='PATCH'
         IF(RUN_TYPE/='NEW' .AND. IC_TYPE(ICV)/='PATCH') &
            IC_DEFINED(ICV) = .FALSE. 

         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IC_X_W(ICV)==UNDEFINED .AND. IC_I_W(ICV)==UNDEFINED_I) THEN 
            IF (NO_I) THEN 
               IC_X_W(ICV) = ZERO 
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_X_w and IC_I_w'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF (IC_X_E(ICV)==UNDEFINED .AND. IC_I_E(ICV)==UNDEFINED_I) THEN 
            IF (NO_I) THEN 
               IC_X_E(ICV) = XLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1100) ICV, 'IC_X_e and IC_I_e'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF (IC_Y_S(ICV)==UNDEFINED .AND. IC_J_S(ICV)==UNDEFINED_I) THEN 
            IF (NO_J) THEN 
               IC_Y_S(ICV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1100) ICV, 'IC_Y_s and IC_J_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF (IC_Y_N(ICV)==UNDEFINED .AND. IC_J_N(ICV)==UNDEFINED_I) THEN 
            IF (NO_J) THEN 
               IC_Y_N(ICV) = YLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1100) ICV, 'IC_Y_n and IC_J_n'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF

         IF (IC_Z_B(ICV)==UNDEFINED .AND. IC_K_B(ICV)==UNDEFINED_I) THEN 
            IF (NO_K) THEN 
               IC_Z_B(ICV) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1100) ICV, 'IC_Z_b and IC_K_b'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF (IC_Z_T(ICV)==UNDEFINED .AND. IC_K_T(ICV)==UNDEFINED_I) THEN 
            IF (NO_K) THEN 
               IC_Z_T(ICV) = ZLENGTH 
            ELSE 
               WRITE(ERR_MSG, 1100) ICV, 'IC_Z_t and IC_K_t'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

      ENDDO   ! end loop over (icv = 1,dimension_ic)

 1100 FORMAT('Error 1100: Initial condition region ',I3,' is ill-',    &
         'defined.',/' > ',A,' are not specified.',/'Please correct ', &
         'the mfix.dat file.')


      DO ICV = 1, DIMENSION_IC

! Skip this check if the IC region is not specified.
         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IC_X_W(ICV)/=UNDEFINED .AND. IC_X_E(ICV)/=UNDEFINED) THEN 
            IF (NO_I) THEN 
               I_W = 1 
               I_E = 1 
            ELSE 
               CALL CALC_CELL (XMIN, IC_X_W(ICV), DX, IMAX, I_W) 
               I_W = I_W + 1 
               CALL CALC_CELL (XMIN, IC_X_E(ICV), DX, IMAX, I_E) 
            ENDIF 
            IF (IC_I_W(ICV)/=UNDEFINED_I .OR. IC_I_E(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_I_W(ICV), I_W, ICV, 'IC - west') 
               CALL LOCATION_CHECK (IC_I_E(ICV), I_E, ICV, 'IC - east') 
            ELSE 
               IC_I_W(ICV) = I_W 
               IC_I_E(ICV) = I_E 
            ENDIF 
         ENDIF 

! Report problems with calculated bounds.
         IF(IC_I_W(ICV) > IC_I_E(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IC_I_E'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_E < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_Z_t and IC_K_t'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IC_Y_S(ICV)/=UNDEFINED .AND. IC_Y_N(ICV)/=UNDEFINED) THEN 
            IF (NO_J) THEN 
               J_S = 1 
               J_N = 1 
            ELSE 
               CALL CALC_CELL (ZERO, IC_Y_S(ICV), DY, JMAX, J_S) 
               J_S = J_S + 1 
               CALL CALC_CELL (ZERO, IC_Y_N(ICV), DY, JMAX, J_N) 
            ENDIF 
            IF (IC_J_S(ICV)/=UNDEFINED_I .OR. IC_J_N(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_J_S(ICV), J_S, ICV, 'IC - south') 
               CALL LOCATION_CHECK (IC_J_N(ICV), J_N, ICV, 'IC - north') 
            ELSE 
               IC_J_S(ICV) = J_S 
               IC_J_N(ICV) = J_N 
            ENDIF 
         ENDIF 

         IF(IC_J_S(ICV) > IC_J_N(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S > IC_J_N'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S >  JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N > JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


         IF (IC_Z_B(ICV)/=UNDEFINED .AND. IC_Z_T(ICV)/=UNDEFINED) THEN 
            IF (NO_K) THEN 
               K_B = 1 
               K_T = 1 
            ELSE 
               CALL CALC_CELL (ZERO, IC_Z_B(ICV), DZ, KMAX, K_B) 
               K_B = K_B + 1 
               CALL CALC_CELL (ZERO, IC_Z_T(ICV), DZ, KMAX, K_T) 
            ENDIF 
            IF (IC_K_B(ICV)/=UNDEFINED_I .OR. IC_K_T(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_K_B(ICV), K_B, ICV, 'IC - bottom') 
               CALL LOCATION_CHECK (IC_K_T(ICV), K_T, ICV, 'IC - top') 
            ELSE 
               IC_K_B(ICV) = K_B 
               IC_K_T(ICV) = K_T 
            ENDIF 
         ENDIF 

         IF(IC_K_B(ICV) > IC_K_T(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > IC_K_T'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


 1101 FORMAT('Error 1101: Initial condition region ',I2,' is ill-',    &
         'defined.',/' > ',A,/'Please correct the mfix.dat file.')


      ENDDO   ! end loop over (icv=1,dimension_ic)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_IC_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GAS_PHASE                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GAS_PHASE(ICV)

      use run, only: RUN_TYPE
      use run, only: ENERGY_EQ
      use run, only: SPECIES_EQ
      use run, only: K_Epsilon

      use physprop, only: RO_G0
      use physprop, only: MU_G0
      use physprop, only: MW_AVG
      use physprop, only: NMAX
      use param, only: DIMENSION_IC
      use scalars, only: NSCALAR
      use ic

      use param1, only: ONE
      use param1, only: ZERO
      use param1, only: UNDEFINED
      use geometry, only: NO_I, NO_J, NO_K

      use error_manager

      implicit none


      INTEGER, INTENT(IN) :: ICV

      INTEGER :: N

      DOUBLE PRECISION :: SUM

      LOGICAL, EXTERNAL :: COMPARE 


      CALL INIT_ERR_MSG("CHECK_IC_GAS_PHASE")



! Check that gas phase velocity components are initialized
      IF(IC_U_G(ICV) == UNDEFINED) THEN 
         IF(NO_I) THEN 
            IC_U_G(ICV) = ZERO 
         ELSE 
            WRITE(ERR_MSG, 1000) trim(iVar('IC_U_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

      IF(IC_V_G(ICV) == UNDEFINED) THEN
         IF(NO_J) THEN 
            IC_V_G(ICV) = ZERO 
         ELSE 
            WRITE(ERR_MSG, 1000) trim(iVar('IC_V_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

      IF(IC_W_G(ICV) == UNDEFINED) THEN
         IF (NO_K) THEN 
            IC_W_G(ICV) = ZERO 
         ELSE 
            WRITE(ERR_MSG, 1000) trim(iVar('IC_W_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

! check that gas phase void fraction is initialized
      IF(IC_EP_G(ICV) == UNDEFINED) THEN 
         WRITE(ERR_MSG, 1000) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

! Check that if the gas phase pressure is initialized and the gas is
! compressible that the gas phase pressure is not zero or negative
      IF(IC_P_G(ICV) /= UNDEFINED) THEN
         IF(RO_G0==UNDEFINED .AND. IC_P_G(ICV)<=ZERO) THEN 
            WRITE(ERR_MSG, 1100) trim(iVar('IC_P_g',ICV)),             &
               iVal(IC_P_G(ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

 1100 FORMAT('Error 1100: Pressure must be greater than 0.0 for ',     &
         'compressible flow',/'Illegal value: ',A,' = ',A,/'Please ',  &
         'correct the mfix.dat file.')


      IF(ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED) THEN
         IF(IC_T_G(ICV)==UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) trim(iVar('IC_T_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

! Gas phase radiation values.
      IF (ENERGY_EQ) THEN 
         IF (IC_GAMA_RG(ICV) < ZERO) THEN 
            WRITE(ERR_MSG, 1001) trim(iVar('IC_GAMA_Rg',ICV)),         &
               iVal(IC_GAMA_RG(ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (IC_GAMA_RG(ICV) > ZERO) THEN 
            IF (IC_T_RG(ICV) == UNDEFINED) THEN 
               WRITE(ERR_MSG, 1000) iVar('IC_T_Rg',ICV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 
      ENDIF 


! First: sum together defiend gas phase species mass fractions.
      SUM = ZERO
      DO N = 1, NMAX(0)
         IF(IC_X_G(ICV,N) /= UNDEFINED) THEN
            SUM = SUM + IC_X_G(ICV,N)
         ELSE
            IC_X_G(ICV,N) = ZERO
         ENDIF
      ENDDO 


! Enforce that the species mass fractions must sum to one.
      IF(.NOT.COMPARE(ONE,SUM)) THEN

         IF(SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1110) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1110 FORMAT('Error 1110: IC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'species equations are solved. Please correct ', &
         'the mfix.dat file.')

         ELSEIF(RO_G0 == UNDEFINED .AND. MW_AVG == UNDEFINED) THEN
            WRITE(ERR_MSG, 1111) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1111 FORMAT('Error 1111: IC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'is compressible and MW_AVG is UNDEFINED.',/     &
         'Please correct the mfix.dat the mfix.dat file.')

         ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
            WRITE(ERR_MSG, 1112) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1112 FORMAT('Error 1112: IC_X_g(',I3,',:) do not sum to ONE or ZERO ',&
         'and they',/'are not needed. Please correct the mfix.dat ',   &
         'the mfix.dat file.')

         ELSE
            IC_X_G(ICV,:) = ZERO
            IC_X_G(ICV,1) = ONE
         ENDIF
      ENDIF 


      DO N = 1, NScalar
         IF(IC_Scalar(ICV,N) == UNDEFINED) IC_Scalar(ICV,N) = ZERO 
      ENDDO  


      IF(K_Epsilon) THEN
         IF (IC_K_Turb_G(ICV) == UNDEFINED) THEN 
            WRITE(ERR_MSG, 1000) iVar('IC_K_Turb_G',ICV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF (IC_E_Turb_G(ICV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) iVar('IC_E_Turb_G',ICV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_GAS_PHASE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_SOLIDS_PHASES(ICV)

!-----------------------------------------------
! Modules
!-----------------------------------------------



      use run, only: ENERGY_EQ
      use run, only: SPECIES_EQ
      use run, only: GRANULAR_ENERGY
      use run, only: SOLVE_ROs
      use physprop, only: SMAX
      use physprop, only: INERT_SPECIES, X_S0
      use physprop, only: NMAX
      use physprop, only: BASE_ROs
      use physprop, only: RO_S0
      use discretelement, only: DES_MMAX
      use param1, only: ONE
      use param1, only: ZERO
      use param1, only: UNDEFINED
      use geometry, only: NO_I, NO_J, NO_K


      USE ic
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICV

      INTEGER :: M, N
      DOUBLE PRECISION SUM, SUM_EP
! Solids phase density in IC region.
      DOUBLE PRECISION :: IC_ROs(1:DIM_M)
! Index of inert species
      INTEGER :: INERT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
! Total number of solids phases (TFM + DEM + MPPIC)
      INTEGER :: MMAX_TOT

!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSS
      LOGICAL, EXTERNAL :: COMPARE 


      CALL INIT_ERR_MSG("CHECK_IC_SOLIDS_PHASES")

! The total number of solids phases (all models).
      MMAX_TOT = SMAX + DES_MMAX

! Calculate EP_s from EP_g if there is only one solids phase.
      IF(MMAX_TOT == 1 .AND. IC_EP_S(ICV,1) == UNDEFINED) THEN
         IC_EP_S(ICV,1) = ONE - IC_EP_g(ICV)
      ENDIF

! Bulk density or solids volume fraction must be explicitly defined 
! if there are more than one solids phase.
      IF(MMAX_TOT > 1 .AND. .NOT.COMPARE(IC_EP_g(ICV),ONE)) THEN
         DO M = 1, SMAX + DES_MMAX
            IF(IC_ROP_S(ICV,M) == UNDEFINED .AND. &
               IC_EP_S(ICV,M) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF

 1400 FORMAT('Error 1400: Insufficient solids phase ',I2,' ',          &
         'information for IC',/'region ',I3,'. ',A,' not specified.',/ &
         'Please correct the mfix.dat file.')

! Determine which solids phases are present.
      DO M = 1, MMAX_TOT
         SKIP(M)=(IC_ROP_S(ICV,M)==UNDEFINED.OR.IC_ROP_S(ICV,M)==ZERO) &
            .AND.(IC_EP_S(ICV,M)==UNDEFINED .OR.IC_EP_S(ICV,M)==ZERO)
      ENDDO

      IF(MMAX_TOT == 1 .AND. IC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

      DO M=1, MMAX_TOT

! check that solids phase m velocity components are initialized
         IF(IC_U_S(ICV,M) == UNDEFINED) THEN 
            IF (SKIP(M) .OR. NO_I) THEN
               IC_U_S(ICV,M) = ZERO 
            ELSE
               WRITE(ERR_MSG, 1000)trim(iVar('IC_U_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(IC_V_S(ICV,M) == UNDEFINED) THEN 
            IF(SKIP(M) .OR. NO_J) THEN 
               IC_V_S(ICV,M) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1000)trim(iVar('IC_V_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(IC_W_S(ICV,M) == UNDEFINED) THEN 
            IF(SKIP(M) .OR. NO_K) THEN
               IC_W_S(ICV,M) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1000)trim(iVar('IC_W_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(ENERGY_EQ .AND. IC_T_S(ICV,M)==UNDEFINED) THEN
            IF(SKIP(M)) THEN
               IC_T_S(ICV,M) = IC_T_G(ICV)
            ELSE 
               WRITE(ERR_MSG, 1000)trim(iVar('IC_T_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

         IF(GRANULAR_ENERGY .AND. IC_THETA_M(ICV,M)==UNDEFINED) THEN
            IF(SKIP(M)) THEN 
               IC_THETA_M(ICV,M) = ZERO 
            ELSE 
               WRITE(ERR_MSG, 1000)trim(iVar('IC_Theta_M',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 
   
         IF(ENERGY_EQ) THEN
            IF(IC_GAMA_RS(ICV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1001)trim(iVar('IC_GAMA_Rs',ICV,M)),     &
                  iVal(IC_GAMA_RS(ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF (IC_GAMA_RS(ICV,M) > ZERO) THEN 
               IF(IC_T_RS(ICV,M) == UNDEFINED) THEN
                  WRITE(ERR_MSG, 1001)trim(iVar('IC_T_Rs',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF 
            ENDIF 
         ENDIF 


! First: sum together defiend species mass fractions.
         SUM = ZERO
         DO N = 1, NMAX(M)
            IF(IC_X_S(ICV,M,N) /= UNDEFINED) THEN
               SUM = SUM + IC_X_S(ICV,M,N)
            ELSE
               IC_X_S(ICV,M,N) = ZERO
            ENDIF
         ENDDO 

! Enforce that the species mass fractions must sum to one.
         IF(.NOT.COMPARE(ONE,SUM)) THEN

            IF(SPECIES_EQ(M)) THEN
               WRITE(ERR_MSG, 1402) ICV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1402 FORMAT('Error 1402: IC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'and the solids phase',/'species equations are solved. ',     &
         'Please correct the mfix.dat file.')

            ELSEIF(SOLVE_ROS(M)) THEN
               WRITE(ERR_MSG, 1403) ICV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1403 FORMAT('Error 1403: IC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'and the solids phase',/'density is calculated. Please ',     &
         'correct the mfix.dat file.')

            ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
                WRITE(ERR_MSG, 1404) ICV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1404 FORMAT('Error 1404: IC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'or ZERO and',/'they are not needed. Please correct the ',    &
         'mfix.dat file.')

            ELSE
               IC_X_S(ICV,M,:) = ZERO
               IC_X_S(ICV,M,1) = ONE
            ENDIF
         ENDIF 

! Set the solids density for the IC region.
         IF(SKIP(M)) THEN
            IC_ROs(M) = merge(BASE_ROs(M), RO_s0(M), SOLVE_ROs(M))

         ELSEIF(.NOT.SOLVE_ROs(M)) THEN
            IC_ROs(M) = RO_s0(M)

         ELSE
! Verify that the species mass fraction for the inert material is not
! zero in the IC region when the solids is present.
            INERT = INERT_SPECIES(M)
            IF(IC_X_S(ICV,M,INERT) == ZERO) THEN
               WRITE(ERR_MSG,1405) M, ICV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1405 FORMAT('Error 1405: No inert species for phase ',I2,' in IC ',   &
         'region',I3,'.',/'Unable to calculate solids phase density. ',&
         'Please refer to the Readme',/' file for required variable ', &
         'soilds density  model input parameters and',/' make the ',   &
         'necessary corrections to the data file.')

! Calculate the solids density.
            IC_ROs(M) = EOSS(BASE_ROs(M), X_s0(M,INERT),               &
               IC_X_S(ICV,M,INERT))
         ENDIF

      ENDDO   ! end loop over (m=1,smax)


! Initialize the sum of the total volume fraction.
      SUM_EP = IC_EP_G(ICV)

      DO M=1, MMAX_TOT

! Clear out both varaibles if this phase is skipped.
         IF(SKIP(M)) THEN
            IC_EP_S(ICV,M)  = ZERO
            IC_ROP_S(ICV,M) = ZERO

! If both input parameters are defined. Make sure they are equivalent.
         ELSEIF(IC_ROP_S(ICV,M) /= UNDEFINED .AND.                     &
            IC_EP_S(ICV,M) /= UNDEFINED) THEN

            IF(.NOT.COMPARE(IC_EP_S(ICV,M)*IC_ROs(M),                  &
               IC_ROP_S(ICV,M))) THEN
               WRITE(ERR_MSG,1406)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1406 FORMAT('Error 1406: Illegal initial condition region : ',I3,/    &
         'IC_EP_s and IC_ROP_s are inconsistent. Please correct the ',/&
         'mfix.dat file.')

! Compute IC_EP_s from IC_ROP_s
         ELSEIF(IC_EP_S(ICV,M) == UNDEFINED) THEN
            IC_EP_S(ICV,M) = IC_ROP_S(ICV,M) / IC_ROs(M)

! Compute IC_ROP_s from IC_EP_s and IC_ROs
         ELSEIF(IC_ROP_S(ICV,M) == UNDEFINED) THEN
            IC_ROP_S(ICV,M) = IC_EP_S(ICV,M) * IC_ROs(M)
! This is a sanity check.
         ELSE

         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + IC_EP_S(ICV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(.NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1410)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1410 FORMAT('Error 1410: Illegal initial condition region : ',I3,/    &
         'Sum of volume fractions does NOT equal ONE. Please correct',/&
         'the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal data in initial condition region ', &
         I3,/A,' defined for unknown solids phase ',I2,'.',/          &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_SOLIDS_PHASES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: INIT_ICBC_FLAGS                                          !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_OVERFLOW(ICV)

      use scalars
      use physprop
      use run
      use ic
      use indices
      use geometry
      use compar
      use discretelement

      use error_manager

      implicit none

! IC region index.
      INTEGER, INTENT(IN) :: ICV

      INTEGER :: M, N
      INTEGER :: MMAX_TOT

      CALL INIT_ERR_MSG("CHECK_IC_OVERFLOW")


! GAS PHASE quantities
! -------------------------------------------->>>                
      IF(IC_U_G(ICV) /= UNDEFINED) THEN 
          WRITE(ERR_MSG, 1010) trim(iVar('IC_U_g',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 
      IF(IC_V_G(ICV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG, 1010) trim(iVar('IC_V_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 
      IF(IC_W_G(ICV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG, 1010) trim(iVar('IC_W_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 
      IF(IC_EP_G(ICV) /= UNDEFINED) THEN 
         WRITE(ERR_MSG, 1010) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 
      IF(IC_T_G(ICV) /= UNDEFINED) THEN 
          WRITE(ERR_MSG, 1010) trim(iVar('IC_T_g',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 
      IF(IC_T_RG(ICV) /= UNDEFINED) THEN 
          WRITE(ERR_MSG, 1010) trim(iVar('IC_T_Rg',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
      ENDIF 

      DO N = 1, DIMENSION_N_G 
         IF(IC_X_G(ICV,N) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_X_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
      ENDDO 
    
      DO N = 1, NScalar 
         IF(IC_Scalar(ICV,N) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_Scalar',ICV)) 
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDDO 
   
      IF(K_Epsilon) THEN
         IF(IC_K_Turb_G(ICV) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_K_Turb_G',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
    
         IF(IC_E_Turb_G(ICV) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_E_Turb_G',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF


! SOLIDS PHASE quantities
! -------------------------------------------->>>
      MMAX_TOT = SMAX + DES_MMAX

      DO M=1, MMAX_TOT

         IF(IC_ROP_S(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_ROP_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
         DO N = 1, DIMENSION_N_S 
            IF(IC_X_S(ICV,M,N) /= UNDEFINED) THEN 
               WRITE(ERR_MSG, 1010) trim(iVar('IC_X_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
            ENDIF 
         ENDDO 
         IF(IC_U_S(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_U_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
         IF(IC_V_S(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_V_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
         IF(IC_W_S(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_W_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
         IF(IC_T_S(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_T_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
         IF(IC_T_RS(ICV,M) /= UNDEFINED) THEN 
            WRITE(ERR_MSG, 1010) trim(iVar('IC_T_Rs',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.) 
         ENDIF 
      ENDDO

      CALL FINL_ERR_MSG
      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined IC region') 


      END SUBROUTINE CHECK_IC_OVERFLOW
