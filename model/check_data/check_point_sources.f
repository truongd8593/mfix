!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_POINT_SOURCES                                     !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_POINT_SOURCES

      use param
      use run
      use ps

      use error_manager

      implicit none

      INTEGER :: PSV

      CALL INIT_ERR_MSG("CHECK_POINT_SOURCES")

      POINT_SOURCE = .FALSE.
      PS_DEFINED = .FALSE.

      DO PSV = 1, DIMENSION_PS
         IF(PS_DEFINED(PSV)) THEN
            CALL GET_PS(PSV)
            CALL CHECK_PS_GAS_PHASE(PSV)
            CALL CHECK_PS_SOLIDS_PHASES(PSV)
         ELSE
            CALL CHECK_PS_OVERFLOW(PSV)
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_POINT_SOURCES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_GEOMETRY                                       !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_GEOMETRY

      use param
      use run, only: SPECIES_EQ
      use physprop, only: MMAX
      use physprop, only: NMAX

      use run
      use rxns
      use ps
      use compar
      use geometry
      use mpi_utility

      use error_manager

      implicit none

      INTEGER :: IJK, I, J, K, M, N

      INTEGER PSV

      CHARACTER*64 eMsg

      DOUBLE PRECISION lSum

      INTEGER :: iErr

! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
      PSV_LP: do PSV = 1, DIMENSION_PS

         IF (PS_X_W(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_X_E(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_Y_S(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_Y_N(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_Z_B(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_Z_T(PSV) /= UNDEFINED)   PS_DEFINED(PSV) = .TRUE. 
         IF (PS_I_W(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         IF (PS_I_E(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         IF (PS_J_S(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         IF (PS_J_N(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         IF (PS_K_B(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         IF (PS_K_T(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 

         IF (.NOT.PS_DEFINED(PSV)) cycle PSV_LP

! Flag that one or more point sources has been detected.
         POINT_SOURCE = .TRUE.

         IF(PS_X_W(PSV)==UNDEFINED .AND. PS_I_W(PSV)==UNDEFINED_I) THEN
            IF(NO_I) THEN
               PS_X_W(PSV) = ZERO 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_X_w and PS_I_w '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)  
            ENDIF
         ENDIF
         IF(PS_X_E(PSV)==UNDEFINED .AND. PS_I_E(PSV)==UNDEFINED_I) THEN
            IF(NO_I) THEN 
               PS_X_E(PSV) = XLENGTH 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_X_e and PS_I_e '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)      
            ENDIF 
         ENDIF 
         IF(PS_Y_S(PSV)==UNDEFINED .AND. PS_J_S(PSV)==UNDEFINED_I) THEN
            IF(NO_J) THEN 
               PS_Y_S(PSV) = ZERO 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_Y_s and PS_J_s '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)      
            ENDIF 
         ENDIF 
         IF(PS_Y_N(PSV)==UNDEFINED .AND. PS_J_N(PSV)==UNDEFINED_I) THEN 
            IF(NO_J) THEN 
               PS_Y_N(PSV) = YLENGTH 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_Y_n and PS_J_n '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 
         IF(PS_Z_B(PSV)==UNDEFINED .AND. PS_K_B(PSV)==UNDEFINED_I) THEN
            IF(NO_K) THEN 
               PS_Z_B(PSV) = ZERO 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_Z_b and PS_K_b '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 
         IF(PS_Z_T(PSV)==UNDEFINED .AND. PS_K_T(PSV)==UNDEFINED_I) THEN
            IF(NO_K) THEN 
               PS_Z_T(PSV) = ZLENGTH 
            ELSE 
               WRITE(ERR_MSG,1101) PSV, 'PS_Z_t and PS_K_t '
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF 
         ENDIF 

 1101 FORMAT('Error 1101: Point source ',I3,' is ill-defined.',/A,     &
         ' are not specified.',/'Please correct the mfix.dat file.')


      ENDDO PSV_LP


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_PS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_GAS_PHASE                                      !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_GAS_PHASE(PSV)

      use param
      use run
      use ps

      use physprop, only: NMAX
      use error_manager

      implicit none

      INTEGER, INTENT(in) :: PSV
      INTEGER :: N
      DOUBLE PRECISION :: SUM

      LOGICAL, EXTERNAL :: COMPARE


      CALL INIT_ERR_MSG("CHECK_PS_GAS_PHASE")



! Mass flow is undefined --> Velocity must also be undefined.
!```````````````````````````````````````````````````````````````````````
      IF(PS_MASSFLOW_G(PSV) == UNDEFINED) THEN
         IF(PS_U_g(PSV) /= UNDEFINED .OR. &
            PS_V_g(PSV) /= UNDEFINED .OR. &
            PS_W_g(PSV) /= UNDEFINED) THEN

            WRITE(ERR_MSG,1100) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Invalid specification for point source ',I3,&
         '.',/A,' is undefined but velocity is given.',/'Please ',    &
         'correct the mfix.dat file.')

         ELSE
            PS_MASSFLOW_G(PSV) = ZERO
            PS_U_g(PSV) = ZERO
            PS_V_g(PSV) = ZERO
            PS_W_g(PSV) = ZERO
         ENDIF

! Mass flow is zero --> Velocity must also be zero.
!```````````````````````````````````````````````````````````````````````
      ELSEIF(PS_MASSFLOW_G(PSV) == ZERO) THEN
         IF(PS_U_g(PSV) /= ZERO .OR. &
            PS_V_g(PSV) /= ZERO .OR. &
            PS_W_g(PSV) /= ZERO) THEN

            WRITE(ERR_MSG,1101) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: Invalid specification for point source ',I3,&
         '.',/A,' is zero but velocity is given.',/'Please correct ', &
         'the mfix.dat file.')

! Mass flow is negative --> ERROR
!```````````````````````````````````````````````````````````````````````
      ELSEIF(PS_MASSFLOW_G(PSV) < ZERO) THEN
         WRITE(ERR_MSG,1102) PSV, trim(iVar('PS_MASSFLOW_G',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1102 FORMAT('Error 1102: Invalid specifications for point source ',I3,&
         '.',/A,' < 0.0. Point sources can only add mass to a system',/&
         'Please correct the mfix.dat file.')


! Mass flow is specified:
! Velocity does not have to be defined (no momentum source). If the
! components are UNDEFINED, zero them out.
!```````````````````````````````````````````````````````````````````````
      ELSE
         IF(PS_U_g(PSV) == UNDEFINED) PS_U_g(PSV) = ZERO
         IF(PS_V_g(PSV) == UNDEFINED) PS_V_g(PSV) = ZERO
         IF(PS_W_g(PSV) == UNDEFINED) PS_W_g(PSV) = ZERO


! Sum together defiend gas phase species mass fractions.
         SUM = ZERO
         DO N = 1, NMAX(0)
            IF(PS_X_G(PSV,N) /= UNDEFINED) THEN
               SUM = SUM + PS_X_G(PSV,N)
            ELSE
               PS_X_G(PSV,N) = ZERO
            ENDIF
         ENDDO 

! Enforce that the species mass fractions must sum to one.
         IF(.NOT.COMPARE(ONE,SUM)) THEN

            IF(SPECIES_EQ(0)) THEN
               WRITE(ERR_MSG, 1110) PSV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1110 FORMAT('Error 1110: PS_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'species equations are solved. Please correct ', &
         'the mfix.dat file.')

            ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
               WRITE(ERR_MSG, 1111) PSV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1111 FORMAT('Error 1111: PS_X_g(',I3,',:) do not sum to ONE or ZERO ',&
         'and they',/'are not needed. Please correct the mfix.dat ',   &
         'the mfix.dat file.')

            ELSE
               PS_X_G(PSV,:) = ZERO
               PS_X_G(PSV,1) = ONE
            ENDIF

         ENDIF 

! Verify that a temperature is provided.
         IF(ENERGY_EQ)THEN
            IF(PS_T_g(PSV) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('PS_T_g',PSV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Verify that a given temperature is physical.
            ELSEIF(PS_T_g(PSV) <= ZERO) THEN
               WRITE(ERR_MSG,1001) PSV, trim(iVar('PS_T_g',PSV)),      &
                  trim(iVal(PS_T_g(PSV)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_PS_GAS_PHASE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_SOLIDS_PHASES                                  !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_SOLIDS_PHASES(PSV)

      use param
      use run
      use ps

      use error_manager
      use discretelement, only: DES_MMAX
      use physprop, only: NMAX
      use physprop, only: SMAX

      use error_manager

      implicit none

      INTEGER, INTENT(in) :: PSV


      INTEGER :: M, N
      INTEGER :: MMAX_TOT

      DOUBLE PRECISION :: SUM

      LOGICAL, EXTERNAL :: COMPARE


      CALL INIT_ERR_MSG("CHECK_PS_SOLIDS_PHASES")

! The total number of solids phases (all models).
      MMAX_TOT = SMAX + DES_MMAX

      DO M=1, MMAX_TOT

! Mass flow is undefined --> Velocity must also be undefined.
!```````````````````````````````````````````````````````````````````````
         IF(PS_MASSFLOW_S(PSV,M) == UNDEFINED) THEN
            IF(PS_U_s(PSV,M) /= UNDEFINED .OR. &
               PS_V_s(PSV,M) /= UNDEFINED .OR. &
               PS_W_s(PSV,M) /= UNDEFINED) THEN

               WRITE(ERR_MSG,1100)PSV, trim(iVar('PS_MASSFLOW_S',PSV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Invalid specification for point source ',I3,&
         '.',/A,' is undefined but velocity is given.',/'Please ',    &
         'correct the mfix.dat file.')

            ELSE
               PS_MASSFLOW_S(PSV,M) = ZERO
               PS_U_s(PSV,M) = ZERO
               PS_V_s(PSV,M) = ZERO
               PS_W_s(PSV,M) = ZERO
            ENDIF

! Mass flow is zero --> Velocity must also be zero.
!```````````````````````````````````````````````````````````````````````
         ELSEIF(PS_MASSFLOW_S(PSV,M) == ZERO) THEN
            IF(PS_U_s(PSV,M) /= ZERO .OR. &
               PS_V_s(PSV,M) /= ZERO .OR. &
               PS_W_s(PSV,M) /= ZERO) THEN

               WRITE(ERR_MSG,1101)PSV, trim(iVar('PS_MASSFLOW_S',PSV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1101 FORMAT('Error 1100: Invalid specification for point source ',I3,&
         '.',/A,' is zero but velocity is given.',/'Please correct ', &
         'the mfix.dat file.')

! Mass flow is negative --> ERROR
!```````````````````````````````````````````````````````````````````````
         ELSEIF(PS_MASSFLOW_S(PSV,M) < ZERO) THEN
            WRITE(ERR_MSG,1102) PSV, trim(iVar('PS_MASSFLOW_S',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1102 FORMAT('Error 1102: Invalid specifications for point source ',I3,&
         '.',/A,' < 0.0. Point sources can only add mass to a system',/&
         'Please correct the mfix.dat file.')


! Mass flow is specified:
!```````````````````````````````````````````````````````````````````````
         ELSE

! Currently, only TFM solids can be used with point sources. However, 
! the could be implemented for PIC solids as well.
            SELECT CASE(SOLIDS_MODEL(M))
            CASE ('DEM','PIC')
               WRITE(ERR_MSG, 1110) PSV, SOLIDS_MODEL(M)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            CASE DEFAULT
            END SELECT

 1110 FORMAT('Error 1110: Invalid specifications for point source ',I3,&
         '.',/'Point sources are not supported for ',A,' solids.',/    &
         'Please correct the mfix.dat file.')

! Velocity does not have to be defined (no momentum source). If the
! components are UNDEFINED, zero them out.
            IF(PS_U_s(PSV,M) == UNDEFINED) PS_U_s(PSV,M) = ZERO
            IF(PS_V_s(PSV,M) == UNDEFINED) PS_V_s(PSV,M) = ZERO
            IF(PS_W_s(PSV,M) == UNDEFINED) PS_W_s(PSV,M) = ZERO

! Sum together defiend gas phase species mass fractions.
            SUM = ZERO
            DO N = 1, NMAX(M)
               IF(PS_X_S(PSV,M,N) /= UNDEFINED) THEN
               SUM = SUM + PS_X_G(PSV,N)
               ELSE
                  PS_X_S(PSV,M,N) = ZERO
               ENDIF
            ENDDO 

! Enforce that the species mass fractions must sum to one.
            IF(.NOT.COMPARE(ONE,SUM)) THEN

               IF(SPECIES_EQ(M)) THEN
                  WRITE(ERR_MSG, 1120) PSV,M
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1120 FORMAT('Error 1120: PS_X_s(',I3,',',I2,',:) do NOT sum to ONE ', &
         'and the solids phase',/'species equations are solved. ',     &
         'Please correct the mfix.dat file.')

               ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
                  WRITE(ERR_MSG, 1121) PSV,M
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1121 FORMAT('Error 1121: PS_X_s(',I3,',',I2,',:) do not sum to ONE ', &
         'or ZERO and they',/'are not needed. Please correct the ',    &
         'mfix.dat the mfix.dat file.')

               ELSE
                  PS_X_S(PSV,M,1)  = ONE
                  PS_X_S(PSV,M,2:) = ZERO
               ENDIF

            ENDIF 

! Verify that a temperature is provided.
            IF(ENERGY_EQ)THEN
               IF(PS_T_s(PSV,M) == UNDEFINED) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('PS_T_s',PSV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Verify that a given temperature is physical.
               ELSEIF(PS_T_s(PSV,M) <= ZERO) THEN
                  WRITE(ERR_MSG,1001) PSV, trim(iVar('PS_T_s',PSV,M)), &
                     trim(iVal(PS_T_s(PSV,M)))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_PS_SOLIDS_PHASES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_PS_OVERFLOW                                       !
!  Author: J. Musser                                  Date: 10-JUN-13  !
!                                                                      !
!  Purpose: Check point source specifications.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_PS_OVERFLOW(PSV)

      use param
      use run
      use ps

      use physprop, only: NMAX
      use error_manager

      implicit none

      INTEGER, INTENT(in) :: PSV
      INTEGER :: M, N
      DOUBLE PRECISION :: SUM

      LOGICAL, EXTERNAL :: COMPARE


      CALL INIT_ERR_MSG("CHECK_PS_OVERFLOW")


! Mass flow is undefined --> Velocity must also be undefined.
!```````````````````````````````````````````````````````````````````````
      IF(PS_MASSFLOW_G(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_MASSFLOW_G',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_U_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_U_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_V_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_V_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_W_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_W_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(PS_T_g(PSV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1010) trim(iVar('PS_T_g',PSV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      DO N = 1, DIM_N_G
         IF(PS_X_G(PSV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_X_G',PSV,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO 

      DO M=1, DIM_M
         IF(PS_MASSFLOW_S(PSV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_MASSFLOW_S',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(PS_U_s(PSV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_U_s',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(PS_V_s(PSV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_V_s',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(PS_W_s(PSV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_W_s',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(PS_T_s(PSV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1010) trim(iVar('PS_T_s',PSV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         DO N = 1, DIM_N_S
            IF(PS_X_S(PSV,M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1010) trim(iVar('PS_X_S',PSV,M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO 
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined PS region.') 

      END SUBROUTINE CHECK_PS_OVERFLOW
