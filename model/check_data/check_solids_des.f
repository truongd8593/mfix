!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_DES_SOLIDS                                        !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Calculate clusters during a DEM simulation
      use discretelement, only: DES_CALC_CLUSTER
      use discretelement, only: CLUSTER_LENGTH_CUTOFF
      use discretelement, only: FACTOR_RLM
      use discretelement, only: MAX_RADIUS

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: CLUSTER_RLM
! NONE

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES")


! Particle-particle collision parameters.
      CALL CHECK_SOLIDS_DES_COLLISION
! DES cohesion model parameters.
      CALL CHECK_SOLIDS_DES_COHESION
! Particle-particle conduction model parameters.
      CALL CHECK_SOLIDS_DES_ENERGY


      IF(DES_CALC_CLUSTER) THEN

! Verify that a cutoff distance was provided.
         IF(CLUSTER_LENGTH_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'CLUSTER_LENGTH_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! Verify that a cutoff distance lands with the neighbor search region.
         CLUSTER_RLM = 1.d0 + CLUSTER_LENGTH_CUTOFF/(2.d0*MAX_RADIUS)
         IF(FACTOR_RLM < CLUSTER_RLM) THEN
            WRITE(ERR_MSG,2000)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2000 FORMAT('Error 2000: CLUSTER_LENGTH_CUTOFF exceeds of neighbor',  &
         ' search',/'distance. Increase FACTOR_RLM.')
      ENDIF



      CALL FINL_ERR_MSG

      RETURN  

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DES




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_DES_ENERGY                                  !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES_ENERGY

      use des_thermo, only: DES_COND_EQ
      use des_thermo, only: DES_MIN_COND_DIST
      use run, only: UNITS

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED

      use error_manager

      IMPLICIT NONE


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES_ENERGY")

! Conduction Equations:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      IF(DES_COND_EQ) THEN

! Set the default value for the minimum distance separating particles'
! surfaces.
         IF(DES_MIN_COND_DIST == UNDEFINED)THEN
            DES_MIN_COND_DIST = 1.0D-04 ! cm
            IF (UNITS == 'SI') DES_MIN_COND_DIST = &
               DES_MIN_COND_DIST/100.0  ! m
         ENDIF
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_DES_ENERGY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_SOLIDS_DES_COHESION                              !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check/set parameters for DES cohesion models.              !
!                                                                      !
!  Comments: Original code moved from CHECK_DES_DATA                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES_COHESION

! Global Variables:
!---------------------------------------------------------------------//

! Runtime Flag: Invoke a cohesion model for DES simulation
      use discretelement, only: USE_COHESION
! Largest discrete particle diameter.
      use discretelement, only: MAX_RADIUS

! Runtime Flag: Invoke Square Well 
      use discretelement, only: SQUARE_WELL
! Runtime Flag: Invoke Van der Waals model.
      use discretelement, only: VAN_DER_WAALS
! User specified parameter for increase neighbor search area.
      use discretelement, only: FACTOR_RLM

! User specified parameters for Van der Waals model.
      use discretelement, only: VDW_INNER_CUTOFF
      use discretelement, only: VDW_OUTER_CUTOFF
      use discretelement, only: HAMAKER_CONSTANT
      use discretelement, only: WALL_VDW_INNER_CUTOFF
      use discretelement, only: WALL_VDW_OUTER_CUTOFF
      use discretelement, only: WALL_HAMAKER_CONSTANT

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED
      use param1, only: ZERO


      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Neighborhood size for Van der Waals force.
      DOUBLE PRECISION :: VDW_NEIGHBORHOOD


! Override the following settings if cohesion not used.
      IF(.NOT.USE_COHESION) THEN
!No more square well in the code 
         SQUARE_WELL = .FALSE. 
         VAN_DER_WAALS = .FALSE.
         WALL_VDW_OUTER_CUTOFF = ZERO
         RETURN
      ENDIF


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES_COHESION")


! Verify that only one cohesion model is specified.
      IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Verify that at a cohesion model is specified.
      ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1001)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Van der Waals model checks.
      IF (VAN_DER_WAALS) THEN

         IF (VDW_INNER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_INNER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(VDW_OUTER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_OUTER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'HAMAKER_CONSTANT'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_INNER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (WALL_VDW_OUTER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_OUTER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'WALL_HAMAKER_CONSTANT'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         VDW_NEIGHBORHOOD = 1.0d0 + (VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS))
         IF (FACTOR_RLM < VDW_NEIGHBORHOOD) THEN
            WRITE(ERR_MSG,1202)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

      ENDIF                


      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT(/1X,70('*')/' From: CHECK_SOLIDS_DES_COHESION',/' Error 1000:',  &
         ' Cannot use SQUARE_WELL and VAN_DER_WAALS cohesion',/        &
         ' models simultaneously.',/1X,70('*')/)            

 1001 FORMAT(/1X,70('*')/' From: CHECK_SOLIDS_DES_COHESION',/' Error 1001:',  &
         ' A cohesion model was not selected. Specify one of the',/    &
         ' available models in the mfix.dat file.',/1X,70('*')/)


!<------------------- Van der Waals model messages. ----------------->!


 1201 FORMAT(/1X,70('*')/' From: CHECK_SOLIDS_DES_COHESION',/' Error 1201:',  &
         ' Missing input data for Van der Waals cohesion model.',/     &
         ' Input parameter ',A,' is UNDEFINED.',/1X,70('*')/)

 1202 FORMAT(/1X,70('*')/' From: CHECK_SOLIDS_DES_COHESION',/' Error 1202: ', &
         ' VDW_OUTER_CUTOFF outside of the neighbor search distance.',/&
         ' Increase FACTOR_RLM to increase the search distance.',/     &
         1X,70('*')/)

      END SUBROUTINE CHECK_SOLIDS_DES_COHESION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DES_COLLISION                              !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DES_COLLISION


! Global Variables:
!---------------------------------------------------------------------//

! User specified collision model
      USE discretelement, only: DES_COLL_MODEL
      USE discretelement, only: DES_COLL_MODEL_ENUM
      USE discretelement, only: LSD
      USE discretelement, only: HERTZIAN
! Number of discrete solids phases
      USE discretelement, only: DES_MMAX
! Particle and wall friction coeff.
      USE discretelement, only: MEW, MEW_W
! Particle and wall normal spring constants
      USE discretelement, only: KN, KN_W
! Particle and wall tangential spring factor := KN/KT
      USE discretelement, only: KT_FAC, KT_W_FAC
! Coefficients of restitution: Normal and Tangential
      USE discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT
      USE discretelement, only: DES_ET_INPUT, DES_ET_WALL_INPUT
! Tangential damping factors := ET/EN
      USE discretelement, only: DES_ETAT_FAC, DES_ETAT_W_FAC
! Particle and wall Young's modulus
      USE discretelement, only: E_YOUNG, EW_YOUNG
! Particle and wall Pooisson ratio
      USE discretelement, only: V_POISSON, VW_POISSON
! Parameter constatns.
      USE param1, only: ONE, ZERO, UNDEFINED

!      USE mpi_utility
      use error_manager

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M
! Number of phase interactions
      INTEGER :: MxM
! Flag to warn user.
      LOGICAL :: FLAG_WARN
! Message for formated output.
      CHARACTER(len=64) :: MSG




! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DES_COLLISION")



! check particle-particle normal restitution coefficient
      MxM = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      DO M = 1, MxM
         IF(DES_EN_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(DES_EN_INPUT(M) > ONE .OR.                             &
            DES_EN_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',M)),          &
               trim(iVal(DES_EN_INPUT(M)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Check particle-wall normal restitution coefficient.
      DO M = 1, DES_MMAX
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(DES_EN_WALL_INPUT(M) > ONE .OR.                        &
            DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)),     &
               trim(iVal(DES_EN_WALL_INPUT(M)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO


! Check coefficient friction
      IF(MEW == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF (MEW < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW', trim(iVal(MEW))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


      IF(MEW_W == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW_W'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(MEW_w < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW_W', trim(iVal(MEW_W))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF




! Check collision model specific model parameters.
      SELECT CASE (trim(DES_COLL_MODEL))

!**********************************************************************!
!*                                                                    *!
!*                     linear spring-dashpot model                    *!
!*                                                                    *!
!**********************************************************************!
      CASE('LSD'); DES_COLL_MODEL_ENUM = LSD

         IF(KN == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) 'KN'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(KN_W == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) 'KN_W'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

! Check for tangential spring constant factors
         IF(KT_FAC == UNDEFINED) THEN
            WRITE (ERR_MSG, 2100) 'KT_FAC'
            CALL FLUSH_ERR_MSG()
         ELSEIF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
            WRITE(ERR_MSG,1001) 'KT_FAC', trim(iVal(KT_FAC))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(KT_W_FAC == UNDEFINED) THEN
            WRITE (ERR_MSG, 2100) 'KT_W_FAC'
            CALL FLUSH_ERR_MSG()
         ELSEIF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
            WRITE(ERR_MSG,1001) 'KT_W_FAC', trim(iVal(KT_W_FAC))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2100 FORMAT('Warning 1102: Tangential spring factor ',A,' not ',      &
         'specified in mfix.dat.',/'It will be defined in cfassign.f ',&
         'as 2/7.')


! Check for tangential damping factor
         IF(DES_ETAT_FAC == UNDEFINED) THEN
            WRITE (ERR_MSG, 2101) 'DES_ETAT_FAC'
            CALL FLUSH_ERR_MSG
         ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
            WRITE(ERR_MSG,1001) 'DES_ETAT_FAC', iVal(DES_ETAT_FAC)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(DES_ETAT_W_FAC == UNDEFINED) THEN
            WRITE (ERR_MSG, 2101) 'DES_ETAT_W_FAC'
            CALL FLUSH_ERR_MSG
         ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
            WRITE(ERR_MSG,1001) 'DES_ETAT_W_FAC', iVal(DES_ETAT_W_FAC)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2101 FORMAT('Warning 2101: Tangential damping factor ',A,' not ',     &
         'specified.',/' The default of 1/2 will be set in cfassign.f.')


! if following are assigned warn user they are discarded
         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
            IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN) THEN
            WRITE(ERR_MSG,2102) 'DES_ET_INPUT'
            CALL FLUSH_ERR_MSG
         ENDIF

         FLAG_WARN = .FALSE.
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
         IF (FLAG_WARN)THEN
            WRITE(ERR_MSG,2102) 'DES_ET_WALL_INPUT'
            CALL FLUSH_ERR_MSG
         ENDIF

 2102 FORMAT('Warning 2102: ',A,' values are not used ',/' with the',  &
         ' linear spring-dashpot collision model.')


!**********************************************************************!
!*                                                                    *!
!*                           Hertzian model                           *!
!*                                                                    *!
!**********************************************************************!
      CASE('HERTZIAN'); DES_COLL_MODEL_ENUM = HERTZIAN

! check young's modulus and poisson ratio
         IF(EW_YOUNG == UNDEFINED ) THEN
            MSG='Wall value for Youngs modulus'
            WRITE(ERR_MSG,1002) 'EW_YOUNG', MSG
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(VW_POISSON == UNDEFINED) THEN
            MSG='Wall value for Poissons ratio'
            WRITE(ERR_MSG,1002) 'VW_POISSON', MSG
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (VW_POISSON > 0.5d0 .OR. VW_POISSON <= -ONE) THEN
            WRITE(ERR_MSG,1001) 'VW_POISSON',iVal(VW_POISSON)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


         DO M = 1, DES_MMAX
            IF(E_YOUNG(M) == UNDEFINED) THEN
               MSG=''; WRITE(MSG,"('Phase ',I2,' Youngs modulus')") M
               WRITE(ERR_MSG,1002) 'E_YOUNG', MSG
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            IF(V_POISSON(M) == UNDEFINED) THEN
               MSG=''; WRITE(MSG,"('Phase ',I2,' Poissons ratio')") M
               WRITE(ERR_MSG,1002) 'V_POISSON', MSG
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(V_POISSON(M) > 0.5d0 .OR.                           &
               V_POISSON(M) <= -ONE) THEN
               WRITE(ERR_MSG,1001) trim(iVar('V_POISSON',M)),          &
                  iVal(V_POISSON(M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO


! check particle-particle tangential restitution coefficient
         MxM = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
         DO M = 1, MxM
            IF(DES_ET_INPUT(M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_ET_INPUT',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

         DO M = 1, MxM
            IF(DES_ET_INPUT(M) > ONE .OR. DES_ET_INPUT(M) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_ET_INPUT',M)),       &
                  iVal(DES_ET_INPUT(M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

! check particle-wall tangential restitution coefficient
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_ET_WALL_INPUT',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
         DO M = 1, DES_MMAX
            IF(DES_ET_WALL_INPUT(M) > ONE .OR.                         &
               DES_ET_WALL_INPUT(M) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_ET_WALL_INPUT',M)),  &
                  trim(iVal(DES_ET_WALL_INPUT(M)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

! If following are assigned warn user they are discarded.
          IF(KN .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'KN'
             CALL FLUSH_ERR_MSG
          ENDIF
          IF(KN_W .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'KN_W'
             CALL FLUSH_ERR_MSG
          ENDIF
          IF(KT_FAC .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'KT_FAC'
             CALL FLUSH_ERR_MSG
          ENDIF
          IF(KT_W_FAC .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'KT_W_FAC'
             CALL FLUSH_ERR_MSG
          ENDIF
          IF(DES_ETAT_FAC .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'DES_ETAT_FAC'
             CALL FLUSH_ERR_MSG
          ENDIF
          IF(DES_ETAT_W_FAC .NE. UNDEFINED) THEN
             WRITE(ERR_MSG, 2200) 'DES_ETAT_W_FAC'
             CALL FLUSH_ERR_MSG
          ENDIF

 2200 FORMAT('Warning 2200: ',A,' values are not used ',/' with the',  &
         ' linear spring-dashpot collision model.')

! Unknown collision model.
      CASE DEFAULT

         WRITE(ERR_MSG,2000) TRIM(DES_COLL_MODEL)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2000 FORMAT('Error 2000: Invalid particle-particle collision model:',&
         A,/'Please correct the mfix.dat file.')

      END SELECT



      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Required input not specified: ',A,/          &
         'Description:',A,/'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DES_COLLISION
