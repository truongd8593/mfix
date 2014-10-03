!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  function: isNaN                                                     !
!  Purpose: check whether argument is NAN                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      LOGICAL FUNCTION isNan(x)

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      double precision x
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      CHARACTER(LEN=80) :: notnumber
!-----------------------------------------------

      isNan = .False.
      WRITE(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contain a letter "N"
! "n" or symbol "?", in which case it is a NaN (Not a Number)

      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        isNan = .TRUE.
         RETURN
      ENDIF

      RETURN
      END FUNCTION isNan


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  function: MAX_VEL_INLET                                             C
!  Purpose: Find maximum velocity at inlets.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 26-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION MAX_VEL_INLET()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE bc
      USE fldvar
      USE geometry
      USE physprop
      USE indices
      USE constant
      USE run
      USE compar
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L, I, J, K, IJK, IJK2, M
!-----------------------------------------------

! initializing
      MAX_VEL_INLET = ZERO

      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
            IF (BC_TYPE(L) == 'MASS_INFLOW' .OR. BC_TYPE(L) == 'P_INFLOW') THEN

               DO K = BC_K_B(L), BC_K_T(L)
                  DO J = BC_J_S(L), BC_J_N(L)
                     DO I = BC_I_W(L), BC_I_E(L)
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)

                        SELECT CASE (BC_PLANE(L))
                        CASE ('S')
                           IJK2 = JM_OF(IJK)
                           IF( ABS(V_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK2))
                        CASE ('N')
                           IF( ABS(V_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK))
                        CASE ('W')
                           IJK2 = IM_OF(IJK)
                           IF( ABS(U_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK2))
                        CASE ('E')
                           IF( ABS(U_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK))
                        CASE ('B')
                           IJK2 = KM_OF(IJK)
                           IF( ABS(W_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK2))
                        CASE ('T')
                           IF( ABS(W_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK))
                        END SELECT

                       IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
                          SELECT CASE (BC_PLANE(L))
                           CASE ('S')
                              IJK2 = JM_OF(IJK)
                              DO M = 1, MMAX
                                 IF( ABS(V_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK2, M))
                              ENDDO
                           CASE ('N')
                              DO M = 1, MMAX
                                IF( ABS(V_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK, M))
                              ENDDO
                           CASE ('W')
                              IJK2 = IM_OF(IJK)
                              DO M = 1, MMAX
                                IF( ABS(U_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK2, M))
                              ENDDO
                           CASE ('E')
                              DO M = 1, MMAX
                                IF( ABS(U_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK, M))
                              ENDDO
                           CASE ('B')
                              IJK2 = KM_OF(IJK)
                              DO M = 1, MMAX
                                IF( ABS(W_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK2, M))
                              ENDDO
                           CASE ('T')
                              DO M = 1, MMAX
                                IF( ABS(W_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK, M))
                              ENDDO
                           END SELECT
                        ENDIF   ! end if (.not.des_continuum_coupled .or. des_continuum_hybrid)

                     ENDDO
                  ENDDO
               ENDDO

           ENDIF
         ENDIF
      ENDDO

      RETURN
      END FUNCTION MAX_VEL_INLET


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: CHECK_VEL_BOUND()                                         C
!  Purpose: Check velocities upper bound to be less than speed of      C
!           sound                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 25-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      LOGICAL FUNCTION CHECK_VEL_BOUND ()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE run
      USE toleranc
      USE compar
      USE mpi_utility
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: M
! Indices
      INTEGER :: IJK
      LOGICAL :: ALL_IS_ERROR
!-----------------------------------------------

!!$omp   parallel do private(IJK)
! initializing
      CHECK_VEL_BOUND = .FALSE.
      ALL_IS_ERROR    = .FALSE.

LOOP_FLUID : DO IJK = IJKSTART3, IJKEND3

         IF (FLUID_AT(IJK)) THEN
            IF(ABS(U_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(V_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(W_G(IJK)) > MAX_INLET_VEL) THEN
               CHECK_VEL_BOUND = .TRUE.
               WRITE(*,1000) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                             EP_g(IJK), U_G(IJK), V_G(IJK), W_G(IJK)
               EXIT LOOP_FLUID
            ENDIF

            IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
               DO M = 1, MMAX
                 IF(ABS(U_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(V_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(W_S(IJK,M)) > MAX_INLET_VEL) THEN
                   CHECK_VEL_BOUND = .TRUE.
                   WRITE(*,1010) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), M, &
                                 EP_s(IJK, M), U_S(IJK,M), V_S(IJK,M), W_S(IJK,M)
                   EXIT LOOP_FLUID
                 ENDIF
               ENDDO
            ENDIF   ! end if(.not.des_continuum_coupled or des_continuum_hybrid)
         ENDIF

      ENDDO LOOP_FLUID

      CALL GLOBAL_ALL_OR(CHECK_VEL_BOUND, ALL_IS_ERROR)
      IF(ALL_IS_ERROR) CHECK_VEL_BOUND = .TRUE.

      RETURN
 1000 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor MAX_INLET_VEL_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Ug = ', G12.5, 'Vg = ', G12.5, 'Wg = ', G12.5)
 1010 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5,/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4,' M = ',I4, /&
            '  ','Eps = ', G12.5,'Us = ', G12.5, 'Vs = ', G12.5, 'Ws = ', G12.5)

      END FUNCTION CHECK_VEL_BOUND

