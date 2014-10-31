!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG                                                !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_U(A_M, B_M, VXF_GS, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use discretelement, only: DES_INTERP_ON
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Flag: TFM and DEM solids exist.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume averaged solids velocity.
      use discretelement, only: DES_U_s
! Number of discrete solids phases
      use discretelement, only: DES_MMAX
! Gas source for coupling gas/solids flows.
      use discretelement, only: VxF_GDS

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays.
      use param, only: DIMENSION_3
! Number of solids phases.
      use param, only: DIMENSION_M

      use functions
      use fun_avg

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Volume times drag coefficient at cell face
      DOUBLE PRECISION, INTENT(IN) :: VXF_GS(DIMENSION_3, DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER


! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IJMK, IJKM, IJMKM
! Solids phase index
      INTEGER :: M
! Face center values of u_sm (i+1/2)
      DOUBLE PRECISION :: USFCM
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION :: tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


! Average the interpoalted drag force from the cell corners to the cell face.
      IF(DES_INTERP_ON)THEN

         AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

!$omp parallel do schedule(guided,50) default(none)                    &
!$omp shared(IJKSTART3, IJKEND3, I_OF, J_OF,   &
!$omp    K_OF, DO_K, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_U)    &
!$omp private(IJK, I, J, K, IJMK, IJKM, IJMKM, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3

            IF(.NOT.FLUID_AT(IJK)) CYCLE

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IJMK = FUNIJK_MAP_C(I, J-1, K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IJMK))
            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,1) + DRAG_BM(IJMK,1))

            IF(DO_K) THEN
               IJKM = FUNIJK_MAP_C(I, J, K-1)
               IJMKM = FUNIJK_MAP_C(I, J-1, K-1)
               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(IJKM) + DRAG_AM(IJMKM))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(IJKM,1) + DRAG_BM(IJMKM,1))
            ENDIF

            A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_U(IJK)
            B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_U(IJK)

         ENDDO
!$omp end parallel do


      ELSE

!$omp parallel do default(none) schedule(guided, 50)                   &
!$omp shared(DES_MMAX, IJKSTART3, IJKEND3, I_OF, DES_CONTINUUM_HYBRID, &
!$omp    DES_U_S, VxF_GDS, VxF_GS, A_M, B_M)                           &
!$omp private(M, IJK, I, USFCM, tmp_A, tmp_B)
         DO M = 1, DES_MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(.NOT.FLUID_AT(IJK)) CYCLE

               I = I_OF(IJK)
               USFCM = AVG_X(DES_U_S(IJK,M),DES_U_S(EAST_OF(IJK),M),I)

               IF (DES_CONTINUUM_HYBRID) THEN
                  tmp_A =  - VXF_GDS(IJK,M)
                  tmp_B =  - VXF_GDS(IJK,M)*USFCM
               ELSE
                  tmp_A =  - VXF_GS(IJK,M)
                  tmp_B =  - VXF_GS(IJK,M)*USFCM
               ENDIF

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
               B_M(IJK,0) = B_M(IJK,0) + tmp_B

            ENDDO
         ENDDO
!$omp end parallel do
      ENDIF

      END SUBROUTINE GAS_DRAG_U


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GAS_DRAG                                                C
!  Purpose: Account for the equal and opposite drag force on the gas   C
!           phase due to particles by introducing the drag as a        C
!           source term.  Face centered.                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GAS_DRAG_V(A_M, B_M, VXF_GS, IER)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use discretelement, only: DES_INTERP_ON
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Flag: TFM and DEM solids exist.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume averaged solids velocity.
      use discretelement, only: DES_V_s
! Number of discrete solids phases
      use discretelement, only: DES_MMAX
! Gas source for coupling gas/solids flows.
      use discretelement, only: VxF_GDS

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays.
      use param, only: DIMENSION_3
! Number of solids phases.
      use param, only: DIMENSION_M

      use functions
      use fun_avg


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Volume times drag coefficient at cell face
      DOUBLE PRECISION, INTENT(IN) :: VXF_GS(DIMENSION_3, DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJKM, IMJKM
! Solids phase index
      INTEGER :: M
! Face center values of v_sm (j+1/2)
      DOUBLE PRECISION :: VSFCM
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!


! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_ON) THEN

         AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

!$omp parallel do schedule (guided,50) default(none)                   &
!$omp shared(IJKSTART3, IJKEND3, I_OF, J_OF,   &
!$omp    K_OF, DO_K, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_V)    &
!$omp private(IJK, I, J, K, IMJK, IJKM, IMJKM, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.FLUID_AT(IJK)) CYCLE

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IMJK = FUNIJK_MAP_C(I-1,J,K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK))
            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,2) + DRAG_BM(IMJK,2))

            IF(DO_K) THEN

               IJKM = FUNIJK(I,J,K-1)
               IMJKM = FUNIJK(I-1,J,K-1)

               tmp_A = tmp_A - AVG_FACTOR*                             &
                  (DRAG_AM(IJKM) + DRAG_AM(IMJKM))
               tmp_B = tmp_B - AVG_FACTOR*                             &
                  (DRAG_BM(IJKM,2) + DRAG_BM(IMJKM,2))
            ENDIF

            A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_V(IJK)
            B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_V(IJK)

         ENDDO
!$omp end parallel do


      ELSE

!$omp parallel do default(none) schedule(guided, 50)                   &
!$omp shared(DES_MMAX, IJKSTART3, IJKEND3, J_OF, DES_CONTINUUM_HYBRID, &
!$omp    DES_V_S, VxF_GDS, VxF_GS, A_M, B_M)                           &
!$omp private(M, IJK, J, VSFCM, tmp_A, tmp_B)
         DO M = 1, DES_MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(.NOT.FLUID_AT(IJK)) CYCLE

               J = J_OF(IJK)
               VSFCM = AVG_Y(DES_V_S(IJK,M),DES_V_S(NORTH_OF(IJK),M),J)

               IF(DES_CONTINUUM_HYBRID) THEN
                   tmp_A =  - VXF_GDS(IJK,M)
                   tmp_B =  - VXF_GDS(IJK,M)*VSFCM
               ELSE
                  tmp_A =  - VXF_GS(IJK,M)
                  tmp_B =  - VXF_GS(IJK,M)*VSFCM
               ENDIF

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
               B_M(IJK,0) = B_M(IJK,0) + tmp_B

            ENDDO
         ENDDO
!$omp end parallel do

      ENDIF

      END SUBROUTINE GAS_DRAG_V


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GAS_DRAG                                                C
!  Purpose: Account for the equal and opposite drag force on the gas   C
!           phase due to particles by introducing the drag as a        C
!           source term.  Face centered.                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GAS_DRAG_W(A_M, B_M, VXF_GS, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use discretelement, only: DES_INTERP_ON
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Flag: TFM and DEM solids exist.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume averaged solids velocity.
      use discretelement, only: DES_W_s
! Number of discrete solids phases
      use discretelement, only: DES_MMAX
! Gas source for coupling gas/solids flows.
      use discretelement, only: VxF_GDS

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays.
      use param, only: DIMENSION_3
! Number of solids phases.
      use param, only: DIMENSION_M

      use functions
      use fun_avg


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Volume times drag coefficient at cell face
      DOUBLE PRECISION, INTENT(IN) :: VXF_GS(DIMENSION_3, DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IMJMK
! Solids phase index
      INTEGER :: M
! Face center values of w_sm (k+1/2)
      DOUBLE PRECISION :: WSFCM
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
! (=0.25 in 3D and =0.5 in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!


! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN


      IF(DES_INTERP_ON) THEN

         AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

!$omp parallel do schedule (guided,50) default(none)                   &
!$omp shared(IJKSTART3, IJKEND3, I_OF, J_OF,   &
!$omp    K_OF, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_W)          &
!$omp private(IJK, I, J, K, IMJK, IJMK, IMJMK, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.FLUID_AT(IJK)) CYCLE

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IMJK = FUNIJK(I-1,J,K)
            IJMK = FUNIJK(I,J-1,K)
            IMJMK = FUNIJK(I-1,J-1,K)

            tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK) +        &
               DRAG_AM(IJMK) + DRAG_AM(IMJMK))

            tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,3) + DRAG_BM(IMJK,3) +    &
               DRAG_BM(IJMK,3) + DRAG_BM(IMJMK,3))

            A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_W(IJK)
            B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_W(IJK)

         ENDDO
!$omp end parallel do

      ELSE

!$omp parallel do default(none) schedule(guided, 50)                   &
!$omp shared(DES_MMAX, IJKSTART3, IJKEND3, K_OF, DES_CONTINUUM_HYBRID, &
!$omp    DES_W_S, VxF_GDS, VxF_GS, A_M, B_M)                           &
!$omp private(M, IJK, K, WSFCM, tmp_A, tmp_B)
         DO M = 1, DES_MMAX
            DO IJK = IJKSTART3, IJKEND3

               IF(FLUID_AT(IJK)) CYCLE

               K = K_OF(IJK)
               WSFCM = AVG_Z(DES_W_S(IJK,M),DES_W_S(TOP_OF(IJK),M),K)

               IF (DES_CONTINUUM_HYBRID) THEN
                  tmp_A =  - VXF_GDS(IJK,M)
                  tmp_B =  - VXF_GDS(IJK,M)*WSFCM
               ELSE
                  tmp_A = - VXF_GS(IJK,M)
                  tmp_B = - VXF_GS(IJK,M)*WSFCM
               ENDIF

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
               B_M(IJK,0) = B_M(IJK,0) + tmp_B

            ENDDO
         ENDDO
!$omp end parallel do

      ENDIF

      RETURN
      END SUBROUTINE GAS_DRAG_W
