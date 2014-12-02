!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIF_PHI_DES(M, A_M, B_M, IER)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE toleranc
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE fun_avg
      USE functions
      USE cutcell

      IMPLICIT NONE

! Phase index
      INTEGER, INTENT(IN) :: M

!  Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

      INTEGER :: IER

! Fluid Cell indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IM, IJKW, IPJK, IJKE
      INTEGER :: IJMK, JM, IJKS, IJPK, IJKN
      INTEGER :: IJKM, KM, IJKB, IJKP, IJKT
!
! Difusion parameter
      DOUBLE PRECISION :: D_f
!
!-----------------------------------------------
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
      DO IJK = ijkstart3, ijkend3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)

         IJKE = EAST_OF(IJK)
         IJKN = NORTH_OF(IJK)

! East face (i+1/2, j, k)
         D_F = ODX_E(I)*AYZ(IJK)
         IF(CUT_TREATMENT_AT(IJK)) THEN
            IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IPJK))) THEN
               D_F = ODX_E(I)*DY(J)*DZ(K)
            ENDIF
         ENDIF

         A_M(IJK,E,M) = D_F
         A_M(IPJK,W,M) = D_F

! West face (i-1/2, j, k)
         IMJK = IM_OF(IJK)
         IF (.NOT.FLUID_AT(IMJK)) THEN
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            D_F = ODX_E(IM)*AYZ(IMJK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IMJK))) THEN
                  D_F = ODX_E(IM)*DY(J)*DZ(K)
               ENDIF
            ENDIF
            A_M(IJK,W,M) = D_F
         ENDIF


! North face (i, j+1/2, k)
         D_F = ODY_N(J)*AXZ(IJK)
         IF(CUT_TREATMENT_AT(IJK)) THEN
            IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJPK))) THEN
               D_F = ODY_N(J)*DX(I)*DZ(K)
            ENDIF
         ENDIF
         A_M(IJK,N,M) = D_F
         A_M(IJPK,S,M) = D_F

! South face (i, j-1/2, k)
         IJMK = JM_OF(IJK)
         IF (.NOT.FLUID_AT(IJMK)) THEN
            JM = JM1(J)
            IJKS = SOUTH_OF(IJK)
            D_F = ODY_N(JM)*AXZ(IJMK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJMK))) THEN
                  D_F = ODY_N(JM)*DX(I)*DZ(K)
               ENDIF
            ENDIF
            A_M(IJK,S,M) = D_F
         ENDIF


! Top face (i, j, k+1/2)
         IF (DO_K) THEN
            IJKP = KP_OF(IJK)
            IJKT = TOP_OF(IJK)
            D_F = OX(I)*ODZ_T(K)*AXY(IJK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJKP))) THEN
                  D_F = OX(I)*ODZ_T(K)*DX(I)*DY(J)
               ENDIF
            ENDIF
            A_M(IJK,T,M) = D_F
            A_M(IJKP,B,M) = D_F


! Bottom face (i, j, k-1/2)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLUID_AT(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               D_F = OX(I)*ODZ_T(KM)*AXY(IJKM)
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJKM))) THEN
                     D_F = OX(I)*ODZ_T(KM)*DX(I)*DY(J)
                  ENDIF
               ENDIF
               A_M(IJK,B,M) = D_F
            ENDIF
         ENDIF
      END DO

      CALL DIF_PHI_IS_DES (A_M, B_M, M, IER)


      RETURN
      END SUBROUTINE DIF_PHI_DES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DIF_PHI_IS_DES(A_m, B_m, M, IER)                       !
!  Author: J.Musser                                   Date: 24-NOV-14  !
!                                                                      !
!  Purpose: Remove diffusive fluxes across internal surfaces.          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIF_PHI_IS_DES(A_M, B_M, M, IER)

      USE param
      USE param1
      USE parallel
      USE matrix
      USE toleranc
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE output
      USE is
      USE fun_avg
      USE function3
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Internal surface
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                       IJKE, IJKN, IJKT, IPJK, IJPK, IJKP
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!-----------------------------------------------
!
! Make user defined internal surfaces non-conducting
!
      DO L = 1, DIMENSION_IS
         IF (IS_DEFINED(L)) THEN
            I1 = IS_I_W(L)
            I2 = IS_I_E(L)
            J1 = IS_J_S(L)
            J2 = IS_J_N(L)
            K1 = IS_K_B(L)
            K2 = IS_K_T(L)

! Limit I1, I2 and all to local processor first ghost layer

            IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

! End of limiting to the first ghost cells of the processor....
            DO K = K1, K2
            DO J = J1, J2
            DO I = I1, I2

               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)

               SELECT CASE (TRIM(IS_PLANE(L)))

               CASE ('E')
                  IJKE = EAST_OF(IJK)
                  IPJK = IP_OF(IJK)

                  D_F = ODX_E(I)*AYZ(IJK)

                  A_M(IJK,E,M) = A_M(IJK,E,M) - D_F
                  A_M(IPJK,W,M) = A_M(IPJK,W,M) - D_F

               CASE ('N')
                  IJKN = NORTH_OF(IJK)
                  IJPK = JP_OF(IJK)

                  D_F = ODY_N(J)*AXZ(IJK)

                  A_M(IJK,N,M) = A_M(IJK,N,M) - D_F
                  A_M(IJPK,S,M) = A_M(IJPK,S,M) - D_F

               CASE ('T')
                  IF (DO_K) THEN
                     IJKT = TOP_OF(IJK)
                     IJKP = KP_OF(IJK)

                     D_F = OX(I)*ODZ_T(K)*AXY(IJK)

                     A_M(IJK,T,M) = A_M(IJK,T,M) - D_F
                     A_M(IJKP,B,M) = A_M(IJKP,B,M) - D_F

                  ENDIF
               CASE DEFAULT
               END SELECT
            END DO
            END DO
            END DO
         ENDIF
      END DO

      RETURN
      END SUBROUTINE DIF_PHI_IS_DES

