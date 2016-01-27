!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS1

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param, only: dimension_3
      USE param1, only: zero, one
      USE fldvar, only: u_s, v_s, w_s, rop_s, ro_s
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE mfix_pic
      USE functions

      use run, only: SOLVE_ROs
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_SIZE
      use particle_filter, only: DES_INTERP_ON
      use physprop, only: mmax
      use sendrecvnode, only: DES_COLLECT_gDATA

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: particles, filter cells, phases
      INTEGER NP, LC, M, MMAX_TOT
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: sVOL(DIMENSION_3,DES_MMAX)
      DOUBLE PRECISION :: sMASS(DIMENSION_3,DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight, and times filter weight
      DOUBLE PRECISION :: VOL_WT, VOLxWEIGHT
! Loop bound for filter
      LOGICAL :: CALC_sMASS
!-----------------------------------------------

      CALC_sMASS = any(SOLVE_ROs)
      MMAX_TOT = DES_MMAX+MMAX

! Initialize arrays

      sVOL(:,:) = ZERO
      if(CALC_sMASS) sMASS(:,:) = ZERO

      IF(MPPIC) THEN
         U_S(:,MMAX+1:MMAX_TOT) = ZERO
         V_S(:,MMAX+1:MMAX_TOT) = ZERO
         IF(DO_K) W_S(:,MMAX+1:MMAX_TOT) = ZERO
      ENDIF

! Calculate the gas phase forces acting on each particle.
!$omp parallel default(none) &
!$omp private(NP, VOL_WT, M, LC, IJK, VOLXWEIGHT) &
!$omp shared(MAX_PIP, PVOL, DES_STAT_WT, PIJK, FILTER_SIZE, MPPIC, &
!$omp       FILTER_WEIGHT, sVOL, U_S, V_S, W_S, DO_K, CALC_sMASS,  &
!$omp       FILTER_CELL, DES_VEL_NEW, sMASS, pMASS, DES_INTERP_ON)
!$omp do
      do NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP)) CYCLE

         VOL_WT = PVOL(NP)
         IF(MPPIC) VOL_WT = VOL_WT*DES_STAT_WT(NP)
! Particle phase for data binning.
         M = PIJK(NP,5)

! Interpolated
!---------------------------------------------------------------------//
         IF(DES_INTERP_ON) THEN

            DO LC=1,FILTER_SIZE
               IJK = FILTER_CELL(LC,NP)
! Particle volume times the weight for this cell.
               VOLxWEIGHT = VOL_WT*FILTER_WEIGHT(LC,NP)


!$omp atomic
               sVOL(IJK,M) = sVOL(IJK,M) + VOLxWEIGHT

! Accumulate total solids momentum (by phase)
               IF(MPPIC) THEN
!$omp atomic
                  U_S(IJK,M) = U_S(IJK,M) + &
                     DES_VEL_NEW(NP,1)*VOLxWEIGHT
!$omp atomic
                  V_S(IJK,M) = V_S(IJK,M) + &
                     DES_VEL_NEW(NP,2)*VOLxWEIGHT

                  IF(DO_K) THEN
!$omp atomic
                     W_S(IJK,M) = W_S(IJK,M) + &
                          DES_VEL_NEW(NP,3)*VOLxWEIGHT
                  ENDIF
               ENDIF

! Accumulate total solids mass (by phase)
               IF(CALC_sMASS) THEN
                  IF(MPPIC) THEN
!$omp atomic
                     sMASS(IJK,M) = sMASS(IJK,M) + pMASS(NP) * &
                        FILTER_WEIGHT(LC,NP) * DES_STAT_WT(NP)
                  ELSE
!$omp atomic
                     sMASS(IJK,M) = sMASS(IJK,M) + pMASS(NP) * &
                        FILTER_WEIGHT(LC,NP)
                  ENDIF
               ENDIF
            ENDDO

! Non-interpolated
!---------------------------------------------------------------------//
         ELSE
            IJK = PIJK(NP,4)
!$omp atomic
            sVOL(IJK,M) = sVOL(IJK,M) + VOL_WT

! Accumulate total solids momentum (by phase)
            IF(MPPIC) THEN
!$omp atomic
               U_S(IJK,M) = U_S(IJK,M) + &
                  DES_VEL_NEW(NP,1)*VOL_WT
!$omp atomic
               V_S(IJK,M) = V_S(IJK,M) + &
                  DES_VEL_NEW(NP,2)*VOL_WT

               IF(DO_K) THEN
!$omp atomic
                  W_S(IJK,M) = W_S(IJK,M) + &
                       DES_VEL_NEW(NP,3)*VOL_WT
               ENDIF
            ENDIF

! Accumulate total solids mass (by phase)
            IF(CALC_sMASS) THEN
               IF(MPPIC) THEN
!$omp atomic
                  sMASS(IJK,M)= sMASS(IJK,M) + pMASS(NP)*DES_STAT_WT(NP)
               ELSE
!$omp atomic
                  sMASS(IJK,M)= sMASS(IJK,M) + pMASS(NP)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!$omp end do
!$omp end parallel


! Summ data interpolted into ghost cells into physical cells
!---------------------------------------------------------------------//
      IF(DES_INTERP_ON) THEN
         CALL DES_COLLECT_gDATA(sVOL(:,MMAX+1:MMAX_TOT))
         IF(CALC_sMASS) CALL DES_COLLECT_gDATA(sMASS(:,MMAX+1:MMAX_TOT))
         IF(MPPIC) THEN
            CALL DES_COLLECT_gDATA(U_s(:,MMAX+1:MMAX_TOT))
            CALL DES_COLLECT_gDATA(V_s(:,MMAX+1:MMAX_TOT))
            CALL DES_COLLECT_gDATA(W_s(:,MMAX+1:MMAX_TOT))
         ENDIF
      ENDIF


! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!---------------------------------------------------------------------//
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)                &
!$omp private(IJK,M,OoSOLVOL)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = MMAX+1, DES_MMAX+MMAX
            IF(sVOL(IJK,M) > ZERO) THEN

! Update the solids density for reacting cases.
               IF(CALC_sMASS) RO_s(IJK,M)=sMASS(IJK,M)/sVOL(IJK,M)

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
               ROP_S(IJK,M) = RO_S(IJK,M)*sVOL(IJK,M)/VOL(IJK)

               IF(MPPIC) THEN
                  OoSOLVOL = ONE/sVOL(IJK,M)
                  U_s(IJK,M) = U_s(IJK,M)*OoSOLVOL
                  V_s(IJK,M) = V_s(IJK,M)*OoSOLVOL
                  IF(DO_K) W_s(IJK,M) = W_s(IJK,M)*OoSOLVOL
               ENDIF
            ELSE
               ROP_S(IJK,M) = ZERO
            ENDIF

         ENDDO
      ENDDO
!$omp end parallel do

      end SUBROUTINE COMP_MEAN_FIELDS1
