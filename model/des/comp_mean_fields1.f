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
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_CELL
      use physprop, only: mmax
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: particles, filter cells, phases
      INTEGER NP, LC, M, MMAX_TOT
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight, and times filter weight
      DOUBLE PRECISION :: VOL_WT, VOLxWEIGHT
! Loop bound for filter
      INTEGER :: LP_BND


!-----------------------------------------------


      SOLVOLINC(:,:) = ZERO
      MMAX_TOT = DES_MMAX+MMAX

      IF(MPPIC) THEN
! initialize only information related to the discrete 'phases' of these 
! continuous variables
         U_S(:,MMAX+1:MMAX_TOT) = ZERO
         V_S(:,MMAX+1:MMAX_TOT) = ZERO
         IF(DO_K) W_S(:,MMAX+1:MMAX_TOT) = ZERO
      ENDIF

! Loop bounds for interpolation.
      LP_BND = merge(27,9,DO_K)

! Calculate the gas phase forces acting on each particle.
!$omp parallel default(none) &
!$omp private(NP, VOL_WT, M, LC, IJK, VOLXWEIGHT) &
!$omp shared(MAX_PIP, PVOL, DES_STAT_WT, PIJK, LP_BND, MPPIC, &
!$omp       FILTER_WEIGHT, SOLVOLINC, U_S, V_S, W_S, DO_K, &
!$omp       FILTER_CELL, DES_VEL_NEW)
!$omp do
      do NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP) .and. .NOT.IS_GHOST(NP)) CYCLE

         VOL_WT = PVOL(NP)
         IF(MPPIC) VOL_WT = VOL_WT*DES_STAT_WT(NP)
! Particle phase for data binning.
         M = PIJK(NP,5)

         DO LC=1,LP_BND
            IJK = FILTER_CELL(LC,NP)
! Particle volume times the weight for this cell.
            VOLxWEIGHT = VOL_WT*FILTER_WEIGHT(LC,NP)
! Accumulate total solids volume (by phase)
!$omp atomic
            SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + VOLxWEIGHT
            IF(MPPIC) THEN
! Accumulate total solids momentum (by phase)
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
         ENDDO
      ENDDO
!$omp end do
!$omp end parallel

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!---------------------------------------------------------------------//
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)                &
!$omp private(IJK,M,OoSOLVOL)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = MMAX+1, DES_MMAX+MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO .AND.MPPIC) THEN
               OoSOLVOL = ONE/SOLVOLINC(IJK,M)
               U_s(IJK,M) = U_s(IJK,M)*OoSOLVOL
               V_s(IJK,M) = V_s(IJK,M)*OoSOLVOL
               IF(DO_K) W_s(IJK,M) = W_s(IJK,M)*OoSOLVOL
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            ROP_S(IJK,M) = RO_S(IJK,M)*SOLVOLINC(IJK,M)/VOL(IJK)

         ENDDO   ! end loop over M=MMAX+1,DES_MMAX+MMAX

      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do


! Halo exchange of solids volume fraction data.
      calL SEND_RECV(ROP_S,2)

      end SUBROUTINE COMP_MEAN_FIELDS1
