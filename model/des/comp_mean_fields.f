!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: COMP_MEAN_FIELDS                                        !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose: Driver routine for calculating continuous field variables  !
!  corresponding to discrete data (ROP_s, u_s, v_s, w_s)               !
!                                                                      !
!  o The diffusion filter is only applied to the the solids bulk       !
!    density because DEM simulations do not utilize the other field    !
!    variables within a time loop.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE COMP_MEAN_FIELDS

! Modules
!---------------------------------------------------------------------//
      use discretelement, only: DES_MMAX
      use fldvar, only: rop_s
! Flag: Diffuse DES field variables.
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use physprop, only: mmax
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//

! Loop counter.
      INTEGER :: M

!......................................................................!

! Calculate field variables from particle data:
      IF(DES_INTERP_MEAN_FIELDS) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL COMP_MEAN_FIELDS_ZERO_ORDER
         CASE(DES_INTERP_GARG) ; CALL COMP_MEAN_FIELDS0
         CASE DEFAULT; CALL COMP_MEAN_FIELDS1
         END SELECT
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

! Apply the diffusion filter.
      IF(DES_DIFFUSE_MEAN_FIELDS) THEN
         DO M=MMAX+1, MMAX+DES_MMAX
            CALL DIFFUSE_MEAN_FIELD(ROP_S(:,M),'ROP_S')
         ENDDO
      ENDIF

! Calculate the gas phase volume fraction from ROP_s.
      CALL CALC_EPG_DES

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

! Modules
!---------------------------------------------------------------------//
      USE param
      USE param1
      USE fldvar
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
      USE physprop, only: MMAX
      USE run, only: solids_model
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Loop counters: partciles, filter cells, phases
      INTEGER NP, M
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight
      DOUBLE PRECISION :: VOL_WT
! total number of 'solids' phases in simulation
      INTEGER :: MMAX_TOT
!......................................................................!
      SOLVOLINC(:,:) = ZERO

      MMAX_TOT = DES_MMAX+MMAX
! initialize only information related to the discrete 'phases' of
! these continuous variables
      U_s(:,MMAX+1:MMAX_TOT) = ZERO
      V_s(:,MMAX+1:MMAX_TOT) = ZERO
      W_s(:,MMAX+1:MMAX_TOT) = ZERO

! Calculate the gas phase forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(IS_NONEXISTENT(NP)) CYCLE
         IF(IS_GHOST(NP) .or. IS_ENTERING_GHOST(NP) .or. &
            IS_EXITING_GHOST(NP)) CYCLE

         VOL_WT = PVOL(NP)
         IF(MPPIC) VOL_WT = VOL_WT*DES_STAT_WT(NP)
! Fluid cell containing the particle
         IJK = PIJK(NP,4)
! Particle phase for data binning.
         M = PIJK(NP,5)
! Accumulate total solids volume (by phase)
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + VOL_WT
! Accumulate total solids momenum-ish (by phase)
         U_S(IJK,M) = U_S(IJK,M) +                             &
            DES_VEL_NEW(NP,1)*VOL_WT
         V_S(IJK,M) = V_S(IJK,M) +                             &
            DES_VEL_NEW(NP,2)*VOL_WT
         IF(DO_K) W_S(IJK,M) = W_S(IJK,M) +                    &
            DES_VEL_NEW(NP,3)*VOL_WT
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!----------------------------------------------------------------//
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = MMAX+1, MMAX_TOT
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OoSOLVOL = ONE/SOLVOLINC(IJK,M)
               U_s(IJK,M) = U_s(IJK,M)*OoSOLVOL
               V_s(IJK,M) = V_s(IJK,M)*OoSOLVOL
               IF(DO_K) W_s(IJK,M) = W_s(IJK,M)*OoSOLVOL
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            ROP_S(IJK,M) = RO_S(IJK,M)*SOLVOLINC(IJK,M)/VOL(IJK)

         ENDDO   ! end loop over M=MMAX+1,MMAX_TOT

      ENDDO     ! end loop over IJK=ijkstart3,ijkend3

! Halo exchange of solids volume fraction data.
      CALL SEND_RECV(ROP_S,2)

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER
