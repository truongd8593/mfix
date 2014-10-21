!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GP                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_gs and calc_des_drag_gs              C
!                                                                      C
!  Comments: The BVK drag model and all drag models with the           C
!            polydisperse correction factor (i.e., suffix _PCF)        C
!            require an average particle diameter. This has been       C
!            loosely defined for discrete particles based on their     C
!            solids phase                                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_DRAG_GP(LL, FLUID_VEL, PARTICLE_VEL)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar
      USE drag
      USE sendrecv
      USE discretelement
      USE ur_facs
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) :: LL
! fluid velocity interpolated to particle position
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: FLUID_VEL
! particle velocity
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: PARTICLE_VEL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices, associated with current particle
      INTEGER :: IJK
! solids phase index, associated with current particle
      INTEGER :: M
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! drag coefficient
      DOUBLE PRECISION :: DgA
! current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
! indices of solids phases (continuous, discrete)
      INTEGER :: CM, DM, L
! temporary shift of total number of solids phases to account for both
! discrete and continuous solids phases used for the hybrid mdoel
      INTEGER :: MAXM
! tmp local variable for the particle diameter of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: DP_loc(2*DIM_M)
! tmp local variable for the solids volume fraction of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: EPs_loc(2*DIM_M)
! tmp local variable for the particle density of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: ROs_loc(2*DIM_M)
! correction factors for implementing polydisperse drag model
! proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION :: F_cor, tmp_sum, tmp_fac
! average particle diameter in polydisperse systems
      DOUBLE PRECISION :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION :: Y_i
! total solids volume fraction
      DOUBLE PRECISION :: phis
! aliases for void fraction, gas density, gas bulk density,
! solids volume fraction, particle diameter, particle density
      DOUBLE PRECISION :: EPG, ROg, ROPg, EP_SM, DPM, ROs
!-----------------------------------------------

! values based on current particle
      IJK = PIJK(LL,4)
! solids phase index of current particle
      M = PIJK(LL,5)

! Assign local variables DP_loc, EPs_loc, and MAXM.  These
! represent arrays for the particle diameter, solids volume
! fraction, and number of particle types (i.e., phases).
      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
         MAXM = DES_MMAX
         DO DM = 1,MAXM
            DP_loc(DM) = DES_D_p0(DM)
            EPs_loc(DM) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
            ROs_loc(DM) = DES_RO_S(DM)
         ENDDO
      ELSE   ! des_continuum_hybrid branch
! For the hybrid model the diameters and solids volume fractions of
! of both discrete and continuous are stored in this single quantity.
! Any loops of solids phases will include all solids phases (discrete
! and continuum)
         MAXM = SMAX + DES_MMAX
! populate DP, EPS starting with discrete phases
         DO DM = 1,DES_MMAX
            DP_loc(DM) = DES_D_p0(DM)
            EPs_loc(DM) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
            ROs_loc(DM) = DES_RO_S(DM)
         ENDDO
         DO CM = 1,SMAX
            L = DES_MMAX + CM
            DP_loc(L) = D_P(IJK,CM)
            EPs_loc(L) = EP_S(IJK,CM)
            ROs_loc(L) = RO_S(IJK,CM)
         ENDDO
      ENDIF   ! end if/else (.not.des_continuum_hybrid)


! magnitude of gas-particle relative velocity
      IF(NO_K)THEN
         VREL = SQRT((FLUID_VEL(1) - PARTICLE_VEL(1))**2 +&
                     (FLUID_VEL(2) - PARTICLE_VEL(2))**2)
      ELSE
         VREL = SQRT((FLUID_VEL(1) - PARTICLE_VEL(1))**2 +&
                     (FLUID_VEL(2) - PARTICLE_VEL(2))**2 +&
                     (FLUID_VEL(3) - PARTICLE_VEL(3))**2)
      ENDIF

! Laminar viscosity at a pressure boundary is given the value of the
! fluid cell next to it. This applies just to the calculation of the
! drag, in other routines the value of viscosity at a pressure boundary
! always has a zero value.
! This will never happen since this subroutine is currently only called
! for fluid_at cells (does not include flow boundaries)
! This points to an inconsitency in calculation of drag between
! continuum and discrete models that is probably not addressed in the
! solution of the gas phase momentum balances
      IF (P_OUTFLOW_AT(IJK)) THEN
         IF( FLUID_AT(EAST_OF(IJK) )) THEN
            Mu = MU_G(EAST_OF(IJK))
         ELSE IF ( FLUID_AT(WEST_OF(IJK)) ) THEN
            Mu = MU_G(WEST_OF(IJK))
         ELSE IF ( FLUID_AT(NORTH_OF(IJK)) ) THEN
            Mu = MU_G(NORTH_OF(IJK))
         ELSE IF ( FLUID_AT(SOUTH_OF(IJK)) ) THEN
            Mu = MU_G(SOUTH_OF(IJK))
         ELSE IF ( FLUID_AT(TOP_OF(IJK)) ) THEN
            Mu = MU_G(TOP_OF(IJK))
         ELSE IF ( FLUID_AT(BOTTOM_OF(IJK)) ) THEN
            Mu = MU_G(BOTTOM_OF(IJK))
         ENDIF
      ELSE
         Mu = MU_G(IJK)
      ENDIF

! calculate the total solids volume fraction
      phis = ZERO
      DO L = 1, MAXM
! this is slightly /= one-ep_g due to round-off
         phis = phis + EPs_loc(L)
      ENDDO

! calculate the average paricle diameter and particle ratio
      DPA = ZERO
      tmp_sum = ZERO
      tmp_fac = ZERO
      DO L = 1, MAXM
         IF (phis .GT. ZERO) THEN
            tmp_fac = EPs_loc(L)/phis
            tmp_sum = tmp_sum + tmp_fac/DP_loc(L)
          ELSE
            tmp_sum = tmp_sum + ONE/DP_loc(L) ! not important, but will avoid NaN's in empty cells
          ENDIF
      ENDDO
      DPA = ONE / tmp_sum
      Y_i = DP_loc(M) * tmp_sum

! assign variables for short dummy arguments
      EPg = EP_G(IJK)
      ROg = RO_G(IJK)
      ROPg = ROP_G(IJK)
      EP_SM = EPs_loc(M)
      DPM = DP_loc(M)
      ROs = ROs_loc(M)

! determine the drag coefficient
      IF (EP_SM <= ZERO) THEN
! this won't happen in DEM case since routine is performed over
! particles not cells as in continuum case
         DgA = ZERO
      ELSEIF (EPg == ZERO) THEN
! this case will already be caught in most drag subroutines whenever
! RE==0 (for correlations in which RE includes EPg). however, this will
! prevent potential divisions by zero in some models by setting it now.
         DgA = ZERO
      ELSE
! determine the drag coefficient
         SELECT CASE(DRAG_TYPE_ENUM)
         CASE (SYAM_OBRIEN)
            CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,DPM)
         CASE (GIDASPOW)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)
         CASE (GIDASPOW_PCF)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (GIDASPOW_BLEND)
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)
         CASE (GIDASPOW_BLEND_PCF)
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (WEN_YU)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM)
         CASE (WEN_YU_PCF)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPA)
         CASE (KOCH_HILL)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)
         CASE (KOCH_HILL_PCF)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (BVK)
            CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (USER_DRAG)
            CALL DRAG_USR(IJK, M, DgA, EPg, Mu, ROg, VREL, DPM, ROs)
         CASE DEFAULT
            CALL START_LOG
            IF(DMP_LOG) WRITE (*, '(A,A)') &
               'Unknown DRAG_TYPE: ', DRAG_TYPE
            WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
            CALL END_LOG
            CALL mfix_exit(myPE)
         END SELECT   ! end selection of drag_type
      ENDIF   ! end if/elseif/else (ep_sm <= zero, ep_g==0)


! Modify drag coefficient to account for possible corrections and
! for differences between Model B and Model A
      IF(DRAG_TYPE_ENUM == GIDASPOW_PCF .OR. &
         DRAG_TYPE_ENUM == GIDASPOW_BLEND_PCF .OR. &
         DRAG_TYPE_ENUM == WEN_YU_PCF .OR. &
         DRAG_TYPE_ENUM == KOCH_HILL_PCF .OR. &
         DRAG_TYPE_ENUM == BVK) THEN
! see erratum by Beetstra et al. (2007) : the correction factor differs
! for model A versus model B.
! application of the correction factor for model A is found from
! the correction factor for model B and neglects the Y_i**3 term
         IF(Model_B) THEN
            IF (M == 1) THEN
               F_cor = (EPg*Y_i + phis*Y_i**2)
            ELSE
               F_cor = (EPg*Y_i + phis*Y_i**2 + &
                  0.064d0*EPg*Y_i**3)
            ENDIF
         ELSE
            F_cor = Y_i
         ENDIF
         DgA = ONE/(Y_i*Y_i) * DgA * F_cor
      ENDIF

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      IF(MODEL_B) THEN
         F_gstmp = DgA * PVOL(LL)/EP_G(IJK)
      ELSE
         F_gstmp = DgA * PVOL(LL)
      ENDIF

! Determine drag force coefficient accounting for any under relaxation
! f_gp() =  single particle drag excluding vector(v_g - v_p)
      F_gp(LL) = (ONE - UR_F_gs) * F_gp(LL) + UR_F_gs * F_gstmp

      RETURN
      END SUBROUTINE DES_DRAG_GP




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Subroutine: DRAG_INTERPOLATION                                       C
!  Purpose: DES - Calculate the fluid velocity interpolated at the      C
!           particle's location and weights. Replace 'interpolator'     C
!                       interface for OpenMP implementation.            C
!                                                                       C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_INTERPOLATION(GSTEN,VSTEN,DESPOS,VELFP,WEIGHTFACTOR)

      use geometry, only: NO_K

        IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: GSTEN
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: VSTEN
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: DESPOS
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: VELFP
        DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(OUT) :: WEIGHTFACTOR
        INTEGER :: II, JJ, KK

        DOUBLE PRECISION, DIMENSION(2) :: XXVAL, YYVAL, ZZVAL
        DOUBLE PRECISION :: DXX, DYY, DZZ
        DOUBLE PRECISION, DIMENSION(3) :: ZETAA

        DXX = GSTEN(2,1,1,1) - GSTEN(1,1,1,1)
        DYY = GSTEN(1,2,1,2) - GSTEN(1,1,1,2)

        ZETAA(1:2) = DESPOS(1:2) - GSTEN(1,1,1,1:2)

        ZETAA(1) = ZETAA(1)/DXX
        ZETAA(2) = ZETAA(2)/DYY

        XXVAL(1)=1-ZETAA(1)
        YYVAL(1)=1-ZETAA(2)
        XXVAL(2)=ZETAA(1)
        YYVAL(2)=ZETAA(2)

        VELFP(:) = 0.D0

        IF(NO_K) THEN
           DO JJ=1,2
              DO II=1,2
                 WEIGHTFACTOR(II,JJ,1) = XXVAL(II)*YYVAL(JJ)
                 VELFP(1:2) = VELFP(1:2) + VSTEN(II,JJ,1,1:2)*WEIGHTFACTOR(II,JJ,1)
              ENDDO
           ENDDO
        ELSE
           DZZ = GSTEN(1,1,2,3) - GSTEN(1,1,1,3)
           ZETAA(3) = DESPOS(3) - GSTEN(1,1,1,3)
           ZETAA(3) = ZETAA(3)/DZZ
           ZZVAL(1)=1-ZETAA(3)
           ZZVAL(2)=ZETAA(3)
           DO KK=1,2
              DO JJ=1,2
                 DO II=1,2
                    WEIGHTFACTOR(II,JJ,KK) = XXVAL(II)*YYVAL(JJ)*ZZVAL(KK)
                    VELFP(1:3) = VELFP(1:3) + VSTEN(II,JJ,KK,1:3)*WEIGHTFACTOR(II,JJ,KK)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

      END SUBROUTINE DRAG_INTERPOLATION
