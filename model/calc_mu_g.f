!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!     Perry, R. H., and Chilton, C. H., Chemical Engineers' Handbook,  C
!        5th Edition, McGraw-Hill Inc., 1973, pp. 248, eqn. 3-133.     C
!     Arnold, J. H., Vapor viscosities and the Sutherland equation,    C
!        Journal of Chemical Physics, 1 (2), 1933, pp. 170-176.        C
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      C
!        Phil. Mag. 5:507-531, 1893.                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MU_G()

! Modules
!-----------------------------------------------
      use compar, only: ijkstart3, ijkend3
      use constant, only: to_si
      use fldvar, only: T_g, epg_ifac
      use functions, only: fluid_at
      use mms, only: use_mms
      use param1, only: undefined, zero
      use physprop, only: mu_g0, mu_g
      use run, only: k_epsilon
      use visc_g, only: l_scale
      use visc_g, only: mu_gt, epmu_gt
      use visc_g, only: lambda_gt, eplambda_gt
      IMPLICIT NONE

! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!-----------------------------------------------
! Cell indices
      INTEGER :: IJK
!-----------------------------------------------

!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's formula with
! Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
! For air  C = 110 (Tb=74.82)
!         mu = 1.71*10-4 poise at T = 273K
            IF (MU_G0 == UNDEFINED) MU_G(IJK) = to_SI*1.7D-4 * &
               (T_G(IJK)/273.0D0)**1.5D0 * (383.D0/(T_G(IJK)+110.D0))

! set viscosity values
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)

! adjust viscosity for tubulence
            IF (K_Epsilon) THEN
               CALL CALC_K_EPSILON_MU(IJK)
            ELSEIF (L_SCALE(IJK) /=ZERO) THEN
               CALL CALC_LSCALE_MU(IJK)
            ENDIF

! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK)= EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)

         ELSE
            MU_G(IJK)  = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
            EPMU_GT(IJK) = ZERO
            EPLAMBDA_GT(IJK) = ZERO
         ENDIF   ! end if (fluid_at(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


! MMS: Force constant gas viscosity at all cells including ghost cells.
      IF (USE_MMS) THEN
         DO IJK = ijkstart3, ijkend3
            MU_G(IJK) = MU_G0
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
         ENDDO
      ENDIF ! end if (USE_MMS)


      RETURN
      END SUBROUTINE CALC_MU_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: compute turbulent eddy viscosity                           C
!  Author: S. Benyahia                                Date: May-13-04  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!   Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent     C
!      flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,      C
!      No. 6, pp. 1203-1228.                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_K_EPSILON_MU(IJK)

! Modules
!-----------------------------------------------
      use constant, only: mu_gmax
      use drag, only: f_gs
      use fldvar, only: k_turb_g, e_turb_g, ro_g
      use fldvar, only: ep_s, ro_s
      use param1, only: one, small_number
      use physprop, only: mu_g
      use run, only: kt_type_enum, ahmadi_1995
      use turb, only: tau_1
      use visc_g, only: mu_gt, lambda_gt
      use visc_s, only: ep_star_array
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Cell indices
      INTEGER, INTENT(IN) :: IJK

! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!-----------------------------------------------
! Solids phase index
      INTEGER :: M
! Constant in turbulent viscosity formulation
      DOUBLE PRECISION :: C_MU
! particle relaxation time
      DOUBLE PRECISION :: Tau_12_st
!-----------------------------------------------

      C_MU = 9D-02

! I'm not very confident about this correction in Peirano paper,
! but it's made available here, uncomment to use it.
! sof@fluent.com --> 02/01/05
!      IF(KT_TYPE_ENUM==SIMONIN_1996 .AND.&
!         F_GS(IJK,1) > SMALL_NUMBER) THEN
! solids phase index used throughout routine...
!         M = 1 ! for solids phase
!         Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
!         X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))
! new definition of C_mu (equation A.12, Peirano et al. (2002),
! Powder tech. 122,69-82)
!         IF( K_12(IJK)/(2.0D0*K_Turb_G(IJK)) < ONE) &
!            C_MU = C_MU/(ONE+ 0.314D0*X_21*Tau_12_st / Tau_1(IJK) * &
!                   (ONE - K_12(IJK)/(2.0D0*K_Turb_G(IJK))) )
!      ENDIF

! Correction in Ahmadi paper (Cao and Ahmadi)
      IF(KT_TYPE_ENUM == AHMADI_1995 .AND.&
         F_GS(IJK,1) > SMALL_NUMBER) THEN
! solids phase index used throughout routine...
         M = 1 ! for solids phase
         Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
         C_MU = C_MU/(ONE+ Tau_12_st/Tau_1(IJK) * &
               (EP_s(IJK,M)/(ONE-EP_star_array(IJK)))**3)
      ENDIF

! Definition of the turbulent viscosity
      MU_GT(IJK) = MU_G(IJK) + RO_G(IJK)*C_MU*&
         K_Turb_G(IJK)**2 / (E_Turb_G(IJK) + SMALL_NUMBER)

      MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
      LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)

      RETURN
      END SUBROUTINE CALC_K_EPSILON_MU


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: compute l_scale0 model                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_LSCALE_MU(IJK)

! Modules
!-----------------------------------------------
      use constant, only: mu_gmax
      use fldvar, only: ro_g
      use physprop, only: mu_g
      use visc_g, only: mu_gt, lambda_gt
      use visc_g, only: l_scale
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Cell indices
      INTEGER, INTENT(IN) :: IJK

! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0

! Local variables
!-----------------------------------------------
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION :: D_g(3,3)
! Gas velocity gradient
      DOUBLE PRECISION :: DelV_g(3,3)
! Second invariant of the deviator of D_g
      DOUBLE PRECISION :: I2_devD_g
!-----------------------------------------------

! Calculate the rate of strain tensor D_g
      CALL CALC_DERIV_VEL_GAS(ijk, DelV_G, D_G)

! Calculate the second invariant of the deviator of D_g
      I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+&
                   (D_G(3,3)-D_G(1,1))**2)/6.D0 + &
                  D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2

      MU_GT(IJK) = MIN(MU_GMAX, &
          MU_G(IJK)+2.0*L_SCALE(IJK)*L_SCALE(IJK)*&
          RO_G(IJK)*SQRT(I2_DEVD_G))
      LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)

      RETURN
      END SUBROUTINE CALC_LSCALE_MU
