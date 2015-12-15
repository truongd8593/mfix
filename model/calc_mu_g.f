!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: MFIX 2.0 mods (previous name CALC_MU_gt)                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI unit                                              C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: compute turbulent eddy viscosity                           C
!  Author: S. Benyahia                                Date: May-13-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!     Perry, R. H., and Chilton, C. H., Chemical Engineers' Handbook,  C
!        5th Edition, McGraw-Hill Inc., 1973, pp. 248, eqn. 3-133.     C
!     Arnold, J. H., Vapor viscosities and the Sutherland equation,    C
!        Journal of Chemical Physics, 1 (2), 1933, pp. 170-176.        C
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      C
!        Phil. Mag. 5:507-531, 1893.                                   C
!                                                                      C
!     Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent   C
!        flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,    C
!        No. 6, pp. 1203-1228.                                         C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_MU_G()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE indices
      USE constant
      USE toleranc
      USE compar
      USE drag
      USE run          !S. Dartevelle
      USE turb
      USE sendrecv
      USE mms
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Cell indices
      INTEGER :: IJK
! Solids phase index
      INTEGER :: M
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION :: D_g(3,3)
! Gas velocity gradient
      DOUBLE PRECISION :: DelV_g(3,3)
! Second invariant of the deviator of D_g
      DOUBLE PRECISION :: I2_devD_g
! Constant in turbulent viscosity formulation
      DOUBLE PRECISION :: C_MU
! particle relaxation time
      DOUBLE PRECISION :: Tau_12_st
!-----------------------------------------------

! JFD: Calling SET_EP_FACTORS here to make sure EPG_IFAC is defined
!      when calc_mu_g is called for the firt time.
!      There may be a better way to do this...
      CALL SET_EP_FACTORS

! solids phase index used throughout routine...
! may be inappropriate for multiple solids phases
      M = 1 ! for solids phase


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

            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK)= EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)

! K_epsilon model
! ---------------------------------------------------------------->>>
            IF (K_Epsilon) THEN

               C_MU = 9D-02

! I'm not very confident about this correction in Peirano paper, but it's made
! available here, uncomment to use it. sof@fluent.com --> 02/01/05
!               IF(SIMONIN .AND. F_GS(IJK,1) > SMALL_NUMBER) THEN
!                  Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
!                  X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))
! new definition of C_mu (equation A.12, Peirano et al. (2002) Powder tech. 122,69-82)
!                  IF( K_12(IJK)/(2.0D0*K_Turb_G(IJK)) < ONE) &
!                     C_MU = C_MU/(ONE+ 0.314D0*X_21*Tau_12_st / Tau_1(IJK) * &
!                            (ONE - K_12(IJK)/(2.0D0*K_Turb_G(IJK))) )
!               ENDIF
! Correction in Ahmadi paper (Cao and Ahmadi)
               IF(AHMADI .AND. F_GS(IJK,1) > SMALL_NUMBER) THEN
                  Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
                  C_MU = C_MU/(ONE+ Tau_12_st/Tau_1(IJK) * &
                         (EP_s(IJK,M)/(ONE-EP_star_array(IJK)))**3)
               ENDIF

! Definition of the turbulent viscosity
               MU_GT(IJK) = MU_G(IJK) + RO_G(IJK)*C_MU*&
                  K_Turb_G(IJK)**2 / (E_Turb_G(IJK) + SMALL_NUMBER)
               MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
               LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
               EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
               EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
            ENDIF
! ----------------------------------------------------------------<<<

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
      END IF ! end if (USE_MMS)



! L_scale0 model
! ---------------------------------------------------------------->>>
!!$omp parallel do &
!!$omp$ schedule(dynamic,chunk_size) &
!!$omp$ private(IJK, DelV_g, D_G, I2_DEVD_G )
      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) .AND. L_SCALE(IJK)/=ZERO) THEN

! Calculate the rate of strain tensor D_g
            CALL CALC_DERIV_VEL_GAS(ijk, DelV_G, D_G)

! Calculate the second invariant of the deviator of D_g
            I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+&
                         (D_G(3,3)-D_G(1,1))**2)/6.D0 + &
                        D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2

            MU_GT(IJK) = MIN(MU_GMAX, MU_G(IJK)+2.0*&
               L_SCALE(IJK)*L_SCALE(IJK)*RO_G(IJK)*SQRT(I2_DEVD_G))
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)

! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
         ENDIF ! end if (fluid_at(ijk) and l_scale(ijk)/=0))
      ENDDO    ! end loop (ijk=ijkstart3,ijkend3)
! end calculations for L_scale0 model
! ----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE CALC_MU_G

