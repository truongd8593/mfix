!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CONVECTION                                         !
!                                                                      !
!  Purpose: This routine is called from the discrete phase. It is used !
!  to determine the rate of heat transfer to the particle from the gas.!
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!  REF: Zhou, Yu, and Zulli, "Particle scale study of heat transfer in !
!       packed and bubbling fluidized beds," AIChE Journal, Vol. 55,   !
!       no 4, pp 868-884, 2009.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_CONVECTION(NP, M, IJK, &
         INTERP_IJK, INTERP_WEIGHTS, FOCUS)

      Use constant
      Use des_thermo
      Use discretelement
      Use fldvar
      Use interpolation
      Use param1
      Use physprop
      Use usr

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Index of particle being looped over
      INTEGER, INTENT(IN) :: NP
! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: M
! IJK value of cell containing particle NP
      INTEGER, INTENT(IN) :: IJK
! IJK indicies of fluid cells involved in interpolation
      INTEGER, INTENT(IN) :: INTERP_IJK(2**DIMN)
! Weights associated with interpolation
      DOUBLE PRECISION, INTENT(IN) :: INTERP_WEIGHTS(2**DIMN)
! Indicates that debugging information for the particle
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!---------------------------------------------------------------------//
! Convective heat transfer coefficient
      DOUBLE PRECISION GAMMA_CP
! Temperature of the gas
      DOUBLE PRECISION Tg
! Weight of specificed IJK value
      DOUBLE PRECISION WEIGHT
! Surface area of particle
      DOUBLE PRECISION Sa
! Convection source
      DOUBLE PRECISION Qcv
      

! Obtain the temperature of the gas. --> Not interpolated.
      Tg = T_g(IJK)
! Obtain the convective heat transfer coefficient of the particle
      CALL DES_CALC_GAMMA(NP, IJK, GAMMA_CP)
! Calculate the surface area of the particle
      Sa = 4.0d0 * Pi * DES_RADIUS(NP) * DES_RADIUS(NP)
! Calculate the rate of heat transfer to the particle
      Qcv = GAMMA_CP * Sa * (Tg - DES_T_s_NEW(NP))
! Store convection source in global energy source array.
      Q_Source(NP) = Q_Source(NP) + Qcv

! Write out the debugging information.
      IF(FOCUS) THEN
         WRITE(*,"(/5X,A)")'From: DES_CONVECTION -'
         WRITE(*,"(8X,A,D12.6)")'Tg: ',Tg
         WRITE(*,"(8X,A,D12.6)")'Tp: ',DES_T_s_NEW(NP)
         WRITE(*,"(8X,A,D12.6)")'GAMMA_CP: ',GAMMA_CP
         WRITE(*,"(8X,A,D12.6)")'Qcv: ',Qcv
         WRITE(*,"(5X,50('-')/)")
      ENDIF

      RETURN

      END SUBROUTINE DES_CONVECTION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_Hgm                                                 !
!                                                                      !
!  Purpose: This routine is called from the continuum phase and        !
!  calculates the source term from the particles to the fluid. This    !
!  routine handles both the interpolated and non-interpolated options. !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Jan-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE  DES_Hgm(S_C, S_P)

      USE compar
      Use constant
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use interpolation
      Use param1
      Use physprop

      IMPLICIT NONE

! Passed Variables
!---------------------------------------------------------------------//
! Source term on LHS.  Must be positive. 
      DOUBLE PRECISION, INTENT(INOUT) :: S_P(DIMENSION_3)
! Source term on RHS 
      DOUBLE PRECISION, INTENT(INOUT) :: S_C(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Index value of particle
      INTEGER lNP, NP
! IJK value of cell containing particle NP
      INTEGER IJK
! Convective heat transfer coefficient
      DOUBLE PRECISION GAMMA_CP
! Surface area of particle
      DOUBLE PRECISION Sa

! Functions
!---------------------------------------------------------------------//
      INCLUDE 'function.inc'

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Interpolation.
!---------------------------------------------------------------------//
! Currently omitted.
!   -> Set the IJK stencil

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)
! Skip indices that do not represent particles
            IF(.NOT.PEA(NP,1)) CYCLE lNP_LP
! Skip indices that represent ghost particles
            IF(PEA(NP,4)) CYCLE lNP_LP

! Obtain the convective heat transfer coefficient of the particle
            CALL DES_CALC_GAMMA(NP, IJK, GAMMA_CP)

! Calculate the surface area of the particle
            Sa = 4.0d0 * Pi * DES_RADIUS(NP)**2

! Calculate the source term components --> Not interpolated.
            S_P(IJK) = S_P(IJK) + (GAMMA_CP * Sa)
            S_C(IJK) = S_C(IJK) + (GAMMA_CP * Sa * DES_T_s_NEW(NP))

         ENDDO lNP_LP ! End loop over particles
      ENDDO IJK_LP ! End loop over fluid cells

      RETURN
      END SUBROUTINE  DES_Hgm


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_CALC_GAMMA                                          !
!                                                                      !
!  Purpose: Calculate the heat transfer coefficient (GAMMA_CP) for     !
!  particle-fluid heat transfer.                                       !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!  REF: Ranz, W.E. and Marshall, W.R., "Friction and transfer          !
!       coefficients for single particles and packed beds," Chemical   !
!       Engineering Science, Vol. 48, No. 5, pp 247-253, 1925.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_CALC_GAMMA(NP, IJK, GAMMA_CP)

      USE compar
      Use constant
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use param1
      Use physprop

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Index value of particle
      INTEGER, INTENT(IN) :: NP
! Index value of fluid cell
      INTEGER, INTENT(IN) :: IJK
! Convective heat transfer coefficient
      DOUBLE PRECISION, INTENT(OUT) :: GAMMA_CP

! Local variables
!---------------------------------------------------------------------//
! Fluid cell indices
      INTEGER IMJK, IJMK, IJKM
! Double precision value for 1/3
      DOUBLE PRECISION, PARAMETER  :: THIRD = (1.0d0/3.0d0)

      DOUBLE PRECISION N_Pr  ! Prandtl Number
      DOUBLE PRECISION N_Re  ! Reynolds Number
      DOUBLE PRECISION N_Nu  ! Nusselt Number

! Magnitude of slip velocity
      DOUBLE PRECISION SLIP
! Fluid velocity
      DOUBLE PRECISION cUg, cVg, cWg
      DOUBLE PRECISION Us, Vs, Ws

! Functions
!---------------------------------------------------------------------//
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

! Initialization
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

      SELECT CASE(TRIM(DES_CONV_CORR))

         CASE ('RANZ_1952') ! (Ranz and Mrshall, 1952)
! Initialize variables
            SLIP = ZERO
            N_Re = ZERO
            N_Nu = ZERO
! Gas velocity in fluid cell IJK
            cUg = AVG_X_E(U_g(IMJK), U_g(IJK), 1)
            cVg = AVG_Y_N(V_g(IJMK), V_g(IJK))
! Particle Velocity
            Us = DES_VEL_NEW(NP,1)
            Vs = DES_VEL_NEW(NP,2)

! Calculate the magnitude of the slip velocity
            IF(DIMN == 2) THEN
               SLIP = SQRT((cUg-Us)**2 + (cVg-Vs)**2)
            ELSE
               cWg = AVG_Z_T(W_g(IJKM), W_g(IJK))
               Ws = DES_VEL_NEW(NP,3)
               SLIP = SQRT((cUg-Us)**2 + (cVg-Vs)**2 + (cWg-Ws)**2)
            ENDIF

! Calculate the Prandtl Number
            IF(K_G(IJK) > ZERO) THEN
               N_Pr = (C_PG(IJK)*MU_G(IJK))/K_G(IJK)
		          ELSE
               N_Pr = LARGE_NUMBER 
        		  ENDIF

! Calculate the particle Reynolds Number
            IF(MU_G(IJK) > ZERO) THEN
               N_Re = (2.0d0*DES_RADIUS(NP)*SLIP*RO_g(IJK)) / MU_g(IJK)
            ELSE
               N_Re = LARGE_NUMBER
            ENDIF

! Calculate the Nusselt Number
            N_Nu = 2.0d0 + 0.6d0 *((N_Re)**HALF * (N_Pr)**THIRD)

! Calculate the convective heat transfer coefficient
            GAMMA_CP = (N_Nu * K_G(IJK))/(2.0d0 * DES_RADIUS(NP))

         CASE DEFAULT
            WRITE(*,*)'INVALID DES CONVECTION MODEL'
            STOP
            CALL MFIX_EXIT
      END SELECT

      RETURN
      END SUBROUTINE DES_CALC_GAMMA
