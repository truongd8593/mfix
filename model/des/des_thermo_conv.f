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
      SUBROUTINE DES_CONVECTION(I, INTERP_IJK, INTERP_WEIGHTS, FOCUS)

      Use constant
      Use des_thermo
      Use discretelement
      Use fldvar
      Use interpolation
      Use param1
      Use physprop

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! Index of particle being looped over
      INTEGER, INTENT(IN) :: I
! IJK indicies of fluid cells involved in interpolation
      INTEGER, INTENT(IN) :: INTERP_IJK(2**DIMN)
! Weights associated with interpolation
      DOUBLE PRECISION, INTENT(IN) :: INTERP_WEIGHTS(2**DIMN)
! Indicates that debugging information for the particle
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!-----------------------------------------------      
! IJK value of cell containing particle I
      INTEGER IJK
! Convective heat transfer coefficient
      DOUBLE PRECISION GAMMA_CP
! Temperature of the gas
      DOUBLE PRECISION Tg
! Weight of specificed IJK value
      DOUBLE PRECISION WEIGHT
! Loop counter
      INTEGER LC
! Surface area of particle
      DOUBLE PRECISION A_S
      
! Initialization
       IJK = PIJK(I,4)

! Obtain the temperature of the gas.
      IF(.NOT.DES_INTERP_ON)THEN
         Tg = T_g(IJK)
      ELSE
         Tg = ZERO
         DO LC=1, 2**DIMN
            IJK = INTERP_IJK(LC)
            WEIGHT = INTERP_WEIGHTS(LC)
            Tg = Tg + (WEIGHT * T_g(IJK))
         ENDDO
      ENDIF

! Obtain the convective heat transfer coefficient of the particle
      CALL DES_CALC_GAMMA(I, GAMMA_CP)

! Calculate the surface area of the particle
      A_S = 4.0d0 * Pi * DES_RADIUS(I)**2

! Calculate the rate of heat transfer to the particle
      Qcv(I) = GAMMA_CP * A_S * (Tg - DES_T_s_NEW(I))

! Write out the debugging information.
      IF(DEBUG_DES .AND. FOCUS) THEN
         WRITE(*,"(//5X,A)")'From: DES_CONVECTION -'
         WRITE(*,"(8X,A,D12.6)")'Tg: ',Tg
         WRITE(*,"(8X,A,D12.6)")'Tp: ',DES_T_s_NEW(I)
         WRITE(*,"(8X,A,D12.6)")'GAMMA_CP: ',GAMMA_CP
         WRITE(*,"(8X,A,D12.6)")'Qcv: ',Qcv
         WRITE(*,"(5X,25('-')/)")
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

! Source term on LHS.  Must be positive. 
      DOUBLE PRECISION, INTENT(INOUT) :: S_P(DIMENSION_3)
! Source term on RHS 
      DOUBLE PRECISION, INTENT(INOUT) :: S_C(DIMENSION_3)

!-----------------------------------------------
! Local variables
!-----------------------------------------------      

! Index value of particle
      INTEGER NP
! IJK value of cell containing particle NP
      INTEGER IJK
! index to track accounted for particles      
      INTEGER PC

! Convective heat transfer coefficient
      DOUBLE PRECISION GAMMA_CP
! Logical used for debugging
      LOGICAL FOCUS
! Surface area of particle
      DOUBLE PRECISION A_S
! IJK values of cells surrounding the particle, including boundary
! conditions (i.e. periodic switches).
      INTEGER INTERP_IJK(2**DIMN)
! Weight of corresponding IJK fluid cell [0,1]
      DOUBLE PRECISION INTERP_WEIGHTS(2**DIMN)

! Weight of specificed IJK value
      DOUBLE PRECISION WEIGHT
! Loop counter
      INTEGER LC


! Initialize the particle counter
      PC = 1
! Loop over the particles in the system
      DO NP = 1, MAX_PIP
! Exit the loop if all the particles in the system have been looped over
         IF(PC .GT. PIP) EXIT
! Cycle the loop if there is no particle associated with this index
         IF(.NOT.PEA(NP,1)) CYCLE

! Set debug flag
         FOCUS = .FALSE.
         IF(NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.

! Obtain the convective heat transfer coefficient of the particle
         CALL DES_CALC_GAMMA(NP, GAMMA_CP)

! Calculate the surface area of the particle
         A_S = 4.0d0 * Pi * DES_RADIUS(NP)**2

! Obtain the temperature of the gas.
         IF(.NOT.DES_INTERP_ON)THEN

! Set the fluid cell index.
            IJK = PIJK(NP,4)
! Calculate the source term components
            S_P(IJK) = S_P(IJK) + ( GAMMA_CP * A_S )
            S_C(IJK) = S_C(IJK) + ( GAMMA_CP * A_S * DES_T_s_NEW(NP) )

         ELSE

! Determine the IJK values for the surrounding fluid cells and calculate
! the cell-centered interpolation weights.
            CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS)
! Distribute the thermal energy to/from the fluid cells
            DO LC = 1, 2**DIMN
! Set variable values
               IJK = INTERP_IJK(LC)
               WEIGHT = INTERP_WEIGHTS(LC)
! Portion the source term to cell IJK
               S_P(IJK) = S_P(IJK) + WEIGHT * (GAMMA_CP * A_S)
               S_C(IJK) = S_C(IJK) + WEIGHT * &
                  (GAMMA_CP * A_S * DES_T_s_NEW(NP))
            ENDDO ! LC-loop
         ENDIF ! End interpolation routine

! Reset the focus particle logical
         FOCUS = .FALSE.
! Index the particle count by one
         PC = PC + 1
      ENDDO

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
      SUBROUTINE DES_CALC_GAMMA(NP, GAMMA_CP)

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

! Index value of particle
      INTEGER, INTENT(IN) :: NP
! Convective heat transfer coefficient
      DOUBLE PRECISION, INTENT(OUT) :: GAMMA_CP

!-----------------------------------------------
! Local variables
!-----------------------------------------------      

! IJK value of cell containing particle NP
      INTEGER IJK
! Fluid cell indices
      DOUBLE PRECISION I, IMJK, IJMK, IJKM
! Double precision value for 1/3
      DOUBLE PRECISION, PARAMETER  :: THIRD = (1.0d0/3.0d0)
! Prandtl Number
      DOUBLE PRECISION Pr
! Particle Reynolds Number
      DOUBLE PRECISION Re_p 
! Nusselt Number
      DOUBLE PRECISION Nu_p
! Slip velocity between the particle and fluid
      DOUBLE PRECISION SLIP_VEL(3)
! Magnitude of slip velocity
      DOUBLE PRECISION SLIP_MAG
! Fluid velocity
      DOUBLE PRECISION UGC, VGC, WGC
!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'



! Initialization
      I =  PIJK(NP,1)
      IJK   = PIJK(NP,4)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

      SELECT CASE(TRIM(DES_CONV_CORR))

         CASE ('RANZ_1952') ! (Ranz and Mrshall, 1952)
! Initialize variables
            SLIP_VEL(:) = ZERO
            SLIP_MAG = ZERO
            Re_p = ZERO
            Nu_p = ZERO
! Gas velocity in fluid cell IJK
            UGC = AVG_X_E(U_g(IMJK), U_g(IJK), 1)
            VGC = AVG_Y_N(V_g(IJMK), V_g(IJK))
            WGC = AVG_Z_T(W_g(IJKM), W_g(IJK))
! Calculate the slip velocity between the particle and fluid
            SLIP_VEL(1) = UGC - DES_VEL_NEW(NP,1)
            SLIP_VEL(2) = VGC - DES_VEL_NEW(NP,2)
            IF(DIMN == 3) SLIP_VEL(3) = WGC - DES_VEL_NEW(NP,3)

! Calculate the magnitude of the slip velocity
            SLIP_MAG = SQRT(DES_DOTPRDCT(SLIP_VEL,SLIP_VEL))

! Calculate the Prandtl Number
            IF(K_G(IJK) > ZERO) THEN
               PR = (C_PG(IJK)*MU_G(IJK))/K_G(IJK)
		          ELSE
               PR = LARGE_NUMBER 
        		  ENDIF

! Calculate the particle Reynolds Number
            IF(MU_G(IJK) > ZERO) THEN
               Re_p = (2.0d0*DES_RADIUS(NP)*SLIP_MAG*RO_g(IJK)) / &
                  MU_g(IJK)
            ELSE
               Re_p = LARGE_NUMBER
            ENDIF

! Calculate the Nusselt Number
            Nu_p = 2 + 0.6d0 * (Re_p)**HALF * (Pr)**THIRD

! Calculate the convective heat transfer coefficient
            GAMMA_CP = (Nu_p * K_G(IJK))/(2.0d0 * DES_RADIUS(NP))

         CASE DEFAULT
            WRITE(*,*)'INVALID DES CONVECTION MODEL'
            STOP
            CALL MFIX_EXIT
      END SELECT

      RETURN
      END SUBROUTINE DES_CALC_GAMMA

