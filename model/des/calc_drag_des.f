!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DRAG_DES                                           !
!                                                                      !
!  Purpose: This subroutine is called from DES routines. It calls      !
!  functions that calcultate the drag force acting on particles. No    !
!  field variables are updated.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DRAG_DES

      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_INTERP_ON
      use discretelement, only: DES_CONTINUUM_HYBRID

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_INTERP_ON) THEN
            CALL DRAG_GS_DES_INTERP0
         ELSE
            CALL DRAG_GS_DES_NONINTERP
         ENDIF
      ENDIF

! Calculate solids-solids drag force on particle. 
      IF(DES_CONTINUUM_HYBRID) CALL DRAG_SS_DEM_NONINTERP


      RETURN
      END SUBROUTINE CALC_DRAG_DES






!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GS                                             C
!                                                                      C
!  Purpose: This subroutine is only called from the CONTINUUM side.    C
!     It performs the following functions:                             C
!     - If non-interpolated, then execution of the code is directed    C
!       to the subroutine drag_gs for the appropriate calculations     C
!       (i.e., calculation of the fluid-solids drag coefficient)       C
!     - If interpolated then, it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. It then determines the            C
!       the contributions of fluid-particle drag to the center         C
!       coefficient of the A matrix and the b (source) vector in the   C
!       matrix equation (A*VEL_FP=b) equation for the fluid phase      C
!       x, y and z momentum balances using F_GP.                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DRAG_DES_2FLUID

      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_INTERP_ON
      use discretelement, only: DES_CONTINUUM_HYBRID

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DES_INTERP_ON) THEN
            CALL DRAG_GS_GAS_INTERP0
         ELSE
            CALL DRAG_GS_GAS_NONINTERP
         ENDIF
      ENDIF

! Calculate solids-solids drag force on particle. 
      IF(DES_CONTINUUM_HYBRID) CALL DRAG_SS_TFM_NONINTERP


      RETURN
      END SUBROUTINE CALC_DRAG_DES_2FLUID
