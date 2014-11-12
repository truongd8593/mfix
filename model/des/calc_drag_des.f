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


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_GS                                             !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  calls the correct routine for calculating the gas drag force.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
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
