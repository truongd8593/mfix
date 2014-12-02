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
      use discretelement, only: DES_CONTINUUM_HYBRID
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE 
      use particle_filter, only: DES_INTERP_GARG

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL DRAG_GS_DES_NONINTERP
         CASE(DES_INTERP_GARG) ; CALL DRAG_GS_DES_INTERP0
         CASE DEFAULT; CALL DRAG_GS_DES_INTERP1
         END SELECT
      ENDIF

! Calculate solids-solids drag force on particle. 
      IF(DES_CONTINUUM_HYBRID) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE DEFAULT; CALL DRAG_SS_DEM_NONINTERP
         END SELECT
      ENDIF

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
      use discretelement, only: DES_CONTINUUM_HYBRID
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE 
      use particle_filter, only: DES_INTERP_GARG

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL DRAG_GS_GAS_NONINTERP
         CASE(DES_INTERP_GARG) ; CALL DRAG_GS_GAS_INTERP0
         CASE DEFAULT; CALL DRAG_GS_GAS_INTERP1
         END SELECT
      ENDIF

! Calculate solids-solids drag force on particle. 
      IF(DES_CONTINUUM_HYBRID) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE DEFAULT; CALL DRAG_SS_TFM_NONINTERP
         END SELECT
      ENDIF

      RETURN
      END SUBROUTINE CALC_DRAG_DES_2FLUID
