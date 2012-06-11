!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  MODULE: COEFF                                                       !
!  Purpose: Contains logic flags that tells the code whether to        !
!           perform the indicated type of calculation when the         !
!           value is true                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE coeff

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use param
      Use param1
!-----------------------------------------------

! Density, particle size, specific heat
! Density is really only invoked for the gas/fluid phase, that is,
! density(1:MMAX) is never used      
      LOGICAL, DIMENSION(:), ALLOCATABLE :: DENSITY, &
                                            PSIZE, SP_HEAT 
              !(0:DIMENSION_M)

! Viscosity, conductivity, diffusivity
      LOGICAL, DIMENSION(:), ALLOCATABLE :: VISC, COND, DIFF
              !(0:DIMENSION_M)

! Granular energy dissipation term
! gran_diss(0) is meaningless              
      LOGICAL, DIMENSION(:), ALLOCATABLE :: GRAN_DISS
              !(0:DIMENSION_M)

! Reaction rates
      LOGICAL :: RRATE

! Exchange functions
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: DRAGCOEF, HEAT_TR
               !(0:DIMENSION_M,0:DIMENSION_M)
      LOGICAL :: WALL_TR
 

      END MODULE coeff
