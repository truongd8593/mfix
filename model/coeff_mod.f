      MODULE coeff


      Use param
      Use param1


!
!                      Flags to tell whether to calculate or not
      LOGICAL, DIMENSION(:), ALLOCATABLE ::  DENSITY, PSIZE,SP_HEAT 
!
!                      Flags to tell whether to calculate or not
      LOGICAL, DIMENSION(:), ALLOCATABLE ::  VISC, COND, DIFF
!
!                      Flag for Reaction rates
      LOGICAL          RRATE
!
!                      Flag for exchange functions
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: DRAGCOEF, HEAT_TR
      LOGICAL          WALL_TR
 
!


      END MODULE coeff                                                                           
