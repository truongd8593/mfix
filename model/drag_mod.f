!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: drag                                                   C
!  Purpose: Common block containing drag arrays                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE drag
 
! Gas-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_gs

! Gas-solids drag coefficient/ep_s
! Needed in computation of some conditional limits of various KT
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: dgA_s

! Solids-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

! Off diagonal friction coefficient in HYS drag relation
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij 

! Temporary storage: Volume x average at momentum cell centers
!      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  VxF_gs
 

      END MODULE drag
