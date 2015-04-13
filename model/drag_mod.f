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

! Solids-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

! Off diagonal friction coefficient in HYS drag relation
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij


      END MODULE drag
