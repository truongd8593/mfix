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
 
! Number specifying the C_Ds function to use.
      INTEGER :: CD_FUNCTION_ENUM

! Enumerated list of C_Ds functions for calculating either the
! 1) Single sphere drag correlation 
! 2) Single sphere drag correlation multipling Reynolds number or
      INTEGER, PARAMETER :: SCHILLER_1933 = 1  ! Schiller and Naumann (1933)
      INTEGER, PARAMETER :: DALLA_1948    = 2  ! Dalla Valle (1948)
      INTEGER, PARAMETER :: DELLINO_2005  = 3  ! Dellino et al. (2005)
      INTEGER, PARAMETER :: TURTON_1986   = 4  ! Turton and Levenspiel (1986)


      END MODULE drag
