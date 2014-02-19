!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: drag.inc                                               C
!  Purpose: Common block containing drag arrays                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE drag
 
! Gas-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_gs

! Solids-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

! Off diagonal friction coefficient in HYS drag relation
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij 

! Temporary storage: Volume x average at momentum cell centers
!      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  VxF_gs
 

! Number specifying the C_Ds function to use.
      INTEGER :: FUN_C_DS_ENUM

! Enumerated list of C_Ds functions for calculating either the
! 1) Single sphere drag correlation 
! 2) Single sphere drag correlation multipling Reynolds number or
      INTEGER, PARAMETER :: C_DS_SN    = 1  ! Schiller and Naumann (1933)
      INTEGER, PARAMETER :: C_DSXRE_DV = 2  ! Dalla Valle (1948)
      INTEGER, PARAMETER :: C_DS_DEL   = 3  ! Dellino et al. (2005)
      INTEGER, PARAMETER :: C_DSxRE_TL = 4  ! Turton and Levenspiel (1986)

!!!HPF$ align F_gs(:, *) with TT(:)
!!!HPF$ align F_ss(:, *) with TT(:)

      END MODULE drag
