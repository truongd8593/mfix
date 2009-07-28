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
 
 
      Use param
      Use param1
 
 
!
!                      Gas-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_gs
!
!                      Solids-solids drag
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

!			Off diagonal friction coefficient in HYS drag relation
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij 

!
!                   temporary storage: Volume x average at momentum cell centers
!      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  VxF_gs
 
!!!HPF$ align F_gs(:, *) with TT(:)
!!!HPF$ align F_ss(:, *) with TT(:)

      END MODULE drag
