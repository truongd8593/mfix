!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_SMASS(SMASS)                                       C
!  Purpose: Determine the weight of solids in the reactor              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, MMAX, ROP_s, DX, DY, DZ, C
!                        X, IMIN1, JMIN1. KMIN1                        C
!  Variables modified: I, J, K, M, IJK                                 C
!                                                                      C
!  Local variables:  None                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_SMASS(SMASS) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE compar 
      USE mpi_utility 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Total Weight of solids in the reactor 
      DOUBLE PRECISION SMASS 
! 
!                      Weight of mth solids phase 
      DOUBLE PRECISION SUM 
! 
!                      Indices 
      INTEGER          I, J, K, IJK 
! 
!                      Solids phase 
      INTEGER          M 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      SMASS = ZERO 
      DO M = 1, MMAX 
         SUM = ZERO 

!!$omp$   parallel do private(IJK) &
!!$omp&   reduction(+:SUM)
         DO IJK = ijkstart3, ijkend3 
	 IF(.NOT.IS_ON_myPE_owns(I_OF(IJK), J_OF(IJK), K_OF(IJK))) cycle
            IF (FLUID_AT(IJK)) SUM = SUM + ROP_S(IJK,M)*VOL(IJK) 
         END DO 
         SMASS = SMASS + SUM 
      END DO 

      call global_all_sum(smass)
      
      RETURN  
      END SUBROUTINE GET_SMASS 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added mpi_utility module and other global reduction (global_all_sum) call

