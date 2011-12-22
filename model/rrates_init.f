!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES_INIT(IER)                                       C
!  Purpose: Initialize reaction rate arrays                            C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES_INIT(IER) 
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE indices
      USE compar       
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      cell index
      INTEGER          IJK
!C
!-----------------------------------------------
      INCLUDE 'function.inc'
      
!!!$omp  parallel do private(ijk)
      DO IJK = IJKSTART3, IJKEND3 
      
         R_gp(IJK, :) = ZERO
         RoX_gc(IJK, :) = ZERO
         R_sp(IJK, :, :) = ZERO
         RoX_sc(IJK, :, :) = ZERO
         SUM_R_G(IJK) = ZERO 
         HOR_G(IJK) = ZERO
         SUM_R_S(IJK, :) = ZERO 
         HOR_S(IJK, :) = ZERO 
	 R_PHASE(IJK, :) = ZERO
	 
      END DO
      RETURN
      END
