!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VAVG_W_g                                               C
!  Purpose: Calculate volume averaged W_g                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-APR-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION VAVG_W_G () 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      Integral of W_g*EP_g for entire volume 
      DOUBLE PRECISION SUM_W_g 
! 
!                      Total volume of computational cells 
      DOUBLE PRECISION SUM_VOL 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Integrate the velocity values for the whole domain,
!
      SUM_W_G = ZERO 
      SUM_VOL = ZERO 

!!$omp   parallel do private(IJK) reduction(+:SUM_VOL,SUM_W_G)

      DO IJK = IJKMIN1, IJKMAX1 
         IF (FLUID_AT(IJK)) THEN 
            SUM_VOL = SUM_VOL + VOL_W(IJK) 
            SUM_W_G = SUM_W_G + W_G(IJK)*EP_G(IJK)*VOL_W(IJK) 
         ENDIF 
      END DO 
      VAVG_W_G = SUM_W_G/SUM_VOL 
!
      RETURN  
      END FUNCTION VAVG_W_G 
