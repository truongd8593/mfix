!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VAVG_U_s(M)                                            C
!  Purpose: Volume average U_s                                         C
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
      DOUBLE PRECISION FUNCTION VAVG_U_S (M) 
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
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER          M
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      Integral of U_s*EP_s for entire volume 
      DOUBLE PRECISION SUM_U_s 
! 
!                      Total volume of computational cells 
      DOUBLE PRECISION SUM_VOL 
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!  Integrate the velocity values for the whole domain,
!
      SUM_U_S = ZERO 
      SUM_VOL = ZERO 

!!$omp   parallel do private(IJK) reduction(+:SUM_VOL,SUM_U_S)

      DO IJK = IJKMIN1, IJKMAX1 
         IF (FLUID_AT(IJK)) THEN 
            SUM_VOL = SUM_VOL + VOL_U(IJK) 
            SUM_U_S = SUM_U_S + U_S(IJK,M)*EP_S(IJK,M)*VOL_U(IJK) 
         ENDIF 
      END DO 
      VAVG_U_S = SUM_U_S/SUM_VOL 
!
      RETURN  
      END FUNCTION VAVG_U_S 
