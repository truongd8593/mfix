!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MW (X_g, DIM, L, NMAX, MW_g)                      C
!  Purpose: Calculate average molecular weight of gas                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-OCT-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:None                                           C
!  Variables modified:None                                             C
!                                                                      C
!  Local variables: SUM, N                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION CALC_MW (X_G, DIM, L, NMAX, MW_G) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Mass fraction array's Ist dimension
      INTEGER          DIM
!
!
!                      Max of X_g array 2nd index and MW_g array index
      INTEGER          NMAX
!
!                      Mass fraction array
!
      DOUBLE PRECISION X_g(DIM, NMAX)
!
!                      Moleculare weight array
!
      DOUBLE PRECISION MW_g(NMAX)
!
!                      Mass fraction array Ist index
      INTEGER          L
!
!  Local variable
!
!                      local sum
      DOUBLE PRECISION SUM
!
!                      local index
      INTEGER          N
!-----------------------------------------------
!
      SUM = ZERO 
      DO N = 1, NMAX 
         SUM = SUM + X_G(L,N)/MW_G(N) 
      END DO 
      CALC_MW = ONE/MAX(SUM,OMW_MAX) 
!
      RETURN  
      END FUNCTION CALC_MW 
      
!// Comments on the modifications for DMP version implementation      
!// 100 EFD Replaced inheritance based dimensioning of arrays, i.e. X_g(*),MW_g(*)
