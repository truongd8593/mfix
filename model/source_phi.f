!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_phi(S_p, S_c, EP, Phi, M, A_m, B_m, IER)        C
!  Purpose: Determine source terms for phi eq. The terms               C
!  appear in the center coefficient and RHS vector.  The center        C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.  S_p must be positive.                   C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE SOURCE_PHI(S_P, S_C, EP, PHI, M, A_M, B_M, IER) 
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
      USE matrix 
      USE scales 
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      Source term on LHS.  Must be positive. 
      DOUBLE PRECISION S_p(DIMENSION_3) 
! 
!                      Source term on RHS 
      DOUBLE PRECISION S_C(DIMENSION_3) 
! 
!                      Phase volume fraction 
      DOUBLE PRECISION EP(DIMENSION_3) 
! 
!                      Phi 
      DOUBLE PRECISION Phi(DIMENSION_3) 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE(2) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!
      DO IJK = 1, IJKMAX2 
!
         IF (FLUID_AT(IJK)) THEN 
!
!             dilute flow
            IF (EP(IJK) <= DIL_EP_S) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
!
               IF (EP(WEST_OF(IJK)) > DIL_EP_S) THEN 
                  A_M(IJK,W,M) = ONE 
               ELSE IF (EP(EAST_OF(IJK)) > DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ONE 
               ELSE 
                  B_M(IJK,M) = -PHI(IJK) 
               ENDIF 
!
!             Normal case
            ELSE 
!
!               Collect the terms
               A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,&
                  S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+S_P(IJK))
!
               B_M(IJK,M) = -S_C(IJK)+B_M(IJK,M)
!			
	    ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_PHI 
