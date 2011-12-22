! Note: This routine is now restricted to Non-Negative scalers when using
! deferred correction. 
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
      USE compar  
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
      INTEGER          IJK, IMJK, IJMK, IJKM 
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
!!!$omp    parallel do private(IJK)
      DO IJK = ijkstart3, ijkend3 
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
! using the average boundary cell values to compute phi (sof, Aug 23 2005)
!
               IMJK = IM_OF(IJK)
	       IJMK = JM_OF(IJK)
	       IJKM = KM_OF(IJK)
               IF (EP(WEST_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE
               IF (EP(EAST_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
               IF (EP(SOUTH_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE
               IF (EP(NORTH_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
               IF(.NOT. NO_K) THEN
	         IF (EP(BOTTOM_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE
                 IF (EP(TOP_OF(IJK)) > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
	       ENDIF 
!               
	       IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
	           A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
	          B_M(IJK,M) = -PHI(IJK)              
	       ELSE
	         A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                  A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
	       ENDIF
!
!             Normal case
            ELSE 
!
!               Collect the terms
               A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,&
                  S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+S_P(IJK))
!
! B_m and A_m are corrected in case deferred corrections computes B_m > S_c
! see CONV_DIF_PHI_DC.
               IF(B_M(IJK,M) < S_C(IJK) .OR. PHI(IJK) == ZERO) THEN
                 B_M(IJK,M) = -S_C(IJK)+B_M(IJK,M)
               ELSE ! disable ELSE statememt if PHI can be negative
                 A_M(IJK,0,M) = A_M(IJK,0,M) - B_M(IJK,M)/PHI(IJK)
                 B_M(IJK,M) = -S_C(IJK)
               ENDIF
!			
	    ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_PHI 

