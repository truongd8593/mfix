!TO DO:
! check the formulation based on MCp.
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_1(IER)                                         C
!  Purpose: Correct the solids volume fraction                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-SEP-96  C
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
      SUBROUTINE CORRECT_1(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE indices
      USE geometry
      USE pscor
      USE ur_facs 
      USE constant
      USE compar 
      USE sendrecv 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      error index
      INTEGER          IER
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      corrected solids volume fraction
      DOUBLE PRECISION EPcor
!
!                      dPodEP_s(EP_s(IJK, M)) * EPp(IJK)
      DOUBLE PRECISION Pp_P

      INTEGER          IJK, IJKE, IJKN, IJKT, M
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'
!
!      DO 200 M = 1, MMAX
      M = MCP 
      IF (CLOSE_PACKED(M)) THEN 
!
!$omp    parallel do &
!$omp&   private( IJK, EPCOR )
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN 
               EPCOR = EP_S(IJK,M) + EPP(IJK) 
               IF (EPCOR>EP_S_CP .AND. EPP(IJK)>ZERO) THEN 
                  EPP(IJK) = UR_FAC(2)*EPP(IJK) 
                  EPCOR = EP_S(IJK,M) + EPP(IJK) 
               ENDIF 
               ROP_S(IJK,M) = MAX(ZERO,RO_S(M)*EPCOR) 
            ENDIF 
         END DO 

!$omp    parallel do &
!$omp&   private( IJK, PP_P, IJKE, IJKN, IJKT )
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN 
!
               PP_P = K_CP(IJK)*EPP(IJK) 
!
               IF (FLOW_AT_E(IJK)) THEN 
                  IJKE = EAST_OF(IJK) 
                  U_S(IJK,M)=U_S(IJK,M)-E_E(IJK)*(K_CP(IJKE)*EPP(IJKE)-PP_P) 
               ENDIF 
!
               IF (FLOW_AT_N(IJK)) THEN 
                  IJKN = NORTH_OF(IJK) 
                  V_S(IJK,M)=V_S(IJK,M)-E_N(IJK)*(K_CP(IJKN)*EPP(IJKN)-PP_P) 
               ENDIF 
!
               IF (DO_K) THEN 
                  IF (FLOW_AT_T(IJK)) THEN 
                     IJKT = TOP_OF(IJK) 
                     W_S(IJK,M) = W_S(IJK,M) - E_T(IJK)*(K_CP(IJKT)*EPP(IJKT)-&
                        PP_P) 
                  ENDIF 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE CORRECT_1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
