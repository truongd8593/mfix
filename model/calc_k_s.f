!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                           C
!  Module name: CALC_K_s(M, IER)                                            C
!  Purpose: Calculate the effective conductivity of solids phases           C
!                                                                           C
!  Author:M. Syamlal                                  Date: 24-APR-96       C
!  Reviewer:                                          Date: dd-mmm-yy       C
!                                                                           C
!  Revision Number: 01                                                      C
!  Purpose: (1) allow to use Bauer & Schlunder's (1978) model in CGS or SI  C
!           (2) If fluid_at(IJK) condition for the Bauer & Schlunder's modelC 
!  Author:  S. Dartevelle                             Date: 10-July-02      C
!  Reviewer:                                          Date: dd-mmm-yy       C
!                                                                           C
!  Literature/Document References:                                          C
!                                                                           C
!  Variables referenced:                                                    C
!  Variables modified:                                                      C
!                                                                           C
!  Local variables:                                                         C
!                                                                           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_K_S(M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop 
      USE fldvar 
      USE geometry 
      USE indices 
      USE constant 
      USE toleranc  
      USE compar 
      USE sendrecv 
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!
!                      Two
      DOUBLE PRECISION, PARAMETER          :: TWO = 2.0d0
      
!              microscopic conductivity of solids in cal/s.cm.K (not modified by the Gas Phase)
      DOUBLE PRECISION Ks_micro
      PARAMETER (Ks_micro = 0.5258E-2)    !(2.2 J/s.m.K, conductivity of Ash) S. Dartevelle
!
!              Constant in conductivity equation
      DOUBLE PRECISION PHI_k
      PARAMETER (PHI_k = 7.26E-3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK
!
!                      Solids phase
      INTEGER          M
!
!                      Quantities in solids conductivity formula
      DOUBLE PRECISION B, R_km, BoR, L_rm
!
!                      Transform K_g(IJK) into the CGS if we work with SI
      DOUBLE PRECISION Kg_micro          !S. Dartevelle, July 9th, 2002
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!1 Cal = 4.183925 J
!
      IF (K_S0 /= UNDEFINED) RETURN  

!!$omp parallel do private(IJK,B,R_km,BoR,L_rm,Kg_micro) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3            
!
        IF (FLUID_AT(IJK)) THEN
          IF (UNITS == 'SI') THEN
            Kg_micro = K_g(IJK)/418.3925  !back to CGS for K_g, Cal/s.cm.K
	  ELSE
            Kg_micro = K_g(IJK)           !it's already within the CGS system
          ENDIF
!
!         All the calculation are in the CGS system
!         Bauer & Schlunder's (1978) theory:
          B = 1.25 * ((ONE - EP_g(IJK))/EP_g(IJK))**(10./9.)
          IF( (ONE - EP_g(IJK)) .GT. DIL_EP_s) THEN
            R_km = Ks_micro/Kg_micro
            BoR  = B/R_km
            L_rm = -(TWO/(ONE - BoR))&
                * ( ((R_km - ONE)/(ONE - BoR)**2)*BoR*LOG(BoR)&
                  + (B - ONE)/(ONE - BoR) + (B + ONE)/TWO)
            !K_s is the macroscopic conductivity, modified by the presence of the gas phase
            K_S(IJK, M) = (Phi_k*R_km + (ONE - Phi_k)*L_rm)*Kg_micro&
			                    /SQRT(ONE - EP_g(IJK))             !Cal/s.cm.K
          ELSE
            K_S(IJK, M) = ZERO
          ENDIF
!
         !An approximate average value for the solids conductivity is 2.5*K_g
!	  K_S(IJK,M) = 2.5*Kg_micro            !in CGS system
!
        ELSE
	  K_S(IJK,M) = ZERO
        ENDIF
!
!to SI, S. Dartevelle:
        IF (UNITS == 'SI') K_s(IJK, M) = 418.3925*K_s(IJK, M)           !J/s.m.K
!
      END DO
!
      CALL send_recv(K_S, 2)
!
      RETURN
      END SUBROUTINE CALC_K_S

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
