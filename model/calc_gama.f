!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_GAMA(HEAT_TR, IER)                               C
!  Purpose: Calculate gas-solids heat transfer coefficients            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-JUL-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for variable grid size capability,           C
!           logic for volume-weighted averaging                        C
!  Author: W. Rogers                                  Date: 20-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Correct the heat transfer coefficient for transpiration    C
!           effect                                                     C
!  Author: M. Syamlal                                 Date: 25-AUG-93  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: MFIX 2.0 mods                                              C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Literature/Document References:                                     C
!    Gunn, D.J.,"Transfer of Heat or Mass to Particles in Fixed and    C
!      Fluidized Beds," Int. J. Heat Mass Transf., 21, 467(1978)       C
!    Bird, R.B., W.E. Stewart, E.N. Lightfoot, Transport Phenomena,    C
!      John Wiley & Sons, New York, 1960.                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: GAMA_gs                             C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_GAMA(HEAT_TR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
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
      USE energy
      USE rxns
      USE indices
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
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I, IJK
!
!                      Solids phase
      INTEGER          M
!
!                      Flag for exchange functions
      LOGICAL          HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M)
!
!                      Cell center value of U_g [Recall U_g (IJK) refers to
!                      U_g at I+1/2, J, K]
      DOUBLE PRECISION UGC
!
!                      Cell center value of U_sm
      DOUBLE PRECISION USCM
!
!                      Cell center value of V_g
      DOUBLE PRECISION VGC
!
!                      Cell center value of V_sm
      DOUBLE PRECISION VSCM
!
!                      Cell center value of W_g
      DOUBLE PRECISION WGC
!
!                      Cell center value of W_sm
      DOUBLE PRECISION WSCM
!
!                      Gas-solids relative velocity
      DOUBLE PRECISION VREL
!
!                      function of Prandtl number, Pr^(1/3)
      DOUBLE PRECISION Pr1o3
!
!                      Reynolds number, Re
      DOUBLE PRECISION Re
!
!                      EP_g^2
      DOUBLE PRECISION EP_g2
!
!                      a factor
      DOUBLE PRECISION FAC
!
      INTEGER          LM
!
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!        this needs to be generalized
      DO M = 1, MMAX 
         IF (HEAT_TR(0,M)) THEN 
!

!!$omp  parallel do &
!!$omp& private(I, IJK, M, EP_g2, Pr1o3, UGC, VGC, WGC, USCM, VSCM, &
!!$omp&  WSCM, VREL, Re, LM, FAC ) &
!!$omp& schedule(dynamic,chunk_size)
         DO IJK = ijkstart3, ijkend3
	 
               IF (FLUIDorP_FLOW_AT(IJK)) THEN 
                  I = I_OF(IJK) 
                  EP_G2 = EP_G(IJK)*EP_G(IJK) 
!
!  Calculate Prandtl number
!
                  if(K_G(IJK) > ZERO) then
                    PR1O3 = (C_PG(IJK)*MU_G(IJK)/K_G(IJK))**(1./3.) 
		  else
                    PR1O3 = LARGE_NUMBER 
		  endif
!
!  Calculate velocity components at the cell center for gas phase
!
                  UGC = HALF*(U_G(IJK)+U_G(IM_OF(IJK))) 
                  VGC = HALF*(V_G(IJK)+V_G(JM_OF(IJK))) 
                  WGC = HALF*(W_G(IJK)+W_G(KM_OF(IJK))) 
!
!  Calculate velocity components at the cell center for solids phase m
!
                  USCM = HALF*(U_S(IJK,M)+U_S(IM_OF(IJK),M)) 
                  VSCM = HALF*(V_S(IJK,M)+V_S(JM_OF(IJK),M)) 
                  WSCM = HALF*(W_S(IJK,M)+W_S(KM_OF(IJK),M)) 
!
!  Calculate the magnitude of gas-solids relative velocity
!
                  VREL=SQRT((UGC-USCM)**2+(VGC-VSCM)**2+(WGC-WSCM)**2) 
                  RE = EP_G(IJK)*D_P(M)*VREL*RO_G(IJK)/MU_G(IJK) 
!
!  Calculate gas-solids heat transfer coefficient (Gunn 1978)
!
                  GAMA_GS(IJK,M) = ((7. - 10.*EP_G(IJK)+5.*EP_G2)*(ONE+0.7*RE**&
                     0.2*PR1O3)+(1.33-2.4*EP_G(IJK)+1.2*EP_G2)*RE**0.7*PR1O3)*(&
                     K_G(IJK)/D_P(M))*(6.*EP_S(IJK,M)/D_P(M)) 
!
!  Correct the heat transfer coefficient for transpiration
!  Bird, Stewart, and Lightfoot (1960, p.663)
!
                  IF (GAMA_GS(IJK,M) > ZERO) THEN 
                     LM = 1 
                     FAC = R_PHASE(IJK,LM)*C_PG(IJK)/GAMA_GS(IJK,M) 
                     IF (ABS(FAC) < 0.1) THEN 
                        GAMA_GS(IJK,M)=GAMA_GS(IJK,M)/(ONE+FAC/2.+FAC*FAC/6.) 
                     ELSE 
                        IF (R_PHASE(IJK,LM) > ZERO) THEN 
                           GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)*EXP((-FAC&
                              ))/(ONE - EXP((-FAC))) 
                        ELSE 
                           GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)/(EXP(FAC)&
                               - ONE) 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

      call send_recv(GAMA_GS,2)

      RETURN  
      END SUBROUTINE CALC_GAMA 


!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
