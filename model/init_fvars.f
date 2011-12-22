!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_FVARS                                             C
!  Purpose: Initialize all field variables as undefined                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-JAN-94  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: ROP_g, EP_g, ROP_s, IJKMAX2, MMAX, U_s, V_s,  C
!                        W_s                                           C
!                                                                      C
!  Variables modified: ROP_go, ROP_so, IJK, M, U_so, V_so, W_so C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE INIT_FVARS 
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
      USE geometry
      USE physprop
      USE indices
      USE scalars
      USE rxns
      USE run
      USE compar
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
!                      Solids phase
      INTEGER          M
!
!                      Species index
      INTEGER          N
!-----------------------------------------------
      IF (IJKMAX2 > 0) THEN 
         EP_G(IJKSTART3:IJKEND3) = UNDEFINED 
         P_G(IJKSTART3:IJKEND3) = UNDEFINED 
         P_STAR(IJKSTART3:IJKEND3) = ZERO 
         RO_G(IJKSTART3:IJKEND3) = UNDEFINED 
         ROP_G(IJKSTART3:IJKEND3) = UNDEFINED 
         T_G(IJKSTART3:IJKEND3) = ZERO 
         U_G(IJKSTART3:IJKEND3) = UNDEFINED 
         V_G(IJKSTART3:IJKEND3) = UNDEFINED 
         W_G(IJKSTART3:IJKEND3) = UNDEFINED 
         IF (NMAX(0) > 0) THEN 
            X_G(IJKSTART3:IJKEND3,:NMAX(0)) = ZERO 
         ENDIF  
      
         IF(Nscalar > 0) Scalar(IJKSTART3:IJKEND3,:Nscalar) = ZERO
         IF(nRR > 0) ReactionRates(IJKSTART3:IJKEND3,:nRR) = ZERO

         IF(K_Epsilon) THEN
           K_Turb_G(IJKSTART3:IJKEND3) = ZERO
           E_Turb_G(IJKSTART3:IJKEND3) = ZERO
         ENDIF
      ENDIF

!!!!$omp parallel do private(M,IJK,N)
      DO M = 1, MMAX 
         IF (IJKMAX2 > 0) THEN 
            ROP_S(IJKSTART3:IJKEND3,M) = UNDEFINED 
            T_S(IJKSTART3:IJKEND3,M) = ZERO 
! add by rong 
            D_P(IJKSTART3:IJKEND3,M)=ZERO
! add by rong
            THETA_M(IJKSTART3:IJKEND3,M) = ZERO 
            P_S(IJKSTART3:IJKEND3,M) = UNDEFINED 
            U_S(IJKSTART3:IJKEND3,M) = UNDEFINED 
            V_S(IJKSTART3:IJKEND3,M) = UNDEFINED 
            W_S(IJKSTART3:IJKEND3,M) = UNDEFINED 
	    P_S(IJKSTART3:IJKEND3,M) = UNDEFINED
	    KTH_S(IJKSTART3:IJKEND3,M) = UNDEFINED
            IF (NMAX(M) > 0) THEN 
               X_S(IJKSTART3:IJKEND3,M,:NMAX(M)) = ZERO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE INIT_FVARS 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 120 Replaced the index for initialization :ijkmax2 --> IJKSTART3:IJKEND3
