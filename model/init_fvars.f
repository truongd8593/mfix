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
      IJK = 1 
      IF (IJKMAX2 > 0) THEN 
!// 200 1010 modified the upper limit from :ijkmax2 --> 0:ijkmax3
!//         EP_G(:IJKMAX2) = UNDEFINED       
         EP_G(0:IJKMAX3) = UNDEFINED 
         P_G(0:IJKMAX3) = UNDEFINED 
         P_STAR(0:IJKMAX3) = ZERO 
         RO_G(0:IJKMAX3) = UNDEFINED 
         ROP_G(0:IJKMAX3) = UNDEFINED 
         T_G(0:IJKMAX3) = ZERO 
         U_G(0:IJKMAX3) = UNDEFINED 
         V_G(0:IJKMAX3) = UNDEFINED 
         W_G(0:IJKMAX3) = UNDEFINED 
         N = 1 
         IF (NMAX(0) > 0) THEN 
            X_G(0:IJKMAX3,:NMAX(0)) = ZERO 
            N = NMAX(0) + 1 
         ENDIF 
!//? what is IJK used for?	 
         IJK = IJKMAX2 + 1 
      ENDIF 

!!$omp parallel do private(M,IJK,N)
      DO M = 1, MMAX 
         IJK = 1 
         IF (IJKMAX2 > 0) THEN 
!// 200 1010 modified the upper limit from :ijkmax2 --> 0:ijkmax3
            ROP_S(0:IJKMAX3,M) = UNDEFINED 
            T_S(0:IJKMAX3,M) = ZERO 
            THETA_M(0:IJKMAX3,M) = ZERO 
            U_S(0:IJKMAX3,M) = UNDEFINED 
            V_S(0:IJKMAX3,M) = UNDEFINED 
            W_S(0:IJKMAX3,M) = UNDEFINED 
            N = 1 
            IF (NMAX(M) > 0) THEN 
               X_S(0:IJKMAX3,M,:NMAX(M)) = ZERO 
               N = NMAX(M) + 1 
            ENDIF 
!//? what is IJK used for?	 	    
            IJK = IJKMAX2 + 1 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE INIT_FVARS 
