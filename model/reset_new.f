!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RESET_NEW                                              C
!  Purpose: Reset the new variables with the stored previous-time-step C
!           values of field variables                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: FEB-6-97   C
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
      SUBROUTINE RESET_NEW 
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
      USE indices
      USE physprop
      USE run
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
      INTEGER :: IJK, M, N 
!-----------------------------------------------
!
!
!
!
!                    Indices
!
      IJK = 1 
      IF (IJKMAX2 > 0) THEN 
         EP_G(:IJKMAX2) = EP_GO(:IJKMAX2) 
         P_G(:IJKMAX2) = P_GO(:IJKMAX2) 
         P_STAR(:IJKMAX2) = P_STARO(:IJKMAX2) 
         RO_G(:IJKMAX2) = RO_GO(:IJKMAX2) 
         ROP_G(:IJKMAX2) = ROP_GO(:IJKMAX2) 
         U_G(:IJKMAX2) = U_GO(:IJKMAX2) 
         V_G(:IJKMAX2) = V_GO(:IJKMAX2) 
         W_G(:IJKMAX2) = W_GO(:IJKMAX2) 
         IF (ENERGY_EQ) T_G(:IJKMAX2) = T_GO(:IJKMAX2) 
         IF (SPECIES_EQ(0)) THEN 
            N = 1 
            IF (NMAX(0) > 0) THEN 
               X_G(:IJKMAX2,:NMAX(0)) = X_GO(:IJKMAX2,:NMAX(0)) 
               N = NMAX(0) + 1 
            ENDIF 
         ENDIF 
         IJK = IJKMAX2 + 1 
      ENDIF 
      DO M = 1, MMAX 


         IJK = 1 
         IF (IJKMAX2 > 0) THEN 
            ROP_S(:IJKMAX2,M) = ROP_SO(:IJKMAX2,M) 
            THETA_M(:IJKMAX2,M) = THETA_MO(:IJKMAX2,M) 
            IF (ENERGY_EQ) T_S(:IJKMAX2,M) = T_SO(:IJKMAX2,M) 
            U_S(:IJKMAX2,M) = U_SO(:IJKMAX2,M) 
            V_S(:IJKMAX2,M) = V_SO(:IJKMAX2,M) 
            W_S(:IJKMAX2,M) = W_SO(:IJKMAX2,M) 
            IF (SPECIES_EQ(M)) THEN 
               N = 1 
               IF (NMAX(M) > 0) THEN 
                  X_S(:IJKMAX2,M,:NMAX(M)) = X_SO(:IJKMAX2,M,:NMAX(M)) 
                  N = NMAX(M) + 1 
               ENDIF 
            ENDIF 
            IJK = IJKMAX2 + 1 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE RESET_NEW 
