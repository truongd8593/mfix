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
!//EFD Nov 11, convert (:IJKMAX2) to just (:)
!// let the compiler determine the array bounds
!
!
         EP_G(:) = EP_GO(:) 
         P_G(:) = P_GO(:) 
         P_STAR(:) = P_STARO(:) 
         RO_G(:) = RO_GO(:) 
         ROP_G(:) = ROP_GO(:) 
         U_G(:) = U_GO(:) 
         V_G(:) = V_GO(:) 
         W_G(:) = W_GO(:) 
         IF (ENERGY_EQ) T_G(:) = T_GO(:) 
         IF (SPECIES_EQ(0)) THEN 
            N = 1 
            IF (NMAX(0) > 0) THEN 
               X_G(:,:NMAX(0)) = X_GO(:,:NMAX(0)) 
               N = NMAX(0) + 1 
            ENDIF 
          ENDIF 
      DO M = 1, MMAX 


            ROP_S(:,M) = ROP_SO(:,M) 
            THETA_M(:,M) = THETA_MO(:,M) 
            IF (ENERGY_EQ) T_S(:,M) = T_SO(:,M) 
            U_S(:,M) = U_SO(:,M) 
            V_S(:,M) = V_SO(:,M) 
            W_S(:,M) = W_SO(:,M) 
            IF (SPECIES_EQ(M)) THEN 
               N = 1 
               IF (NMAX(M) > 0) THEN 
                  X_S(:,M,:NMAX(M)) = X_SO(:,M,:NMAX(M)) 
                  N = NMAX(M) + 1 
               ENDIF 
            ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE RESET_NEW 
