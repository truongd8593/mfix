!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_CONSTProp                                          C
!  Purpose: This module sets all the constant physical properties      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-MAY-97  C
!                                                                      C
!  Local variables: IJK
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_CONSTPROP 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE visc_s
      USE visc_g
      USE energy
      USE geometry
      USE indices
      USE physprop
      USE constant
      USE run
      USE funits 
      USE drag
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
      INTEGER :: IJK, M, N 
!-----------------------------------------------
      INCLUDE 'function.inc'

! Initialize transport coefficients to zero everywhere

      MU_gt(:) = ZERO
      LAMBDA_GT(:) = ZERO
      MU_s(:, :) = ZERO
      LAMBDA_s_c(:, :) = ZERO
      LAMBDA_s(:, :) = ZERO
      K_g(:) = ZERO
      K_s(:, :) = ZERO
      DIF_g(:, :) = ZERO
      DIF_S(:, :, :) = ZERO 
!
      IF (ENERGY_EQ .OR. L_SCALE0/=ZERO) THEN 
         RECALC_VISC_G = .TRUE. 
      ELSE 
         RECALC_VISC_G = .FALSE. 
      ENDIF 
!
!
!  Set specified constant physical properties values
!

      DO IJK = ijkstart3, ijkend3
      
         IF (WALL_AT(IJK)) THEN 
            RO_G(IJK) = ZERO 
            MU_G(IJK) = ZERO 
            K_G(IJK) = ZERO 
            C_PG(IJK) = ZERO 
            MW_MIX_G(IJK) = ZERO 
         ELSE 
            IF (RO_G0 /= UNDEFINED) RO_G(IJK) = RO_G0 
            IF (MU_G0 /= UNDEFINED) MU_G(IJK) = MU_G0 
            IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0 
            IF (C_PG0 /= UNDEFINED) C_PG(IJK) = C_PG0 
            IF (MW_AVG /= UNDEFINED) MW_MIX_G(IJK) = MW_AVG 
         ENDIF 
         HOR_G(IJK) = ZERO 
	 
	 DO M = 1, MMAX
           HOR_S(IJK,M) = ZERO 
	 END DO
      END DO 
      
      IF (DIF_G0 /= UNDEFINED) THEN 
         N = 1 
         IF (NMAX(0) > 0) THEN 
            IJK = 1 
            IF (IJKMAX2 > 0) THEN 
               DIF_G(IJKSTART3:IJKEND3,:NMAX(0)) = DIF_G0 
               IJK = IJKMAX2 + 1 
            ENDIF 
            N = NMAX(0) + 1 
         ENDIF 
      ENDIF 
!
      DO M = 1, MMAX 
         DO IJK = ijkstart3, ijkend3
            IF (WALL_AT(IJK)) THEN 
               P_S(IJK,M) = ZERO 
               MU_S(IJK,M) = ZERO 
               LAMBDA_S(IJK,M) = ZERO 
               ALPHA_S(IJK,M) = ZERO 
               K_S(IJK,M) = ZERO 
               C_PS(IJK,M) = ZERO 
            ELSE 
               IF (MU_S0 /= UNDEFINED) THEN 
                  P_S(IJK,M) = ZERO 
                  MU_S(IJK,M) = MU_S0 
                  LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M) 
                  ALPHA_S(IJK,M) = ZERO 
               ENDIF 
               IF (K_S0 /= UNDEFINED) K_S(IJK,M) = K_S0 
               IF (C_PS0 /= UNDEFINED) C_PS(IJK,M) = C_PS0 
            ENDIF 
         END DO 
	 
         IF (DIF_S0 /= UNDEFINED) THEN 
            N = 1 
            IF (NMAX(M) > 0) THEN 
               IJK = 1 
               IF (IJKMAX2 > 0) THEN 
                  DIF_S(IJKSTART3:IJKEND3,M,:NMAX(M)) = DIF_S0 
                  IJK = IJKMAX2 + 1 
               ENDIF 
               N = NMAX(M) + 1 
            ENDIF 
         ENDIF 
      END DO 
      
      IF (RO_G0 == ZERO) THEN 
         M = 1 
         IF (MMAX > 0) THEN 
            IJK = 1 
            IF (IJKMAX2 > 0) THEN 
               F_GS(IJKSTART3:IJKEND3,:MMAX) = ZERO 
               IJK = IJKMAX2 + 1 
            ENDIF 
            M = MMAX + 1 
         ENDIF 
      ENDIF 
!

      RETURN  
      END SUBROUTINE SET_CONSTPROP 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 120 Replaced the index for initialization: e.g. in DIF_G :ijkmax2 --> IJKSTART3:IJKEND3
