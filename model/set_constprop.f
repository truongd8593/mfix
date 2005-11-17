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
      INTEGER :: IJK, M, N, I, J
      DOUBLE PRECISION SUM, SUM_EP, old_value, DP_TMP(MMAX)
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
      F_GS(:,:) = ZERO                 !S. Dartevelle, LANL, 2004
      F_SS(:,:) = ZERO                 !S. Dartevelle, LANL, 2004
!
      IF (ENERGY_EQ .OR. L_SCALE0/=ZERO .OR. K_EPSILON) THEN 
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
            IF (C_PG0 /= UNDEFINED) C_PG(IJK) = C_PG0 
            IF (MW_AVG /= UNDEFINED) MW_MIX_G(IJK) = MW_AVG 
	    IF(FLUID_AT(IJK)) THEN
              IF (MU_G0 /= UNDEFINED) MU_G(IJK) = MU_G0 
              IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0 
              IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = DIF_G0 
	    ENDIF
         ENDIF 
      END DO 
!
      DO M = 1, MMAX 
         DO IJK = ijkstart3, ijkend3
            IF (WALL_AT(IJK)) THEN 
               P_S(IJK,M) = ZERO 
               MU_S(IJK,M) = ZERO 
! add by rong
               D_p(IJK,M)=D_P0(M)
             
! add by rong
               LAMBDA_S(IJK,M) = ZERO 
               ALPHA_S(IJK,M) = ZERO 
               K_S(IJK,M) = ZERO 
               C_PS(IJK,M) = ZERO 
            ELSE 
               IF (C_PS0 /= UNDEFINED) C_PS(IJK,M) = C_PS0 
	       IF(FLUID_AT(IJK))THEN
                 IF (MU_S0 /= UNDEFINED) THEN 
                    P_S(IJK,M) = ZERO 
                    MU_S(IJK,M) = MU_S0 
                    LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M) 
                    ALPHA_S(IJK,M) = ZERO 
                 ENDIF 
                 IF (K_S0 /= UNDEFINED) K_S(IJK,M) = K_S0 
                 IF (DIF_S0 /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = DIF_S0
	       ENDIF
! add by rong; modified by sof (June 17 2005)
                 IF (D_P0(M)/=UNDEFINED)  D_P(IJK,M)=D_P0(M)
! add by rong 
            ENDIF 
! set ep_star_array to user input ep_star in all cells. sof--> Nov-17-05
            EP_star_array(ijk) = ep_star
	    IF(EP_S_MAX(M) == UNDEFINED) EP_S_MAX(M) = ONE-EP_STAR
	 END DO 
	 
      END DO 
      
      IF (RO_G0 == ZERO) THEN 
         IF (MMAX > 0) THEN 
            IF (IJKMAX2 > 0) THEN 
               F_GS(IJKSTART3:IJKEND3,:MMAX) = ZERO 
            ENDIF 
         ENDIF 
      ENDIF 
!
!
! start sof modifications: 05/04-2005
!
! initializing the new indexing system
! this doesn't need to be done if no correlation is used to compute ep_star
!
      IF(YU_STANDISH .OR. FEDORS_LANDEL) THEN
	IF(.NOT. CALL_DQMOM) THEN
         DO I = 1, MMAX
           IF(D_P0(I) == UNDEFINED) RETURN
	   DP_TMP(I) = D_P0(I)
           M_MAX(I) = I
         END DO
!
! rearrange the indices from coarsest particles to finest to be used in CALC_ep_star
! I did this here because it may need to be done for auto_restart
         DO I = 1, MMAX	 
	   DO J = I, MMAX
	    
	     IF(DP_TMP(I) < DP_TMP(J)) THEN
	       old_value = DP_TMP(I)
	       DP_TMP(I) = DP_TMP(J)
	       DP_TMP(J) = old_value
	     ENDIF
	   
	   END DO
         END DO
!
	 DO I = 1, MMAX	 
	   DO J = 1, MMAX
	     
	     IF(DP_TMP(I) == D_P0(J) .AND. D_P0(I) .NE. D_P0(J)) THEN
	       M_MAX(I) = J 
	     ENDIF
	   
	   END DO
         END DO
        ENDIF ! for .not. call_dqmom
      ENDIF ! for Yu-standish or Fedors-Landel
!
! end of sof modifications: 05/04-2005
!

      RETURN  
      END SUBROUTINE SET_CONSTPROP 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 120 Replaced the index for initialization: e.g. in DIF_G :ijkmax2 --> IJKSTART3:IJKEND3
