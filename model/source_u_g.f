!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_g(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
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
      SUBROUTINE SOURCE_U_G(A_M, B_M, IER) 
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
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
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
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I, IJK, IJKE, IPJK, IJKM, IPJKM 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPGA 
! 
!                      Average density 
      DOUBLE PRECISION ROPGA, ROGA 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Average viscosity 
      DOUBLE PRECISION MUGA 
! 
!                      Average viscosity 
      DOUBLE PRECISION EPMU_gte, EPMU_gbe, EPMUGA 
! 
!                      Average W_g 
      DOUBLE PRECISION Wge 
! 
!                      Average dW/Xdz 
      DOUBLE PRECISION dWoXdz 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcf, Vtza 
! 
!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_CHEM and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      M = 0 
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN  
!
!$omp    parallel do private(I, IJK, IJKE, IJKM, IPJK, IPJKM,     &
!$omp&                  ISV, Sdp, V0, Vpm, Vmt, Vbf,              &
!$omp&                  Vcf, EPMUGA, VTZA, WGE, PGE, ROGA,        &
!$omp&                  MUGA, ROPGA, EPGA )
!
!!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
      IF (CALL_CHEM .or. CALL_ISAT) THEN
         SUM_R_G_temp = SUM_R_G
         SUM_R_G = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
         IJKE = EAST_OF(IJK) 
         IJKM = KM_OF(IJK) 
         IPJK = IP_OF(IJK) 
         IPJKM = IP_OF(IJKM) 
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
         IF (IP_AT_E(IJK)) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
!
!       dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
!
            IF (EP_G(WEST_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,W,M) = ONE 
            ELSE IF (EP_G(EAST_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,E,M) = ONE 
            ELSE 
               B_M(IJK,M) = -U_G(IJK) 
            ENDIF 
         ELSE 
!
!       Surface forces
!
!         Pressure term
            PGE = P_G(IJKE) 
            IF (CYCLIC_X_PD) THEN 
               IF (IMAP(I_OF(IJK)).EQ.IMAX1) PGE = P_G(IJKE) - DELP_X 
            ENDIF 
            IF (MODEL_B) THEN 
               SDP = -P_SCALE*(PGE - P_G(IJK))*AYZ(IJK) 
!
            ELSE 
               SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK) 
!
            ENDIF 
!
!       Volumetric forces
            ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
            ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK) 
!
!         Previous time step
            V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK) 
!
!         pressure drop through porous media
            IF (SIP_AT_E(IJK)) THEN 
               ISV = IS_ID_AT_E(IJK) 
               MUGA = AVG_X(MU_G(IJK),MU_G(IJKE),I) 
               VPM = MUGA/IS_PC(ISV,1) 
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + HALF*IS_PC(ISV,2)*ROPGA*ABS(&
                  U_G(IJK)) 
            ELSE 
               VPM = ZERO 
            ENDIF 
!
!         Interphase mass transfer
            VMT = HALF * (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/VOL_U(IJK) 
!
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFX_G(IJK) 
!
            ELSE                                 !Model A 
               VBF = ROPGA*BFX_G(IJK) 
!
            ENDIF 
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           centrifugal force
               WGE = AVG_X(HALF*(W_G(IJK)+W_G(IJKM)),HALF*(W_G(IPJK)+W_G(IPJKM)&
                  ),I) 
               VCF = ROPGA*WGE**2*OX_E(I) 
!
!           -(2mu/x)*(u/x) part of Tau_zz/X
               EPMUGA = AVG_X(MU_GT(IJK),MU_GT(IJKE),I) 
               VTZA = 2.*EPMUGA*OX_E(I)*OX_E(I) 
            ELSE 
               VCF = ZERO 
               VTZA = ZERO 
            ENDIF 
!
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT)+VTZA)*VOL_U(IJK)) 
            B_M(IJK,M) = -(SDP + TAU_U_G(IJK)+((V0+ZMAX((-VMT)))*U_GO(IJK)+VBF+&
               VCF)*VOL_U(IJK))+B_M(IJK,M)
	ENDIF 
      END DO 
      CALL SOURCE_U_G_BC (A_M, B_M, IER) 
!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_CHEM .or. CALL_ISAT) THEN
         SUM_R_G = SUM_R_G_temp
      END IF  
!     CHEM & ISAT end (nan xie)
!

      RETURN  
      END SUBROUTINE SOURCE_U_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_g_BC(A_m, B_m, IER)                           C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
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
      SUBROUTINE SOURCE_U_G_BC(A_M, B_M, IER) 
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
      USE constant
      USE physprop
      USE fldvar
      USE visc_g 
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_g 
      USE bc
      USE output
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
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IP, IPJK 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Turbulent shear stress
      DOUBLE PRECISION  W_F_Slip
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      M = 0 
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
               IJK = FUNIJK(I1,J1,K1) 	       
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
	 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = -ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
!
      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
  	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = -ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      J1 = JMAX2 
      DO K1 = kmin3, kmax3      
         DO I1 = imin3, imax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = -ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = -ONE 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           A_M(IJK,T,M) = -ONE 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           A_M(IJK,B,M) = -ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
			   IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(L)-ODY_N(J)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,S,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,S,M) = -(HALF*BC_HW_G(L)-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN  
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,T,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,0,M)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I)) 
                              A_M(IJK,T,M)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN   
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,B,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
                                 )) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
                                 )) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
! wall functions for U-momentum are specify in this section of the code
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL'   .OR.  &
	             BC_TYPE(L) == 'NO_SLIP_WALL'    .OR.  &
		     BC_TYPE(L) == 'FREE_SLIP_WALL'  .AND. &
		     K_Epsilon                            )THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN  
			     CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                             A_M(IJK,N,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
			     CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                             A_M(IJK,S,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
			     CALL Wall_Function(IJK,TOP_OF(IJK),ODZ_T(K)*OX_E(I),W_F_Slip)
                             A_M(IJK,T,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN
			     CALL Wall_Function(IJK,BOTTOM_OF(IJK),ODZ_T(KM)*OX_E(I),W_F_Slip)
                             A_M(IJK,B,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO
! end of wall functions
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IM = IM1(I) 
                           IMJK = IM_OF(IJK) 
                           A_M(IMJK,E,M) = ZERO 
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM)) 
                           A_M(IMJK,N,M) = ZERO 
                           A_M(IMJK,S,M) = ZERO 
                           A_M(IMJK,T,M) = ZERO 
                           A_M(IMJK,B,M) = ZERO 
                           A_M(IMJK,0,M) = -ONE 
                           B_M(IMJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ELSE IF (BC_PLANE(L) == 'E') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
!
                           IP = IP1(I) 
                           IPJK = IP_OF(IJK) 
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I) 
                           A_M(IPJK,W,M) = ZERO 
                           A_M(IPJK,N,M) = ZERO 
                           A_M(IPJK,S,M) = ZERO 
                           A_M(IPJK,T,M) = ZERO 
                           A_M(IPJK,B,M) = ZERO 
                           A_M(IPJK,0,M) = -ONE 
                           B_M(IPJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -U_G(IJK) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_G(IJKW) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_U_G_BC 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Wall_Function(IJK1,IJK2,ODX_WF,W_F_Slip)               C
!  Purpose: Calculate Slip velocity using wall functions               C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE Wall_Function(IJK1,IJK2,ODX_WF,W_F_Slip)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop 
      USE fldvar
      USE visc_g  
      USE geometry 
      USE indices 
      USE bc
      USE compar 
      USE turb        
      USE mpi_utility 
      IMPLICIT NONE
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!                      IJK indices for wall cell and fluid cell
      INTEGER          IJK1, IJK2

!                      ODX_WF: 1/dx, and W_F_Slip: value of turb. shear stress at walls
      DOUBLE PRECISION ODX_WF, W_F_Slip

!                      C_mu and Kappa are constants in turb. viscosity and Von Karmen const.
      DOUBLE PRECISION C_mu, Kappa
!-----------------------------------------------
!
!
	C_mu = 0.09D+0
	Kappa = 0.42D+0
	
		W_F_Slip = (ONE - ONE/ODX_WF* RO_g(IJK2)*C_mu**0.25   &
			   *SQRT(K_Turb_G(IJK2))/MU_gT(IJK2)	      &
			   *Kappa/LOG(9.81D+0/(ODX_WF*2.D+0)*         &
			    RO_g(IJK2)*C_mu**0.25*                    &
			   SQRT(K_Turb_G(IJK2))/MU_g(IJK2)))
!
      RETURN  
      END SUBROUTINE Wall_Function

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
