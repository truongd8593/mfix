!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_W_g(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
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
      SUBROUTINE SOURCE_W_G(A_M, B_M, IER) 
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
      USE compar        !//d
      USE sendrecv      !// 400      
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
      INTEGER          I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,& 
                       IJKE, IJKW, IJKTE, IJKTW, IM, IPJK 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at top cell 
      DOUBLE PRECISION PgT 
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
!                      Average coefficients 
      DOUBLE PRECISION Cte, Ctw, EPMUoX 
! 
!                      Average U_g 
      DOUBLE PRECISION Ugt 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sxzb 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcoa, Vcob, Vxza, Vxzb 
! 
!                      error message 
      CHARACTER*80     LINE 
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
      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN  
!
!// 350 1223 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

!$omp  parallel do private( I, J, K, IJK, IJKT, ISV, Sdp, V0, Vpm, Vmt, Vbf, &
!$omp&  PGT, ROGA, IMJK, IJKP, IMJKP, IJKW, IJKTE, IJKTW, IM, IPJK,  &
!$omp&  CTE, CTW, SXZB, EPMUOX, VXZA, VXZB, UGT, VCOA, VCOB, IJKE,&
!$omp&  MUGA, ROPGA, EPGA, LINE)  
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K) 
         IF (IP_AT_T(IJK)) THEN 
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
            IF (EP_G(BOTTOM_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,B,M) = ONE 
            ELSE IF (EP_G(TOP_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,T,M) = ONE 
            ELSE 
               B_M(IJK,M) = -W_G(IJK) 
            ENDIF 
         ELSE 
!
!       Surface forces
!
!         Pressure term
            PGT = P_G(IJKT) 
            IF (CYCLIC_Z_PD) THEN 
               IF (CYCLIC_AT_T(IJK)) PGT = P_G(IJKT) - DELP_Z 
            ENDIF 
            IF (MODEL_B) THEN 
               SDP = -P_SCALE*(PGT - P_G(IJK))*AXY(IJK) 
!
            ELSE 
               SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK) 
!
            ENDIF 
!
!       Volumetric forces
            ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K) 
            ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K) 
!
!         Previous time step
            V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT 
!
!         pressure drop through porous media
            IF (SIP_AT_T(IJK)) THEN 
               ISV = IS_ID_AT_T(IJK) 
               MUGA = AVG_Z(MU_G(IJK),MU_G(IJKT),K) 
               VPM = MUGA/IS_PC(ISV,1) 
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + HALF*IS_PC(ISV,2)*ROPGA*ABS(&
                  W_G(IJK)) 
            ELSE 
               VPM = ZERO 
            ENDIF 
!
!         Interphase mass transfer
            VMT = AVG_Z(SUM_R_G(IJK),SUM_R_G(IJKT),K) 
!
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFZ_G(IJK) 
!
            ELSE 
               VBF = ROPGA*BFZ_G(IJK) 
!
            ENDIF 
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           Coriolis force
               IMJK = IM_OF(IJK) 
               IJKP = KP_OF(IJK) 
               IMJKP = KP_OF(IMJK) 
               UGT = AVG_Z(HALF*(U_G(IJK)+U_G(IMJK)),HALF*(U_G(IJKP)+U_G(IMJKP)&
                  ),K) 
               IF (UGT > ZERO) THEN 
                  VCOA = ROPGA*UGT*OX(I) 
                  VCOB = ZERO 
               ELSE 
                  VCOA = ZERO 
                  VCOB = -ROPGA*UGT*W_G(IJK)*OX(I) 
               ENDIF 
!
!           Term from tau_xz
               IJKE = EAST_OF(IJK) 
               IJKW = WEST_OF(IJK) 
               IJKTE = TOP_OF(IJKE) 
               IJKTW = TOP_OF(IJKW) 
               IM = IM1(I) 
               IPJK = IP_OF(IJK) 
!
               CTE = HALF*AVG_Z_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(&
                  MU_GT(IJKT),MU_GT(IJKTE),I),K)*OX_E(I)*AYZ_W(IJK) 
!
               CTW = HALF*AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(&
                  MU_GT(IJKTW),MU_GT(IJKT),IM),K)*DY(J)*(HALF*(DZ(K)+DZ(KP1(K))&
                  )) 
!
!
               IF (CTE > CTW) THEN 
                  SXZB = ZERO 
                  A_M(IJK,E,M) = A_M(IJK,E,M) + CTE 
                  A_M(IJK,W,M) = A_M(IJK,W,M) - CTW 
               ELSE 
                  SXZB = W_G(IJK)*(CTW - CTE) + W_G(IPJK)*CTE - W_G(IMJK)*CTW 
               ENDIF 
!
!           part of tau_xz/x
               EPMUOX = AVG_Z(MU_GT(IJK),MU_GT(IJKT),K)*OX(I) 
               VXZA = EPMUOX*OX(I) 
               IF (ODX_E(I) > ODX_E(IM)) THEN 
                  VXZB = ZERO 
                  A_M(IJK,E,M) = A_M(IJK,E,M) + HALF*EPMUOX*ODX_E(I)*VOL_W(IJK) 
                  A_M(IJK,W,M)=A_M(IJK,W,M)-HALF*EPMUOX*ODX_E(IM)*VOL_W(IJK) 
               ELSE 
                  VXZB = W_G(IJK)*HALF*EPMUOX*(ODX_E(IM)-ODX_E(I)) + W_G(IPJK)*&
                     HALF*EPMUOX*ODX_E(I) - W_G(IMJK)*HALF*EPMUOX*ODX_E(IM) 
               ENDIF 
            ELSE 
               VCOA = ZERO 
               VCOB = ZERO 
               SXZB = ZERO 
               VXZA = ZERO 
               VXZB = ZERO 
            ENDIF 
!
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT)+VCOA+VXZA)*VOL_W(&
               IJK)) 
            B_M(IJK,M) = -(SDP + TAU_W_G(IJK)+SXZB+((V0+ZMAX((-VMT)))*W_GO(IJK)&
               +VBF+VCOB+VXZB)*VOL_W(IJK)) + B_M(IJK, M) 
         ENDIF 
      END DO 
      CALL SOURCE_W_G_BC (A_M, B_M, IER) 
!
!//? probably need to communicate A_M and B_M here or in solve_vel_star in order
!//? collect the COMMs
!!!      CALL SEND_RECV(A_M, 2)
!!!      CALL SEND_RECV(B_M, 2)

      RETURN  
      END SUBROUTINE SOURCE_W_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_W_g_BC(A_m, B_m, IER)                           C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
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
      SUBROUTINE SOURCE_W_G_BC(A_M, B_M, IER) 
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
      USE compar        !//d
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
      INTEGER          I,  J, K, KM, I1, I2, J1, J2, K1, K2, IJK,& 
                       IM, JM, IJKB, IJKM, IJKP 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
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
      J1 = 1 
      DO K1 = kmin3,kmax3 
         DO I1 = imin3,imax3 
!// 360 1224 Check if current i,j,k resides on this PE
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
!// 360 1223 Check if current i,j,k resides on this PE
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
      I1 = 1 
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3 
!// 360 1223 Check if current i,j,k resides on this PE
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = -ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      I1 = IMAX2 
      DO K1 = kmin3,kmax3 
         DO J1 = jmin3,jmax3 
!// 360 1223 Check if current i,j,k resides on this PE
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = -ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
!// 360 1223 Check if current i,j,k resides on this PE
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE	    	    	 		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           A_M(IJK,E,M) = -ONE 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           A_M(IJK,W,M) = -ONE 
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = -ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = -ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2
!// 360 1223 Check if current i,j,k resides on this PE
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		      
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
!// 360 1223 Check if current i,j,k resides on this PE
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                        IJK = FUNIJK(I,J,K) 
                        IM = IM1(I) 
                        JM = JM1(J) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,E,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              IF (CYLINDRICAL) THEN 
                                 A_M(IJK,0,M) = -(HALF*(BC_HW_G(L)+OX_E(I))+&
                                    ODX_E(I)) 
                                 A_M(IJK,E,M) = -(HALF*(BC_HW_G(L)+OX_E(I))-&
                                    ODX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(I)) 
                                 A_M(IJK,E,M) = -(HALF*BC_HW_G(L)-ODX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ENDIF 
                           ENDIF 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,W,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              IF (CYLINDRICAL) THEN 
                                 A_M(IJK,W,M) = -(HALF*(BC_HW_G(L)-OX_E(IM))-&
                                    ODX_E(IM)) 
                                 A_M(IJK,0,M) = -(HALF*(BC_HW_G(L)-OX_E(IM))+&
                                    ODX_E(IM)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ELSE 
                                 A_M(IJK,W,M) = -(HALF*BC_HW_G(L)-ODX_E(IM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(IM)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ENDIF 
                           ENDIF 
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(L)-ODY_N(J)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,S,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              A_M(IJK,S,M) = -(HALF*BC_HW_G(L)-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                           ENDIF 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'B') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2
!// 360 1223 Check if current i,j,k resides on this PE
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			 
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ONE 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'B') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2
!// 360 1223 Check if current i,j,k resides on this PE
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			 
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ONE 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IJKM = KM_OF(IJK) 
                           A_M(IJKM,E,M) = ZERO 
                           A_M(IJKM,W,M) = ZERO 
                           A_M(IJKM,N,M) = ZERO 
                           A_M(IJKM,S,M) = ZERO 
                           A_M(IJKM,T,M) = ZERO 
                           A_M(IJKM,B,M) = ONE 
                           A_M(IJKM,0,M) = -ONE 
                           B_M(IJKM,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ELSE IF (BC_PLANE(L) == 'T') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1223 Check if current i,j,k resides on this PE
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			
                           IJK = FUNIJK(I,J,K) 
!
                           IJKP = KP_OF(IJK) 
                           A_M(IJKP,E,M) = ZERO 
                           A_M(IJKP,W,M) = ZERO 
                           A_M(IJKP,N,M) = ZERO 
                           A_M(IJKP,S,M) = ZERO 
                           A_M(IJKP,T,M) = ONE 
                           A_M(IJKP,B,M) = ZERO 
                           A_M(IJKP,0,M) = -ONE 
                           B_M(IJKP,M) = ZERO 
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
!// 360 1223 Check if current i,j,k resides on this PE
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -W_G(IJK) 
                        IF (BC_PLANE(L) == 'B') THEN 
                           IJKB = BOTTOM_OF(IJK) 
                           A_M(IJKB,E,M) = ZERO 
                           A_M(IJKB,W,M) = ZERO 
                           A_M(IJKB,N,M) = ZERO 
                           A_M(IJKB,S,M) = ZERO 
                           A_M(IJKB,T,M) = ZERO 
                           A_M(IJKB,B,M) = ZERO 
                           A_M(IJKB,0,M) = -ONE 
                           B_M(IJKB,M) = -W_G(IJKB) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_W_G_BC 
