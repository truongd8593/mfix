!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_W_s(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for W_s momentum eq. The terms      C
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
!  Revision Number:                                                    C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE SOURCE_W_S(A_M, B_M, IER) 
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
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
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
      INTEGER          I, J, K, IJK, IJKT, IMJK, IJMK, IJKM, IJKP, IMJKP 
      INTEGER          IJKE, IJKW, IJKTE, IJKTW, IM, IPJK 
! 
!                      Phase index 
      INTEGER          M,MM 
      DOUBLE PRECISION   SUM_EPS_CP 
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
      DOUBLE PRECISION ROPGA 
! 
!                      Average density differences 
      DOUBLE PRECISION dro1, dro2, droa 
! 
!                      Average quantities 
      DOUBLE PRECISION ugt, Cte, Ctw, EPMUoX 
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
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sdps, Sxzb, Vxza, Vxzb 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vmt, Vbf, Vcoa, Vcob 
! 
!                      error message 
      CHARACTER*80     LINE
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
! 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Z_EQ(M)) THEN 
!
!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
            IF (CALL_DI .or. CALL_ISAT) THEN
               SUM_R_S_temp = SUM_R_S
               SUM_R_S = ZERO
            END IF
!     CHEM & ISAT end (nan xie)
!
!
!$omp  parallel do &
!$omp& private(IJK, I, J, K, IJKT, EPGA, ISV, &
!$omp& PGT,SDP,SDPS,   ROPGA,V0,VMT,   DRO1,DRO2,DROA,VBF, &
!$omp& IMJK,IJKP,IMJKP,  UGT,VCOA,VCOB, &
!$omp& IJKE,IJKW,IJKTE,IJKTW,IM,IPJK, &
!$omp& CTE,CTW,SXZB,  EPMUOX,VXZA,VXZB )
            DO IJK = ijkstart3, ijkend3 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK)
               IMJK = IM_OF(IJK)
	       IJMK = JM_OF(IJK)
	       IJKM = KM_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               EPGA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K) 
               IF (IP_AT_T(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (SIP_AT_T(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_T(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 
!
!           dilute flow
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
! using the average boundary cell values to compute U_s (sof, Aug 23 2005)
!
                  IF (EP_S(WEST_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE 
                  IF (EP_S(EAST_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
                  IF (EP_S(SOUTH_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE 
                  IF (EP_S(NORTH_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF (EP_S(BOTTOM_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE 
                  IF (EP_S(TOP_OF(IJK),M) > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
!               
	          IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
	              A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
	             B_M(IJK,M) = -W_S(IJK,M)        
	          ELSE
	            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
	          ENDIF
!
               ELSE 
!
!           Surface forces
!
!             Pressure term
                  PGT = P_G(IJKT) 
                  IF (CYCLIC_Z_PD) THEN 
                     IF (CYCLIC_AT_T(IJK)) PGT = P_G(IJKT) - DELP_Z 
                  ENDIF 
                  IF (MODEL_B) THEN 
                     SDP = ZERO 
!
                  ELSE 
                     SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK) 
!
                  ENDIF 
!
                  IF (CLOSE_PACKED(M)) THEN 
		     SUM_EPS_CP=0.0 
		     DO MM=1,MMAX
		       IF (CLOSE_PACKED(MM))&
			     SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
		     END DO
		     SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                     SDPS = -((P_S(IJKT,M)-P_S(IJK,M))+(EP_S(IJK,M)/SUM_EPS_CP)*&
		     (P_STAR(IJKT)-P_STAR(IJK&
                        )))*AXY(IJK) 
                  ELSE 
                     SDPS = -(P_S(IJKT,M)-P_S(IJK,M))*AXY(IJK) 
                  ENDIF 
!
!             Volumetric forces
                  ROPGA = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K) 
!
!             Previous time step
                  V0 = AVG_Z(ROP_SO(IJK,M),ROP_SO(IJKT,M),K)*ODT 
!
!             Interphase mass transfer
                  VMT = AVG_Z(SUM_R_S(IJK,M),SUM_R_S(IJKT,M),K) 
!
!             Body force
                  IF (MODEL_B) THEN 
                     DRO1 = (RO_S(M)-RO_G(IJK))*EP_S(IJK,M) 
                     DRO2 = (RO_S(M)-RO_G(IJKT))*EP_S(IJKT,M) 
                     DROA = AVG_Z(DRO1,DRO2,K) 
!
                     VBF = DROA*BFZ_S(IJK,M) 
!
!
                  ELSE 
                     VBF = ROPGA*BFZ_S(IJK,M) 
!
                  ENDIF 
!
!
!             Special terms for cylindrical coordinates
               VCOA = ZERO 
               VCOB = ZERO 
               SXZB = ZERO 
               VXZA = ZERO 
               VXZB = ZERO 
	       CTE  = ZERO
	       CTW  = ZERO
                  IF (CYLINDRICAL) THEN 
!
!               Coriolis force
                     IMJK = IM_OF(IJK) 
                     IJKP = KP_OF(IJK) 
                     IMJKP = KP_OF(IMJK) 
                     UGT = AVG_Z(HALF*(U_S(IJK,M)+U_S(IMJK,M)),HALF*(U_S(IJKP,M&
                        )+U_S(IMJKP,M)),K) 
                     IF (UGT > ZERO) THEN 
                        VCOA = ROPGA*UGT*OX(I) 
                        VCOB = ZERO 
                     ELSE 
                        VCOA = ZERO 
                        VCOB = -ROPGA*UGT*W_S(IJK,M)*OX(I) 
                     ENDIF 
!
!           Term from tau_xz: intergral of (1/x)*(d/dx)(x*mu*(-w/x))
                     IJKE = EAST_OF(IJK) 
                     IJKW = WEST_OF(IJK) 
                     IJKTE = TOP_OF(IJKE) 
                     IJKTW = TOP_OF(IJKW) 
                     IM = IM1(I) 
                     IPJK = IP_OF(IJK) 
!
                     CTE = HALF*AVG_Z_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                        AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*OX_E(I)*AYZ_W(&
                        IJK) 
!
                     CTW = HALF*AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                        AVG_X_H(MU_S(IJKTW,M),MU_S(IJKT,M),IM),K)*DY(J)*(HALF*(&
                        DZ(K)+DZ(KP1(K)))) !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity
!
!           (mu/x)*(dw/dx) part of tau_xz/x
                     EPMUOX = AVG_Z(MU_S(IJK,M),MU_S(IJKT,M),K)*OX(I) 
                     VXZB = ZERO 
                     A_M(IJK,E,M) = A_M(IJK,E,M) + HALF*EPMUOX*ODX_E(I)*&
                           VOL_W(IJK) 
                     A_M(IJK,W,M) = A_M(IJK,W,M) - HALF*EPMUOX*ODX_E(IM)*&
                           VOL_W(IJK) 
!
!           -(mu/x)*(w/x) part of tau_xz/x
                     VXZA = EPMUOX*OX(I) 
                  ELSE 
                     VCOA = ZERO 
                     VCOB = ZERO 
                     SXZB = ZERO 
                     VXZA = ZERO 
                     VXZB = ZERO 
                  ENDIF 
!
!             Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(&
                     IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+ZMAX(VMT)+VCOA+VXZA&
                     )*VOL_W(IJK)+ CTE - CTW) 
                  A_M(IJK,E,M) = A_M(IJK,E,M) - CTE
                  A_M(IJK,W,M) = A_M(IJK,W,M) + CTW 
                  B_M(IJK,M) = -(SDP + SDPS + TAU_W_S(IJK,M)+SXZB+((V0+ZMAX((-&
                     VMT)))*W_SO(IJK,M)+VBF+VCOB+VXZB)*VOL_W(IJK)) + B_m(IJK, M)
               ENDIF 
            END DO 
            CALL SOURCE_W_S_BC (A_M, B_M, M, IER) 
!
!     CHEM & ISAT begin (nan xie)
            IF (CALL_DI .or. CALL_ISAT) THEN
               SUM_R_S = SUM_R_S_temp
            END IF
!     CHEM & ISAT end (nan xie)
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE SOURCE_W_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_W_s_BC(A_m, B_m, M, IER)                        C
!  Purpose: Determine source terms for W_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
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
      SUBROUTINE SOURCE_W_S_BC(A_M, B_M, M, IER) 
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
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
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
! 
!                      coefficient for granular bc 
      DOUBLE PRECISION hw 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
!
!  Set the default boundary conditions
!
      J1 = 1 
      DO K1 = kmin3,kmax3 
         DO I1 = imin3,imax3 
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
      I1 = 1 
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3 
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
!
               IF (BC_JJ_PS(L) == 0) THEN 
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
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
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
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE	    	    	 		     

                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
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
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,E,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_WW_S(L,M) 
                              ELSE 
                                 IF (CYLINDRICAL) THEN 
                                    A_M(IJK,0,M) = -(HALF*(BC_HW_S(L,M)+OX_E(I)&
                                       )+ODX_E(I)) 
                                    A_M(IJK,E,M) = -(HALF*(BC_HW_S(L,M)+OX_E(I)&
                                       )-ODX_E(I)) 
                                 ELSE 
                                    A_M(IJK,0,M)=-(HALF*BC_HW_S(L,M)+ODX_E(I)) 
                                    A_M(IJK,E,M)=-(HALF*BC_HW_S(L,M)-ODX_E(I)) 
                                 ENDIF 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,W,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_WW_S(L,M) 
                              ELSE 
                                 IF (CYLINDRICAL) THEN 
                                    A_M(IJK,W,M) = -(HALF*(BC_HW_S(L,M)-OX_E(IM&
                                       ))-ODX_E(IM)) 
                                    A_M(IJK,0,M) = -(HALF*(BC_HW_S(L,M)-OX_E(IM&
                                       ))+ODX_E(IM)) 
                                 ELSE 
                                    A_M(IJK,W,M) = -(HALF*BC_HW_S(L,M)-ODX_E(IM&
                                       )) 
                                    A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(IM&
                                       )) 
                                 ENDIF 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,N,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_WW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J)) 
                                 A_M(IJK,N,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,S,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_WW_S(L,M) 
                              ELSE 
                                 A_M(IJK,S,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
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
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE	    	    	 		     
		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -W_S(IJK,M) 
                        IF (BC_PLANE(L) == 'B') THEN 
                           IJKB = BOTTOM_OF(IJK) 
                           A_M(IJKB,E,M) = ZERO 
                           A_M(IJKB,W,M) = ZERO 
                           A_M(IJKB,N,M) = ZERO 
                           A_M(IJKB,S,M) = ZERO 
                           A_M(IJKB,T,M) = ZERO 
                           A_M(IJKB,B,M) = ZERO 
                           A_M(IJKB,0,M) = -ONE 
                           B_M(IJKB,M) = -W_S(IJKB,M) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_W_S_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_W_s(I1, I2, J1, J2, K1, K2, L, M, A_m, b_m)      C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal, A. Srivastava,                 Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!  Modified: S. Benyahia, Fluent Inc.                 Date: 02-FEB-05  C
!      Added the argument L to calc_grbdry                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE JJ_BC_W_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
      USE visc_s 
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
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
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                       JM, IM
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
!                      coefficient for granular bc
      DOUBLE PRECISION hw, gw, cw
      
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
	    
               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
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
!
               IF (FLUID_AT(EAST_OF(IJK))) THEN 
                  IF (EP_S(EAST_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,E,M) = -ONE 
                  ELSE 
! start anuj 04/20
                     IF (FRICTION .AND. (ONE-EP_G(EAST_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, EAST_OF(IJK), 'E', 'W', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, EAST_OF(IJK), 'E', 'W', M, L,&
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_WW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_WW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     IF (CYLINDRICAL) THEN 
                        A_M(IJK,E,M) = -(HALF*(HW + OX_E(I)*GW)-ODX_E(I)*GW) 
                        A_M(IJK,0,M) = -(HALF*(HW + OX_E(I)*GW)+ODX_E(I)*GW) 
                     ELSE 
                        A_M(IJK,E,M) = -(HALF*HW - ODX_E(I)*GW) 
                        A_M(IJK,0,M) = -(HALF*HW + ODX_E(I)*GW) 
                     ENDIF 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                  IF (EP_S(WEST_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,W,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(WEST_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, WEST_OF(IJK), 'W', 'W', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, WEST_OF(IJK), 'W', 'W', M, L,&
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_WW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_WW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     IF (CYLINDRICAL) THEN 
                        A_M(IJK,W,M) = -(HALF*(HW - OX_E(IM)*GW)-ODX_E(IM)*GW) 
                        A_M(IJK,0,M) = -(HALF*(HW - OX_E(IM)*GW)+ODX_E(IM)*GW) 
                     ELSE 
                        A_M(IJK,W,M) = -(HALF*HW - ODX_E(IM)*GW) 
                        A_M(IJK,0,M) = -(HALF*HW + ODX_E(IM)*GW) 
                     ENDIF 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,N,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(NORTH_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, NORTH_OF(IJK), 'N', 'W', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'W', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_WW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_WW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,N,M) = -(HALF*HW - ODY_N(J)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(J)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,S,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(SOUTH_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, SOUTH_OF(IJK), 'S', 'W', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'W', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_WW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_WW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,S,M) = -(HALF*HW - ODY_N(JM)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(JM)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE JJ_BC_W_S 


!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
