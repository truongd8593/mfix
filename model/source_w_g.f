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
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
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
      USE compar  
      USE sendrecv  
      USE ghdtheory
      USE drag  
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
      USE quadric
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

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
      INTEGER          M, L, MM
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
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop

!			Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag

! 
!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: JM,IP,JP,IJMK,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKM,KM,KP
      INTEGER :: IJKTN,IJKWT,IJKST
      DOUBLE PRECISION :: We,Ww,Wn,Ws,Wt,Wb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_GT_E,MU_GT_W,MU_GT_N,MU_GT_S,MU_GT_T,MU_GT_B,MU_GT_CUT
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

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
!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_G_temp = SUM_R_G
         SUM_R_G = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
! 
!
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
               IF (KMAP(K_OF(IJK)).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z 
            ENDIF 
            IF (MODEL_B) THEN 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*(PGT - P_G(IJK))*AXY(IJK) 
               ELSE
                   SDP = -P_SCALE*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
! Original terms
!               SDP = -P_SCALE*(PGT - P_G(IJK))*AXY(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK)
               ELSE
                   SDP = -P_SCALE*EPGA*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
! Original terms
!               SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
            ENDIF 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
!       Volumetric forces
               ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K) 
               ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K) 
!         Previous time step
               V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT 
            ELSE
!       Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IJKT)*ROP_G(IJKT))/(VOL(IJK) + VOL(IJKT))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + VOL(IPJK)*RO_G(IJKT) )/(VOL(IJK) + VOL(IJKT))
!         Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKT)*ROP_GO(IJKT))*ODT/(VOL(IJK) + VOL(IJKT))  
            ENDIF


! Original terms
!       Volumetric forces
!            ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K) 
!            ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K) 
!         Previous time step
!            V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
               VMT = AVG_Z(SUM_R_G(IJK),SUM_R_G(IJKT),K) 
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IJKT)*SUM_R_G(IJKT))/(VOL(IJK) + VOL(IJKT))  
            ENDIF
! Original terms
!            VMT = AVG_Z(SUM_R_G(IJK),SUM_R_G(IJKT),K) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFZ_G(IJK) 
!
            ELSE 
               VBF = ROPGA*BFZ_G(IJK) 
!
            ENDIF  

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
		  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
		    DO L = 1,SMAX
		      avgRop = AVG_Z(ROP_S(IJK,L),ROP_S(IJKT,L),K)
		      if(avgRop > ZERO) Ghd_drag = Ghd_drag +&
		           AVG_Z(F_GS(IJK,L),F_GS(IJKT,L),K) * JoiZ(IJK,L) / avgRop
		    ENDDO
		  ENDIF
! end of modifications for GHD theory

! Additional force for HYS drag force
		 avgDrag = ZERO
                 HYS_drag = ZERO
		 IF (TRIM(DRAG_TYPE) .EQ. 'HYS') THEN
		     DO MM=1,MMAX
		        DO L = 1,MMAX
		           IF (L /= MM) THEN
		              avgDrag = AVG_Z(beta_ij(IJK,MM,L),beta_ij(IJKT,MM,L),K)
		              HYS_drag = HYS_drag + avgDrag * (W_g(ijk) - W_s(IJK,L))
		           ENDIF
		        ENDDO
		   ENDDO
		 ENDIF
! end of modifications for HYS drag
!
!         Special terms for cylindrical coordinates
               VCOA = ZERO 
               VCOB = ZERO 
               SXZB = ZERO 
               VXZA = ZERO 
               VXZB = ZERO 
	       CTE  = ZERO
	       CTW  = ZERO
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
!           Term from tau_xz: intergral of (1/x)*(d/dx)(x*mu*(-w/x))
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
                  )) !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity
!
!           (mu/x)*(dw/dx) part of tau_xz/x
               EPMUOX = AVG_Z(MU_GT(IJK),MU_GT(IJKT),K)*OX(I) 
               VXZB = ZERO 
               A_M(IJK,E,M) = A_M(IJK,E,M) + HALF*EPMUOX*ODX_E(I)*VOL_W(IJK) 
               A_M(IJK,W,M) = A_M(IJK,W,M) - HALF*EPMUOX*ODX_E(IM)*VOL_W(IJK) 
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
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT)+VCOA + VXZA)*VOL_W(&
               IJK) + CTE - CTW) 
            A_M(IJK,E,M) = A_M(IJK,E,M) - CTE
            A_M(IJK,W,M) = A_M(IJK,W,M) + CTW 
	       
            B_M(IJK,M) = -(SDP + TAU_W_G(IJK)+SXZB+((V0+ZMAX((-VMT)))*W_GO(IJK)&
               +VBF+VCOB+VXZB+Ghd_drag+HYS_drag)*VOL_W(IJK)) + B_M(IJK, M) 
         ENDIF 
      END DO 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G(A_M, B_M, IER)

      CALL SOURCE_W_G_BC (A_M, B_M, IER) 

      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G_BC(A_M, B_M, IER)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
!     CHEM & ISAT begin (nan xie)
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_G = SUM_R_G_temp
      END IF  
!     CHEM & ISAT end (nan xie)
!
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
!                      Turb. Shear at walls
      DOUBLE PRECISION W_F_Slip 
!
!                      C_mu and Kappa are constants in turb. viscosity and Von Karmen const.
      DOUBLE PRECISION C_mu, Kappa
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
      C_mu = 0.09D+0
      Kappa = 0.42D+0
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
                                 A_M(IJK,0,M) = -(HALF*(BC_HW_G(L)-OX_E(I))+&
                                    ODX_E(I)) 
                                 A_M(IJK,E,M) = -(HALF*(BC_HW_G(L)-OX_E(I))-&
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
! wall functions for V-momentum are specify in this section of the code
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
			   IF (CYLINDRICAL) THEN
			       W_F_Slip = (ONE/(ODX_E(I)+HALF*OX_E(I)))*	    &
			       (ODX_E(I) -OX_E(I)-				    &
			       RO_g(EAST_OF(IJK))*C_mu**0.25			    &
			       *SQRT(K_Turb_G((EAST_OF(IJK))))/MU_gT(EAST_OF(IJK))  &
			       *Kappa/LOG(9.81D+0	/ODX_E(I)		    &
			       /(2.D+0)*RO_g(EAST_OF(IJK))*C_mu**0.25*		    &
			       SQRT(K_Turb_G((EAST_OF(IJK))))/MU_g(EAST_OF(IJK))))
			   ELSE
			       CALL Wall_Function(IJK,EAST_OF(IJK),ODX_E(I),W_F_Slip)
			   ENDIF
                           A_M(IJK,E,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN
			   IF (CYLINDRICAL) THEN
			       W_F_Slip =  (ONE/(ONE*ODX_E(IM)+HALF*OX_E(IM)))*	    &
			       (ONE*ODX_E(IM) -OX_E(IM)-			    &
			       RO_g(WEST_OF(IJK))*C_mu**0.25			    &
			       *SQRT(K_Turb_G((WEST_OF(IJK))))/MU_gT(WEST_OF(IJK))  &
			       *Kappa/LOG(9.81D+0	/ODX_E(IM)		    &
			       /(2.D+0)*RO_g(WEST_OF(IJK))*C_mu**0.25*		    &
			       SQRT(K_Turb_G((WEST_OF(IJK))))/MU_g(WEST_OF(IJK))))
			   ELSE
			       CALL Wall_Function(IJK,WEST_OF(IJK),ODX_E(IM),W_F_Slip)
			   ENDIF
                           A_M(IJK,W,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN
			   CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                           A_M(IJK,N,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
			   CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                           A_M(IJK,S,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
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

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
