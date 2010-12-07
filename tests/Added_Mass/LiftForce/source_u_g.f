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
      INTEGER          I, J, K, IJK, IJKE, IPJK, IJKM, IPJKM, IMJK, IJMK, IPJMK, IJPK, IJKP
! 
!                      Phase index 
      INTEGER          M, L, IM
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
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop

!			Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag

!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Use, Usw, Vsw, Vse, Usn, Uss, Wsb, Wst, Wse, Usb, Ust
!			Lift force
      DOUBLE PRECISION F_Lift, Vgw, Vge, Vgc, U_slip, gradVg
! 
!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!      INTEGER :: J,K,IMJK 
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
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN  
!
!
!!     CHEM & ISAT begin (nan xie)
! Set the source terms zero
      IF (CALL_DI .or. CALL_ISAT) THEN
         SUM_R_G_temp = SUM_R_G
         SUM_R_G = ZERO
      END IF
!     CHEM & ISAT end (nan xie)
!
!$omp    parallel do private(I, IJK, IJKE, IJKM, IPJK, IPJKM,     &
!$omp&                  ISV, Sdp, V0, Vpm, Vmt, Vbf,              &
!$omp&                  Vcf, EPMUGA, VTZA, WGE, PGE, ROGA,        &
!$omp&                  MUGA, ROPGA, EPGA )
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
	 J = J_OF(IJK)
	 K = K_OF(IJK)
         IJKE = EAST_OF(IJK) 
         IJKM = KM_OF(IJK) 
         IPJK = IP_OF(IJK)  
         IMJK = IM_OF(IJK) 
         IPJKM = IP_OF(IJKM) 
	 IJMK = JM_OF(IJK)
	 IPJMK = IP_OF(IJMK)
	 IJPK = JP_OF(IJK)
	 IJKP = KP_OF(IJK)
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*(PGE - P_G(IJK))*AYZ(IJK) 
               ELSE
                   SDP = -P_SCALE*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
               ENDIF
! Original terms
!                   SDP = -P_SCALE*(PGE - P_G(IJK))*AYZ(IJK)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK) 
               ELSE
                   SDP = -P_SCALE*EPGA*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
               ENDIF
! Original terms
!                   SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            ENDIF 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
!       Volumetric forces
               ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
               ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK) 
!         Previous time step
               V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK)  
!         Added mass implicit transient term {Cv eps rop_g dU/dt}
               IF(Added_Mass) THEN
	         ROP_MA = AVG_X(ROP_g(IJK)*EP_s(IJK,M_AM),ROP_g(IJKE)*EP_s(IJKE,M_AM),I)
		 V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ELSE
!       Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IPJK)*ROP_G(IJKE))/(VOL(IJK) + VOL(IPJK))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + VOL(IPJK)*RO_G(IJKE) )/(VOL(IJK) + VOL(IPJK))
!         Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*ODT/(VOL(IJK) + VOL(IPJK))  
!         Added mass implicit transient term {Cv eps rop_g dU/dt}
               IF(Added_Mass) THEN
                 ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM) + VOL(IPJK)*ROP_g(IJKE)*EP_s(IJKE,M_AM) )/(VOL(IJK) + VOL(IPJK))
		 V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ENDIF
!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term dvg/dt - dVs/dt to virtual mass term			    
	    F_vir = ZERO
	    IF(Added_Mass.AND.(.NOT.CUT_U_TREATMENT_AT(IJK))) THEN        
	      F_vir = ( (U_s(IJK,M_AM) - U_sO(IJK,M_AM)) )*ODT*VOL_U(IJK)
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)    
	      Usw = AVG_X_E(U_S(IMJK,M_AM),U_s(IJK,M_AM),I)
	      Use = AVG_X_E(U_s(IJK,M_AM),U_s(IPJK,M_AM),IP1(I))
	      
	      Vsw = AVG_Y_N(V_s(IJMK,M_AM),V_s(IJK,M_AM))  
	      Vse = AVG_Y_N(V_s(IPJMK,M_AM),V_s(IPJK,M_AM)) 
	      Uss = AVG_Y(U_s(IJMK,M_AM),U_s(IJK,M_AM),JM1(J))
	      Usn = AVG_Y(U_s(IJK,M_AM),U_s(IJPK,M_AM),J)

	      IF(DO_K) THEN
	         Wsb = AVG_Z_T(W_s(IJKM,M_AM),W_s(IJK,M_AM))  
	         Wst = AVG_Z_T(W_s(IPJKM,M_AM),W_s(IPJK,M_AM)) 
		 Wse = AVG_X(Wsb,Wst,I)
	         Usb = AVG_Z(U_s(IJKM,M_AM),U_s(IJK,M_AM),KM1(K))
	         Ust = AVG_Z(U_s(IJK,M_AM),U_s(IJKP,M_AM),K)
	         F_vir = F_vir + Wse*OX_E(I) * (Ust - Usb) *AXY(IJK)
		 IF (CYLINDRICAL) F_vir = F_vir - Wse**2*OX_E(I) ! centrifugal force
	      ENDIF
!
! adding convective terms (U dU/dx + V dU/dy + W dU/dz) to virtual mass
	      F_vir = F_vir + U_s(IJK,M_AM)*(Use - Usw)*AYZ(IJK) + &
	         AVG_X(Vsw,Vse,I) * (Usn - Uss)*AXZ(IJK)
	         
	    
	      F_vir = F_vir * Cv * ROP_MA
!
! Lift Force implemented for a special case, i.e. this is not a general implementation.
	      Vgw = AVG_Y_N(V_g(IJMK),V_g(IJK))  
	      Vge = AVG_Y_N(V_g(IPJMK),V_g(IPJK)) 
	      Vgc = AVG_X(Vgw,Vge,I) 
	      
	      U_slip = (Vgc - AVG_X(Vsw,Vse,I))
	      gradVg = (Vge - Vgw)*AYZ(IJK)
	      
	      F_lift = U_slip * gradVg
	      F_lift = F_lift * 0.288d0*ROP_MA ! Lift coefficient Cl = 0.288
!
! end of Lift Force
	    ENDIF
!
!!! END VIRTUAL MASS SECTION

! Original terms
!       Volumetric forces
!            ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
!            ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK) 
!         Previous time step
!            V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK)  
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
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
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
               VMT = HALF * (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/VOL_U(IJK)  
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/(VOL(IJK) + VOL(IPJK))  
            ENDIF
! Original terms
!            VMT = HALF * (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/VOL_U(IJK)  
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFX_G(IJK) 
!
            ELSE                                 !Model A 
               VBF = ROPGA*BFX_G(IJK) 
!
            ENDIF 

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
		  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
		    DO L = 1,SMAX
		      avgRop = AVG_X(ROP_S(IJK,L),ROP_S(IJKE,L),I)
		      if(avgRop > ZERO) Ghd_drag = Ghd_drag +&
		           AVG_X(F_GS(IJK,L),F_GS(IJKE,L),I) * JoiX(IJK,L) / avgRop
		    ENDDO
		  ENDIF
! end of modifications for GHD theory

! Additional force for HYS drag force, do not use with mixture GHD theory
                avgDrag = ZERO
                  HYS_drag = ZERO
		  IF (TRIM(DRAG_TYPE) .EQ. 'HYS' .AND. TRIM(KT_TYPE) /= 'GHD') THEN
		     DO IM=1,MMAX
                        DO L = 1,MMAX
		           IF (L /= IM) THEN
		              avgDrag = AVG_X(beta_ij(IJK,IM,L),beta_ij(IJKE,IM,L),I)
		              HYS_drag = HYS_drag + avgDrag * (U_g(ijk) - U_s(IJK,L))
		           ENDIF
		        ENDDO
		     ENDDO
		  ENDIF
! end of modifications for HYS drag
!

!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           centrifugal force
               WGE = AVG_X(HALF*(W_G(IJK)+W_G(IJKM)),HALF*(W_G(IPJK)+W_G(IPJKM)&
                  ),I) 
               VCF = ROPGA*WGE**2*OX_E(I) 
	       IF(Added_Mass) VCF = VCF + Cv*ROP_MA*WGE**2*OX_E(I) ! virtual mass contribution.
!
!           -(2mu/x)*(u/x) part of Tau_zz/X
               EPMUGA = AVG_X(MU_GT(IJK),MU_GT(IJKE),I) 
               VTZA = 2.d0*EPMUGA*OX_E(I)*OX_E(I) 
            ELSE 
               VCF = ZERO 
               VTZA = ZERO 
            ENDIF 
!
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT)+VTZA)*VOL_U(IJK)) 
            B_M(IJK,M) = -(SDP + TAU_U_G(IJK)+((V0+ZMAX((-VMT)))*U_GO(IJK)+VBF+&
               VCF+Ghd_drag+HYS_drag)*VOL_U(IJK))+B_M(IJK,M) 

            B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit part of virtual mass force
            B_M(IJK,M) = B_M(IJK,M) + F_Lift ! explicit term of Lift force

	ENDIF 
      END DO 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G(A_M, B_M, IER)

      CALL SOURCE_U_G_BC (A_M, B_M, IER) 

      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G_BC(A_M, B_M, IER)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
!     CHEM & ISAT begin (nan xie)
!
      IF (CALL_DI .or. CALL_ISAT) THEN
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
