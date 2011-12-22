!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_g(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for V_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96  C
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
      SUBROUTINE SOURCE_V_G(A_M, B_M, IER) 
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
      USE vshear
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
      INTEGER          I, J, K, IJK, IJKN
! 
!                      Phase index 
      INTEGER          M, L, IM
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at north cell 
      DOUBLE PRECISION PgN 
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
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf 
!
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop

!			Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag

!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Vsn, Vss, U_se, Usw, Vse, Vsw, Wst, Wsb, Vst, Vsb

! 
! loezos 
      DOUBLE PRECISION VSH_n,VSH_s,VSH_e,VSH_w,VSH_p,Source_conv
      DOUBLE PRECISION SRT
! loezos

!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER ::          JM,IP,JP,KM,KP
      INTEGER ::          IMJK,IPJK,IJMK,IJPK,IJKP,IJKM,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK
      INTEGER ::          IJKT,IJKTN,IJKB,IJKBN, IJPKM
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
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

!
      M = 0 
      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN  
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
!!!$omp  parallel do private( I, J, K, IJK, IJKN, ISV, Sdp, V0, Vpm, Vmt, Vbf, &
!!!$omp&  PGN, ROGA, MUGA, ROPGA, EPGA,VSH_n,VSH_s,VSH_e,VSH_w,&
!!!$omp&  VSH_p,Source_conv, SRT ) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 
         IJKN = NORTH_OF(IJK)  
         IMJK = IM_OF(IJK) 
         IPJK = IP_OF(IJK)   
         IJMK = JM_OF(IJK) 
         IJPK = JP_OF(IJK) 
	 IMJPK = IM_OF(IJPK)
	 IJKM = KM_OF(IJK) 
	 IJPKM = KM_OF(IJPK) 
	 IJKP = KP_OF(IJK) 
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
         IF (IP_AT_N(IJK)) THEN 
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
            IF (EP_G(SOUTH_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,S,M) = ONE 
            ELSE IF (EP_G(NORTH_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,N,M) = ONE 
            ELSE 
               B_M(IJK,M) = -V_G(IJK) 
            ENDIF 
!
         ELSEIF (BLOCKED_V_CELL_AT(IJK)) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
         ELSE
!
!       Surface forces
!
!         Pressure term
            PGN = P_G(IJKN) 
            IF (CYCLIC_Y_PD) THEN 
               IF (JMAP(J_OF(IJK)).EQ.JMAX1)PGN = P_G(IJKN) - DELP_Y 
            ENDIF 
            IF (MODEL_B) THEN 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                  SDP = -P_SCALE*(PGN - P_G(IJK))*AXZ(IJK) 
               ELSE
                  SDP = -P_SCALE*(PGN * A_VPG_N(IJK)  - P_G(IJK) * A_VPG_S(IJK) )
               ENDIF
! Original terms
!                  SDP = -P_SCALE*(PGN - P_G(IJK))*AXZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
               IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                  SDP = -P_SCALE*EPGA*(PGN - P_G(IJK))*AXZ(IJK) 
               ELSE
                  SDP = -P_SCALE*EPGA*(PGN * A_VPG_N(IJK)  - P_G(IJK) * A_VPG_S(IJK) )
               ENDIF

! Original terms
!                  SDP = -P_SCALE*EPGA*(PGN - P_G(IJK))*AXZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!
            ENDIF 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
!       Volumetric forces
               ROPGA = AVG_Y(ROP_G(IJK),ROP_G(IJKN),J) 
               ROGA = AVG_Y(RO_G(IJK),RO_G(IJKN),J) 
!         Previous time step
               V0 = AVG_Y(ROP_GO(IJK),ROP_GO(IJKN),J)*ODT 
!         Added mass implicit transient term {Cv eps rop_g dV/dt}
	       IF(Added_Mass) THEN
                 ROP_MA = AVG_Y(ROP_g(IJK)*EP_s(IJK,M_AM),ROP_g(IJKN)*EP_s(IJKN,M_AM),J)
	         V0 = V0 + Cv * ROP_MA * ODT
	       ENDIF
            ELSE
!       Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IJKN)*ROP_G(IJKN))/(VOL(IJK) + VOL(IJKN))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + VOL(IJKN)*RO_G(IJKN) )/(VOL(IJK) + VOL(IJKN))
!         Previous time step
            V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKN)*ROP_GO(IJKN))*ODT/(VOL(IJK) + VOL(IJKN))  
!         Added mass implicit transient term {Cv eps rop_g dV/dt}
	       IF(Added_Mass) THEN
                 ROP_MA  = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM)  + VOL(IJKN)*ROP_g(IJKN)*EP_s(IJKN,M_AM) )/(VOL(IJK) + VOL(IJKN))
	         V0 = V0 + Cv * ROP_MA * ODT
	       ENDIF
            ENDIF
!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term dVs/dt to virtual mass term		    
	    F_vir = ZERO
	    IF(Added_Mass.AND.(.NOT.CUT_V_TREATMENT_AT(IJK))) THEN        
	      F_vir = ( (V_s(IJK,M_AM) - V_sO(IJK,M_AM)) )*ODT*VOL_V(IJK)
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)     
	      Vss = AVG_Y_N(V_S(IJMK,M_AM),V_s(IJK,M_AM))
	      Vsn = AVG_Y_N(V_s(IJK,M_AM),V_s(IJPK,M_AM))  
	      
	      U_se = AVG_Y(U_s(IJK,M_AM),U_s(IJPK,M_AM),J)
	      Usw = AVG_Y(U_s(IMJK,M_AM),U_s(IMJPK,M_AM),J)
	      Vse = AVG_X(V_s(IJK,M_AM),V_s(IPJK,M_AM),IP1(I))
	      Vsw = AVG_X(V_s(IMJK,M_AM),V_s(IJK,M_AM),I)
	      
	      IF(DO_K) THEN
	         Wst = AVG_Y(W_s(IJK,M_AM),W_s(IJPK,M_AM),J)
	         Wsb = AVG_Y(W_s(IJKM,M_AM),W_s(IJPKM,M_AM),J)
	         Vst = AVG_Z(V_s(IJK,M_AM),V_s(IJKP,M_AM),KP1(K))
	         Vsb = AVG_Z(V_s(IJKM,M_AM),V_s(IJK,M_AM),K)
	         F_vir = F_vir + AVG_Z_T(Wsb,Wst)*OX(I) * (Vst - Vsb)*AXY(IJK)
	      ENDIF
!
! adding convective terms (U dV/dx + V dV/dy) to virtual mass; W dV/dz added above.
	      F_vir = F_vir + V_s(IJK,M_AM)*(Vsn - Vss)*AXZ(IJK) + &
	              AVG_X_E(Usw,U_se,IP1(I))*(Vse - Vsw)*AYZ(IJK)
	    
	      F_vir = F_vir * Cv * ROP_MA
	    ENDIF
!
!!! END VIRTUAL MASS SECTION

! Original terms
!       Volumetric forces
!            ROPGA = AVG_Y(ROP_G(IJK),ROP_G(IJKN),J) 
!            ROGA = AVG_Y(RO_G(IJK),RO_G(IJKN),J) 
!         Previous time step
!            V0 = AVG_Y(ROP_GO(IJK),ROP_GO(IJKN),J)*ODT 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!         pressure drop through porous media
            IF (SIP_AT_N(IJK)) THEN 
               ISV = IS_ID_AT_N(IJK) 
               MUGA = AVG_Y(MU_G(IJK),MU_G(IJKN),J) 
               VPM = MUGA/IS_PC(ISV,1) 
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + HALF*IS_PC(ISV,2)*ROPGA*ABS(&
                  V_G(IJK)) 
            ELSE 
               VPM = ZERO 
            ENDIF 
!
!         Interphase mass transfer
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN 
              VMT = AVG_Y(SUM_R_G(IJK),SUM_R_G(IJKN),J) 
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IJKN)*SUM_R_G(IJKN))/(VOL(IJK) + VOL(IJKN))  
            ENDIF

! Original terms
!            VMT = AVG_Y(SUM_R_G(IJK),SUM_R_G(IJKN),J) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFY_G(IJK) 
!
            ELSE                                 !Model A 
               VBF = ROPGA*BFY_G(IJK) 
!
            ENDIF 

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
		  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
		    DO L = 1,SMAX
		      avgRop = AVG_Y(ROP_S(IJK,L),ROP_S(IJKN,L),J)
		      if(avgRop > ZERO) Ghd_drag = Ghd_drag +&
		           AVG_Y(F_GS(IJK,L),F_GS(IJKN,L),J) * JoiY(IJK,L) / avgRop
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
		              avgDrag = AVG_Y(beta_ij(IJK,IM,L),beta_ij(IJKN,IM,L),J)
		              HYS_drag = HYS_drag + avgDrag * (V_g(ijk) - V_s(IJK,L))
		           ENDIF
		        ENDDO
                     ENDDO        
		  ENDIF
! end of modifications for HYS drag

! loezos	 Source terms from convective mom. flux
	         IF (SHEAR) THEN
		SRT=(2d0*V_sh/XLENGTH)

		VSH_p=VSH(IJK)

		VSH_n=VSH_p
		VSH_s=VSH_p		

		VSH_e=VSH(IJK)+SRT*1d0/oDX_E(I)
		VSH_w=VSH(IJK)-SRT*1d0/oDX_E(IM1(I))


		Source_conv=A_M(IJK,N,m)*VSH_n+A_M(IJK,S,m)*VSH_s&
		+A_M(IJK,W,m)*VSH_w+A_M(IJK,E,m)*VSH_e&
		-(A_M(IJK,N,m)+A_M(IJK,S,m)+A_M(IJK,W,m)+A_M(IJK,E,m))&
		*VSH_p


		ELSE 
		Source_conv=0d0
		END IF
	

!
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT))*VOL_V(IJK)) 
            B_M(IJK,M) = -(SDP + TAU_V_G(IJK)&
	       +Source_conv+((V0+ZMAX((-VMT)))*V_GO(IJK)+VBF+Ghd_drag+HYS_drag)&
               *VOL_V(IJK))+B_M(IJK,M) 
            B_M(IJK,M) = B_M(IJK,M) - F_vir ! adding explicit-part of virtual mass force.
         ENDIF 
      END DO 
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IF(CARTESIAN_GRID) CALL CG_SOURCE_V_G(A_M, B_M, IER)

      CALL SOURCE_V_G_BC(A_M, B_M, IER)

      IF(CARTESIAN_GRID) CALL CG_SOURCE_V_G_BC(A_M, B_M, IER)
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
      END SUBROUTINE SOURCE_V_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_g_BC(A_m, B_m, IER)                           C
!  Purpose: Determine source terms for V_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
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
      SUBROUTINE SOURCE_V_G_BC(A_M, B_M, IER) 
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
      INTEGER          I,  J, K, JM, I1, I2, J1, J2, K1, K2, IJK,& 
                       IM, KM, IJKS, IJMK, IJPK 
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
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
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
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           A_M(IJK,W,M) = ONE 
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
                        IM = IM1(I) 
                        KM = KM1(K) 
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
                              B_M(IJK,M) = -BC_VW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(I)) 
                              A_M(IJK,E,M) = -(HALF*BC_HW_G(L)-ODX_E(I)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN  
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,W,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_VW_G(L) 
                           ELSE 
                              A_M(IJK,W,M) = -(HALF*BC_HW_G(L)-ODX_E(IM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(IM)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN   
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,T,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_VW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(K)*OX(I)) 
                              A_M(IJK,T,M) = -(HALF*BC_HW_G(L)-ODZ_T(K)*OX(I)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN    
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,B,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_VW_G(L) 
                           ELSE 
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX(I)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX(I)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_VW_G(L) 
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
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
			     CALL Wall_Function(IJK,EAST_OF(IJK),ODX_E(I),W_F_Slip)
                             A_M(IJK,E,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_VW_G(L)
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN
			     CALL Wall_Function(IJK,WEST_OF(IJK),ODX_E(IM),W_F_Slip)
                             A_M(IJK,W,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_VW_G(L) 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN
			     CALL Wall_Function(IJK,TOP_OF(IJK),ODZ_T(K)*OX(I),W_F_Slip)
                             A_M(IJK,T,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_VW_G(L) 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN
			     CALL Wall_Function(IJK,BOTTOM_OF(IJK),ODZ_T(KM)*OX(I),W_F_Slip)
                             A_M(IJK,B,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_VW_G(L) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO
! end of wall functions 
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'S') THEN 
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
                           A_M(IJK,S,M) = ONE 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'S') THEN 
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
                           A_M(IJK,S,M) = ONE 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IJMK = JM_OF(IJK) 
                           A_M(IJMK,E,M) = ZERO 
                           A_M(IJMK,W,M) = ZERO 
                           A_M(IJMK,N,M) = ZERO 
                           A_M(IJMK,S,M) = ONE 
                           A_M(IJMK,T,M) = ZERO 
                           A_M(IJMK,B,M) = ZERO 
                           A_M(IJMK,0,M) = -ONE 
                           B_M(IJMK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ELSE IF (BC_PLANE(L) == 'N') THEN 
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
                           IJPK = JP_OF(IJK) 
                           A_M(IJPK,E,M) = ZERO 
                           A_M(IJPK,W,M) = ZERO 
                           A_M(IJPK,N,M) = ONE 
                           A_M(IJPK,S,M) = ZERO 
                           A_M(IJPK,T,M) = ZERO 
                           A_M(IJPK,B,M) = ZERO 
                           A_M(IJPK,0,M) = -ONE 
                           B_M(IJPK,M) = ZERO 
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
!              write(*,*) K1
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
                        B_M(IJK,M) = -V_G(IJK) 
                        IF (BC_PLANE(L) == 'S') THEN 
                           IJKS = SOUTH_OF(IJK) 
                           A_M(IJKS,E,M) = ZERO 
                           A_M(IJKS,W,M) = ZERO 
                           A_M(IJKS,N,M) = ZERO 
                           A_M(IJKS,S,M) = ZERO 
                           A_M(IJKS,T,M) = ZERO 
                           A_M(IJKS,B,M) = ZERO 
                           A_M(IJKS,0,M) = -ONE 
                           B_M(IJKS,M) = -V_G(IJKS) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_V_G_BC 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
