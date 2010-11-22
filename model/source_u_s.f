!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for U_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.  The center        C
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
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
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
      SUBROUTINE SOURCE_U_S(A_M, B_M, IER) 
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
      use kintheory
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
      INTEGER          I, IJK,IMJK, IJMK, IJKE, IJKM, IPJK, IPJKM 
! 
!                      Phase index 
      INTEGER          M, MM, L 
      DOUBLE PRECISION   SUM_EPS_CP
! 
!                      Internal surface number 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
! 
!                      average volume fraction 
      DOUBLE PRECISION EPSA, EPStmp, epse, epsw, epsn, epss, &
                       epst, epsb, epsMix, epsMixE
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
! 
!                      Average density difference 
      DOUBLE PRECISION dro1, dro2, droa 
! 
!                      Average quantities 
      DOUBLE PRECISION wse, EPMUSA 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sdps 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vmt, Vbf, Vcf, Vtza, Vmttmp
!
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop

!			Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag

!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Uge, Ugw, Vgw, Vge, Ugn, Ugs, Wgb, Wgt, Wge, Ugb, Ugt
! 
!                      error message 
      CHARACTER*80     LINE(2) 
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: J,K,IM,JM,IP,JP,KM,KP
      INTEGER :: IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKP,IJKT,IJKTE,IJKB,IJKBE
      DOUBLE PRECISION :: Ue,Uw,Un,Us,Ut,Ub
      DOUBLE PRECISION :: P_CUT,z_plane,DH,Nx,Ny,Nz,B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------


      DO M = 1, MMAX 
        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
          IF (MOMENTUM_X_EQ(M)) THEN 

! CHEM & ISAT (nan xie)
! Set the source terms zero
            IF (CALL_DI .or. CALL_ISAT) THEN
              SUM_R_S_temp = SUM_R_S
              SUM_R_S = ZERO
            ENDIF


!$omp  parallel do private( IJK, IJKE, ISV, Sdp, Sdps, V0, Vmt, Vbf, &
!$omp&  I,PGE,DRO1,DRO2,DROA, IJKM,IPJK,IPJKM,  WSE,VCF,EPMUSA,VTZA, &
!$omp&  EPSA, EPStmp, ROPSA, LINE,SUM_EPS_CP,MM) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3 

! Wall or impermeable internal surface
                I = I_OF(IJK) 
		J = J_OF(IJK)
		K = K_OF(IJK)
                IJKE = EAST_OF(IJK)
                IMJK = IM_OF(IJK) 
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK) 
                IPJK = IP_OF(IJK) 
                IJPK = JP_OF(IJK) 
		IPJMK = IP_OF(IJMK)
                IPJKM = IP_OF(IJKM)
		IJKP = KP_OF(IJK)
                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixE= ZERO          
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I) 
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixE to be used for modelB
                    epsMixE = epsMixE + EP_S(IJKE,L)
		    IF(IP_AT_E(IJK)) THEN
		       U_S(IJK,L) = ZERO
                    ELSEIF(SIP_AT_E(IJK)) THEN 
		       ISV = IS_ID_AT_E(IJK)
		       U_S(IJK,L) = IS_VEL_S(ISV,L)
		    ENDIF
                  ENDDO                        
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
                ENDIF                  
                IF (IP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
                ELSEIF (SIP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_E(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 

! Dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                      EPSw = ZERO
                      EPSe = ZERO
                      EPSn = ZERO
                      EPSs = ZERO
                      EPSt = ZERO
                      EPSb = ZERO
                      DO L = 1, SMAX
                        EPSw = EPSw + EP_S(WEST_OF(IJK),L)
                        EPSe = EPSe + EP_S(EAST_OF(IJK),L)
                        EPSn = EPSn + EP_S(NORTH_OF(IJK),L)
                        EPSs = EPSs + EP_S(SOUTH_OF(IJK),L)
                        IF(.NOT. NO_K) THEN
                          EPSt = EPSt + EP_S(TOP_OF(IJK),L)
                          EPSb = EPSb + EP_S(BOTTOM_OF(IJK),L)
                        ENDIF
                      ENDDO
                  ELSE
                      EPSw = EP_S(WEST_OF(IJK),M)
                      EPSe = EP_S(EAST_OF(IJK),M)
                      EPSn = EP_S(NORTH_OF(IJK),M)
                      EPSs = EP_S(SOUTH_OF(IJK),M)
                      IF(.NOT. NO_K) THEN
                        EPSt = EP_S(TOP_OF(IJK),M)
                        EPSb = EP_S(BOTTOM_OF(IJK),M)
                      ENDIF
                  ENDIF
! using the average boundary cell values to compute U_s (sof, Aug 23 2005)
                  IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE 
                  IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
                  IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE 
                  IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF(.NOT. NO_K) THEN
                    IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE 
                    IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
                  ENDIF
 
                  IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
                    A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
                    B_M(IJK,M) = -U_S(IJK,M)          
                  ELSE
                    A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
                  ENDIF

! Normal case
                ELSE 

! Surface forces
!   Pressure terms
                  PGE = P_G(IJKE) 
                  IF (CYCLIC_X_PD) THEN 
                    IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X 
                  ENDIF 
                  IF (MODEL_B) THEN 
                    SDP = ZERO 
                  ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        SDP = -P_SCALE*EPSA*(PGE - P_G(IJK))*AYZ(IJK) 
                     ELSE
                        SDP = -P_SCALE*EPSA*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
                     ENDIF
! Original terms
!                     SDP = -P_SCALE*EPSA*(PGE - P_G(IJK))*AYZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

                  ENDIF 

                  IF (CLOSE_PACKED(M)) THEN 
                    IF(SMAX > 1 .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                      SUM_EPS_CP=0.0 
                      DO MM=1,SMAX
                        IF (CLOSE_PACKED(MM))&
                          SUM_EPS_CP=SUM_EPS_CP+AVG_X(EP_S(IJK,MM),EP_S(IJKE,MM),I)
                      ENDDO
                      SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                      SDPS = -( (P_S(IJKE,M)-P_S(IJK,M))+(EPSA/SUM_EPS_CP)*&
                        (P_STAR(IJKE)-P_STAR(IJK)) )*AYZ(IJK) 
                    ELSE
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        SDPS =-((P_S(IJKE,M)-P_S(IJK,M))+(P_STAR(IJKE)-P_STAR(IJK)))*AYZ(IJK)  
                     ELSE
                        SDPS =-((P_S(IJKE,M)* A_UPG_E(IJK)-P_S(IJK,M)* A_UPG_W(IJK)) & 
                              +(P_STAR(IJKE)* A_UPG_E(IJK)-P_STAR(IJK)* A_UPG_W(IJK)))
                     ENDIF
! Original terms
!                     SDPS =-((P_S(IJKE,M)-P_S(IJK,M))+(P_STAR(IJKE)-P_STAR(IJK)))*AYZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                    ENDIF
                  ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        SDPS = -(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK) 
                     ELSE
                        SDPS = - (P_S(IJKE,M) * A_UPG_E(IJK) - P_S(IJK,M) * A_UPG_W(IJK))
                     ENDIF
! Original terms
!                     SDPS = -(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

                  ENDIF 

! Shear stress terms

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
! Volumetric forces
                     ROPSA = HALF * (VOL(IJK)*ROP_S(IJK,M) + VOL(IPJK)*ROP_S(IJKE,M))/VOL_U(IJK) 
! Previous time step
                     V0 = HALF * (VOL(IJK)*ROP_SO(IJK,M) + VOL(IPJK)*ROP_SO(IJKE,M))*ODT/VOL_U(IJK)  
!         Added mass implicit transient term {Cv eps rop_g dU/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
		       ROP_MA = AVG_X(ROP_g(IJK)*EP_s(IJK,M),ROP_g(IJKE)*EP_s(IJKE,M),I)
		       V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ELSE
! Volumetric forces
                     ROPSA =  (VOL(IJK)*ROP_S(IJK,M) + VOL(IPJK)*ROP_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
! Previous time step
                     V0 = (VOL(IJK)*ROP_SO(IJK,M) + VOL(IPJK)*ROP_SO(IJKE,M))*ODT/(VOL(IJK) + VOL(IPJK))
!         Added mass implicit transient term {Cv eps rop_g dU/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
                       ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + VOL(IPJK)*ROP_g(IJKE)*EP_s(IJKE,M))/(VOL(IJK) + VOL(IPJK))
		       V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ENDIF
!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term dvg/dt - dVs/dt to virtual mass term    
	    F_vir = ZERO
	    IF(Added_Mass .AND. M==M_AM .AND.(.NOT.CUT_U_TREATMENT_AT(IJK))) THEN      
	      F_vir = ( (U_G(IJK) - U_GO(IJK)) )*ODT*VOL_U(IJK)
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)    
	      Ugw = AVG_X_E(U_G(IMJK),U_G(IJK),I)
	      Uge = AVG_X_E(U_G(IJK),U_G(IPJK),IP1(I))  
	      Vgw = AVG_Y_N(V_G(IJMK),V_G(IJK))  
	      
	      Vge = AVG_Y_N(V_G(IPJMK),V_G(IPJK)) 
	      Ugs = AVG_Y(U_G(IJMK),U_G(IJK),JM1(J))
	      Ugn = AVG_Y(U_G(IJK),U_G(IJPK),J)
	      
	      IF(DO_K) THEN
	         Wgb = AVG_Z_T(W_g(IJKM),W_g(IJK))  
	         Wgt = AVG_Z_T(W_g(IPJKM),W_g(IPJK)) 
		 Wge = AVG_X(Wgb,Wgt,I)
	         Ugb = AVG_Z(U_g(IJKM),U_g(IJK),KM1(K))
	         Ugt = AVG_Z(U_g(IJK),U_g(IJKP),K)
	         F_vir = F_vir + Wge*OX_E(I) * (Ugt*AXY(IJKP) - Ugb*(AXY(IJK)))
		 IF (CYLINDRICAL) F_vir = F_vir - Wge**2*OX_E(I) ! centrifugal force
	      ENDIF
!
! adding convective terms (U dU/dx + V dU/dy + W dU/dz) to virtual mass
	      F_vir = F_vir + U_g(IJK)*(Uge*AYZ(IPJK) - Ugw*AYZ(IJK)) + &
	         AVG_X(Vgw,Vge,I) * (Ugn*AXZ(IJPK) - Ugs*(AXZ(IJK)))
	                      
	    
	      F_vir = F_vir * Cv * ROP_MA
	    ENDIF
!
!!! END VIRTUAL MASS SECTION

! Original terms
! Volumetric forces
!                  ROPSA = HALF * (VOL(IJK)*ROP_S(IJK,M) + VOL(IPJK)*ROP_S(IJKE,M))/VOL_U(IJK) 
!             Previous time step
!                  V0 = HALF * (VOL(IJK)*ROP_SO(IJK,M) + VOL(IPJK)*ROP_SO(IJKE,M))*ODT/VOL_U(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! Interphase mass transfer
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                    VMTtmp = ZERO
                    DO L = 1,SMAX
                      VMTtmp = VMTtmp + HALF*(VOL(IJK)*SUM_R_S(IJK,L) + &
                             VOL(IPJK)*SUM_R_S(IJKE,L))/VOL_U(IJK)
                    ENDDO
                    VMT = VMTtmp
                  ELSE
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                     VMT = HALF * (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IPJK)*SUM_R_S(IJKE,M))/VOL_U(IJK)  
                  ELSE
                     VMT = (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IPJK)*SUM_R_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                  ENDIF
! Original terms
!                  VMT = HALF * (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IPJK)*SUM_R_S(IJKE,M))/VOL_U(IJK)  
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  ENDIF

! Body force
                  IF (MODEL_B) THEN 
                    IF (TRIM(KT_TYPE) /= 'GHD') THEN
                      DRO1 = (RO_S(M)-RO_G(IJK))*EP_S(IJK,M) 
                      DRO2 = (RO_S(M)-RO_G(IJKE))*EP_S(IJKE,M) 
                      DROA = AVG_X(DRO1,DRO2,I) 
                      VBF = DROA*BFX_S(IJK,M) 
                    ELSE ! GHD and M = MMAX
                      DRO1 = ROP_S(IJK,M)  - RO_G(IJK) *epsMix
                      DRO2 = ROP_S(IJKE,M) - RO_G(IJKE)*epsMixE
                      DROA = AVG_X(DRO1,DRO2,I) 
                      VBF = DROA*BFX_S(IJK,M) 
                    ENDIF
                  ELSE ! model A
                    VBF = ROPSA*BFX_S(IJK,M) 
                  ENDIF 

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                    DO L = 1,SMAX
                      avgRop = AVG_X(ROP_S(IJK,L),ROP_S(IJKE,L),I)
                      if(avgRop > ZERO) Ghd_drag = Ghd_drag - &
                           AVG_X(F_GS(IJK,L),F_GS(IJKE,L),I) * JoiX(IJK,L) / avgRop
                    ENDDO
                  ENDIF
! end of modifications for GHD theory

! Additional force for HYS drag force, do not use with mixture GHD theory
                  HYS_drag = ZERO
                  IF (TRIM(DRAG_TYPE) .EQ. 'HYS' .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                     DO L = 1,MMAX
                        IF (L /= M) THEN
                           avgDrag = AVG_X(beta_ij(IJK,M,L),beta_ij(IJKE,M,L),I)
                           HYS_drag = HYS_drag - avgDrag * (U_g(ijk) - U_s(IJK,L))
                        ENDIF
                     ENDDO
                  ENDIF
! end of modifications for HYS drag


! Special terms for cylindrical coordinates
                  IF (CYLINDRICAL) THEN 

!   centrifugal force
                    WSE = AVG_X(HALF*(W_S(IJK,M)+W_S(IJKM,M)),HALF*(W_S(IPJK,M&
                        )+W_S(IPJKM,M)),I) 
                    VCF = ROPSA*WSE**2*OX_E(I) 
		    IF(Added_Mass .AND. M==M_AM) VCF = VCF + Cv*ROP_MA*WSE**2*OX_E(I) ! virtual mass contribution.

!   -(2mu/x)*(u/x) part of Tau_zz/X
                    EPMUSA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                    VTZA = 2.d0*EPMUSA*OX_E(I)*OX_E(I) 
                  ELSE 
                    VCF = ZERO 
                    VTZA = ZERO 
                  ENDIF 

! Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(&
                     IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+ZMAX(VMT)+VTZA)*&
                     VOL_U(IJK)) 

                  IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN 
                    B_M(IJK,M) = -(SDP + KTMOM_U_S(IJK,M) + SDPS + TAU_U_S(IJK,M)+&
                    ((V0+ZMAX((-VMT)))*U_SO(IJK,M)+VBF+VCF+HYS_drag)*VOL_U(IJK))+B_M(IJK,M) 
                  ELSE
                    B_M(IJK,M) = -(SDP + SDPS + TAU_U_S(IJK,M)+&
                    ((V0+ZMAX((-VMT)))*U_SO(IJK,M)+VBF+VCF+Ghd_drag+HYS_drag)*VOL_U(IJK))+B_M(IJK,M)
                  ENDIF

                  B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit term of virtual mass force
                ENDIF ! end if sip or ip or dilute flow branch
            ENDDO 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CARTESIAN_GRID) CALL CG_SOURCE_U_S (A_M, B_M, M, IER)

            CALL SOURCE_U_S_BC (A_M, B_M, M, IER) 

            IF(CARTESIAN_GRID) CALL CG_SOURCE_U_S_BC (A_M, B_M, M, IER)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================


! CHEM & ISAT (nan xie)
            IF (CALL_DI .or. CALL_ISAT) THEN
                SUM_R_S = SUM_R_S_temp
            ENDIF
 
          ENDIF  
        ENDIF ! for GHD Theory

      ENDDO 
      
      RETURN  
      END SUBROUTINE SOURCE_U_S 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s_BC(A_m, B_m, M, IER)                        C
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
      SUBROUTINE SOURCE_U_S_BC(A_M, B_M, M, IER) 
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
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IPJK, IP 
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
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,N,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J)) 
                                 A_M(IJK,N,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,S,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,S,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,T,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*&
                                    OX_E(I)) 
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,B,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX_E(I)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
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
                        B_M(IJK,M) = -U_S(IJK,M) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_S(IJKW,M) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
	       
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_U_S_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_U_s(I1, I2, J1, J2, K1, K2, L, M, A_m, b_m)      C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava,                Date: 14-APR-98  C
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
      SUBROUTINE JJ_BC_U_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                       JM, KM, IPJK
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
!                      coefficients for granular bc
      DOUBLE PRECISION hw, gw, cw
!      
 !-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!
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
!
               IF (FLUID_AT(NORTH_OF(IJK))) THEN 
	          IPJK = IP_OF(NORTH_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
		  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,N,M) = ONE 
                  ELSE 
! start anuj 4/20
                     IF (FRICTION .AND. (ONE-EP_G(NORTH_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, NORTH_OF(IJK), 'N', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
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
	          IPJK = IP_OF(SOUTH_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,S,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(SOUTH_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, SOUTH_OF(IJK), 'S', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
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
!
               ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
	          IPJK = IP_OF(TOP_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,T,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(TOP_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, TOP_OF(IJK), 'T', 'U', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'U', M, L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
	          IPJK = IP_OF(BOTTOM_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,B,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(BOTTOM_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, BOTTOM_OF(IJK), 'B', 'U', L&
                           , M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
     
      RETURN  
      END SUBROUTINE JJ_BC_U_S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
