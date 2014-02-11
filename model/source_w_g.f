!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g                                              C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_G(A_M, B_M, IER) 

!-----------------------------------------------
! Modules
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
      USE cutcell
      USE quadric
      USE mms
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices 
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,& 
                 IJKE, IJKW, IJKTE, IJKTW, IM, IPJK,   &
                 IJKM, IJMK, IJMKP, IJPK
! Phase index 
      INTEGER :: M, L, MM
! Internal surface 
      INTEGER :: ISV 
! Pressure at top cell 
      DOUBLE PRECISION :: PgT 
! Average volume fraction 
      DOUBLE PRECISION :: EPGA 
! Average density 
      DOUBLE PRECISION :: ROPGA, ROGA 
! Average viscosity 
      DOUBLE PRECISION :: MUGA 
! Average coefficients 
      DOUBLE PRECISION :: Cte, Ctw, EPMUoX 
! Average U_g 
      DOUBLE PRECISION Ugt 
! Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sxzb 
! Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcoa, Vcob, Vxza, Vxzb 
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, U_se, Usw, Ust, Vsb, Vst, &
                          Wse, Wsw, Wsn, Wss, Wst, Wsb
      DOUBLE PRECISION F_vir
! error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------

! Set reference phase to gas
      M = 0 

      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN  

!
! trailing commas added to the first two lines in the 
! private list
! Charles Crosby
! CHPC, 26 September 2013
!

!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IJKT, IJKM, IJKP, IMJK, IPJK, IJMK,     &
!$omp          IMJKP, IJPK, IJMKP, IJKTE, IJKTW, IM, IJKW, IJKE,     &
!$omp          EPGA, PGT, SDP, ROPGA, ROGA, V0, ISV, MUGA, Vpm,      &
!$omp          Vmt, Vbf, F_vir, Ghd_drag, avgRop, HYS_drag, avgDrag, &
!$omp          MM, L, VXZA, VXZB, VCOA, VCOB, CTE, CTW, UGT,         &
!$omp          SXZB, EPMUOX) 
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IMJKP = KP_OF(IMJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IJMKP = KP_OF(IJMK)

         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K) 

! Impermeable internal surface
         IF (IP_AT_T(IJK)) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
! set velocity equal to that of bottom or top cell if solids are 
! present in those cells else set velocity equal to known value 
            IF (EP_G(BOTTOM_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,B,M) = ONE 
            ELSE IF (EP_G(TOP_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,T,M) = ONE 
            ELSE 
               B_M(IJK,M) = -W_G(IJK) 
            ENDIF 

! Cartesian grid implementation
         ELSEIF (BLOCKED_W_CELL_AT(IJK)) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 

! Normal case
         ELSE 

! Surface forces

! Pressure term
            PGT = P_G(IJKT) 
            IF (CYCLIC_Z_PD) THEN 
               IF (KMAP(K_OF(IJK)).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z 
            ENDIF 
            IF (MODEL_B) THEN 
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*(PGT - P_G(IJK))*AXY(IJK) 
               ELSE
                   SDP = -P_SCALE*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
            ELSE 
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK)
               ELSE
                   SDP = -P_SCALE*EPGA*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
            ENDIF 

            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K) 
               ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K) 

! Previous time step
               V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT 

! Added mass implicit transient term {Cv eps rop_g dW/dt}
               IF(Added_Mass) THEN
                  ROP_MA = AVG_Z(ROP_g(IJK)*EP_s(IJK,M_AM),&
                     ROP_g(IJKT)*EP_s(IJKT,M_AM),K)
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IJKT)*ROP_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
               ROGA  = (VOL(IJK)*RO_G(IJK) + VOL(IJKT)*RO_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKT)*ROP_GO(IJKT))*&
                  ODT/(VOL(IJK) + VOL(IJKT))  
! Added mass implicit transient term {Cv eps rop_g dW/dt}
               IF(Added_Mass) THEN
                  ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM) + &
                     VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M_AM))/&
                     (VOL(IJK) + VOL(IJKT))
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dWs/dt to virtual mass term		    
            F_vir = ZERO
            IF(Added_Mass.AND.(.NOT.CUT_W_TREATMENT_AT(IJK))) THEN
               F_vir = ( (W_s(IJK,M_AM) - W_sO(IJK,M_AM)) )*ODT*VOL_W(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
               Wsb = AVG_Z_T(W_S(IJKM,M_AM),W_s(IJK,M_AM))
               Wst = AVG_Z_T(W_s(IJK,M_AM),W_s(IJKP,M_AM))  
               U_se = AVG_Z(U_s(IJK,M_AM),U_s(IJKP,M_AM),K)
               Usw = AVG_Z(U_s(IMJK,M_AM),U_s(IMJKP,M_AM),K)
               Ust = AVG_X_E(Usw,U_se,IP1(I))
               Wse = AVG_X(W_s(IJK,M_AM),W_s(IPJK,M_AM),IP1(I))
               Wsw = AVG_X(W_s(IMJK,M_AM),W_s(IJK,M_AM),I)
               Vsb = AVG_Y_N(V_s(IJMK,M_AM),V_s(IJK,M_AM))  
               Vst = AVG_Y_N(V_s(IJMKP,M_AM),V_s(IJKP,M_AM)) 
               Wss = AVG_Y(W_s(IJMK,M_AM),W_s(IJK,M_AM),JM1(J))
               Wsn = AVG_Y(W_s(IJK,M_AM),W_s(IJPK,M_AM),J)

! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass.
               F_vir = F_vir + W_s(IJK,M_AM)*OX(I)* &
                  (Wst - Wsb)*AXY(IJK) + Ust*(Wse - Wsw)*AYZ(IJK) + &
                  AVG_Z(Vsb,Vst,K)*(Wsn - Wss)*AXZ(IJK)

! Coriolis force
               IF (CYLINDRICAL) F_vir = F_vir + &
                  Ust*W_s(IJK,M_AM)*OX(I)  
               F_vir = F_vir * Cv * ROP_MA
            ENDIF

! pressure drop through porous media
            IF (SIP_AT_T(IJK)) THEN 
               ISV = IS_ID_AT_T(IJK) 
               MUGA = AVG_Z(MU_G(IJK),MU_G(IJKT),K) 
               VPM = MUGA/IS_PC(ISV,1) 
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + &
                  HALF*IS_PC(ISV,2)*ROPGA*ABS(W_G(IJK)) 
            ELSE 
               VPM = ZERO 
            ENDIF 

! Interphase mass transfer
            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
               VMT = AVG_Z(SUM_R_G(IJK),SUM_R_G(IJKT),K) 
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IJKT)*SUM_R_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))  
            ENDIF

! Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFZ_G(IJK) 
            ELSE 
               VBF = ROPGA*BFZ_G(IJK) 
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

! Additional force for HYS drag force, do not use with mixture GHD theory
            avgDrag = ZERO
            HYS_drag = ZERO
            IF (TRIM(DRAG_TYPE) .EQ. 'HYS' .AND. TRIM(KT_TYPE) /= 'GHD') THEN
               DO MM=1,MMAX
                  DO L = 1,MMAX
                     IF (L /= MM) THEN
                        avgDrag = AVG_Z(beta_ij(IJK,MM,L),beta_ij(IJKT,MM,L),K)
                        HYS_drag = HYS_drag + avgDrag * (W_g(IJK) - W_s(IJK,L))
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

! Special terms for cylindrical coordinates
            VCOA = ZERO 
            VCOB = ZERO 
            SXZB = ZERO 
            VXZA = ZERO 
            VXZB = ZERO 
            CTE  = ZERO
            CTW  = ZERO
            IF (CYLINDRICAL) THEN 
! Coriolis force
               IMJK = IM_OF(IJK) 
               IJKP = KP_OF(IJK) 
               IMJKP = KP_OF(IMJK) 
               UGT = AVG_Z(HALF*(U_G(IJK)+U_G(IMJK)),&
                  HALF*(U_G(IJKP)+U_G(IMJKP)),K) 
               IF (UGT > ZERO) THEN 
                  VCOA = ROPGA*UGT*OX(I) 
                  VCOB = ZERO 
! virtual mass contribution                  
                  IF(Added_Mass) VCOA = VCOA + Cv*ROP_MA*UGT*OX(I) 
               ELSE 
                  VCOA = ZERO 
                  VCOB = -ROPGA*UGT*W_G(IJK)*OX(I) 
! virtual mass contribution                  
                  IF(Added_Mass) VCOB = VCOB - Cv*ROP_MA*UGT*W_G(IJK)*OX(I)
               ENDIF 

! Term from tau_xz: integral of (1/x)*(d/dx)(x*mu*(-w/x))
               IJKE = EAST_OF(IJK) 
               IJKW = WEST_OF(IJK) 
               IJKTE = TOP_OF(IJKE) 
               IJKTW = TOP_OF(IJKW) 
               IM = IM1(I) 
               IPJK = IP_OF(IJK) 
               CTE = HALF*AVG_Z_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),&
                                  AVG_X_H(MU_GT(IJKT),MU_GT(IJKTE),I),K)*&
                                  OX_E(I)*AYZ_W(IJK) 
               CTW = HALF*AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),&
                                  AVG_X_H(MU_GT(IJKTW),MU_GT(IJKT),IM),K)*&
                                  DY(J)*(HALF*(DZ(K)+DZ(KP1(K)))) 
                    ! same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! (mu/x)*(dw/dx) part of tau_xz/x
               EPMUOX = AVG_Z(MU_GT(IJK),MU_GT(IJKT),K)*OX(I) 
               VXZB = ZERO 
               A_M(IJK,E,M) = A_M(IJK,E,M) + HALF*EPMUOX*ODX_E(I)*VOL_W(IJK) 
               A_M(IJK,W,M) = A_M(IJK,W,M) - HALF*EPMUOX*ODX_E(IM)*VOL_W(IJK) 

! -(mu/x)*(w/x) part of tau_xz/x
               VXZA = EPMUOX*OX(I) 
            ELSE 
               VCOA = ZERO 
               VCOB = ZERO 
               SXZB = ZERO 
               VXZA = ZERO 
               VXZB = ZERO 
            ENDIF 

! Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
               A_M(IJK,N,M)+A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+&
               (V0+VPM+ZMAX(VMT)+VCOA + VXZA)*VOL_W(IJK) + CTE - CTW) 
            A_M(IJK,E,M) = A_M(IJK,E,M) - CTE
            A_M(IJK,W,M) = A_M(IJK,W,M) + CTW 
            B_M(IJK,M) = B_M(IJK,M) - ( SDP + TAU_W_G(IJK) + SXZB + &
               ( (V0+ZMAX((-VMT)))*W_GO(IJK) + VBF + VCOB + VXZB + &
               Ghd_drag+HYS_drag)*VOL_W(IJK) ) 
! adding explicit term of virtual mass force
            B_M(IJK,M) = B_M(IJK,M) - F_vir 
! MMS Source term.
            IF(USE_MMS) B_M(IJK,M) = &
               B_M(IJK,M) - MMS_W_G_SRC(IJK)*VOL_W(IJK)
         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
!$omp end parallel do

! modifications for cartesian grid implementation 
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G(A_M, B_M, IER)
! modifications for bc
      CALL SOURCE_W_G_BC (A_M, B_M, IER) 
! modifications for cartesian grid implementation 
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G_BC(A_M, B_M, IER)

      RETURN  
      END SUBROUTINE SOURCE_W_G 



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g_BC                                           C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_G_BC(A_M, B_M, IER) 

!-----------------------------------------------
! Modules
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
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! C_mu is constant in turbulent viscosity 
      DOUBLE PRECISION, PARAMETER :: C_mu = 0.09D0
! Kappa is Von Karmen constant
      DOUBLE PRECISION, PARAMETER :: Kappa = 0.42D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Boundary condition 
      INTEGER :: L 
! Indices 
      INTEGER :: I, J, K, KM, I1, I2, J1, J2, K1, K2, IJK,& 
                 IM, JM, IJKB, IJKM, IJKP 
! Phase index
      INTEGER :: M 
! Turbulent shear at walls
      DOUBLE PRECISION W_F_Slip 

!-----------------------------------------------
! Include statment functions
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------

! Set reference phase to gas
      M = 0 


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default 
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the top and bottom xy planes do not have to be explictly addressed for
! the w-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a 
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_t branch in the above source routine).
! ---------------------------------------------------------------->>>

! south xz plane
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
            ELSEIF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO

! north xz plane
      J1 = JMAX2 
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3 
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSEIF (FS_WALL_AT(IJK)) THEN 
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO 

! west yz plane
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
            ELSEIF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO 

! east yz plane
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
            ELSEIF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO 
! End setting the default boundary conditions
! ----------------------------------------------------------------<<<

! Setting user specified boundary conditions

      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
! ---------------------------------------------------------------->>>
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
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
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
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN 
                           A_M(IJK,W,M) = -ONE 
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = -ONE 
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = -ONE 
                        ENDIF 
                     ENDDO 
                  ENDDO 
               ENDDO

            ELSEIF (BC_TYPE(L) == 'FREE_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN
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
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = ONE 
                        ENDIF 
                     ENDDO 
                  ENDDO 
               ENDDO 

            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN
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
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
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
                                 A_M(IJK,0,M) = -( HALF*(BC_HW_G(L)-&
                                    OX_E(I)) + ODX_E(I) ) 
                                 A_M(IJK,E,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(I)) - ODX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(I)) 
                                 A_M(IJK,E,M) = -(HALF*BC_HW_G(L)-ODX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ENDIF 
                           ENDIF 
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN  
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,W,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              IF (CYLINDRICAL) THEN 
                                 A_M(IJK,W,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(IM)) - ODX_E(IM)) 
                                 A_M(IJK,0,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(IM)) + ODX_E(IM)) 
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ELSE 
                                 A_M(IJK,W,M) = -(HALF*BC_HW_G(L)-ODX_E(IM))
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(IM))
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                              ENDIF 
                           ENDIF 
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_WW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L) 
                           ENDIF 
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
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
                     ENDDO 
                  ENDDO 
               ENDDO

! Setting wall boundary conditions when K_EPSILON
! wall functions for V-momentum are specify in this section of the code
            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL'   .OR.  &
                    BC_TYPE(L) == 'NO_SLIP_WALL'    .OR.  &
                    BC_TYPE(L) == 'FREE_SLIP_WALL'  .AND. &
                    K_Epsilon )THEN
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
                              W_F_Slip = ( ONE/&
                                 (ODX_E(I)+HALF*OX_E(I)) )*          &
                                 ( ODX_E(I) - OX_E(I) -              &
                                   RO_g(EAST_OF(IJK))*C_mu**0.25*    &
                                   SQRT(K_Turb_G((EAST_OF(IJK))))/   &
                                   MU_gT(EAST_OF(IJK))*Kappa/        &
                                   LOG(9.81D0/ODX_E(I)/(2.D0)*       &
                                   RO_g(EAST_OF(IJK))*C_mu**0.25*    &
                                   SQRT(K_Turb_G((EAST_OF(IJK))))/   &
                                   MU_g(EAST_OF(IJK))) )
                           ELSE
                              CALL Wall_Function(IJK,EAST_OF(IJK),&
                                 ODX_E(I),W_F_Slip)
                           ENDIF
                           A_M(IJK,E,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (CYLINDRICAL) THEN
                              W_F_Slip =  (ONE/&
                                 (ONE*ODX_E(IM) + HALF*OX_E(IM)))*    &
                                 ( ONE*ODX_E(IM) - OX_E(IM) -         &
                                   RO_g(WEST_OF(IJK))*C_mu**0.25*     &
                                   SQRT(K_Turb_G((WEST_OF(IJK))))/    &
                                   MU_gT(WEST_OF(IJK))*Kappa/         &
                                   LOG(9.81D0/ODX_E(IM)/(2.D0)*       &
                                   RO_g(WEST_OF(IJK))*C_mu**0.25*     &
                                   SQRT(K_Turb_G((WEST_OF(IJK))))/    &
                                   MU_g(WEST_OF(IJK))) )
                           ELSE
                              CALL Wall_Function(IJK,WEST_OF(IJK),&
                                 ODX_E(IM),W_F_Slip)
                           ENDIF
                           A_M(IJK,W,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                           A_M(IJK,N,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                           A_M(IJK,S,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE 
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_WW_G(L)
                        ENDIF 
                     ENDDO 
                  ENDDO 
               ENDDO

! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'B') THEN 
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the boundary cell equal to the velocity of 
! the adjacent fluid cell                       
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
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDIF
! end setting of p_inflow or p_otuflow flow boundary conditions
! ----------------------------------------------------------------<<<
               
! Setting outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L) == 'OUTFLOW') THEN 
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
                           IJKM = KM_OF(IJK) 
                           A_M(IJKM,E,M) = ZERO 
                           A_M(IJKM,W,M) = ZERO 
                           A_M(IJKM,N,M) = ZERO 
                           A_M(IJKM,S,M) = ZERO 
                           A_M(IJKM,T,M) = ZERO 
                           A_M(IJKM,B,M) = ONE 
                           A_M(IJKM,0,M) = -ONE 
                           B_M(IJKM,M) = ZERO 
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ELSEIF (BC_PLANE(L) == 'T') THEN 
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
                           IJKP = KP_OF(IJK) 
                           A_M(IJKP,E,M) = ZERO 
                           A_M(IJKP,W,M) = ZERO 
                           A_M(IJKP,N,M) = ZERO 
                           A_M(IJKP,S,M) = ZERO 
                           A_M(IJKP,T,M) = ONE 
                           A_M(IJKP,B,M) = ZERO 
                           A_M(IJKP,0,M) = -ONE 
                           B_M(IJKP,M) = ZERO 
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDIF
! end setting of outflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting bc that are defined but not nsw, fsw, psw, p_inflow,
! p_outflow, or outflow (at this time, this section addresses 
! mass_inflow and mass_outflow type boundaries)
! ---------------------------------------------------------------->>>
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
! setting the velocity in the boundary cell equal to what is known
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -W_G(IJK) 
                        IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
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
                     ENDDO 
                  ENDDO 
               ENDDO 
            ENDIF   ! end if/else (bc_type)
                    ! ns, fs, psw; else
                    ! p_inflow, p_outflow, or outflow; else
! end setting of 'else' flow boundary conditions
! (mass_inflow/mass_outflow)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN  
      END SUBROUTINE SOURCE_W_G_BC  



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_W_G                                        C
!  Purpose: Adds point sources to the gas phase W-Momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_W_G(A_M, B_M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar    
      use constant
      use geometry
      use indices
      use physprop
      use ps
      use run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M) 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!----------------------------------------------- 
! Local variables
!----------------------------------------------- 
! Indices 
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lKT, lKB
! terms of bm expression
      DOUBLE PRECISION :: pSource
!----------------------------------------------- 
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'function.inc'
!----------------------------------------------- 

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each 
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_W_g(PSV)) < small_number) cycle PS_LP

         if(PS_W_g(PSV) < ZERO) then
            lKB = PS_K_B(PSV)-1
            lKT = PS_K_T(PSV)-1
         else
            lKB = PS_K_B(PSV)
            lKT = PS_K_T(PSV)
         endif

         do k = lKB, lKT
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK,M) = B_M(IJK,M) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
