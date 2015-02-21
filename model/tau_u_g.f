!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_U_g(TAU_U_g, IER)                             C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
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
      SUBROUTINE CALC_TAU_U_G(TAU_U_G)
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
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE bc
      USE quadric
      USE cutcell
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
!                      TAU_U_g
      DOUBLE PRECISION TAU_U_g(DIMENSION_3)
!
!                      Indices
      INTEGER          I, J, JM, K, KM, IJK, IJKE, IPJK, IP, IMJK, IJKN,&
                       IJKNE, IJKS, IJKSE, IPJMK, IJMK, IJKT, IJKTE,&
                       IJKB, IJKBE, IJKM, IPJKM
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Average viscosity
      DOUBLE PRECISION EPMU_gte, EPMU_gbe, EPMUGA
!
!                      Average dW/Xdz
      DOUBLE PRECISION dWoXdz
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz
!
!                      Source terms (Volumetric)
      DOUBLE PRECISION Vtzb
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: V_NODE_AT_NE,V_NODE_AT_NW,V_NODE_AT_SE,V_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TE,W_NODE_AT_TW,W_NODE_AT_BE,W_NODE_AT_BW
      DOUBLE PRECISION :: dvdx_at_N,dvdx_at_S
      DOUBLE PRECISION :: dwdx_at_T,dwdx_at_B
      DOUBLE PRECISION :: Xi,Yi,Zi,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: MU_GT_CUT,SSY_CUT,SSZ_CUT
      DOUBLE PRECISION :: UW_g,VW_g,WW_g
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!$omp  parallel do default(none) &
!$omp              private( IJK, I, IJKE, EPGA, IP, J, JM, K, KM,  &
!$omp              IPJK,IMJK,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJMK,IJKT,IJKTE,  &
!$omp              IJKB,IJKBE,IJKM,IPJKM, &
!$omp              SBV,  SSX,SSY,  EPMU_GTE,EPMU_GBE, SSZ, &
!$omp              EPMUGA,DWOXDZ,VTZB,BCV,NOC_UG, uw_g, vw_g, ww_g, cut_tau_ug, bct,mu_gt_cut,del_h,nx,ny,nz,v_node_at_ne, v_node_at_sw,v_node_at_se, v_node_at_nw, w_node_at_bw, w_node_at_be,w_node_at_tw,w_node_at_te,xi, vi,yi,zi,sx,sy,sz,dvdx_at_S,dvdx_at_N,ssy_cut,dwdx_at_B,dwdx_at_T,wi,ssz_cut) &
!$omp  shared(ijkstart3, ijkend3, do_k, bc_type, cut_u_cell_at, bc_u_id, tau_u_g, bc_uw_g, bc_vw_g, bc_ww_g, x_u, y_u, z_u, vol, x_v, y_v,z_v, wall_v_at, area_u_cut,oneodx_e_u,oneodx_e_v,oneodx_e_w,x_w,y_w,z_w,bc_hw_g,blocked_v_cell_at,blocked_w_cell_at,cylindrical,ox_e,mu_g,v_g,w_g,odx_e,axy_u,ox,vol_u,odx,odz,axz_u,wall_w_at,jm1,k_of,cartesian_grid,lambda_gt,trd_g,ayz,u_g,ayz_u,i_of,ip1,j_of,km1,cg_safe_mode,mu_gt,ep_g)
      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
         IF ( .NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN
            IP = IP1(I)
            J = J_OF(IJK)
            JM = JM1(J)
            K = K_OF(IJK)
            KM = KM1(K)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKNE = EAST_OF(IJKN)
            IJKS = SOUTH_OF(IJK)
            IJKSE = EAST_OF(IJKS)
            IPJMK = JM_OF(IPJK)
            IJMK = JM_OF(IJK)
            IJKT = TOP_OF(IJK)
            IJKTE = EAST_OF(IJKT)
            IJKB = BOTTOM_OF(IJK)
            IJKBE = EAST_OF(IJKB)
            IJKM = KM_OF(IJK)
            IPJKM = IP_OF(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(3)==1)) THEN
!
!       Surface forces
!
!         bulk viscosity term
               SBV = (LAMBDA_GT(IJKE)*TRD_G(IJKE)-LAMBDA_GT(IJK)*TRD_G(IJK))*AYZ(&
                  IJK)
!
!         shear stress terms
               SSX = MU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ODX(IP)*AYZ_U(IJK) - MU_GT(&
                  IJK)*(U_G(IJK)-U_G(IMJK))*ODX(I)*AYZ_U(IMJK)
               SSY = AVG_X_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE)&
                  ,MU_GT(IJKNE),J),I)*(V_G(IPJK)-V_G(IJK))*ODX_E(I)*AXZ_U(IJK) - &
                  AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(IJKSE),&
                  MU_GT(IJKE),JM),I)*(V_G(IPJMK)-V_G(IJMK))*ODX_E(I)*AXZ_U(IJMK)
               EPMU_GTE = AVG_X_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)
               EPMU_GBE = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT&
                  (IJKBE),MU_GT(IJKE),KM),I)
               SSZ = EPMU_GTE*(W_G(IPJK)-W_G(IJK))*ODX_E(I)*AXY_U(IJK) - EPMU_GBE*&
                  (W_G(IPJKM)-W_G(IJKM))*ODX_E(I)*AXY_U(IJKM)

            ELSE ! CARTESIAN GRID CASE

!         bulk viscosity term
               SBV =  (LAMBDA_GT(IJKE)*TRD_G(IJKE)) * AYZ_U(IJK) &
                     -(LAMBDA_GT(IJK) *TRD_G(IJK) ) * AYZ_U(IMJK)

!         shear stress terms
               IF(.NOT.CUT_U_CELL_AT(IJK))   THEN

                  SSX = MU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ONEoDX_E_U(IJK)*AYZ_U(IJK) - MU_GT(&
                  IJK)*(U_G(IJK)-U_G(IMJK))*ONEoDX_E_U(IMJK)*AYZ_U(IMJK)

                  SSY = AVG_X_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE)&
                  ,MU_GT(IJKNE),J),I)*(V_G(IPJK)-V_G(IJK))*ONEoDX_E_V(IJK)*AXZ_U(IJK) - &
                  AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(IJKSE),&
                  MU_GT(IJKE),JM),I)*(V_G(IPJMK)-V_G(IJMK))*ONEoDX_E_V(IJMK)*AXZ_U(IJMK)

                  IF(DO_K) THEN
                     EPMU_GTE = AVG_X_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                     IJKE),MU_GT(IJKTE),K),I)
                     EPMU_GBE = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT&
                     (IJKBE),MU_GT(IJKE),KM),I)

                     SSZ = EPMU_GTE*(W_G(IPJK)-W_G(IJK))*ONEoDX_E_W(IJK)*AXY_U(IJK) - EPMU_GBE*&
                     (W_G(IPJKM)-W_G(IJKM))*ONEoDX_E_W(IJKM)*AXY_U(IJKM)

                  ELSE
                     SSZ = ZERO
                  ENDIF

               ELSE

                  BCV = BC_U_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)

                     CASE ('CG_NSW')
                        CUT_TAU_UG = .TRUE.
                        NOC_UG     = .TRUE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     CASE ('CG_FSW')
                        CUT_TAU_UG = .FALSE.
                        NOC_UG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     CASE('CG_PSW')
                        IF(BC_HW_G(BC_U_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                           CUT_TAU_UG = .TRUE.
                           NOC_UG     = .TRUE.
                           UW_g = BC_UW_G(BCV)
                           VW_g = BC_VW_G(BCV)
                           WW_g = BC_WW_G(BCV)
                        ELSEIF(BC_HW_G(BC_U_ID(IJK))==ZERO) THEN   ! same as FSW
                           CUT_TAU_UG = .FALSE.
                           NOC_UG     = .FALSE.
                        ELSE                              ! partial slip
                           CUT_TAU_UG = .FALSE.
                           NOC_UG     = .FALSE.
                           UW_g = ZERO
                           VW_g = ZERO
                           WW_g = ZERO
                        ENDIF
                     CASE ('NONE')
                        TAU_U_G(IJK) = ZERO
                        CYCLE
                  END SELECT

                  IF(CUT_TAU_UG) THEN
                     MU_GT_CUT =  (VOL(IJK)*MU_GT(IJK) + VOL(IPJK)*MU_GT(IJKE))/(VOL(IJK) + VOL(IPJK))
                  ELSE
                     MU_GT_CUT = ZERO
                  ENDIF

!           SSX:

                  CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U(IJK),Y_U(IJK),Z_U(IJK),Del_H,Nx,Ny,Nz)

                  SSX = MU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ONEoDX_E_U(IJK)*AYZ_U(IJK) - MU_GT(&
                        IJK)*(U_G(IJK)-U_G(IMJK))*ONEoDX_E_U(IMJK)*AYZ_U(IMJK) &
                      - MU_GT_CUT * (U_g(IJK) - UW_g) / DEL_H * (Nx**2) * Area_U_CUT(IJK)

!           SSY:

                  V_NODE_AT_NE = ((.NOT.BLOCKED_V_CELL_AT(IPJK)).AND.(.NOT.WALL_V_AT(IPJK)))
                  V_NODE_AT_NW = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.(.NOT.WALL_V_AT(IJK)))
                  V_NODE_AT_SE = ((.NOT.BLOCKED_V_CELL_AT(IPJMK)).AND.(.NOT.WALL_V_AT(IPJMK)))
                  V_NODE_AT_SW = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

                  IF(V_NODE_AT_NE.AND.V_NODE_AT_NW) THEN

                     Vi = HALF * (V_G(IPJK) + V_G(IJK))
                     Xi = HALF * (X_V(IPJK) + X_V(IJK))
                     Yi = HALF * (Y_V(IPJK) + Y_V(IJK))
                     Zi = HALF * (Z_V(IPJK) + Z_V(IJK))
                     Sx = X_V(IPJK) - X_V(IJK)
                     Sy = Y_V(IPJK) - Y_V(IJK)
                     Sz = Z_V(IPJK) - Z_V(IJK)

                     CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     dvdx_at_N =  (V_G(IPJK) - V_G(IJK)) * ONEoDX_E_V(IJK)

                     IF(NOC_UG) dvdx_at_N = dvdx_at_N - ((Vi-VW_g)*ONEoDX_E_V(IJK)/DEL_H*(Sy*Ny+Sz*Nz))

                  ELSE
                     dvdx_at_N =  ZERO
                  ENDIF

                  IF(V_NODE_AT_SE.AND.V_NODE_AT_SW) THEN

                     Vi = HALF * (V_G(IPJMK) + V_G(IJMK))
                     Xi = HALF * (X_V(IPJMK) + X_V(IJMK))
                     Yi = HALF * (Y_V(IPJMK) + Y_V(IJMK))
                     Zi = HALF * (Z_V(IPJMK) + Z_V(IJMK))
                     Sx = X_V(IPJMK) - X_V(IJMK)
                     Sy = Y_V(IPJMK) - Y_V(IJMK)
                     Sz = Z_V(IPJMK) - Z_V(IJMK)
!
                     CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     dvdx_at_S =  (V_G(IPJMK) - V_G(IJMK)) * ONEoDX_E_V(IJMK)

                     IF(NOC_UG) dvdx_at_S = dvdx_at_S - ((Vi-VW_g)*ONEoDX_E_V(IJMK)/DEL_H*(Sy*Ny+Sz*Nz))

                  ELSE
                     dvdx_at_S =  ZERO
                  ENDIF

                    IF(V_NODE_AT_NW) THEN
                       CALL GET_DEL_H(IJK,'U_MOMENTUM',X_V(IJK),Y_V(IJK),Z_V(IJK),Del_H,Nx,Ny,Nz)
                       SSY_CUT = - MU_GT_CUT * (V_G(IJK) - VW_g) / DEL_H * (Nx*Ny) * Area_U_CUT(IJK)
                    ELSE
                       SSY_CUT =  ZERO
                    ENDIF


                  SSY = AVG_X_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE)&
                     ,MU_GT(IJKNE),J),I)*dvdx_at_N*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(IJKSE),&
                     MU_GT(IJKE),JM),I)*dvdx_at_S*AXZ_U(IJMK) &
                    + SSY_CUT

!           SSZ:

                  IF(DO_K) THEN

                     W_NODE_AT_TE = ((.NOT.BLOCKED_W_CELL_AT(IPJK)).AND.(.NOT.WALL_W_AT(IPJK)))
                     W_NODE_AT_TW = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))
                     W_NODE_AT_BE = ((.NOT.BLOCKED_W_CELL_AT(IPJKM)).AND.(.NOT.WALL_W_AT(IPJKM)))
                     W_NODE_AT_BW = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                     IF(W_NODE_AT_TE.AND.W_NODE_AT_TW) THEN

                        Wi = HALF * (W_G(IPJK) + W_G(IJK))
                        Xi = HALF * (X_W(IPJK) + X_W(IJK))
                        Yi = HALF * (Y_W(IPJK) + Y_W(IJK))
                        Zi = HALF * (Z_W(IPJK) + Z_W(IJK))
                        Sx = X_W(IPJK) - X_W(IJK)
                        Sy = Y_W(IPJK) - Y_W(IJK)
                        Sz = Z_W(IPJK) - Z_W(IJK)

                        CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                        dwdx_at_T =  (W_G(IPJK) - W_G(IJK)) * ONEoDX_E_W(IJK)

                        IF(NOC_UG) dwdx_at_T = dwdx_at_T - ((Wi-WW_g) * ONEoDX_E_W(IJK)/DEL_H*(Sy*Ny+Sz*Nz))

                     ELSE
                        dwdx_at_T =  ZERO
                     ENDIF


                     IF(W_NODE_AT_BE.AND.W_NODE_AT_BW) THEN

                        Wi = HALF * (W_G(IPJKM) + W_G(IJKM))
                        Xi = HALF * (X_W(IPJKM) + X_W(IJKM))
                        Yi = HALF * (Y_W(IPJKM) + Y_W(IJKM))
                        Zi = HALF * (Z_W(IPJKM) + Z_W(IJKM))
                        Sx = X_W(IPJKM) - X_W(IJKM)
                        Sy = Y_W(IPJKM) - Y_W(IJKM)
                        Sz = Z_W(IPJKM) - Z_W(IJKM)

                        CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                        dwdx_at_B =  (W_G(IPJKM) - W_G(IJKM)) * ONEoDX_E_W(IJKM)

                        IF(NOC_UG) dwdx_at_B = dwdx_at_B - ((Wi-WW_g) * ONEoDX_E_W(IJKM)/DEL_H*(Sy*Ny+Sz*Nz))

                     ELSE
                        dwdx_at_B =  ZERO
                     ENDIF

                     IF(W_NODE_AT_TW) THEN
                        CALL GET_DEL_H(IJK,'U_MOMENTUM',X_W(IJK),Y_W(IJK),Z_W(IJK),Del_H,Nx,Ny,Nz)
                        SSZ_CUT = - MU_GT_CUT * (W_G(IJK) - WW_g) / DEL_H * (Nx*Nz) * Area_U_CUT(IJK)
                     ELSE
                        SSZ_CUT =  ZERO
                     ENDIF

                     EPMU_GTE = AVG_X_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                        IJKE),MU_GT(IJKTE),K),I)
                     EPMU_GBE = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT&
                        (IJKBE),MU_GT(IJKE),KM),I)
                     SSZ =   EPMU_GTE*dwdx_at_T*AXY_U(IJK)  &
                           - EPMU_GBE*dwdx_at_B*AXY_U(IJKM) &
                           + SSZ_CUT

                  ELSE

                    SSZ = ZERO

                  ENDIF ! DO_K

               ENDIF  ! CUT_CELL

!  Original terms
!            SSX = MU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ODX(IP)*AYZ_U(IJK) - MU_GT(&
!               IJK)*(U_G(IJK)-U_G(IMJK))*ODX(I)*AYZ_U(IMJK)

!            SSY = AVG_X_H(AVG_Y_H(MU_GT(IJK),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE)&
!               ,MU_GT(IJKNE),J),I)*(V_G(IPJK)-V_G(IJK))*ODX_E(I)*AXZ_U(IJK) - &
!               AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJK),JM),AVG_Y_H(MU_GT(IJKSE),&
!               MU_GT(IJKE),JM),I)*(V_G(IPJMK)-V_G(IJMK))*ODX_E(I)*AXZ_U(IJMK)


!            EPMU_GTE = AVG_X_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
!               IJKE),MU_GT(IJKTE),K),I)
!            EPMU_GBE = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT&
!               (IJKBE),MU_GT(IJKE),KM),I)
!            SSZ = EPMU_GTE*(W_G(IPJK)-W_G(IJK))*ODX_E(I)*AXY_U(IJK) - EPMU_GBE*&
!               (W_G(IPJKM)-W_G(IJKM))*ODX_E(I)*AXY_U(IJKM)

            ENDIF ! CARTESIAN GRID
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN
!
!           modify Ssz: integral of (1/x) (d/dz) (mu*(-w/x))
               SSZ = SSZ - (EPMU_GTE*(HALF*(W_G(IPJK)+W_G(IJK))*OX_E(I))*AXY_U(&
                  IJK)-EPMU_GBE*(HALF*(W_G(IPJKM)+W_G(IJKM))*OX_E(I))*AXY_U(&
                  IJKM))
!
!           -(2mu/x)*(1/x)*dw/dz part of Tau_zz/X
               EPMUGA = AVG_X(MU_G(IJK),MU_G(IJKE),I)
               DWOXDZ = HALF*((W_G(IJK)-W_G(IJKM))*OX(I)*ODZ(K)+(W_G(IPJK)-W_G(&
                  IPJKM))*OX(IP)*ODZ(K))
               VTZB = -2.d0*EPMUGA*OX_E(I)*DWOXDZ
            ELSE
               VTZB = ZERO
            ENDIF
!
!         Add the terms
            TAU_U_G(IJK) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK)
         ELSE
            TAU_U_G(IJK) = ZERO
         ENDIF
      END DO

      call send_recv(tau_u_g,2)

      RETURN

      END SUBROUTINE CALC_TAU_U_G

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
