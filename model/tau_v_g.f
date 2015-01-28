!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_V_g(TAU_V_g, IER)                             C
!  Purpose: Cross terms in the gradient of stress in V_g momentum      c
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
      SUBROUTINE CALC_TAU_V_G(TAU_V_G)
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
      USE sendrecv
      USE compar
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
!                      TAU_V_g
      DOUBLE PRECISION TAU_V_g(DIMENSION_3)
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKN, JP, IM,  KM, IJPK, IJMK,&
                       IJKE, IJKNE, IJKW, IJKNW, IMJPK, IMJK, IJKT,&
                       IJKTN, IJKB, IJKBN, IJKM, IJPKM
!
!                      Average volume fraction
      DOUBLE PRECISION EPGA
!
!                      Source terms (Surface)
      DOUBLE PRECISION Sbv, Ssx, Ssy, Ssz
!
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: U_NODE_AT_NE,U_NODE_AT_NW,U_NODE_AT_SE,U_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TN,W_NODE_AT_TS,W_NODE_AT_BN,W_NODE_AT_BS
      DOUBLE PRECISION :: U_SUM,X_SUM,Y_SUM,Z_SUM
      DOUBLE PRECISION :: dudy_at_E,dudy_at_W
      DOUBLE PRECISION :: dwdy_at_T,dwdy_at_B
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: MU_GT_CUT,SSX_CUT,SSZ_CUT
      DOUBLE PRECISION :: UW_g,VW_g,WW_g
      INTEGER :: N_SUM
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!!!!$omp  parallel do private( IJK, I, IJKE, EPGA,  J,  K, KM,  &
!!!!$omp&  IMJK,IJKN,IJKNE,IJMK,IJKT,  &
!!!!$omp&  JP,IM,IJPK,IJKW,IJKNW,IMJPK,IJKTN,IJKBN,IJPKM,IJKB,IJKM, &
!!!!$omp&  SBV,  SSX,SSY,   SSZ)  &
!!!!$omp&  schedule(static)
      DO IJK = IJKSTART3, IJKEND3
         J = J_OF(IJK)
         IJKN = NORTH_OF(IJK)
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J)
         IF ( .NOT.IP_AT_N(IJK) .AND. EPGA>DIL_EP_S) THEN
            JP = JP1(J)
            I = I_OF(IJK)
            IM = IM1(I)
            K = K_OF(IJK)
            KM = KM1(K)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKNE = EAST_OF(IJKN)
            IJKW = WEST_OF(IJK)
            IJKNW = NORTH_OF(IJKW)
            IMJPK = IM_OF(IJPK)
            IMJK = IM_OF(IJK)
            IJKT = TOP_OF(IJK)
            IJKTN = NORTH_OF(IJKT)
            IJKB = BOTTOM_OF(IJK)
            IJKBN = NORTH_OF(IJKB)
            IJKM = KM_OF(IJK)
            IJPKM = JP_OF(IJKM)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(4)==1)) THEN
!
!       Surface forces
!
!         bulk viscosity term
               SBV = (LAMBDA_GT(IJKN)*TRD_G(IJKN)-LAMBDA_GT(IJK)*TRD_G(IJK))*AXZ(&
                  IJK)
!
!         shear stress terms
               SSX = AVG_Y_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN)&
                  ,MU_GT(IJKNE),I),J)*(U_G(IJPK)-U_G(IJK))*ODY_N(J)*AYZ_V(IJK) - &
                  AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(IJKNW),&
                  MU_GT(IJKN),IM),J)*(U_G(IMJPK)-U_G(IMJK))*ODY_N(J)*AYZ_V(IMJK)
               SSY = MU_GT(IJKN)*(V_G(IJPK)-V_G(IJK))*ODY(JP)*AXZ_V(IJK) - MU_GT(&
                  IJK)*(V_G(IJK)-V_G(IJMK))*ODY(J)*AXZ_V(IJMK)
               SSZ = AVG_Y_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(IJKN)&
                  ,MU_GT(IJKTN),K),J)*(W_G(IJPK)-W_G(IJK))*ODY_N(J)*AXY_V(IJK) - &
                  AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT(IJKBN),&
                  MU_GT(IJKN),KM),J)*(W_G(IJPKM)-W_G(IJKM))*ODY_N(J)*AXY_V(IJKM)

            ELSE  ! CARTESIAN GRID CASE
!
!       Surface forces
!
!         bulk viscosity term
               SBV =  (LAMBDA_GT(IJKN)*TRD_G(IJKN)) * AXZ_V(IJK) &
                     -(LAMBDA_GT(IJK) *TRD_G(IJK) ) * AXZ_V(IJMK)
!
!         shear stress terms

               IF(.NOT.CUT_V_CELL_AT(IJK))   THEN

                  SSX = AVG_Y_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN)&
                     ,MU_GT(IJKNE),I),J)*(U_G(IJPK)-U_G(IJK))*ONEoDY_N_U(IJK)*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(IJKNW),&
                     MU_GT(IJKN),IM),J)*(U_G(IMJPK)-U_G(IMJK))*ONEoDY_N_U(IMJK)*AYZ_V(IMJK)

                  SSY = MU_GT(IJKN)*(V_G(IJPK)-V_G(IJK))*ONEoDY_N_V(IJK)*AXZ_V(IJK) - MU_GT(&
                     IJK)*(V_G(IJK)-V_G(IJMK))*ONEoDY_N_V(IJMK)*AXZ_V(IJMK)



                  IF(DO_K) THEN

                     SSZ = AVG_Y_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(IJKN)&
                        ,MU_GT(IJKTN),K),J)*(W_G(IJPK)-W_G(IJK))*ONEoDY_N_W(IJK)*AXY_V(IJK) - &
                        AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT(IJKBN),&
                        MU_GT(IJKN),KM),J)*(W_G(IJPKM)-W_G(IJKM))*ONEoDY_N_W(IJKM)*AXY_V(IJKM)

                  ELSE
                     SSZ = ZERO
                  ENDIF

               ELSE  ! CUT CELL

                  BCV = BC_V_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)
                     CASE ('CG_NSW')
                        CUT_TAU_VG = .TRUE.
                        NOC_VG     = .TRUE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     CASE ('CG_FSW')
                        CUT_TAU_VG = .FALSE.
                        NOC_VG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     CASE('CG_PSW')
                        IF(BC_HW_G(BC_V_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                           CUT_TAU_VG = .TRUE.
                           NOC_VG     = .TRUE.
                           UW_g = BC_UW_G(BCV)
                           VW_g = BC_VW_G(BCV)
                           WW_g = BC_WW_G(BCV)
                        ELSEIF(BC_HW_G(BC_V_ID(IJK))==ZERO) THEN   ! same as FSW
                           CUT_TAU_VG = .FALSE.
                           NOC_VG     = .FALSE.
                           UW_g = ZERO
                           VW_g = ZERO
                           WW_g = ZERO
                        ELSE                              ! partial slip
                           CUT_TAU_VG = .FALSE.
                           NOC_VG     = .FALSE.
                        ENDIF
                     CASE ('NONE')
                        TAU_V_G(IJK) = ZERO
                        CYCLE
                  END SELECT


                  IF(CUT_TAU_VG) THEN
                     MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJPK)*MU_GT(IJKN))/(VOL(IJK) + VOL(IJPK))
                  ELSE
                     MU_GT_CUT = ZERO
                  ENDIF

!           SSX:

                  U_NODE_AT_NE = ((.NOT.BLOCKED_U_CELL_AT(IJPK)).AND.(.NOT.WALL_U_AT(IJPK)))
                  U_NODE_AT_SE = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.(.NOT.WALL_U_AT(IJK)))
                  U_NODE_AT_NW = ((.NOT.BLOCKED_U_CELL_AT(IMJPK)).AND.(.NOT.WALL_U_AT(IMJPK)))
                  U_NODE_AT_SW = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))


                  IF(U_NODE_AT_NE.AND.U_NODE_AT_SE) THEN

                     Ui = HALF * (U_G(IJPK) + U_G(IJK))
                     Xi = HALF * (X_U(IJPK) + X_U(IJK))
                     Yi = HALF * (Y_U(IJPK) + Y_U(IJK))
                     Zi = HALF * (Z_U(IJPK) + Z_U(IJK))
                     Sx = X_U(IJPK) - X_U(IJK)
                     Sy = Y_U(IJPK) - Y_U(IJK)
                     Sz = Z_U(IJPK) - Z_U(IJK)


                     CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     dudy_at_E =  (U_G(IJPK) - U_G(IJK)) * ONEoDY_N_U(IJK)

                     IF(NOC_VG) dudy_at_E = dudy_at_E - ((Ui-UW_g) * ONEoDY_N_U(IJK)/DEL_H*(Sx*Nx+Sz*Nz))

                  ELSE
                     dudy_at_E =  ZERO
                  ENDIF

                  IF(U_NODE_AT_NW.AND.U_NODE_AT_SW) THEN

                     Ui = HALF * (U_G(IMJPK) + U_G(IMJK))
                     Xi = HALF * (X_U(IMJPK) + X_U(IMJK))
                     Yi = HALF * (Y_U(IMJPK) + Y_U(IMJK))
                     Zi = HALF * (Z_U(IMJPK) + Z_U(IMJK))
                     Sx = X_U(IMJPK) - X_U(IMJK)
                     Sy = Y_U(IMJPK) - Y_U(IMJK)
                     Sz = Z_U(IMJPK) - Z_U(IMJK)

                     CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     dudy_at_W =  (U_G(IMJPK) - U_G(IMJK)) * ONEoDY_N_U(IMJK)

                     IF(NOC_VG) dudy_at_W = dudy_at_W - ((Ui-UW_g) * ONEoDY_N_U(IMJK)/DEL_H*(Sx*Nx+Sz*Nz))

                  ELSE
                     dudy_at_W =  ZERO
                  ENDIF

                  U_SUM = ZERO
                  X_SUM = ZERO
                  Y_SUM = ZERO
                  Z_SUM = ZERO
                  N_SUM = 0

                 IF(U_NODE_AT_SE) THEN
                    CALL GET_DEL_H(IJK,'V_MOMENTUM',X_U(IJK),Y_U(IJK),Z_U(IJK),Del_H,Nx,Ny,Nz)
                    SSX_CUT = - MU_GT_CUT * (U_G(IJK) - UW_g) / DEL_H * (Ny*Nx) * Area_V_CUT(IJK)
                 ELSE
                    SSX_CUT =  ZERO
                 ENDIF

                  SSX = AVG_Y_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN)&
                     ,MU_GT(IJKNE),I),J)*dudy_at_E*AYZ_V(IJK) - &
                     AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(IJKNW),&
                     MU_GT(IJKN),IM),J)*dudy_at_W*AYZ_V(IMJK) &
                    + SSX_CUT

!           SSY:
                  CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V(IJK),Y_V(IJK),Z_V(IJK),Del_H,Nx,Ny,Nz)

                  SSY = MU_GT(IJKN)*(V_G(IJPK)-V_G(IJK))*ONEoDY_N_V(IJK)*AXZ_V(IJK) - MU_GT(&
                        IJK)*(V_G(IJK)-V_G(IJMK))*ONEoDY_N_V(IJMK)*AXZ_V(IJMK) &
                      - MU_GT_CUT * (V_g(IJK) - VW_g) / DEL_H * (Ny**2) * Area_V_CUT(IJK)

!           SSZ:
                  IF(DO_K) THEN

                     W_NODE_AT_TN = ((.NOT.BLOCKED_W_CELL_AT(IJPK)).AND.(.NOT.WALL_W_AT(IJPK)))
                     W_NODE_AT_TS = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))
                     W_NODE_AT_BN = ((.NOT.BLOCKED_W_CELL_AT(IJPKM)).AND.(.NOT.WALL_W_AT(IJPKM)))
                     W_NODE_AT_BS = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                     IF(W_NODE_AT_TN.AND.W_NODE_AT_TS) THEN

                        Wi = HALF * (W_G(IJPK) + W_G(IJK))
                        Xi = HALF * (X_W(IJPK) + X_W(IJK))
                        Yi = HALF * (Y_W(IJPK) + Y_W(IJK))
                        Zi = HALF * (Z_W(IJPK) + Z_W(IJK))
                        Sx = X_W(IJPK) - X_W(IJK)
                        Sy = Y_W(IJPK) - Y_W(IJK)
                        Sz = Z_W(IJPK) - Z_W(IJK)

                        CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                        dwdy_at_T =  (W_G(IJPK) - W_G(IJK)) * ONEoDY_N_W(IJK)

                        IF(NOC_VG) dwdy_at_T = dwdy_at_T - ((Wi-WW_g) * ONEoDY_N_W(IJK)/DEL_H*(Sx*Nx+Sz*Nz))

                     ELSE
                        dwdy_at_T =  ZERO
                     ENDIF

                     IF(W_NODE_AT_BN.AND.W_NODE_AT_BS) THEN

                        Wi = HALF * (W_G(IJPKM) + W_G(IJKM))
                        Xi = HALF * (X_W(IJPKM) + X_W(IJKM))
                        Yi = HALF * (Y_W(IJPKM) + Y_W(IJKM))
                        Zi = HALF * (Z_W(IJPKM) + Z_W(IJKM))
                        Sx = X_W(IJPKM) - X_W(IJKM)
                        Sy = Y_W(IJPKM) - Y_W(IJKM)
                        Sz = Z_W(IJPKM) - Z_W(IJKM)

                        CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                        dwdy_at_B =  (W_G(IJPKM) - W_G(IJKM)) * ONEoDY_N_W(IJKM)

                        IF(NOC_VG) dwdy_at_B = dwdy_at_B - ((Wi-WW_g) * ONEoDY_N_W(IJKM)/DEL_H*(Sx*Nx+Sz*Nz))

                     ELSE
                        dwdy_at_B =  ZERO
                     ENDIF

                     IF(W_NODE_AT_TS) THEN
                        CALL GET_DEL_H(IJK,'V_MOMENTUM',X_W(IJK),Y_W(IJK),Z_W(IJK),Del_H,Nx,Ny,Nz)
                        SSZ_CUT = - MU_GT_CUT * (W_G(IJK) - WW_g) / DEL_H * (Ny*Nz) * Area_V_CUT(IJK)
                     ELSE
                        SSZ_CUT =  ZERO
                     ENDIF

                     SSZ = AVG_Y_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(IJKN)&
                        ,MU_GT(IJKTN),K),J)*dwdy_at_T*AXY_V(IJK) - &
                        AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT(IJKBN),&
                        MU_GT(IJKN),KM),J)*dwdy_at_B*AXY_V(IJKM) &
                       + SSZ_CUT

                  ELSE

                     SSZ = ZERO

                  ENDIF  ! DO_K


               END IF  ! CUT CELL


! Original terms
!            SSX = AVG_Y_H(AVG_X_H(MU_GT(IJK),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN)&
!              ,MU_GT(IJKNE),I),J)*(U_G(IJPK)-U_G(IJK))*ODY_N(J)*AYZ_V(IJK) - &
!               AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJK),IM),AVG_X_H(MU_GT(IJKNW),&
!               MU_GT(IJKN),IM),J)*(U_G(IMJPK)-U_G(IMJK))*ODY_N(J)*AYZ_V(IMJK)

!            SSY = MU_GT(IJKN)*(V_G(IJPK)-V_G(IJK))*ODY(JP)*AXZ_V(IJK) - MU_GT(&
!               IJK)*(V_G(IJK)-V_G(IJMK))*ODY(J)*AXZ_V(IJMK)

!            SSZ = AVG_Y_H(AVG_Z_H(MU_GT(IJK),MU_GT(IJKT),K),AVG_Z_H(MU_GT(IJKN)&
!              ,MU_GT(IJKTN),K),J)*(W_G(IJPK)-W_G(IJK))*ODY_N(J)*AXY_V(IJK) - &
!               AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJK),KM),AVG_Z_H(MU_GT(IJKBN),&
!               MU_GT(IJKN),KM),J)*(W_G(IJPKM)-W_G(IJKM))*ODY_N(J)*AXY_V(IJKM)

            ENDIF  ! CARTESIAN GRID
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!
!
!         Add the terms
            TAU_V_G(IJK) = SBV + SSX + SSY + SSZ
         ELSE
            TAU_V_G(IJK) = ZERO
         ENDIF
      END DO
      call send_recv(tau_v_g,2)
      RETURN
      END SUBROUTINE CALC_TAU_V_G

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
