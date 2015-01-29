!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_W_s                                            C
!  Purpose: Cross terms in the gradient of stress in W_s momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_TAU_W_S(TAU_W_S)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: zero, half, one
      USE fun_avg

! number of solids phases
      USE physprop, only: smax, mmax

! x,y,z-components of solids velocity
      USE fldvar, only: u_s, v_s, w_s
! solids phase particle bulk density and material density
      USE fldvar, only: rop_s, ro_s, ep_s

! viscous solids transport coefficients
      USE visc_s, only: mu_s, lambda_s
! trace of rate of strain tensor
      USE visc_s, only: trd_s

! dilute threshold
      USE toleranc, only: dil_ep_s

! kinetic theories
      USE run, only: kt_type_enum, ghd_2007
! runtime flag for treating the system as if shearing
      USE run, only: shear

! primarily needed for function.inc
      USE compar
      USE geometry
      USE indices

! for sendrecv calls
      USE sendrecv

! for cutcell:
! wall velocity for slip conditions
      USE bc, only: bc_hw_s, bc_uw_s, bc_vw_s, bc_ww_s
      USE bc, only: bc_type
      USE quadric
      USE cutcell

      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! TAU_W_s
      DOUBLE PRECISION, INTENT(OUT) :: TAU_W_s(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, J, I, IM, IJKP, IMJK, IJKN, IJKNT, IJKS,&
                 IJKST, IJMKP, IJMK, IJKE, IJKTE, IJKW, IJKTW,&
                 IMJKP, K, IJKT, JM, KP, IJKM
! Phase index
      INTEGER :: M, L
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp
! Average velocity gradients
      DOUBLE PRECISION :: duodz
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Source terms (Volumetric)
      DOUBLE PRECISION :: Vxz

! for cartesian grid implementation:
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: U_NODE_AT_ET,U_NODE_AT_EB,U_NODE_AT_WT,U_NODE_AT_WB
      LOGICAL :: V_NODE_AT_NT,V_NODE_AT_NB,V_NODE_AT_ST,V_NODE_AT_SB

      DOUBLE PRECISION :: dudz_at_E,dudz_at_W
      DOUBLE PRECISION :: dvdz_at_N,dvdz_at_S
      DOUBLE PRECISION :: MU_S_CUT,SSX_CUT,SSY_CUT
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_s,VW_s,WW_s
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!-----------------------------------------------

      DO M = 1, MMAX
        IF(KT_TYPE_ENUM == GHD_2007 .AND. M /= MMAX) CYCLE


!!$omp  parallel do private( IJK, I, IJKE, EPSA, EPStmp  J,  K,   &
!!$omp&  JM,IJKP,IJKNT,IJKS,IJKST,IJMKP,IJKTE,IJKTW,IMJKP,KP, &
!!$omp&  DUODZ,VXZ, &
!!$omp&  IM,IJKW, IPJK,&
!!$omp&  IMJK,IJKN,IJMK,IJKT,  &
!!$omp&  IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ) &
!!$omp&  schedule(static)

        DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
          IF(WALL_AT(IJK)) cycle

          K = K_OF(IJK)
          IJKT = TOP_OF(IJK)

          IF (KT_TYPE_ENUM == GHD_2007) THEN
             EPStmp = ZERO
             DO L = 1, SMAX
                EPStmp = EPStmp + AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
             ENDDO
             EPSA = EPStmp
          ELSE
             EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
          ENDIF

          IF ( .NOT.SIP_AT_T(IJK) .AND. EPSA>DIL_EP_S) THEN
            J = J_OF(IJK)
            I = I_OF(IJK)
            IM = IM1(I)
            JM = JM1(J)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKNT = TOP_OF(IJKN)
            IJKS = SOUTH_OF(IJK)
            IJKST = TOP_OF(IJKS)
            IJMKP = JM_OF(IJKP)
            IJMK = JM_OF(IJK)
            IJKE = EAST_OF(IJK)
            IJKTE = EAST_OF(IJKT)
            IJKW = WEST_OF(IJK)
            IJKTW = WEST_OF(IJKT)
            IMJKP = KP_OF(IMJK)
            KP = KP1(K)
            IJKM = KM_OF(IJK)

            IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(5)==1)) THEN
! NON CARTESIAN GRID CASE
! ---------------------------------------------------------------->>>
! Surface forces

! bulk viscosity term
              SBV = (LAMBDA_S(IJKT,M)*TRD_S(IJKT,M)-&
                     LAMBDA_S(IJK,M)*TRD_S(IJK,M))*AXY(IJK)

! shear stress terms
              SSX = AVG_Z_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                            AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                    (U_S(IJKP,M)-U_S(IJK,M))*OX_E(I)*ODZ_T(K)*AYZ_W(IJK) - &
                    AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                            AVG_X_H(MU_S(IJKTW,M),MU_S(IJKT,M),IM),K)*&
                    (U_S(IMJKP,M)-U_S(IMJK,M))*ODZ_T(K)*DY(J)*(HALF*(DZ(K)+DZ(KP)))
                    !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity
              SSY = AVG_Z_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),&
                            AVG_Y_H(MU_S(IJKT,M),MU_S(IJKNT,M),J),K)*&
                    (V_S(IJKP,M)-V_S(IJK,M))*OX(I)*ODZ_T(K)*AXZ_W(IJK) - &
                    AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
                            AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                    (V_S(IJMKP,M)-V_S(IJMK,M))*OX(I)*ODZ_T(K)*AXZ_W(IJMK)
              SSZ = MU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
                    OX(I)*ODZ(KP)*AXY_W(IJK) - &
                    MU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
                    OX(I)*ODZ(K)* AXY_W(IJKM)
! ----------------------------------------------------------------<<<

            ELSE

! CARTESIAN GRID CASE
! ---------------------------------------------------------------->>>
! Surface forces

! bulk viscosity term
              SBV =  (LAMBDA_S(IJKT,M)*TRD_S(IJKT,M)) * AXY_W(IJK) - &
                     (LAMBDA_S(IJK,M) *TRD_S(IJK,M) ) * AXY_W(IJKM)

! shear stress terms
              IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
                SSX = AVG_Z_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                      (U_S(IJKP,M)-U_S(IJK,M))*ONEoDZ_T_U(IJK)*AYZ_W(IJK) - &
                      AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                              AVG_X_H(MU_S(IJKTW,M),MU_S(IJKT,M),IM),K)*&
                      (U_S(IMJKP,M)-U_S(IMJK,M))*ONEoDZ_T_U(IMJK)*AYZ_W(IMJK)
                SSY = AVG_Z_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),&
                              AVG_Y_H(MU_S(IJKT,M),MU_S(IJKNT,M),J),K)*&
                      (V_S(IJKP,M)-V_S(IJK,M))*ONEoDZ_T_V(IJK)*AXZ_W(IJK) - &
                      AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
                              AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                      (V_S(IJMKP,M)-V_S(IJMK,M))*ONEoDZ_T_V(IJMK)*AXZ_W(IJMK)
                SSZ = MU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
                      ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
                      MU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
                      ONEoDZ_T_W(IJKM)*AXY_W(IJKM)

              ELSE   ! CUT_CELL
                BCV = BC_W_ID(IJK)
                IF(BCV > 0 ) THEN
                   BCT = BC_TYPE(BCV)
                ELSE
                   BCT = 'NONE'
                ENDIF

                SELECT CASE (BCT)
                  CASE ('CG_NSW')
                     CUT_TAU_WS = .TRUE.
                     NOC_WS     = .TRUE.
                     UW_s = ZERO
                     VW_s = ZERO
                     WW_s = ZERO
                  CASE ('CG_FSW')
                     CUT_TAU_WS = .FALSE.
                     NOC_WS     = .FALSE.
                     UW_s = ZERO
                     VW_s = ZERO
                     WW_s = ZERO
                  CASE('CG_PSW')
                     IF(BC_HW_S(BC_W_ID(IJK),M)==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_WS = .TRUE.
                        NOC_WS     = .TRUE.
                        UW_s = BC_UW_S(BCV,M)
                        VW_s = BC_VW_S(BCV,M)
                        WW_s = BC_WW_S(BCV,M)
                     ELSEIF(BC_HW_S(BC_W_ID(IJK),M)==ZERO) THEN   ! same as FSW
                        CUT_TAU_WS = .FALSE.
                        NOC_WS     = .FALSE.
                        UW_s = ZERO
                        VW_s = ZERO
                        WW_s = ZERO
                     ELSE                              ! partial slip
                        CUT_TAU_WS = .FALSE.
                        NOC_WS     = .FALSE.
                     ENDIF
                  CASE ('NONE')
                     TAU_W_S(IJK,M) = ZERO
                     CYCLE
                END SELECT

                IF(CUT_TAU_WS) THEN
                   MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + &
                      VOL(IJKP)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKP))
                ELSE
                   MU_S_CUT = ZERO
                ENDIF

! SSX:
                U_NODE_AT_ET = ((.NOT.BLOCKED_U_CELL_AT(IJKP)).AND.(.NOT.WALL_U_AT(IJKP)))
                U_NODE_AT_EB = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.(.NOT.WALL_U_AT(IJK)))
                U_NODE_AT_WT = ((.NOT.BLOCKED_U_CELL_AT(IMJKP)).AND.(.NOT.WALL_U_AT(IMJKP)))
                U_NODE_AT_WB = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))

                IF(U_NODE_AT_ET.AND.U_NODE_AT_EB) THEN
                   Ui = HALF * (U_S(IJKP,M) + U_S(IJK,M))
                   Xi = HALF * (X_U(IJKP) + X_U(IJK))
                   Yi = HALF * (Y_U(IJKP) + Y_U(IJK))
                   Zi = HALF * (Z_U(IJKP) + Z_U(IJK))
                   Sx = X_U(IJKP) - X_U(IJK)
                   Sy = Y_U(IJKP) - Y_U(IJK)
                   Sz = Z_U(IJKP) - Z_U(IJK)
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dudz_at_E =  (U_S(IJKP,M) - U_S(IJK,M)) * ONEoDZ_T_U(IJK)
                   IF(NOC_WS) dudz_at_E = dudz_at_E - ((Ui-UW_s) * &
                                 ONEoDZ_T_U(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
                ELSE
                   dudz_at_E =  ZERO
                ENDIF

                IF(U_NODE_AT_WT.AND.U_NODE_AT_WB) THEN
                   Ui = HALF * (U_S(IMJKP,M) + U_S(IMJK,M))
                   Xi = HALF * (X_U(IMJKP) + X_U(IMJK))
                   Yi = HALF * (Y_U(IMJKP) + Y_U(IMJK))
                   Zi = HALF * (Z_U(IMJKP) + Z_U(IMJK))
                   Sx = X_U(IMJKP) - X_U(IMJK)
                   Sy = Y_U(IMJKP) - Y_U(IMJK)
                   Sz = Z_U(IMJKP) - Z_U(IMJK)
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dudz_at_W =  (U_S(IMJKP,M) - U_S(IMJK,M)) * ONEoDZ_T_U(IMJK)
                   IF(NOC_WS) dudz_at_W = dudz_at_W - ((Ui-UW_s) * &
                                 ONEoDZ_T_U(IMJK)/DEL_H*(Sx*Nx+Sy*Ny))
                ELSE
                   dudz_at_W =  ZERO
                ENDIF

                IF(U_NODE_AT_EB) THEN
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',X_U(IJK),Y_U(IJK),&
                                  Z_U(IJK),Del_H,Nx,Ny,Nz)
                   SSX_CUT = - MU_S_CUT * (U_S(IJK,M) - UW_s) / DEL_H * &
                             (Nz*Nx) * Area_W_CUT(IJK)
                ELSE
                   SSX_CUT =  ZERO
                ENDIF

                SSX = AVG_Z_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)*&
                      dudz_at_E*AYZ_W(IJK) - &
                      AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                              AVG_X_H(MU_S(IJKTW,M),MU_S(IJKT,M),IM),K)*&
                      dudz_at_W*AYZ_W(IMJK) + SSX_CUT

! SSY:
                V_NODE_AT_NT = ((.NOT.BLOCKED_V_CELL_AT(IJKP)).AND.(.NOT.WALL_V_AT(IJKP)))
                V_NODE_AT_NB = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.(.NOT.WALL_V_AT(IJK)))
                V_NODE_AT_ST = ((.NOT.BLOCKED_V_CELL_AT(IJMKP)).AND.(.NOT.WALL_V_AT(IJMKP)))
                V_NODE_AT_SB = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

                IF(V_NODE_AT_NT.AND.V_NODE_AT_NB) THEN
                   Vi = HALF * (V_S(IJKP,M) + V_S(IJK,M))
                   Xi = HALF * (X_V(IJKP) + X_V(IJK))
                   Yi = HALF * (Y_V(IJKP) + Y_V(IJK))
                   Zi = HALF * (Z_V(IJKP) + Z_V(IJK))
                   Sx = X_V(IJKP) - X_V(IJK)
                   Sy = Y_V(IJKP) - Y_V(IJK)
                   Sz = Z_V(IJKP) - Z_V(IJK)
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dvdz_at_N =  (V_S(IJKP,M) - V_S(IJK,M)) * ONEoDZ_T_V(IJK)
                   IF(NOC_WS) dvdz_at_N = dvdz_at_N - ((Vi-VW_s) * &
                                 ONEoDZ_T_V(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
                ELSE
                   dvdz_at_N =  ZERO
                ENDIF

                IF(V_NODE_AT_ST.AND.V_NODE_AT_SB) THEN
                   Vi = HALF * (V_S(IJMKP,M) + V_S(IJMK,M))
                   Xi = HALF * (X_V(IJMKP) + X_V(IJMK))
                   Yi = HALF * (Y_V(IJMKP) + Y_V(IJMK))
                   Zi = HALF * (Z_V(IJMKP) + Z_V(IJMK))
                   Sx = X_V(IJMKP) - X_V(IJMK)
                   Sy = Y_V(IJMKP) - Y_V(IJMK)
                   Sz = Z_V(IJMKP) - Z_V(IJMK)
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dvdz_at_S =  (V_S(IJMKP,M) - V_S(IJMK,M)) * ONEoDZ_T_V(IJMK)
                   IF(NOC_WS) dvdz_at_S = dvdz_at_S - ((Vi-VW_s) * &
                                 ONEoDZ_T_V(IJMK)/DEL_H*(Sx*Nx+Sy*Ny))
                ELSE
                   dvdz_at_S =  ZERO
                ENDIF

                IF(V_NODE_AT_NB) THEN
                   CALL GET_DEL_H(IJK,'W_MOMENTUM',X_V(IJK),Y_V(IJK),&
                                  Z_V(IJK),Del_H,Nx,Ny,Nz)
                   SSY_CUT = - MU_S_CUT * (V_S(IJK,M) - VW_s) / DEL_H * &
                             (Nz*Ny) * Area_W_CUT(IJK)
                ELSE
                   SSY_CUT =  ZERO
                ENDIF

                SSY = AVG_Z_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),&
                              AVG_Y_H(MU_S(IJKT,M),MU_S(IJKNT,M),J),K)*&
                      dvdz_at_N*AXZ_W(IJK) - &
                      AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
                              AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)*&
                      dvdz_at_S*AXZ_W(IJMK) + SSY_CUT

! SSZ:
                CALL GET_DEL_H(IJK,'W_MOMENTUM',X_W(IJK),Y_W(IJK),&
                               Z_W(IJK),Del_H,Nx,Ny,Nz)
                SSZ = MU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
                      ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
                      MU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
                      OX(I)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM) &
                      - MU_S_CUT * (W_S(IJK,M) - WW_s) / DEL_H * &
                      (Nz**2) * Area_W_CUT(IJK)

              ENDIF  ! end if/else cut_cell
! ----------------------------------------------------------------<<<
            ENDIF   ! end if/else cartesian grid


! Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN
! modify Szz to include integral of (1/x)*(d/dz)(2*u/x)
               SSZ = SSZ + &
                  MU_S(IJKT,M)*(U_S(IJKP,M)+U_S(IMJKP,M))*OX(I)*AXY_W(IJK) - &
                  MU_S(IJK,M)*(U_S(IJK,M)+U_S(IMJK,M))*OX(I)*AXY_W(IJKM)

! (mu_s/x)* part of tau_xz/X
               IF (OX_E(IM) /= UNDEFINED) THEN
                  DUODZ = (U_S(IMJKP,M)-U_S(IMJK,M))*OX_E(IM)*ODZ_T(K)
               ELSE
                  DUODZ = ZERO
               ENDIF
               VXZ = AVG_Z(MU_S(IJK,M),MU_S(IJKT,M),K)*OX(I)*HALF*&
                  ((U_S(IJKP,M)-U_S(IJK,M))*OX_E(I)*ODZ_T(K)+DUODZ)
            ELSE
               VXZ = ZERO
            ENDIF

! Add the terms
              TAU_W_S(IJK,M) = SBV + SSX + SSY + SSZ + VXZ*VOL_W(IJK)
          ELSE
            TAU_W_S(IJK,M) = ZERO
          ENDIF   ! end if/else .not.sip_at_t .and. epsa>dil_ep_s

        ENDDO   ! end do ijk
      ENDDO

      call send_recv(tau_w_s,2)
      RETURN
      END SUBROUTINE CALC_TAU_W_S
