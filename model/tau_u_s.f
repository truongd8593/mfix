!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_Tau_U_s(TAU_U_s, IER)                             C
!  Purpose: Cross terms in the gradient of stress in U_s momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
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
      SUBROUTINE CALC_TAU_U_S(TAU_U_S, IER) 
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
      USE vshear
      USE sendrecv  
      USE compar   
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
! 
!                      Error index 
      INTEGER          IER 
! 
!                      TAU_U_s 
      DOUBLE PRECISION TAU_U_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Indices 
      INTEGER          I, J, JM, K, KM, IJK, IJKE, IPJK, IP, IMJK, IJKN,& 
                       IJKNE, IJKS, IJKSE, IPJMK, IJMK, IJKT, IJKTE,& 
                       IJKB, IJKBE, IJKM, IPJKM 
! 
!                      Phase index 
      INTEGER          M, L
!
! loezos
       INTEGER I1,J1    
! loezos
 
!                      Average volume fraction 
      DOUBLE PRECISION EPSA, EPStmp
! 
!                      Average EP_s*viscosity 
      DOUBLE PRECISION EPMU_ste, EPMU_sbe, EPMUSA 
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
!                      error message 
      CHARACTER*80     LINE 
! 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: IM,JP,KP
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: V_NODE_AT_NE,V_NODE_AT_NW,V_NODE_AT_SE,V_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TE,W_NODE_AT_TW,W_NODE_AT_BE,W_NODE_AT_BW
      DOUBLE PRECISION :: V_SUM,W_SUM,X_SUM,Y_SUM,Z_SUM,Vc,Wc
      DOUBLE PRECISION :: Xvc,Yvc,Zvc,Xwc,Ywc,Zwc,Nxv,Nyv,Nzv,Nxw,Nyw,Nzw
      DOUBLE PRECISION :: dvdx_at_N,dvdx_at_S
      DOUBLE PRECISION :: dwdx_at_T,dwdx_at_B
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: MU_S_CUT,SSY_CUT,SSZ_CUT
      INTEGER :: N_SUM
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT  
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!


      DO M = 1, MMAX 
        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
! loezos
! update to true velocity
        IF (SHEAR) THEN        
!$omp  parallel do private(IJK)
          DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN 
              V_S(IJK,m)=V_S(IJK,m)+VSH(IJK)
            ENDIF
          ENDDO     
        ENDIF
! loezos


!!$omp  parallel do private( IJK, I, IJKE, EPSA, EPStmp, IP, J, JM, K, KM,  &
!!$omp&  IPJK,IMJK,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJMK,IJKT,IJKTE,  &
!!$omp&  IJKB,IJKBE,IJKM,IPJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ, EPMU_STE,EPMU_SBE, &
!!$omp&  EPMUSA,DWOXDZ,VTZB ) &
!!$omp&  schedule(static)

        DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
              EPStmp = ZERO                
              DO L = 1, SMAX
                EPStmp = EPStmp + AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I) 
              ENDDO
              EPSA = EPStmp
            ELSE                  
              EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
            ENDIF 

            IF ( .NOT.SIP_AT_E(IJK) .AND. EPSA>DIL_EP_S) THEN 
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
! Surface forces

! bulk viscosity term
                 SBV = (LAMBDA_S(IJKE,M)*TRD_S(IJKE,M)-LAMBDA_S(IJK,M)*TRD_S(IJK,&
                     M))*AYZ(IJK) 

! shear stress terms
                 SSX = MU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*ODX(IP)*AYZ_U(IJK)&
                      - MU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*ODX(I)*AYZ_U(IMJK) 
                 SSY = AVG_X_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
                     IJKE,M),MU_S(IJKNE,M),J),I)*(V_S(IPJK,M)-V_S(IJK,M))*ODX_E(I)&
                     *AXZ_U(IJK) - AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
                     AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*(V_S(IPJMK,M)-V_S(&
                     IJMK,M))*ODX_E(I)*AXZ_U(IJMK) 
                 EPMU_STE = AVG_X_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(&
                     MU_S(IJKE,M),MU_S(IJKTE,M),K),I) 
                 EPMU_SBE = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),AVG_Z_H(&
                     MU_S(IJKBE,M),MU_S(IJKE,M),KM),I) 
                 SSZ = EPMU_STE*(W_S(IPJK,M)-W_S(IJK,M))*ODX_E(I)*AXY_U(IJK) - &
                     EPMU_SBE*(W_S(IPJKM,M)-W_S(IJKM,M))*ODX_E(I)*AXY_U(IJKM) 

              ELSE ! CARTESIAN GRID CASE

! Surface forces

! bulk viscosity term
                 SBV =  (LAMBDA_S(IJKE,M)*TRD_S(IJKE,M)) * AYZ_U(IJK) &
                       -(LAMBDA_S(IJK,M) *TRD_S(IJK,M) ) * AYZ_U(IMJK) 

!         shear stress terms

                 IF(.NOT.CUT_U_CELL_AT(IJK))   THEN            

                    SSX = MU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*ONEoDX_E_U(IJK)*AYZ_U(IJK) - MU_S(&
                       IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*ONEoDX_E_U(IMJK)*AYZ_U(IMJK) 

                    SSY = AVG_X_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(IJKE,M)&
                       ,MU_S(IJKNE,M),J),I)*(V_S(IPJK,M)-V_S(IJK,M))*ONEoDX_E_V(IJK)*AXZ_U(IJK) - &
                       AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),AVG_Y_H(MU_S(IJKSE,M),&
                       MU_S(IJKE,M),JM),I)*(V_S(IPJMK,M)-V_S(IJMK,M))*ONEoDX_E_V(IJMK)*AXZ_U(IJMK) 


                    IF(DO_K) THEN
                       EPMU_STE = AVG_X_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                          IJKE,M),MU_S(IJKTE,M),K),I) 
                       EPMU_SBE = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),AVG_Z_H(MU_S&
                          (IJKBE,M),MU_S(IJKE,M),KM),I) 

                       SSZ = EPMU_STE*(W_S(IPJK,M)-W_S(IJK,M))*ONEoDX_E_W(IJK)*AXY_U(IJK) - EPMU_SBE*&
                          (W_S(IPJKM,M)-W_S(IJKM,M))*ONEoDX_E_W(IJKM)*AXY_U(IJKM) 

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
                          CUT_TAU_US = .TRUE.
                          NOC_US     = .TRUE.
                       CASE ('CG_FSW')
                          CUT_TAU_US = .FALSE.
                          NOC_US     = .FALSE.
                       CASE('CG_PSW')
                          IF(BC_HW_S(BC_U_ID(IJK),M)==UNDEFINED) THEN   ! same as NSW
                             CUT_TAU_US = .TRUE.
                             NOC_US     = .TRUE.
                          ELSEIF(BC_HW_S(BC_U_ID(IJK),M)==ZERO) THEN   ! same as FSW
                             CUT_TAU_US = .FALSE.
                             NOC_US     = .FALSE.
                          ELSE                              ! partial slip
                             CUT_TAU_US = .FALSE.
                             NOC_US     = .FALSE.
                          ENDIF

                       CASE ('NONE')
                          TAU_U_S(IJK,M) = ZERO 
                          RETURN  
                    END SELECT 

                    IF(CUT_TAU_US) THEN
                       MU_S_CUT =  (VOL(IJK)*MU_S(IJK,M) + VOL(IPJK)*MU_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                    ELSE
                       MU_S_CUT = ZERO
                    ENDIF

!           SSX:

                    CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U(IJK),Y_U(IJK),Z_U(IJK),Del_H,Nx,Ny,Nz)

                    SSX = MU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*ONEoDX_E_U(IJK)*AYZ_U(IJK) - MU_S(&
                          IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*ONEoDX_E_U(IMJK)*AYZ_U(IMJK) &
                        - MU_S_CUT * (U_S(IJK,M) - ZERO) / DEL_H * (Nx**2) * Area_U_CUT(IJK)         

!           SSY:

                    V_NODE_AT_NE = ((.NOT.BLOCKED_V_CELL_AT(IPJK)).AND.(.NOT.WALL_V_AT(IPJK)))
                    V_NODE_AT_NW = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.(.NOT.WALL_V_AT(IJK)))
                    V_NODE_AT_SE = ((.NOT.BLOCKED_V_CELL_AT(IPJMK)).AND.(.NOT.WALL_V_AT(IPJMK)))
                    V_NODE_AT_SW = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

                    IF(V_NODE_AT_NE.AND.V_NODE_AT_NW) THEN

                       Vi = HALF * (V_S(IPJK,M) + V_S(IJK,M))
                       Xi = HALF * (X_V(IPJK) + X_V(IJK))
                       Yi = HALF * (Y_V(IPJK) + Y_V(IJK))
                       Zi = HALF * (Z_V(IPJK) + Z_V(IJK))
                       Sx = X_V(IPJK) - X_V(IJK)
                       Sy = Y_V(IPJK) - Y_V(IJK)
                       Sz = Z_V(IPJK) - Z_V(IJK)

                       CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                       dvdx_at_N =  (V_S(IPJK,M) - V_S(IJK,M)) * ONEoDX_E_V(IJK)

                       IF(NOC_US) dvdx_at_N = dvdx_at_N - (Vi * ONEoDX_E_V(IJK) /DEL_H*(Sy*Ny+Sz*Nz))            

                    ELSE
                       dvdx_at_N =  ZERO
                    ENDIF


                    IF(V_NODE_AT_SE.AND.V_NODE_AT_SW) THEN

                       Vi = HALF * (V_S(IPJMK,M) + V_S(IJMK,M))
                       Xi = HALF * (X_V(IPJMK) + X_V(IJMK))
                       Yi = HALF * (Y_V(IPJMK) + Y_V(IJMK))
                       Zi = HALF * (Z_V(IPJMK) + Z_V(IJMK))
                       Sx = X_V(IPJMK) - X_V(IJMK)
                       Sy = Y_V(IPJMK) - Y_V(IJMK)
                       Sz = Z_V(IPJMK) - Z_V(IJMK)

                       CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                       dvdx_at_S =  (V_S(IPJMK,M) - V_S(IJMK,M)) * ONEoDX_E_V(IJMK)

                       IF(NOC_US) dvdx_at_S = dvdx_at_S - (Vi * ONEoDX_E_V(IJMK)/DEL_H*(Sy*Ny+Sz*Nz))        

                    ELSE
                       dvdx_at_S =  ZERO
                    ENDIF

                    IF(V_NODE_AT_NW) THEN
                       CALL GET_DEL_H(IJK,'U_MOMENTUM',X_V(IJK),Y_V(IJK),Z_V(IJK),Del_H,Nx,Ny,Nz)
                       SSY_CUT = - MU_S_CUT * (V_S(IJK,M) - ZERO) / DEL_H * (Nx*Ny) * Area_U_CUT(IJK)        
                    ELSE
                       SSY_CUT =  ZERO     
                    ENDIF

                    SSY = AVG_X_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(IJKE,M)&
                       ,MU_S(IJKNE,M),J),I)*dvdx_at_N*AXZ_U(IJK) - &
                       AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),AVG_Y_H(MU_S(IJKSE,M),&
                       MU_S(IJKE,M),JM),I)*dvdx_at_S*AXZ_U(IJMK) &
                      + SSY_CUT

!           SSZ:

                    IF(DO_K) THEN  

                       W_NODE_AT_TE = ((.NOT.BLOCKED_W_CELL_AT(IPJK)).AND.(.NOT.WALL_W_AT(IPJK)))
                       W_NODE_AT_TW = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))
                       W_NODE_AT_BE = ((.NOT.BLOCKED_W_CELL_AT(IPJKM)).AND.(.NOT.WALL_W_AT(IPJKM)))
                       W_NODE_AT_BW = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                       IF(W_NODE_AT_TE.AND.W_NODE_AT_TW) THEN

                          Wi = HALF * (W_S(IPJK,M) + W_S(IJK,M))
                          Xi = HALF * (X_W(IPJK) + X_W(IJK))
                          Yi = HALF * (Y_W(IPJK) + Y_W(IJK))
                          Zi = HALF * (Z_W(IPJK) + Z_W(IJK))
                          Sx = X_W(IPJK) - X_W(IJK)
                          Sy = Y_W(IPJK) - Y_W(IJK)
                          Sz = Z_W(IPJK) - Z_W(IJK)


                          CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                          dwdx_at_T =  (W_S(IPJK,M) - W_S(IJK,M)) * ONEoDX_E_W(IJK)  
      
                          IF(NOC_US) dwdx_at_T = dwdx_at_T - (Wi * ONEoDX_E_W(IJK)/DEL_H*(Sy*Ny+Sz*Nz))    

                       ELSE
                          dwdx_at_T =  ZERO
                       ENDIF


                       IF(W_NODE_AT_BE.AND.W_NODE_AT_BW) THEN

                          Wi = HALF * (W_S(IPJKM,M) + W_S(IJKM,M))
                          Xi = HALF * (X_W(IPJKM) + X_W(IJKM))
                          Yi = HALF * (Y_W(IPJKM) + Y_W(IJKM))
                          Zi = HALF * (Z_W(IPJKM) + Z_W(IJKM))
                          Sx = X_W(IPJKM) - X_W(IJKM)
                          Sy = Y_W(IPJKM) - Y_W(IJKM)
                          Sz = Z_W(IPJKM) - Z_W(IJKM)

                          CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                          dwdx_at_B =  (W_S(IPJKM,M) - W_S(IJKM,M)) * ONEoDX_E_W(IJKM) 

                          IF(NOC_US) dwdx_at_B = dwdx_at_B  - (Wi * ONEoDX_E_W(IJKM)/DEL_H*(Sy*Ny+Sz*Nz))      

                       ELSE
                          dwdx_at_B =  ZERO
                       ENDIF

                       IF(W_NODE_AT_TW) THEN
                          CALL GET_DEL_H(IJK,'U_MOMENTUM',X_W(IJK),Y_W(IJK),Z_W(IJK),Del_H,Nx,Ny,Nz)
                          SSZ_CUT = - MU_S_CUT * (W_S(IJK,M) - ZERO) / DEL_H * (Nx*Nz) * Area_U_CUT(IJK) 
                       ELSE
                          SSZ_CUT =  ZERO
                       ENDIF

                       EPMU_STE = AVG_X_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                          IJKE,M),MU_S(IJKTE,M),K),I) 
                       EPMU_SBE = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),AVG_Z_H(MU_S&
                          (IJKBE,M),MU_S(IJKE,M),KM),I) 
                       SSZ =   EPMU_STE*dwdx_at_T*AXY_U(IJK)  &
                             - EPMU_SBE*dwdx_at_B*AXY_U(IJKM) &
                             + SSZ_CUT
                    ELSE

                      SSZ = ZERO

                    ENDIF  ! DO_K

                 END IF  ! CUT CELL

!  Original terms

!               SSX = MU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*ODX(IP)*AYZ_U(IJK)&
!                   - MU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*ODX(I)*AYZ_U(IMJK) 
!               SSY = AVG_X_H(AVG_Y_H(MU_S(IJK,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
!                  IJKE,M),MU_S(IJKNE,M),J),I)*(V_S(IPJK,M)-V_S(IJK,M))*ODX_E(I)&
!                  *AXZ_U(IJK) - AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJK,M),JM),&
!                  AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*(V_S(IPJMK,M)-V_S(&
!                  IJMK,M))*ODX_E(I)*AXZ_U(IJMK) 
!               EPMU_STE = AVG_X_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),AVG_Z_H(&
!                  MU_S(IJKE,M),MU_S(IJKTE,M),K),I) 
!               EPMU_SBE = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),AVG_Z_H(&
!                  MU_S(IJKBE,M),MU_S(IJKE,M),KM),I) 
!               SSZ = EPMU_STE*(W_S(IPJK,M)-W_S(IJK,M))*ODX_E(I)*AXY_U(IJK) - &
!                  EPMU_SBE*(W_S(IPJKM,M)-W_S(IJKM,M))*ODX_E(I)*AXY_U(IJKM) 

              ENDIF  ! CARTESIAN GRID
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! Special terms for cylindrical coordinates
              IF (CYLINDRICAL) THEN 
!   modify Ssz: integral of (1/x) (d/dz) (mu*(-w/x))
                SSZ = SSZ - (EPMU_STE*(HALF*(W_S(IPJK,M)+W_S(IJK,M))*OX_E(I))&
                     *AXY_U(IJK)-EPMU_SBE*(HALF*(W_S(IPJKM,M)+W_S(IJKM,M))*OX_E&
                     (I))*AXY_U(IJKM)) 

!  -(2mu/x)*(1/x)*dw/dz part of Tau_zz/X
                EPMUSA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                DWOXDZ = HALF*((W_S(IJK,M)-W_S(IJKM,M))*OX(I)*ODZ(K)+(W_S(&
                     IPJK,M)-W_S(IPJKM,M))*OX(IP)*ODZ(K)) 
                VTZB = -2.d0*EPMUSA*OX_E(I)*DWOXDZ 
              ELSE 
                VTZB = ZERO 
              ENDIF 

! Add the terms
              TAU_U_S(IJK,M) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK) 
            ELSE 
              TAU_U_S(IJK,M) = ZERO 
            ENDIF 
        ENDDO 

! loezos 
        IF (SHEAR) THEN
!$omp  parallel do private(IJK) 
          DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN   
              V_S(IJK,m)=V_S(IJK,m)-VSH(IJK)
            ENDIF
          ENDDO
        ENDIF
! loezos
       
        ENDIF !  end if GHD theory
      ENDDO 
      call send_recv(tau_u_s,2)
      RETURN  
      END SUBROUTINE CALC_TAU_U_S 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
