!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_s(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_s momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
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
      SUBROUTINE CG_SOURCE_W_S(A_M, B_M, M, IER) 
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
      INTEGER          I, J, K, IJK, IJKT, IMJK, IJMK, IJKM, IJKP, IMJKP 
      INTEGER          IJKE, IJKW, IJKTE, IJKTW, IM, IPJK 
! 
!                      Phase index 
      INTEGER          M, MM, L
      DOUBLE PRECISION   SUM_EPS_CP 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at top cell 
      DOUBLE PRECISION PgT 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPSA, EPStmp, epse, epsw, epsn, epss, &
                       epst, epsb, epsMix, epsMixT
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
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
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sdps, Sxzb, Vxza, Vxzb 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vmt, Vbf, Vcoa, Vcob, Vmttmp
!
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop
! 
!                      error message 
      CHARACTER*80     LINE
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: JM,IP,JP,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,KM,KP,IJMKP
      INTEGER :: IJKTN,IJKWT,IJKST
      DOUBLE PRECISION :: We,Ww,Wn,Ws,Wt,Wb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Uge, Ugw, Vgb, Vgt, Wge, Wgw, Wgn, Wgs, Wgt, Wgb,Ugc,Vgc,Vgn,Vgs
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


      IF(CG_SAFE_MODE(5)==1) RETURN


        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
          IF (MOMENTUM_Z_EQ(M)) THEN 

!$omp  parallel do &
!$omp& private(IJK, I, J, K, IJKT, EPSA, EPStmp, ISV, &
!$omp& PGT,SDP,SDPS,   ROPSA,V0,VMT,   DRO1,DRO2,DROA,VBF, &
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
                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN ! with ghd theory, m = mmax
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixT= ZERO   
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)  
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixT to be used for modelB
                    epsMixT = epsMixT + EP_S(IJKT,L)
                  ENDDO
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K) 
                ENDIF                 
                IF (IP_AT_T(IJK)) THEN 
!
!        do nothing
!     
                ELSEIF (SIP_AT_T(IJK)) THEN 
!
!        do nothing
!     

! dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN 
!
!        do nothing
!     
! Normal case
                ELSE 

                  BCV = BC_W_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)
                     CASE ('CG_NSW')
                        NOC_WS = .TRUE.
                        MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKT)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                        A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_W_CUT(IJK)/DELH_W(IJK) 
                     CASE ('CG_FSW')
                        NOC_WS = .FALSE.
                     CASE('CG_PSW')
                        IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_WS = .TRUE.
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKT)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_W_CUT(IJK)/DELH_W(IJK) 
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                           NOC_WS = .FALSE.
                        ELSE                              ! partial slip
                           NOC_WS = .FALSE.
                        ENDIF
                     CASE ('NONE')
                        NOC_WS = .FALSE.
                  END SELECT 

                  IF(NOC_WS) THEN

                     J = J_OF(IJK) 
                     K = K_OF(IJK)

                     IM = I - 1 
                     JM = J - 1 
                     KM = K - 1

                     IP = I + 1 
                     JP = J + 1 
                     KP = K + 1
    
                     IMJK = FUNIJK(IM,J,K)
                     IJMK = FUNIJK(I,JM,K)
                     IPJK = FUNIJK(IP,J,K)
                     IJPK = FUNIJK(I,JP,K)
                     IJKP = FUNIJK(I,J,KP)
                     IJKM = FUNIJK(I,J,KM)

                     We = Theta_We_bar(IJK)  * W_S(IJK,M)  + Theta_We(IJK)  * W_S(IPJK,M)
                     Ww = Theta_We_bar(IMJK) * W_S(IMJK,M) + Theta_We(IMJK) * W_S(IJK,M)

                     Wn = Theta_Wn_bar(IJK)  * W_S(IJK,M)  + Theta_Wn(IJK)  * W_S(IJPK,M)
                     Ws = Theta_Wn_bar(IJMK) * W_S(IJMK,M) + Theta_Wn(IJMK) * W_S(IJK,M)

                     Wt = Theta_Wt_bar(IJK)  * W_S(IJK,M)  + Theta_Wt(IJK)  * W_S(IJKP,M)
                     Wb = Theta_Wt_bar(IJKM) * W_S(IJKM,M) + Theta_Wt(IJKM) * W_S(IJK,M)
      
                     IPJK = IP_OF(IJK) 
                     IJPK = JP_OF(IJK) 
                     IJKE = EAST_OF(IJK) 

                     ijkt = top_of(ijk)

                     IF (WALL_AT(IJK)) THEN 
                        IJKC = IJKT
                     ELSE 
                        IJKC = IJK 
                     ENDIF 
                     IP = IP1(I) 
                     IJKN = NORTH_OF(IJK) 
                     IJKNE = EAST_OF(IJKN)

                     JM = JM1(J) 
                     IPJMK = IP_OF(IJMK) 
                     IJKS = SOUTH_OF(IJK) 
                     IJKSE = EAST_OF(IJKS) 

                     KP = KP1(K) 
                     IJKT = TOP_OF(IJK) 
                     IJKE = EAST_OF(IJK) 
                     IJKP = KP_OF(IJK) 
                     IJKTN = NORTH_OF(IJKT) 
                     IJKTE = EAST_OF(IJKT) 
                     IJKW = WEST_OF(IJK) 
                     IJKWT = TOP_OF(IJKW) 
                     IJKS = SOUTH_OF(IJK) 
                     IJKST = TOP_OF(IJKS) 

                     MU_S_E = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)

                     MU_S_W = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                              AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)

                     MU_S_N = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                              AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)

                     MU_S_S = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                              AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)

                     MU_S_T = MU_S(IJKT,M)
                     MU_S_B = MU_S(IJKC,M)

                     B_NOC =     MU_S_E * Ayz_W(IJK)  * We * NOC_W_E(IJK)  &
                             -   MU_S_W * Ayz_W(IMJK) * Ww * NOC_W_E(IMJK) &
                             +   MU_S_N * Axz_W(IJK)  * Wn * NOC_W_N(IJK)  &
                             -   MU_S_S * Axz_W(IJMK) * Ws * NOC_W_N(IJMK) &
                             +   MU_S_T * Axy_W(IJK)  * Wt * NOC_W_T(IJK)  *2.0d0&
                             -   MU_S_B * Axy_W(IJKM) * Wb * NOC_W_T(IJKM) *2.0D0

                     B_M(IJK,M) = B_M(IJK,M)   +  B_NOC
                  ENDIF

                  IF(CUT_W_TREATMENT_AT(IJK)) THEN
!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term  dWg/dt to virtual mass term			    
                     F_vir = ZERO
	             IF(Added_Mass.AND. M == M_AM ) THEN 

                        F_vir = ( (W_g(IJK) - W_gO(IJK)) )*ODT*VOL_W(IJK)

                        I = I_OF(IJK) 
                        J = J_OF(IJK) 
                        K = K_OF(IJK)
   
                        IM = I - 1 
                        JM = J - 1 
                        KM = K - 1

                        IP = I + 1 
                        JP = J + 1 
                        KP = K + 1

                        IMJK = FUNIJK(IM,J,K)
                        IJMK = FUNIJK(I,JM,K)
                        IPJK = FUNIJK(IP,J,K)
                        IJPK = FUNIJK(I,JP,K)
                        IJKP = FUNIJK(I,J,KP)
                        IJKM = FUNIJK(I,J,KM)

                        IMJKP = KP_OF(IMJK)
                        IJMKP = KP_OF(IJMK)

                        IJKE = EAST_OF(IJK) 
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)    


                        Wge = Theta_We_bar(IJK) * W_g(IJK) + Theta_We(IJK) * W_g(IPJK)
                        Wgw = Theta_We_bar(IMJK) * W_g(IMJK) + Theta_We(IMJK) * W_g(IJK)

                        Uge = Theta_W_te(IJK) * U_g(IJK) + Theta_W_be(IJK) * U_g(IJKP)
                        Ugw = Theta_W_te(IMJK) * U_g(IMJK) + Theta_W_be(IMJK) * U_g(IMJKP)

                        Ugc = (DELX_we(IJK) * Ugw + DELX_ww(IJK) * Uge) / (DELX_we(IJK) + DELX_ww(IJK))


                        Wgn = Theta_Wn_bar(IJK) * W_g(IJK) + Theta_Wn(IJK) * W_g(IJPK)
                        Wgs = Theta_Wn_bar(IJMK) * W_g(IJMK) + Theta_Wn(IJMK) * W_g(IJK)

                        Vgn =  Theta_W_tn(IJK)  * V_g(IJK)  + Theta_W_bn(IJK)  * V_g(IJKP)
                        Vgs =  Theta_W_tn(IJMK) * V_g(IJMK) + Theta_W_bn(IJMK) * V_g(IJMKP)

                        Vgc = (DELY_wn(IJK) * Vgs + DELY_ws(IJK) * Vgn) / (DELY_wn(IJK) + DELY_ws(IJK))


                        Wgt = Theta_Wt_bar(IJK)  * W_g(IJK)  + Theta_Wt(IJK)  * W_g(IJKP)
                        Wgb = Theta_Wt_bar(IMJK) * W_g(IMJK) + Theta_Wt(IMJK) * W_g(IMJKP)

!
! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass

                        F_vir = F_vir +  Ugc * (Wge*AYZ(IPJK) - Wgw*AYZ(IJK))    + &
                                         Vgc * (Wgn*AXZ(IJPK) - Wgs*(AXZ(IJK)))  + &
                                         W_g(IJK)*(Wgt*AYZ(IJKP) - Wgb*AYZ(IJK))

	         
                        ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M))/(VOL(IJK) + VOL(IJKT))

	                F_vir = F_vir * Cv * ROP_MA

                        B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit part of virtual mass force

                     ENDIF
!
!!! END VIRTUAL MASS SECTION

                  ENDIF


                ENDIF   ! end if sip or ip or dilute flow branch
            ENDDO 

          ENDIF 
        ENDIF ! for GHD Theory

      RETURN  
      END SUBROUTINE CG_SOURCE_W_S 


!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_s_BC(A_m, B_m, M, IER)                     C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_s momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
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
      SUBROUTINE CG_SOURCE_W_S_BC(A_M, B_M, M, IER) 
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

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: Del_H,Nx,Ny,Nz,Um,Vm,Wm,VdotN
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT

!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'

      IF(CG_SAFE_MODE(5)==1) RETURN

      DO IJK = ijkstart3, ijkend3

         BCV = BC_W_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT)  
   
            CASE ('CG_NSW')
    
               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 

                  B_M(IJK,M) = ZERO

               ENDIF 

            CASE ('CG_FSW')

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 

!                  B_M(IJK,M) = - W_s(W_MASTER_OF(IJK),M)  ! Velocity of master node

                  B_M(IJK,M) = ZERO 

                  IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN

                     IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                        A_M(IJK,E,M) = ONE 
                     ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                        A_M(IJK,W,M) = ONE 
                     ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                        A_M(IJK,N,M) = ONE 
                     ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                        A_M(IJK,S,M) = ONE 
                     ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                        A_M(IJK,T,M) = ONE 
                     ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                        A_M(IJK,B,M) = ONE 
                     ENDIF 

                  ENDIF

               ENDIF  

            CASE ('CG_PSW')

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 


                  IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = ZERO
                  ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO 

                     IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN   

                        IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 

                     ENDIF

                  ELSE                              ! partial slip

                  ENDIF

               ENDIF

            CASE ('CG_MI')

               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 

               B_M(IJK,M) = - BC_W_s(BC_W_ID(IJK),M)

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,E,M) = ZERO 
                  A_M(IJKB,W,M) = ZERO 
                  A_M(IJKB,N,M) = ZERO 
                  A_M(IJKB,S,M) = ZERO 
                  A_M(IJKB,T,M) = ZERO 
                  A_M(IJKB,B,M) = ZERO 
                  A_M(IJKB,0,M) = -ONE 
                  B_M(IJKB,M) = - BC_W_s(BC_W_ID(IJK),M)  

               ENDIF

            CASE ('CG_PO')

               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 

               ENDIF

         END SELECT 

         BCV = BC_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT)  

            CASE ('CG_MI')

               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO    
               A_M(IJK,0,M) = -ONE 

               B_M(IJK,M) = - BC_W_s(BC_ID(IJK),M)  


               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,E,M) = ZERO 
                  A_M(IJKB,W,M) = ZERO 
                  A_M(IJKB,N,M) = ZERO 
                  A_M(IJKB,S,M) = ZERO 
                  A_M(IJKB,T,M) = ZERO 
                  A_M(IJKB,B,M) = ZERO 
                  A_M(IJKB,0,M) = -ONE 
                  B_M(IJKB,M) = - BC_W_s(BC_ID(IJK),M)  

               ENDIF


            CASE ('CG_PO')

               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 

               ENDIF

         END SELECT 

      ENDDO

      RETURN

!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END SUBROUTINE CG_SOURCE_W_S_BC 
