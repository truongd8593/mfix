!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_V_s(A_m, B_m, M, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_s momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_S(A_M, B_M, M, IER) 
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
      USE vshear
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
      INTEGER          I, J, K, IJK,IMJK, IJMK, IJKM, IJKN
! 
!                      Phase index 
      INTEGER          M, MM, L
      DOUBLE PRECISION   SUM_EPS_CP 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at north cell 
      DOUBLE PRECISION PgN 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPSA, EPStmp, epse, epsw, epsn, epss, &
                       epst, epsb, epsMix, epsMixN
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
! 
!                      Average density difference 
      DOUBLE PRECISION dro1, dro2, droa 
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
      DOUBLE PRECISION V0, Vmt, Vbf, Vmttmp 
!
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop
!
! loezos
      DOUBLE PRECISION VSH_n,VSH_s,VSH_e,VSH_w,VSH_p,Source_conv
      DOUBLE PRECISION SRT
! loezos
 
!                      error message 
      CHARACTER*80     LINE(2) 
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER ::          IM,JM,IP,JP,KM,KP
      INTEGER ::          IPJK,IJPK,IJKP,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK,IJPKM
      INTEGER ::          IJKT,IJKTN,IJKB,IJKBN
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      DOUBLE PRECISION :: VW_s
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Vgn, Vgs, Uge, Ugw, Ugc,Vge, Vgw, Wgt, Wgb, Wgc,Vgt, Vgb
      DOUBLE PRECISION :: ep_star_avg,TH_avg,EPs_avg,F_2,H_W
!                      radial distribution function 
      DOUBLE PRECISION , EXTERNAL :: G_0CS
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

      IF(CG_SAFE_MODE(4)==1) RETURN


        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN

          IF (MOMENTUM_Y_EQ(M)) THEN 

!$omp  parallel do private( I, J, K, IJK, IJKN, ISV, Sdp, Sdps, V0, Vmt, &
!$omp&  PGN,DRO1,DRO2,DROA, Vbf, ROPSA, EPSA, EPStmp, VSH_n,VSH_s,VSH_e,&
!$omp&  VSH_w,VSH_p,Source_conv, SRT,SUM_EPS_CP,MM) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3
                I = I_OF(IJK) 
                J = J_OF(IJK) 
                K = K_OF(IJK)
                IMJK = IM_OF(IJK)
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK) 
                IJKN = NORTH_OF(IJK)
                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixN= ZERO  
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J) 
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixN to be used for modelB
                    epsMixN = epsMixN + EP_S(IJKN,L)
                  ENDDO                        
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J) 
                ENDIF 
                IF (IP_AT_N(IJK)) THEN 
!
!        do nothing
!     
                ELSEIF (SIP_AT_N(IJK)) THEN 
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


                  BCV = BC_V_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)
                     CASE ('CG_NSW')
                        NOC_VS = .TRUE.
                        VW_s = ZERO
                        MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                        A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_S_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
                     CASE ('CG_FSW')
                        NOC_VS = .FALSE.
                        VW_s = ZERO
                     CASE('CG_PSW')
                        IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
                           NOC_VS = .FALSE.
                           VW_s = BC_VW_S(BCV,M)
                           CALL CG_CALC_GRBDRY(IJK, 'V_MOMENTUM', M, BCV, F_2)
                           A_M(IJK,0,M) = A_M(IJK,0,M) - Area_V_CUT(IJK)*F_2
                           B_M(IJK,M) = B_M(IJK,M) - Area_V_CUT(IJK)*F_2*VW_s
                        ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_VS = .TRUE.
                           VW_s = BC_VW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                           A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_S_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * VW_s * Area_V_CUT(IJK)/DELH_V(IJK) 
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                           NOC_VS = .FALSE.
                           VW_s = ZERO
                        ELSE                              ! partial slip
                           NOC_VS = .FALSE.
                           VW_s = BC_VW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_V_CUT(IJK)*(BC_HW_S(BCV,M)) 
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * VW_s * Area_V_CUT(IJK)*(BC_HW_S(BCV,M))
                        ENDIF

                     CASE ('NONE')
                        NOC_VS = .FALSE.

                  END SELECT 

                  IF(NOC_VS) THEN


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

                     Vn = Theta_Vn_bar(IJK)  * V_S(IJK,M)  + Theta_Vn(IJK)  * V_S(IJPK,M)
                     Vs = Theta_Vn_bar(IJMK) * V_S(IJMK,M) + Theta_Vn(IJMK) * V_S(IJK,M)

                     Ve = Theta_Ve_bar(IJK)  * V_S(IJK,M)  + Theta_Ve(IJK)  * V_S(IPJK,M)
                     Vw = Theta_Ve_bar(IMJK) * V_S(IMJK,M) + Theta_Ve(IMJK) * V_S(IJK,M)

                     IPJK = IP_OF(IJK)
                     IMJK = IM_OF(IJK)
                     IJPK = JP_OF(IJK)
                     IJMK = JM_OF(IJK)
                     IJKP = KP_OF(IJK)
                     IJKM = KM_OF(IJK)
                     I = I_OF(IJK) 
                     J = J_OF(IJK) 
                     K = K_OF(IJK)     
                     IJKN = NORTH_OF(IJK) 
                     IF (WALL_AT(IJK)) THEN 
                        IJKC = IJKN 
                     ELSE 
                        IJKC = IJK 
                     ENDIF 
                     JP = JP1(J) 
                     IJKE = EAST_OF(IJK) 
                     IJKNE = EAST_OF(IJKN) 

                     IM = IM1(I) 
                     IJKW = WEST_OF(IJK) 
                     IJKWN = NORTH_OF(IJKW) 
                     IMJPK = JP_OF(IMJK) 

                     IJKT = TOP_OF(IJK) 
                     IJKTN = NORTH_OF(IJKT) 

                     IJKB = BOTTOM_OF(IJK) 
                     IJKBN = NORTH_OF(IJKB) 

                     MU_S_E = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)

                     MU_S_W = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                              AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)

                     MU_S_N = MU_S(IJKN,M)

                     MU_S_S = MU_S(IJKC,M)

                     B_NOC =   MU_S_N * Axz_V(IJK)  * (Vn-VW_s) * NOC_V_N(IJK)   *2.0d0&
                             - MU_S_S * Axz_V(IJMK) * (Vs-VW_s) * NOC_V_N(IJMK)  *2.0d0&
                             + MU_S_E * Ayz_V(IJK)  * (Ve-VW_s) * NOC_V_E(IJK)   &
                             - MU_S_W * Ayz_V(IMJK) * (Vw-VW_s) * NOC_V_E(IMJK)  

                     IF(DO_K) THEN

                        Vt = Theta_Vt_bar(IJK)  * V_S(IJK,M)  + Theta_Vt(IJK)  * V_S(IJKP,M)
                        Vb = Theta_Vt_bar(IJKM) * V_S(IJKM,M) + Theta_Vt(IJKM) * V_S(IJK,M)

                        MU_S_T = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                 AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)

                        MU_S_B = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                 AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)


                        B_NOC = B_NOC + MU_S_T * Axy_V(IJK)  * (Vt-VW_s) * NOC_V_T(IJK)   &
                                      - MU_S_B * Axy_V(IJKM) * (Vb-VW_s) * NOC_V_T(IJKM) 

                     ENDIF

                     B_M(IJK,M) = B_M(IJK,M) + B_NOC

                  ENDIF


                  IF(CUT_V_TREATMENT_AT(IJK)) THEN

!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term dVg/dt to virtual mass term		    
	             F_vir = ZERO
                     IF(Added_Mass.AND. M==M_AM ) THEN         
                        F_vir = ( (V_g(IJK) - V_gO(IJK)) )*ODT*VOL_V(IJK)

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

                        IMJPK = IM_OF(IJPK)

                        IJKN = NORTH_OF(IJK)  
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)     

                        Vge = Theta_Ve_bar(IJK)  * V_g(IJK)  + Theta_Ve(IJK)  * V_g(IPJK)
                        Vgw = Theta_Ve_bar(IMJK) * V_g(IMJK) + Theta_Ve(IMJK) * V_g(IJK)

                        Uge =  Theta_V_ne(IJK)  * U_g(IJK)  + Theta_V_se(IJK)  * U_g(IJPK)
                        Ugw =  Theta_V_ne(IMJK) * U_g(IMJK) + Theta_V_se(IMJK) * U_g(IMJPK)

                        Ugc = (DELX_ve(IJK) * Ugw + DELX_vw(IJK) * Uge) / (DELX_ve(IJK) + DELX_vw(IJK))

                        Vgn = Theta_Vn_bar(IJK)  * V_g(IJK)  + Theta_Vn(IJK)  * V_g(IJPK)
                        Vgs = Theta_Vn_bar(IJMK) * V_g(IJMK) + Theta_Vn(IJMK) * V_g(IJK)
	      
	                IF(DO_K) THEN

                           IJPKM = KM_OF(IJPK) 

                           Vgt = Theta_Vt_bar(IJK)  * V_g(IJK)  + Theta_Vt(IJK)  * V_g(IJKP)
                           Vgb = Theta_Vt_bar(IJKM) * V_g(IJKM) + Theta_Vt(IJKM) * V_g(IJK)

                           Wgt = Theta_V_nt(IJK)  * W_g(IJK)  + Theta_V_st(IJK)  * W_g(IJPK)
                           Wgb = Theta_V_nt(IJKM) * W_g(IJKM) + Theta_V_st(IJKM) * W_g(IJPKM)
                           Wgc = (DELZ_vt(IJK) * Wgb + DELZ_vb(IJK) * Wgt) / (DELZ_vt(IJK) + DELZ_vb(IJK))

                           F_vir = F_vir +  Wgc* (Vgt - Vgb)*AXY(IJK)

                        ENDIF
!
! adding convective terms (U dV/dx + V dV/dy) to virtual mass; W dV/dz added above.
                        F_vir = F_vir + V_g(IJK)*(Vgn - Vgs)*AXZ(IJK) + &
                        Ugc*(Vge - Vgw)*AYZ(IJK)
	    
                        ROP_MA =  (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + VOL(IJKN)*ROP_g(IJKN)*EP_s(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                        F_vir = F_vir * Cv * ROP_MA

                        B_M(IJK,M) = B_M(IJK,M) - F_vir ! adding explicit-part of virtual mass force.
                     ENDIF
!
!!! END VIRTUAL MASS SECTION
                  ENDIF
                  
                ENDIF   ! end if sip or ip or dilute flow branch
            ENDDO

    
          ENDIF  
        ENDIF ! for GHD Theory
   
      RETURN  
      END SUBROUTINE CG_SOURCE_V_S 

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_V_s_BC(A_m, B_m, M, IER)                     C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_s momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_S_BC(A_M, B_M, M, IER) 
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

      IF(CG_SAFE_MODE(4)==1) RETURN


      DO IJK = ijkstart3, ijkend3

         BCV = BC_V_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT) 
   
            CASE ('CG_NSW')
    
               IF(WALL_V_AT(IJK)) THEN

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

               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 

!                  B_M(IJK,M) = - V_s(V_MASTER_OF(IJK),M)    ! Velocity of master node

                  B_M(IJK,M) = ZERO 

                  IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                     IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                        A_M(IJK,E,M) = ONE 
                     ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                        A_M(IJK,W,M) = ONE 
                     ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                        A_M(IJK,N,M) = ONE 
                     ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                        A_M(IJK,S,M) = ONE 
                     ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                        A_M(IJK,T,M) = ONE 
                     ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                        A_M(IJK,B,M) = ONE 
                     ENDIF 

                  ENDIF

               ENDIF  

            CASE ('CG_PSW')

               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 


                  IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_VW_S(BCV,M)
                  ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW


                     B_M(IJK,M) = ZERO 

                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 

                     ENDIF

                  ELSE                              ! partial slip  (WARNING:currently same as FSW)

                     B_M(IJK,M) = ZERO 

                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 

                     ENDIF

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

               B_M(IJK,M) = - BC_V_s(BC_V_ID(IJK),M)  

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO 
                  A_M(IJKS,W,M) = ZERO 
                  A_M(IJKS,N,M) = ZERO 
                  A_M(IJKS,S,M) = ZERO 
                  A_M(IJKS,T,M) = ZERO 
                  A_M(IJKS,B,M) = ZERO 
                  A_M(IJKS,0,M) = -ONE 
                  B_M(IJKS,M) = - BC_V_s(BC_V_ID(IJK),M)  

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

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJK,S,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO

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

               B_M(IJK,M) = - BC_V_s(BC_ID(IJK),M)  

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO 
                  A_M(IJKS,W,M) = ZERO 
                  A_M(IJKS,N,M) = ZERO 
                  A_M(IJKS,S,M) = ZERO 
                  A_M(IJKS,T,M) = ZERO 
                  A_M(IJKS,B,M) = ZERO 
                  A_M(IJKS,0,M) = -ONE 
                  B_M(IJKS,M) = - BC_V_s(BC_ID(IJK),M)  

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

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJK,S,M) = ONE 
                  A_M(IJK,0,M) = -ONE 

               ENDIF

         END SELECT 
    


      ENDDO

      RETURN 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
 
      END SUBROUTINE CG_SOURCE_V_S_BC 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor

