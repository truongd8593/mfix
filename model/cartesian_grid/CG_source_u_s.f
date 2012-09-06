!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_s(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_s momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_S(A_M, B_M, M, IER) 
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
      DOUBLE PRECISION :: UW_s
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Uge, Ugw, Vgw, Vge, Vgn, Vgs,Vgc, Ugn, Ugs, Wgb, Wgt, Wge, Wgc, Ugb, Ugt

      DOUBLE PRECISION :: ep_star_avg,TH_avg,EPs_avg,F_2
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

      IF(CG_SAFE_MODE(3)==1) RETURN


        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
          IF (MOMENTUM_X_EQ(M)) THEN 

!!!$omp  parallel do private( IJK, IJKE, ISV, Sdp, Sdps, V0, Vmt, Vbf, &
!!!$omp&  I,PGE,DRO1,DRO2,DROA, IJKM,IPJK,IPJKM,  WSE,VCF,EPMUSA,VTZA, &
!!!$omp&  EPSA, EPStmp, ROPSA, LINE,SUM_EPS_CP,MM) &
!!!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3 

! Wall or impermeable internal surface
                I = I_OF(IJK) 
                IJKE = EAST_OF(IJK)
                IMJK = IM_OF(IJK) 
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK) 
                IPJK = IP_OF(IJK) 
                IPJKM = IP_OF(IJKM)
                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixE= ZERO          
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I) 
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixE to be used for modelB
                    epsMixE = epsMixE + EP_S(IJKE,L)
                  ENDDO                        
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
                ENDIF                  
                IF (IP_AT_E(IJK)) THEN 
!
!                  do nothing
!
                ELSEIF (SIP_AT_E(IJK)) THEN 
!
!                  do nothing
!
! Dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN 
!
!                  do nothing
!
                ELSE 

                  BCV = BC_U_ID(IJK)

                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE(BCV)
                  ELSE
                     BCT = 'NONE'
                  ENDIF

                  SELECT CASE (BCT)
                     CASE ('CG_NSW')
                        NOC_US = .TRUE.
                        UW_s = ZERO
                        MU_S_CUT =  (VOL(IJK)*MU_S(IJK,M) + VOL(IPJK)*MU_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                        A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_U_CUT(IJK)/DELH_U(IJK)  
                     CASE ('CG_FSW')
                        NOC_US = .FALSE.
                        UW_s = ZERO
                     CASE('CG_PSW')
                        IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
                           NOC_US = .FALSE.
                           UW_s = BC_UW_S(BCV,M)
                           CALL CG_CALC_GRBDRY(IJK, 'U_MOMENTUM', M, BCV, F_2)
                           A_M(IJK,0,M) = A_M(IJK,0,M) - Area_U_CUT(IJK)*F_2
                           B_M(IJK,M) = B_M(IJK,M) - Area_U_CUT(IJK)*F_2*UW_s
                        ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_US = .TRUE.
                           UW_s = BC_UW_S(BCV,M)
                           MU_S_CUT =  (VOL(IJK)*MU_S(IJK,M) + VOL(IPJK)*MU_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_U_CUT(IJK)/DELH_U(IJK)  
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * UW_s * Area_U_CUT(IJK)/DELH_U(IJK) 
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                           NOC_US = .FALSE.
                           UW_s = ZERO
                        ELSE                              ! partial slip
                           NOC_US = .FALSE.
                           UW_s = BC_UW_S(BCV,M)
                           MU_S_CUT =  (VOL(IJK)*MU_S(IJK,M) + VOL(IPJK)*MU_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_U_CUT(IJK)*(BC_HW_S(BCV,M)) 
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * UW_s * Area_U_CUT(IJK)*(BC_HW_S(BCV,M))
                        ENDIF



                     CASE ('NONE', 'CG_MI')
                        NOC_US = .FALSE.


                  END SELECT 

                  IF(NOC_US)THEN

                     I = I_OF(IJK) 
                     J = J_OF(IJK) 
                     K = K_OF(IJK)

                     IMJK = IM_OF(IJK)
                     IJMK = JM_OF(IJK)
                     IJKM = KM_OF(IJK)
                     IPJK = IP_OF(IJK)
                     IJPK = JP_OF(IJK)
                     IJKP = KP_OF(IJK)

                     Ue = Theta_Ue_bar(IJK)  * U_S(IJK,M)  + Theta_Ue(IJK)  * U_S(IPJK,M)
                     Uw = Theta_Ue_bar(IMJK) * U_S(IMJK,M) + Theta_Ue(IMJK) * U_S(IJK,M)
 
                     Un = Theta_Un_bar(IJK)  * U_S(IJK,M)  + Theta_Un(IJK)  * U_S(IJPK,M)
                     Us = Theta_Un_bar(IJMK) * U_S(IJMK,M) + Theta_Un(IJMK) * U_S(IJK,M)
      
                     IJKE = EAST_OF(IJK) 

                     IF (WALL_AT(IJK)) THEN 
                        IJKC = IJKE 
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

                     KM = KM1(K) 
                     IPJMK = IP_OF(IJMK) 
                     IJKS = SOUTH_OF(IJK) 
                     IJKSE = EAST_OF(IJKS) 

                     IJKT = TOP_OF(IJK) 
                     IJKTE = EAST_OF(IJKT) 

                     IJKB = BOTTOM_OF(IJK) 
                     IJKBE = EAST_OF(IJKB) 


                     MU_S_E = MU_S(IJKE,M)
                     MU_S_W = MU_S(IJKC,M)


                     MU_S_N = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                              AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)

                     MU_S_S = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                               AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)


                     B_NOC =   MU_S_E * Ayz_U(IJK)  * (Ue-UW_s) * NOC_U_E(IJK)  *2.0d0&
                           -   MU_S_W * Ayz_U(IMJK) * (Uw-UW_s) * NOC_U_E(IMJK) *2.0d0&
                           +   MU_S_N * Axz_U(IJK)  * (Un-UW_s) * NOC_U_N(IJK)  &
                           -   MU_S_S * Axz_U(IJMK) * (Us-UW_s) * NOC_U_N(IJMK) 



                     IF(DO_K) THEN

                        Ut = Theta_Ut_bar(IJK)  * U_S(IJK,M)  + Theta_Ut(IJK)  * U_S(IJKP,M)
                        Ub = Theta_Ut_bar(IJKM) * U_S(IJKM,M) + Theta_Ut(IJKM) * U_S(IJK,M)

                        MU_S_T = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                 AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)

                        MU_S_B = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                 AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)

                        B_NOC = B_NOC  +   MU_S_T * Axy_U(IJK)  * (Ut-UW_s) * NOC_U_T(IJK)  &
                                       -   MU_S_B * Axy_U(IJKM) * (Ub-UW_s) * NOC_U_T(IJKM) 
                     ENDIF
 
                     B_M(IJK,M) = B_M(IJK,M)   +  B_NOC

                  ENDIF


                  IF(CUT_U_TREATMENT_AT(IJK)) THEN
!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term  dUg/dt to virtual mass term			    
                     F_vir = ZERO
                    
                     IF(Added_Mass.AND. M==M_AM ) THEN        

                        F_vir = ( U_g(IJK) - U_gO(IJK) )*ODT*VOL_U(IJK)

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

                        IPJMK = IP_OF(IJMK)

                        IJKE = EAST_OF(IJK) 
!
! defining gas-particles velocity at momentum cell faces (or scalar cell center)    
                        Uge = Theta_Ue_bar(IJK)  * U_g(IJK)  + Theta_Ue(IJK)  * U_g(IPJK)
                        Ugw = Theta_Ue_bar(IMJK) * U_g(IMJK) + Theta_Ue(IMJK) * U_g(IJK)

                        Ugn = Theta_Un_bar(IJK)  * U_g(IJK)  + Theta_Un(IJK)  * U_g(IJPK)
                        Ugs = Theta_Un_bar(IJMK) * U_g(IJMK) + Theta_Un(IJMK) * U_g(IJK)

                        Vgn =  Theta_U_ne(IJK)  * V_g(IJK)  + Theta_U_nw(IJK)  * V_g(IPJK)
                        Vgs =  Theta_U_ne(IJMK) * V_g(IJMK) + Theta_U_nw(IJMK) * V_g(IPJMK)

                        Vgc = (DELY_un(IJK) * Vgs + DELY_us(IJK) * Vgn) / (DELY_un(IJK) + DELY_us(IJK))

                        IF(DO_K) THEN

                           IPJKM = IP_OF(IJKM) 

                           Ugt = Theta_Ut_bar(IJK)  * U_g(IJK)  + Theta_Ut(IJK)  * U_g(IJKP)
                           Ugb = Theta_Ut_bar(IJKM) * U_g(IJKM) + Theta_Ut(IJKM) * U_g(IJK)

                           Wgt = Theta_U_te(IJK)  * W_g(IJK)  + Theta_U_tw(IJK)  * W_g(IPJK)
                           Wgb = Theta_U_te(IJKM) * W_g(IJKM) + Theta_U_tw(IJKM) * W_g(IPJKM)
	                   Wgc = (DELZ_ut(IJK) * Wgb + DELZ_ub(IJK) * Wgt) / (DELZ_ut(IJK) + DELZ_ub(IJK))

                           F_vir = F_vir +  Wgc* (Ugt - Ugb)*AXY(IJK)

                        ENDIF
!
! adding convective terms (U dU/dx + V dU/dy + W dU/dz) to virtual mass
                        F_vir = F_vir + U_g(IJK)*(Uge - Ugw)*AYZ(IJK) + &
	                                Vgc * (Ugn - Ugs)*AXZ(IJK)

                        ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + VOL(IPJK)*ROP_g(IJKE)*EP_s(IJKE,M))/(VOL(IJK) + VOL(IPJK))

	                F_vir = F_vir * Cv * ROP_MA

                        B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit part of virtual mass force

                     ENDIF
!
!!! END VIRTUAL MASS SECTION

                   ENDIF

                ENDIF ! end if sip or ip or dilute flow branch
            ENDDO 
          ENDIF  
        ENDIF ! for GHD Theory

      
      RETURN  
      END SUBROUTINE CG_SOURCE_U_S
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_s_BC(A_m, B_m, M, IER)                     C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_s momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_S_BC(A_M, B_M, M, IER) 
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


      IF(CG_SAFE_MODE(3)==1) RETURN


      DO IJK = ijkstart3, ijkend3

         BCV = BC_U_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE(BCV)
         ELSE
            BCT = 'NONE'
         ENDIF

         SELECT CASE (BCT)  

            CASE ('CG_NSW')

               IF(WALL_U_AT(IJK)) THEN

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

               IF(WALL_U_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 

!                  B_M(IJK,M) = - U_s(U_MASTER_OF(IJK),M)  ! Velocity of master node

                  B_M(IJK,M) = ZERO 

                  IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN

                     IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                        A_M(IJK,E,M) = ONE 
                     ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                        A_M(IJK,W,M) = ONE 
                     ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                        A_M(IJK,N,M) = ONE 
                     ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                        A_M(IJK,S,M) = ONE 
                     ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                        A_M(IJK,T,M) = ONE 
                     ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                        A_M(IJK,B,M) = ONE 
                     ENDIF 

                  ENDIF

               ENDIF

            CASE ('CG_PSW')

               IF(WALL_U_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 


                  IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_UW_S(BCV,M)
                  ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO 

                     IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN
 
                        IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 

                     ENDIF

                  ELSE                              ! partial slip  (WARNING:currently same as FSW)

                     B_M(IJK,M) = ZERO 

                     IF(DABS(NORMAL_U(IJK,1))/=ONE) THEN
 
                        IF (U_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                           A_M(IJK,E,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                           A_M(IJK,W,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == TOP_OF(IJK)) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSEIF (U_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN 
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

               IF(BC_U_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_U_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_U(IJK,1)  
               ENDIF

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO 
                  A_M(IJKW,W,M) = ZERO 
                  A_M(IJKW,N,M) = ZERO 
                  A_M(IJKW,S,M) = ZERO 
                  A_M(IJKW,T,M) = ZERO 
                  A_M(IJKW,B,M) = ZERO 
                  A_M(IJKW,0,M) = -ONE 

                  IF(BC_U_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKW,M) = - BC_U_s(BCV,M)
                  ELSE
                     B_M(IJKW,M) = - BC_VELMAG_s(BCV,M)*NORMAL_U(IJK,1)  
                  ENDIF

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

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJK,W,M) = ONE
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

               IF(BC_U_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_U_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,1)  
               ENDIF 


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO 
                  A_M(IJKW,W,M) = ZERO 
                  A_M(IJKW,N,M) = ZERO 
                  A_M(IJKW,S,M) = ZERO 
                  A_M(IJKW,T,M) = ZERO 
                  A_M(IJKW,B,M) = ZERO 
                  A_M(IJKW,0,M) = -ONE 

                  IF(BC_U_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKW,M) = - BC_U_s(BCV,M)
                  ELSE
                     B_M(IJKW,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,1)
                  ENDIF

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


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJK,W,M) = ONE
                  A_M(IJK,0,M) = -ONE 

               ENDIF

         END SELECT 


      ENDDO

      RETURN  
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      END SUBROUTINE CG_SOURCE_U_S_BC 
