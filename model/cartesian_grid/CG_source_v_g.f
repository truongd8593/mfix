!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_g(A_m, B_m, IER)                              C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_G(A_M, B_M, IER) 
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
      INTEGER          M, L
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
      INTEGER ::          IM,JM,IP,JP,KM,KP
      INTEGER ::          IMJK,IPJK,IJMK,IJPK,IJKP,IJKM,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK,IJPKM
      INTEGER ::          IJKT,IJKTN,IJKB,IJKBN
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_GT_E,MU_GT_W,MU_GT_N,MU_GT_S,MU_GT_T,MU_GT_B,MU_GT_CUT
      DOUBLE PRECISION :: VW_g
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!			virtual (added) mass
      DOUBLE PRECISION F_vir, ROP_MA, Vsn, Vss, U_se, Usc, Usw, Vse, Vsw, Wst, Wsb, Wsc, Vst, Vsb
! Wall function
      DOUBLE PRECISION :: W_F_Slip

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

      IF(CG_SAFE_MODE(4)==1) RETURN

!
      M = 0 
      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN  
!
!$omp  parallel do private( I, J, K, IJK, IJKN, ISV, Sdp, V0, Vpm, Vmt, Vbf, &
!$omp&  PGN, ROGA, MUGA, ROPGA, EPGA,VSH_n,VSH_s,VSH_e,VSH_w,&
!$omp&  VSH_p,Source_conv, SRT ) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
         EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
         IF (IP_AT_N(IJK)) THEN 
!
!        do nothing
!     
!       dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN 
!
!        do nothing
!     
         ELSE 
!
            BCV = BC_V_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE(BCV)
            ELSE
               BCT = 'NONE'
            ENDIF

            SELECT CASE (BCT)
               CASE ('CG_NSW')
                  NOC_VG = .TRUE.
                  VW_g = ZERO
                  MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJKN)*MU_GT(IJKN))/(VOL(IJK) + VOL(IJKN))  

                  IF(.NOT.K_EPSILON) THEN
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_GT_CUT * Area_V_CUT(IJK)/DELH_V(IJK)         
                  ELSE
                     CALL Wall_Function(IJK,IJK,ONE/DELH_V(IJK),W_F_Slip)
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_GT_CUT * Area_V_CUT(IJK)*(ONE-W_F_Slip)/DELH_V(IJK)         
                  ENDIF
               CASE ('CG_FSW')
                  NOC_VG = .FALSE.
                  VW_g = ZERO
               CASE('CG_PSW')
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_VG = .TRUE.
                     VW_g = BC_VW_G(BCV)
                     MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJKN)*MU_GT(IJKN))/(VOL(IJK) + VOL(IJKN))  
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_GT_CUT * Area_V_CUT(IJK)/DELH_V(IJK)  
                     B_M(IJK,M) = B_M(IJK,M) - MU_GT_CUT * VW_g * Area_V_CUT(IJK)/DELH_V(IJK)        
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_VG = .FALSE.
                     VW_g = ZERO
                  ELSE                              ! partial slip
                     NOC_VG = .FALSE.
                     VW_g = BC_VW_G(BCV)
                     MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJKN)*MU_GT(IJKN))/(VOL(IJK) + VOL(IJKN))  
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_GT_CUT * Area_V_CUT(IJK)*(BC_HW_G(BCV))
                     B_M(IJK,M) = B_M(IJK,M) - MU_GT_CUT * VW_g * Area_V_CUT(IJK)*(BC_HW_G(BCV))
                  ENDIF
               CASE ('NONE')
                  NOC_VG = .FALSE.
            END SELECT 


            IF(NOC_VG) THEN

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


               Vn = Theta_Vn_bar(IJK)  * V_g(IJK)  + Theta_Vn(IJK)  * V_g(IJPK)
               Vs = Theta_Vn_bar(IJMK) * V_g(IJMK) + Theta_Vn(IJMK) * V_g(IJK)

               Ve = Theta_Ve_bar(IJK)  * V_g(IJK)  + Theta_Ve(IJK)  * V_g(IPJK)
               Vw = Theta_Ve_bar(IMJK) * V_g(IMJK) + Theta_Ve(IMJK) * V_g(IJK)


               IPJK = IP_OF(IJK)
               IMJK = IM_OF(IJK)
               IJPK = JP_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJk)     
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


               MU_GT_E = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),&
                         AVG_X_H(MU_GT(IJKN),MU_GT(IJKNE),I),J)

               MU_GT_W = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),&
                         AVG_X_H(MU_GT(IJKWN),MU_GT(IJKN),IM),J)

               MU_GT_N = MU_GT(IJKN)

               MU_GT_S = MU_GT(IJKC)


               B_NOC =     MU_GT_N * Axz_V(IJK)  * (Vn-VW_g) * NOC_V_N(IJK)   *2.0d0&
                       -   MU_GT_S * Axz_V(IJMK) * (Vs-VW_g) * NOC_V_N(IJMK)  *2.0d0&
                       +   MU_GT_E * Ayz_V(IJK)  * (Ve-VW_g) * NOC_V_E(IJK)   &
                       -   MU_GT_W * Ayz_V(IMJK) * (Vw-VW_g) * NOC_V_E(IMJK)  

               IF(DO_K) THEN
  
                  Vt = Theta_Vt_bar(IJK)  * V_g(IJK)  + Theta_Vt(IJK)  * V_g(IJKP)
                  Vb = Theta_Vt_bar(IJKM) * V_g(IJKM) + Theta_Vt(IJKM) * V_g(IJK)

                  MU_GT_T = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),&
                            AVG_Z_H(MU_GT(IJKN),MU_GT(IJKTN),K),J)

                  MU_GT_B = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),&
                            AVG_Z_H(MU_GT(IJKBN),MU_GT(IJKN),KM),J)

                  B_NOC = B_NOC + MU_GT_T * Axy_V(IJK)  * (Vt-VW_g) * NOC_V_T(IJK)   &
                                - MU_GT_B * Axy_V(IJKM) * (Vb-VW_g) * NOC_V_T(IJKM) 

               ENDIF

                  B_M(IJK,M) = B_M(IJK,M) + B_NOC

            ENDIF


            IF(CUT_V_TREATMENT_AT(IJK)) THEN

!
!!! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term dVs/dt to virtual mass term		    
	       F_vir = ZERO
	       IF(Added_Mass) THEN        
	          F_vir = ( (V_s(IJK,M_AM) - V_sO(IJK,M_AM)) )*ODT*VOL_V(IJK)

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


                  Vse = Theta_Ve_bar(IJK)  * V_s(IJK,M_AM)  + Theta_Ve(IJK)  * V_s(IPJK,M_AM)
                  Vsw = Theta_Ve_bar(IMJK) * V_s(IMJK,M_AM) + Theta_Ve(IMJK) * V_s(IJK,M_AM)

                  U_se =  Theta_V_ne(IJK)  * U_s(IJK,M_AM)  + Theta_V_se(IJK)  * U_s(IJPK,M_AM)
                  Usw =  Theta_V_ne(IMJK) * U_s(IMJK,M_AM) + Theta_V_se(IMJK) * U_s(IMJPK,M_AM)

                  Usc = (DELX_ve(IJK) * Usw + DELX_vw(IJK) * U_se) / (DELX_ve(IJK) + DELX_vw(IJK))


                  Vsn = Theta_Vn_bar(IJK)  * V_s(IJK,M_AM)  + Theta_Vn(IJK)  * V_s(IJPK,M_AM)
                  Vss = Theta_Vn_bar(IJMK) * V_s(IJMK,M_AM) + Theta_Vn(IJMK) * V_s(IJK,M_AM)


	      
	          IF(DO_K) THEN

                     IJPKM = KM_OF(IJPK) 

                     Vst = Theta_Vt_bar(IJK)  * V_s(IJK,M_AM)  + Theta_Vt(IJK)  * V_s(IJKP,M_AM)
                     Vsb = Theta_Vt_bar(IJKM) * V_s(IJKM,M_AM) + Theta_Vt(IJKM) * V_s(IJK,M_AM)

                     Wst = Theta_V_nt(IJK)  * W_s(IJK,M_AM)  + Theta_V_st(IJK)  * W_s(IJPK,M_AM)
                     Wsb = Theta_V_nt(IJKM) * W_s(IJKM,M_AM) + Theta_V_st(IJKM) * W_s(IJPKM,M_AM)
		     Wsc = (DELZ_vt(IJK) * Wsb + DELZ_vb(IJK) * Wst) / (DELZ_vt(IJK) + DELZ_vb(IJK))

	             F_vir = F_vir +  Wsc* (Vst - Vsb)*AXY(IJK)

	          ENDIF
!
! adding convective terms (U dV/dx + V dV/dy) to virtual mass; W dV/dz added above.
	          F_vir = F_vir + V_s(IJK,M_AM)*(Vsn - Vss)*AXZ(IJK) + &
	              Usc*(Vse - Vsw)*AYZ(IJK)
	    

                  ROP_MA  = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM)  + VOL(IJKN)*ROP_g(IJKN)*EP_s(IJKN,M_AM) )/(VOL(IJK) + VOL(IJKN))
	          F_vir = F_vir * Cv * ROP_MA

                  B_M(IJK,M) = B_M(IJK,M) - F_vir ! adding explicit-part of virtual mass force.
	       ENDIF
!
!!! END VIRTUAL MASS SECTION
            ENDIF

         ENDIF 
      END DO 
!
      RETURN  
      END SUBROUTINE CG_SOURCE_V_G 


!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_V_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for V_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_V_G_BC(A_M, B_M, IER) 
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
                       IM, KM, IJKS, IJMK, IJPK ,IMJK
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

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: Del_H,Nx,Ny,Nz,Um,Vm,Wm,VdotN
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT

      INTEGER   IJKDB1,IJKDB2 

      LOGICAL :: PRINT_FLAG
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'

      IF(CG_SAFE_MODE(4)==1) RETURN

!
      M = 0 

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

                  IF(K_EPSILON) THEN

                     I = I_OF(IJK) 
                     J = J_OF(IJK) 
                     K = K_OF(IJK)

                     IM = I - 1 
                     JM = J - 1 
                     KM = K - 1

                  IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                     IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
                        CALL Wall_Function(IJK,EAST_OF(IJK),ODX_E(I),W_F_Slip)
                        A_M(IJK,E,M) = W_F_Slip
                        A_M(IJK,0,M) = -ONE 
                     ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
                        CALL Wall_Function(IJK,WEST_OF(IJK),ODX_E(IM),W_F_Slip)
                        A_M(IJK,W,M) = W_F_Slip
                        A_M(IJK,0,M) = -ONE 
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

            CASE ('CG_FSW')

               IF(WALL_V_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 

!                  B_M(IJK,M) = - V_g(V_MASTER_OF(IJK))    ! Velocity of master node

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


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_VW_G(BCV)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
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

                  ELSE                              ! partial slip

                     B_M(IJK,M) = ZERO 

                     I = I_OF(IJK) 
                     J = J_OF(IJK) 
                     K = K_OF(IJK)

                     IM = I - 1 
                     JM = J - 1 
                     KM = K - 1

                     IMJK = FUNIJK(IM,J,K)

                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN

                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN 
!                          A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(I)) 
!                          A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ODX_E(I)) 
!                          B_M(IJK,M) = -BC_HW_G(BCV)*BC_VW_G(BCV) 

                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_V(IJK)) 
                           A_M(IJK,E,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_V(IJK)) 
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_VW_G(BCV) 

                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN 
!                          A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ODX_E(IM)) 
!                          A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODX_E(IM)) 
!                          B_M(IJK,M) = -BC_HW_G(BCV)*BC_VW_G(BCV) 

                           A_M(IJK,W,M) = -(HALF*BC_HW_G(BCV)-ONEoDX_E_V(IMJK)) 
                           A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ONEoDX_E_V(IMJK)) 
                           B_M(IJK,M) = -BC_HW_G(BCV)*BC_VW_G(BCV) 


                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN 

                              A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(BCV)-ODY_N(J)) 
                              B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV) 

!                           print*,'vg master at north'
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN 

                              A_M(IJK,S,M) = -(HALF*BC_HW_G(BCV)-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(BCV)+ODY_N(JM)) 
                              B_M(IJK,M) = -BC_HW_G(BCV)*BC_UW_G(BCV) 

!                           print*,'vg master at south'
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

               B_M(IJK,M) = - BC_V_g(BC_V_ID(IJK))  

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO 
                  A_M(IJKS,W,M) = ZERO 
                  A_M(IJKS,N,M) = ZERO 
                  A_M(IJKS,S,M) = ZERO 
                  A_M(IJKS,T,M) = ZERO 
                  A_M(IJKS,B,M) = ZERO 
                  A_M(IJKS,0,M) = -ONE 
                  B_M(IJKS,M) = - BC_V_g(BC_V_ID(IJK))  

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

               B_M(IJK,M) = - BC_V_g(BC_ID(IJK))  


               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJKS,E,M) = ZERO 
                  A_M(IJKS,W,M) = ZERO 
                  A_M(IJKS,N,M) = ZERO 
                  A_M(IJKS,S,M) = ZERO 
                  A_M(IJKS,T,M) = ZERO 
                  A_M(IJKS,B,M) = ZERO 
                  A_M(IJKS,0,M) = -ONE    
                  B_M(IJKS,M) = - BC_V_g(BC_ID(IJK))  

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
 
      END SUBROUTINE CG_SOURCE_V_G_BC 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
