!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_g(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_G(A_M, B_M, IER) 
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
      USE compar    
      USE sendrecv  

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
      INTEGER          I, IJK, IJKE, IPJK, IJKM, IPJKM 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
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
!                      Average viscosity 
      DOUBLE PRECISION EPMU_gte, EPMU_gbe, EPMUGA 
! 
!                      Average W_g 
      DOUBLE PRECISION Wge 
! 
!                      Average dW/Xdz 
      DOUBLE PRECISION dWoXdz 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcf, Vtza 
! 
!                      error message 
      CHARACTER*80     LINE 
!
!     FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_G_temp(DIMENSION_3)

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER :: J,K,IM,JM,IP,JP,KM,KP
      INTEGER :: IMJK,IJMK,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKP
      INTEGER :: IJKT,IJKTE,IJKB,IJKBE
      DOUBLE PRECISION :: Ue,Uw,Un,Us,Ut,Ub
      DOUBLE PRECISION :: P_CUT,z_plane,DH,Nx,Ny,Nz,B_NOC
      DOUBLE PRECISION :: x_circle,y_circle, Angle
      DOUBLE PRECISION :: MU_GT_E,MU_GT_W,MU_GT_N,MU_GT_S,MU_GT_T,MU_GT_B,MU_GT_CUT
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT

      INTEGER   IJKDB1,IJKDB2 

      LOGICAL :: PRINT_FLAG
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
!
      M = 0 
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN  
!
!$omp    parallel do private(I, IJK, IJKE, IJKM, IPJK, IPJKM,     &
!$omp&                  ISV, Sdp, V0, Vpm, Vmt, Vbf,              &
!$omp&                  Vcf, EPMUGA, VTZA, WGE, PGE, ROGA,        &
!$omp&                  MUGA, ROPGA, EPGA )
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
         IJKE = EAST_OF(IJK) 
         IJKM = KM_OF(IJK) 
         IPJK = IP_OF(IJK) 
         IPJKM = IP_OF(IJKM) 
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
         IF (IP_AT_E(IJK)) THEN 
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
            BCV = BC_U_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE(BCV)
            ELSE
               BCT = 'NONE'
            ENDIF

            SELECT CASE (BCT)
               CASE ('CG_NSW')
                  NOC_UG = .TRUE.
                  MU_GT_CUT =  (VOL(IJK)*MU_GT(IJK) + VOL(IPJK)*MU_GT(IJKE))/(VOL(IJK) + VOL(IPJK))
                  A_M(IJK,0,M) = A_M(IJK,0,M) - MU_GT_CUT * Area_U_CUT(IJK)/DELH_U(IJK) 
               CASE ('CG_FSW')
                  NOC_UG = .FALSE.
               CASE('CG_PSW')
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_UG = .TRUE.
                     MU_GT_CUT =  (VOL(IJK)*MU_GT(IJK) + VOL(IPJK)*MU_GT(IJKE))/(VOL(IJK) + VOL(IPJK))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_GT_CUT * Area_U_CUT(IJK)/DELH_U(IJK) 
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_UG = .FALSE.
                  ELSE                              ! partial slip
                     NOC_UG = .FALSE.
                  ENDIF
               CASE ('NONE')
                  NOC_UG = .FALSE.
            END SELECT 

            IF(NOC_UG) THEN

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


               Ue = Theta_Ue_bar(IJK)  * U_g(IJK)  + Theta_Ue(IJK)  * U_g(IPJK)
               Uw = Theta_Ue_bar(IMJK) * U_g(IMJK) + Theta_Ue(IMJK) * U_g(IJK)

               Un = Theta_Un_bar(IJK)  * U_g(IJK)  + Theta_Un(IJK)  * U_g(IJPK)
               Us = Theta_Un_bar(IJMK) * U_g(IJMK) + Theta_Un(IJMK) * U_g(IJK)

     
               IPJK = IP_OF(IJK) 
               IJPK = JP_OF(IJK) 
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

               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 

               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB)

               MU_GT_E = MU_GT(IJKE)
               MU_GT_W = MU_GT(IJKC)

               MU_GT_N = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),&
                                 AVG_Y_H(MU_GT(IJKE),MU_GT(IJKNE),J),I)

               MU_GT_S = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),&
                                 AVG_Y_H(MU_GT(IJKSE),MU_GT(IJKE),JM),I)

               B_NOC =     MU_GT_E * Ayz_U(IJK)  * Ue * NOC_U_E(IJK)  *2.0d0&
                       -   MU_GT_W * Ayz_U(IMJK) * Uw * NOC_U_E(IMJK) *2.0d0&
                       +   MU_GT_N * Axz_U(IJK)  * Un * NOC_U_N(IJK)  &
                       -   MU_GT_S * Axz_U(IJMK) * Us * NOC_U_N(IJMK) 

               IF(DO_K) THEN

                  Ut = Theta_Ut_bar(IJK)  * U_g(IJK)  + Theta_Ut(IJK)  * U_g(IJKP)
                  Ub = Theta_Ut_bar(IJKM) * U_g(IJKM) + Theta_Ut(IJKM) * U_g(IJK)

                  MU_GT_T = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),&
                            AVG_Z_H(MU_GT(IJKE),MU_GT(IJKTE),K),I)

                  MU_GT_B = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),&
                            AVG_Z_H(MU_GT(IJKBE),MU_GT(IJKE),KM),I)

                  B_NOC = B_NOC  +   MU_GT_T * Axy_U(IJK)  * Ut * NOC_U_T(IJK)  &
                                 -   MU_GT_B * Axy_U(IJKM) * Ub * NOC_U_T(IJKM) 
               ENDIF

               B_M(IJK,M) = B_M(IJK,M)   +  B_NOC


            ENDIF

         ENDIF 

      END DO 


      RETURN  
      END SUBROUTINE CG_SOURCE_U_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_U_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for U_g momentum eq.                                                C
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
      SUBROUTINE CG_SOURCE_U_G_BC(A_M, B_M, IER) 
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
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IP, IPJK 
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
!
      M = 0 

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
 
!                  B_M(IJK,M) = - U_g(U_MASTER_OF(IJK))  ! Velocity of master node

                  B_M(IJK,M) = ZERO 

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


            CASE ('CG_PSW')

               IF(WALL_U_AT(IJK)) THEN

                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = ZERO
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO 
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

               B_M(IJK,M) = - BC_U_g(BC_U_ID(IJK))


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO 
                  A_M(IJKW,W,M) = ZERO 
                  A_M(IJKW,N,M) = ZERO 
                  A_M(IJKW,S,M) = ZERO 
                  A_M(IJKW,T,M) = ZERO 
                  A_M(IJKW,B,M) = ZERO 
                  A_M(IJKW,0,M) = -ONE 
                  B_M(IJKW,M) = - BC_U_g(BC_U_ID(IJK))  

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

               B_M(IJK,M) = - BC_U_g(BC_ID(IJK))  


               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN

                  A_M(IJKW,E,M) = ZERO 
                  A_M(IJKW,W,M) = ZERO 
                  A_M(IJKW,N,M) = ZERO 
                  A_M(IJKW,S,M) = ZERO 
                  A_M(IJKW,T,M) = ZERO 
                  A_M(IJKW,B,M) = ZERO 
                  A_M(IJKW,0,M) = -ONE 
                  B_M(IJKW,M) = - BC_U_g(BC_ID(IJK))  

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

      END SUBROUTINE CG_SOURCE_U_G_BC 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
