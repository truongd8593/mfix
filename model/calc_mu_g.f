!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MU_g(IER)                                         C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: MFIX 2.0 mods (previous name CALC_MU_gt)                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE CALC_MU_G(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE indices
      USE constant
      USE compar    !//d
      USE sendrecv  !// 400
      USE funits !// AIKEPARDBG
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2./3. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                      IM, JM, KM
      INTEGER          IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK, IJMKP, &
                      IJMKM, IJPKM

!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION D_g(3,3)
!
!                      U_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION U_g_N
!
!                      U_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION U_g_S
!
!                      U_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION U_g_T
!
!                      U_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION U_g_B
!
!                      U_g at the center of the THETA cell-(i, j, k)
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION U_g_C
!
!                      V_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION V_g_E
!
!                      V_g at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION V_g_W
!
!                      V_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION V_g_T
!
!                      V_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION V_g_B
!
!                      W_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION W_g_E
!
!                      W_g at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION W_g_W
!
!                      W_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION W_g_N
!
!                      W_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION W_g_S
!
!                      W_g at the center of the THETA cell-(i, j, k).
!                      Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION W_g_C
!
!                      Second invariant of the deviator of D_g
      DOUBLE PRECISION I2_devD_g
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!
!
!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
!// 350 1112 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3 
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
!
!  Molecular viscosity
!
            IF (MU_G0 == UNDEFINED) MU_G(IJK) = 1.7D-4*(T_G(IJK)/273.0)**1.5*(&
               383./(T_G(IJK)+110.)) 
            MU_GT(IJK) = MU_G(IJK) 
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK) 
         ENDIF 
      END DO 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): after 1st DO ijk loop in CALC_MU_G ')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP

!!$omp parallel do &
!!$omp$ schedule(dynamic,chunk_size) &
!!$omp$ private(IJK, I,J,K,IM,JM,KM, &
!!$omp& IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJMK,IMJKP, &
!!$omp& IMJKM,IPJKM,IPJMK,IJMKP,IJMKM,IJPKM, &
!!$omp& U_G_N,U_G_S,U_G_T,U_G_B,V_G_E,V_G_W,V_G_T,V_G_B, &
!!$omp$ W_G_N,W_G_S,W_G_E,W_G_W,  U_G_C,W_G_C, D_G,I2_DEVD_G )


!//? 1112 need equivalents of ijkmin1, ijkmax1 based on each processors subdomain
! Following defns assumes DO_K = TRUE see set_max2.f for details
!      ijmax2_L = (iend3-istart3+1)*(jend3-jstart3+1)
!      ijkmin1_L = ijmax3 + 1
!      ijkmax1_L = ((iend3-istart3+1)*(jend3-jstart3+1)*(kend3-kstart3+1))- ijmax3

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): IJKMIN1= ',I4,' IJKMAX1= ',I4,' IJMAX2 = ',I4)") &
!           myPE,IJKMIN1,IJKMAX1,IJMAX2    !//AIKEPARDBG
!    write(*,"('(PE ',I2,'): IJMIN1_L= ',I4,' IJKMAX1_L= ',I4,'  IJMAX2_L= ',I4)") &
!           myPE, ((iend3-istart3+1)*(jend3-jstart3+1)+1),   &  
!	  (((iend3-istart3+1)*(jend3-jstart3+1)*(kend3-kstart3+1))- (iend3-istart3+1)*(jend3-jstart3+1)), &
!          ((iend3-istart3+1)*(jend3-jstart3+1)) 

!//? 1112 Make sure all the variables used in following DO loop are up to date
!//       at the boundaries as it has references to K-1. K-1 on PE 1 at the
!//       bottom face refers to the top face layer of PE 0!!!

!// 350 changed do loop limits ijkmin1, ijkmax1 ==> ijkstart3,ijkend3     
!      DO IJK = IJKMIN1, IJKMAX1 
      DO IJK = ijkstart3, ijkend3
!//? replace the crowd once ijmax_L, ijmin1_L, ijmax1_L are defined & set globally
!      DO IJK = ((iend3-istart3+1)*(jend3-jstart3+1)+1), &
!       ((iend3-istart3+1)*(jend3-jstart3+1)*(kend3-kstart3+1)) - (iend3-istart3+1)*(jend3-jstart3+1)


!//SP
         IF ( .NOT.WALL_AT(IJK) .AND. L_SCALE(IJK)/=ZERO) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
!// 360 1025 Check if current i,j,k resides on this PE	    
!            IF (.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
!/AIKERPARDBG
!            write(UNIT_LOG,"('(PE ',I2,'): from calc_mu_g ',&                    
!                & /,9X,'K = ',I6,'  KM = ',I6, &   
!                & /,9X,'J = ',I6,'  JM = ',I6,&    
!		& /,9X,'I = ',I6,'  IM = ',I6)") &  !//AIKEPARDBG
!                & myPE,K,KM,J,JM,I,IM !//AIKEPARDBG
	    
            IM = IM1(I) 
            JM = JM1(J) 
            KM = KM1(K) 

            IMJK = IM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK) 
            IMJPK = IM_OF(IJPK) 
            IMJMK = IM_OF(IJMK) 
            IMJKP = IM_OF(IJKP) 
            IMJKM = IM_OF(IJKM) 
            IPJKM = IP_OF(IJKM) 
            IPJMK = IP_OF(IJMK) 
            IJMKP = JM_OF(IJKP) 
            IJMKM = JM_OF(IJKM) 
            IJPKM = JP_OF(IJKM) 	    
!
!         Find fluid velocity values at faces of the cell
            U_G_N = AVG_Y(AVG_X_E(U_G(IMJK),U_G(IJK),I),AVG_X_E(U_G(IMJPK),U_G(&
               IJPK),I),J)                       !i, j+1/2, k 
            U_G_S = AVG_Y(AVG_X_E(U_G(IMJMK),U_G(IJMK),I),AVG_X_E(U_G(IMJK),U_G&
               (IJK),I),JM)                      !i, j-1/2, k 
            U_G_T = AVG_Z(AVG_X_E(U_G(IMJK),U_G(IJK),I),AVG_X_E(U_G(IMJKP),U_G(&
               IJKP),I),K)                       !i, j, k+1/2 
            U_G_B = AVG_Z(AVG_X_E(U_G(IMJKM),U_G(IJKM),I),AVG_X_E(U_G(IMJK),U_G&
               (IJK),I),KM)                      !i, j, k-1/2 
            V_G_E = AVG_X(AVG_Y_N(V_G(IJMK),V_G(IJK)),AVG_Y_N(V_G(IPJMK),V_G(&
               IPJK)),I)                         !i+1/2, j, k 
            V_G_W = AVG_X(AVG_Y_N(V_G(IMJMK),V_G(IMJK)),AVG_Y_N(V_G(IJMK),V_G(&
               IJK)),IM)                         !i-1/2, j, k 
            V_G_T = AVG_Z(AVG_Y_N(V_G(IJMK),V_G(IJK)),AVG_Y_N(V_G(IJMKP),V_G(&
               IJKP)),K)                         !i, j, k+1/2 
            V_G_B = AVG_Z(AVG_Y_N(V_G(IJMKM),V_G(IJKM)),AVG_Y_N(V_G(IJMK),V_G(&
               IJK)),KM)                         !i, j, k-1/2 
            W_G_N = AVG_Y(AVG_Z_T(W_G(IJKM),W_G(IJK)),AVG_Z_T(W_G(IJPKM),W_G(&
               IJPK)),J)                         !i, j+1/2, k 
            W_G_S = AVG_Y(AVG_Z_T(W_G(IJMKM),W_G(IJMK)),AVG_Z_T(W_G(IJKM),W_G(&
               IJK)),JM)                         !i, j-1/2, k 
            W_G_E = AVG_X(AVG_Z_T(W_G(IJKM),W_G(IJK)),AVG_Z_T(W_G(IPJKM),W_G(&
               IPJK)),I)                         !i+1/2, j, k 
            W_G_W = AVG_X(AVG_Z_T(W_G(IMJKM),W_G(IMJK)),AVG_Z_T(W_G(IJKM),W_G(&
               IJK)),IM)                         !i-1/2, j, k 
!/AIKERPARDBG
!            write(UNIT_LOG,"('U_G_N= ',E12.4 &                    
!                & ,2X,'U_G_S= ',E12.4,2X,'U_G_T= ',E12.4,2X,'U_G_B= ',E12.4 &   
!!                & ,/,'V_G_E= ',E12.4,2X,'V_G_W= ',E12.4,2X,'V_G_T= ',E12.4 &    
!!		& ,2X,'V_G_B= ',E12.4
!                & ,/,'W_G_N= ',E12.4,2X,'W_G_S= ',E12.4,2X,'W_G_E= ',E12.4 &    
!		& ,2X,'W_G_W= ',E12.4)") &    !//AIKEPARDBG
!                  U_G_N,U_G_S,U_G_T,U_G_B, &
!!		  V_G_E,V_G_W,V_G_T,V_G_B, &
!		  W_G_N,W_G_S,W_G_E,W_G_W     !//AIKEPARDBG

!
            IF (CYLINDRICAL) THEN 
!                                                !i, j, k
               U_G_C = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
!                                                !i, j, k
               W_G_C = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            ELSE 
               U_G_C = ZERO 
               W_G_C = ZERO 
            ENDIF 

!
!         Find components of fluid phase strain rate
!         tensor, D_g, at center of the cell - (i, j, k)
            D_G(1,1) = (U_G(IJK)-U_G(IMJK))*ODX(I) 
            D_G(1,2) = HALF*((U_G_N - U_G_S)*ODY(J)+(V_G_E-V_G_W)*ODX(I)) 
            D_G(1,3) = HALF*((W_G_E - W_G_W)*ODX(I)+(U_G_T-U_G_B)*(OX(I)*ODZ(K)&
               )-W_G_C*OX(I)) 
            D_G(2,1) = D_G(1,2) 
            D_G(2,2) = (V_G(IJK)-V_G(IJMK))*ODY(J) 
            D_G(2,3)=HALF*((V_G_T-V_G_B)*(OX(I)*ODZ(K))+(W_G_N-W_G_S)*ODY(J)) 
            D_G(3,1) = D_G(1,3) 
            D_G(3,2) = D_G(2,3) 
            D_G(3,3) = (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) + U_G_C*OX(I) 
!
!         Calculate the second invariant of the deviator of D_g
!
            I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+(D_G(3,3&
               )-D_G(1,1))**2)/6. + D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2 
!
            MU_GT(IJK) = MIN(MU_GMAX,MU_GT(IJK)+2.0*L_SCALE(IJK)*L_SCALE(IJK)*&
               RO_G(IJK)*SQRT(I2_DEVD_G)) 
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK) 
	    
         ENDIF 
      END DO 

!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): before COMM in CALC_MU_G ')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP

!//S 1113 try to move this COMM to the end of transport_prop to do all COMMs
!//       at certain locations, provided that no data dependency in between.

!//? 1112 need to update the boundaries for variables : MU_GT,LAMBDA_GT
!// 400 1112 update the boundaries for recently calculated field vars
!!!      call send_recv(mu_gt,2)
!!!      call send_recv(lambda_gt,2) 
           
!//AIKEPARDBG
!      write(*,"('(PE ',I2,'): eof CALC_MU_G ')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP
      
      RETURN  
      END SUBROUTINE CALC_MU_G 
