!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_e(A_m, VxF_gs, d_e, IER)                        C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- East                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
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
      SUBROUTINE CALC_D_E(A_M, VXF_GS, D_E, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales 
      USE compar     !//d
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
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA, EPSA
!
!                      F/(a0 + F), F/(a1 + F)
      DOUBLE PRECISION FoA0pF,     FoA1pF
!
!                      Indices
      INTEGER          I, IJK, IJKE
!
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (MOMENTUM_X_EQ(0) .AND. MOMENTUM_X_EQ(1)) THEN 
!
!//I? check any data dependency for I direction decomposition in following loop
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp  parallel do private( I,IJK, IJKE, EPGA, EPSA, FoA0pF, FoA1pF), &
!$omp&  schedule(static)

         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN 
               D_E(IJK,0) = ZERO 
               D_E(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
               EPSA = AVG_X(EP_S(IJK,1),EP_S(IJKE,1),I) 
               IF (MODEL_B) THEN 
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_E(IJK,0) = P_SCALE*AYZ(IJK)/((-A_M(IJK,0,0))-A_M(IJK,0,1&
                        )*FOA1PF) 
                     D_E(IJK,1) = FOA1PF*D_E(IJK,0) 
                  ELSE 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_E(IJK,0) = P_SCALE*AYZ(IJK)/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_E(IJK,0) = ZERO 
                     ENDIF 
                     D_E(IJK,1) = ZERO 
                  ENDIF 
!
               ELSE                              !Model A 
!                 MFIX convention: center coeff is negative
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
                     FOA0PF = VXF_GS(IJK,1)/((-A_M(IJK,0,0))+VXF_GS(IJK,1)) 
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_E(IJK,0) = P_SCALE*AYZ(IJK)*(EPGA + EPSA*FOA1PF)/((-A_M(&
                        IJK,0,0))-A_M(IJK,0,1)*FOA1PF) 
                     D_E(IJK,1) = P_SCALE*AYZ(IJK)*(EPSA + EPGA*FOA0PF)/((-A_M(&
                        IJK,0,1))-A_M(IJK,0,0)*FOA0PF) 
                  ELSE 
                     IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                        D_E(IJK,0) = P_SCALE*AYZ(IJK)*EPGA/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_E(IJK,0) = ZERO 
                     ENDIF 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_E(IJK,1) = P_SCALE*AYZ(IJK)*EPSA/(-A_M(IJK,0,1)) 
                     ELSE 
                        D_E(IJK,1) = ZERO 
                     ENDIF 
!
                  ENDIF 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_X_EQ(0)) THEN 

!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp    parallel do &
!$omp&   private( IJK, I, IJKE, EPGA )

         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN 
               D_E(IJK,0) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_E(IJK,0) = P_SCALE*AYZ(IJK)/((-A_M(IJK,0,0))+VXF_GS(IJK,&
                        1)) 
                  ELSE 
                     D_E(IJK,0) = P_SCALE*AYZ(IJK)*EPGA/((-A_M(IJK,0,0))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_E(IJK,0) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_X_EQ(1)) THEN 

!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
!$omp    parallel do &
!$omp&   private( IJK, I, IJKE, EPSA )

         DO IJK = ijkstart3, ijkend3 
            IF (IP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN 
               D_E(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               EPSA = AVG_X(EP_S(IJK,1),EP_S(IJKE,1),I) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_E(IJK,1) = ZERO 
                  ELSE 
                     D_E(IJK,1) = P_SCALE*AYZ(IJK)*EPSA/((-A_M(IJK,0,1))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_E(IJK,1) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE CALC_D_E 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_n(A_m, VxF_gs, d_n, IER)                        C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- North                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
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
      SUBROUTINE CALC_D_N(A_M, VXF_GS, D_N, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales 
      USE compar   !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA, EPSA
!
!                      F/(a0 + F), F/(a1 + F)
      DOUBLE PRECISION FoA0pF,     FoA1pF
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKN
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (MOMENTUM_Y_EQ(0) .AND. MOMENTUM_Y_EQ(1)) THEN 
!
!//I? check any data dependency for I direction decomposition in following loop
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp  parallel do private( J, K, IJK, IJKN, EPGA, EPSA, FoA0pF, FoA1pF), &
!$omp&  schedule(static)
         DO IJK = ijkstart3, ijkend3 
            IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_N(IJK,0) = ZERO 
               D_N(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
               EPSA = AVG_Y(EP_S(IJK,1),EP_S(IJKN,1),J) 
               IF (MODEL_B) THEN 
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
!
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_N(IJK,0) = P_SCALE*AXZ(IJK)/((-A_M(IJK,0,0))-A_M(IJK,0,1&
                        )*FOA1PF) 
                     D_N(IJK,1) = FOA1PF*D_N(IJK,0) 
!
                  ELSE 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_N(IJK,0) = P_SCALE*AXZ(IJK)/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_N(IJK,0) = ZERO 
                     ENDIF 
                     D_N(IJK,1) = ZERO 
                  ENDIF 
!
               ELSE                              !Model A 
!                 MFIX convention: center coeff is negative
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
                     FOA0PF = VXF_GS(IJK,1)/((-A_M(IJK,0,0))+VXF_GS(IJK,1)) 
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_N(IJK,0) = P_SCALE*AXZ(IJK)*(EPGA + EPSA*FOA1PF)/((-A_M(&
                        IJK,0,0))-A_M(IJK,0,1)*FOA1PF) 
                     D_N(IJK,1) = P_SCALE*AXZ(IJK)*(EPSA + EPGA*FOA0PF)/((-A_M(&
                        IJK,0,1))-A_M(IJK,0,0)*FOA0PF) 
                  ELSE 
                     IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                        D_N(IJK,0) = P_SCALE*AXZ(IJK)*EPGA/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_N(IJK,0) = ZERO 
                     ENDIF 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_N(IJK,1) = P_SCALE*AXZ(IJK)*EPSA/(-A_M(IJK,0,1)) 
                     ELSE 
                        D_N(IJK,1) = ZERO 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_Y_EQ(0)) THEN 
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
      
!$omp    parallel do &
!$omp&   private( IJK, I,J,K,  IJKN, EPGA )
         DO IJK = ijkstart3, ijkend3 
            IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_N(IJK,0) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_N(IJK,0) = P_SCALE*AXZ(IJK)/((-A_M(IJK,0,0))+VXF_GS(IJK,&
                        1)) 
                  ELSE 
                     D_N(IJK,0) = P_SCALE*AXZ(IJK)*EPGA/((-A_M(IJK,0,0))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_N(IJK,0) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_Y_EQ(1)) THEN 

!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp    parallel do &
!$omp&   private( IJK, I,J,K, IJKN, EPGA,EPSA )
         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_N(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               EPGA = AVG_Y(EP_G(IJK),EP_G(IJKN),J) 
               EPSA = AVG_Y(EP_S(IJK,1),EP_S(IJKN,1),J) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_N(IJK,1) = ZERO 
                  ELSE 
                     D_N(IJK,1) = P_SCALE*AXZ(IJK)*EPSA/((-A_M(IJK,0,1))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_N(IJK,1) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE CALC_D_N 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_d_t(A_m, VxF_gs, d_t, IER)                        C
!  Purpose: calculte coefficients linking velocity correction to       C
!           pressure correction -- Top                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
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
      SUBROUTINE CALC_D_T(A_M, VXF_GS, D_T, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE scales 
      USE compar    !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Volume x average at momentum cell centers
      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M)
!                      
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M)
!
!                      Average volume fraction at momentum cell centers
      DOUBLE PRECISION EPGA, EPSA
!
!                      F/(a0 + F), F/(a1 + F)
      DOUBLE PRECISION FoA0pF,     FoA1pF
!
!                      Indices
      INTEGER          I, J, K, IJK, IJKT
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
      IF (MOMENTUM_Z_EQ(0) .AND. MOMENTUM_Z_EQ(1)) THEN 
!
!//I? check any data dependency for I direction decomposition in following loop
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp  parallel do private( J, K, IJK, IJKT, EPGA, EPSA, FoA0pF, FoA1pF), &
!$omp&  schedule(static)
         DO IJK = ijkstart3, ijkend3 
            IF (IP_AT_T(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_T(IJK,0) = ZERO 
               D_T(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K) 
               EPSA = AVG_Z(EP_S(IJK,1),EP_S(IJKT,1),K) 
               IF (MODEL_B) THEN 
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
!
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_T(IJK,0) = P_SCALE*AXY(IJK)/((-A_M(IJK,0,0))-A_M(IJK,0,1&
                        )*FOA1PF) 
                     D_T(IJK,1) = FOA1PF*D_T(IJK,0) 
!
                  ELSE 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_T(IJK,0) = P_SCALE*AXY(IJK)/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_T(IJK,0) = ZERO 
                     ENDIF 
                     D_T(IJK,1) = ZERO 
                  ENDIF 
!
               ELSE                              !Model A 
!                 MFIX convention: center coeff is negative
                  IF (VXF_GS(IJK,1) > SMALL_NUMBER) THEN 
                     FOA0PF = VXF_GS(IJK,1)/((-A_M(IJK,0,0))+VXF_GS(IJK,1)) 
                     FOA1PF = VXF_GS(IJK,1)/((-A_M(IJK,0,1))+VXF_GS(IJK,1)) 
                     D_T(IJK,0) = P_SCALE*AXY(IJK)*(EPGA + EPSA*FOA1PF)/((-A_M(&
                        IJK,0,0))-A_M(IJK,0,1)*FOA1PF) 
                     D_T(IJK,1) = P_SCALE*AXY(IJK)*(EPSA + EPGA*FOA0PF)/((-A_M(&
                        IJK,0,1))-A_M(IJK,0,0)*FOA0PF) 
                  ELSE 
                     IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                        D_T(IJK,0) = P_SCALE*AXY(IJK)*EPGA/(-A_M(IJK,0,0)) 
                     ELSE 
                        D_T(IJK,0) = ZERO 
                     ENDIF 
                     IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                        D_T(IJK,1) = P_SCALE*AXY(IJK)*EPSA/(-A_M(IJK,0,1)) 
                     ELSE 
                        D_T(IJK,1) = ZERO 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_Z_EQ(0)) THEN 
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    

!$omp    parallel do &
!$omp&   private( IJK, I,J,K,IJKT, EPGA )
         DO IJK = ijkstart3, ijkend3
            IF (IP_AT_T(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_T(IJK,0) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,0)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_T(IJK,0) = P_SCALE*AXY(IJK)/((-A_M(IJK,0,0))+VXF_GS(IJK,&
                        1)) 
                  ELSE 
                     D_T(IJK,0) = P_SCALE*AXY(IJK)*EPGA/((-A_M(IJK,0,0))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_T(IJK,0) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ELSE IF (MOMENTUM_Z_EQ(1)) THEN
!// 350 1225 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3    
       
!$omp    parallel do &
!$omp&   private( IJK, I,J,K, IJKT, EPSA )
         DO IJK = ijkstart3, ijkend3  
            IF (IP_AT_T(IJK) .OR. MFLOW_AT_N(IJK)) THEN 
               D_T(IJK,1) = ZERO 
            ELSE 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               EPSA = AVG_Z(EP_S(IJK,1),EP_S(IJKT,1),K) 
!                 MFIX convention: center coeff is negative
               IF ((-A_M(IJK,0,1)) > SMALL_NUMBER) THEN 
                  IF (MODEL_B) THEN 
                     D_T(IJK,1) = ZERO 
                  ELSE 
                     D_T(IJK,1) = P_SCALE*AXY(IJK)*EPSA/((-A_M(IJK,0,1))+VXF_GS&
                        (IJK,1)) 
                  ENDIF 
               ELSE 
                  D_T(IJK,1) = ZERO 
               ENDIF 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE CALC_D_T 
