!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BC_phi(BC_phif, BC_Phiw, BC_hw_Phi, BC_C_Phi, M,       C
!                              A_m, B_m, IER)                          C
!  Purpose: Set up the phi boundary conditions                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
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
      SUBROUTINE BC_PHI(BC_PHIF,BC_PHIW,BC_HW_PHI,BC_C_PHI,M,A_M,B_M,IER) 
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
      USE toleranc 
      USE run
      USE physprop
      USE fldvar
      USE visc_s
      USE geometry
      USE output
      USE indices
      USE bc
      USE compar    
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
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                      IM, JM, KM
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
!                      Boundary conditions
      DOUBLE PRECISION BC_phif(DIMENSION_BC), BC_Phiw(DIMENSION_BC), &
                      BC_hw_Phi(DIMENSION_BC), BC_C_Phi(DIMENSION_BC)
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Set up the default walls as non-conducting.
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3, jmax3 
            DO I1 = imin3, imax3 
!//SP ----> Not very efficient - can check only for K
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1) 
               IF (DEFAULT_WALL_AT(IJK)) THEN 
!
                  A_M(KP_OF(IJK),B,M) = ZERO 
!
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
         K1 = KMAX2 
         DO J1 = jmin3, jmax3 
            DO I1 = imin3, imax3 
!//SP ----> Not very efficient - can check only for K
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1) 
               IF (DEFAULT_WALL_AT(IJK)) THEN 
!
                  A_M(KM_OF(IJK),T,M) = ZERO 
!
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
!
      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3 
!//SP ----> Not very efficient - can check only for J
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(JP_OF(IJK),S,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      J1 = JMAX2 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3 
!//SP ----> Not very efficient - can check only for J
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(JM_OF(IJK),N,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 

      I1 = imin2 
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
!//SP ----> Not very efficient - can check only for I
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(IP_OF(IJK),W,M) = ZERO 
!
               A_M(IJK,E,M) = ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 

      I1 = IMAX2 
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
!//SP ----> Not very efficient - can check only for I
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               A_M(IM_OF(IJK),E,M) = ZERO 
!
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      !first set the bc for walls then overwrite where ever inflow/outflows are
      !defined so that the order in which the bcs are defined in the data file
      !does not matter.  Here set wall bcs . . .
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L)=='NO_SLIP_WALL' .OR. BC_TYPE(L)=='FREE_SLIP_WALL'&
                .OR. BC_TYPE(L)=='PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 		     
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IM = IM1(I) 
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(EAST_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,E,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(I)) 
                              A_M(IJK,E,M) = -(HALF*BC_HW_PHI(L)-ODX_E(I)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,W,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,W,M) = -(HALF*BC_HW_PHI(L)-ODX_E(IM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(IM)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_PHI(L)-ODY_N(J)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,S,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,S,M) = -(HALF*BC_HW_PHI(L)-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(JM)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,T,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(K)) 
                              A_M(IJK,T,M)=-(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(K)) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN 
                              A_M(IJK,B,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_PHIW(L) 
                           ELSE 
                              A_M(IJK,B,M) = -(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(KM&
                                 )) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(KM&
                                 )) 
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L&
                                 )) 
                           ENDIF 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      
      
      !. . . then set bcs for non-wall cells
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L)=='NO_SLIP_WALL' .OR. BC_TYPE(L)=='FREE_SLIP_WALL'&
                .OR. BC_TYPE(L)=='PAR_SLIP_WALL') THEN 
               !Dummy statement to do nothing.  The bcs were set in the previous loop
		I1 = 1  
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW' .OR. &
                  BC_TYPE(L)=='MASS_OUTFLOW' .OR. BC_TYPE(L)=='OUTFLOW') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
 	       	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
!
                        SELECT CASE (BC_PLANE(L))  
                        CASE ('E')  
                           A_M(IJK,E,M) = ONE 
                        CASE ('W')  
                           A_M(IJK,W,M) = ONE 
                        CASE ('N')  
                           A_M(IJK,N,M) = ONE 
                        CASE ('S')  
                           A_M(IJK,S,M) = ONE 
                        CASE ('T')  
                           A_M(IJK,T,M) = ONE 
                        CASE ('B')  
                           A_M(IJK,B,M) = ONE 
                        END SELECT 
!
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                     END DO 
                  END DO 
               END DO 
            ELSE 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
 	                IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -BC_PHIF(L) 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 


      RETURN  
      END SUBROUTINE BC_PHI 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 1206 change do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
