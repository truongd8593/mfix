!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_BC0                                                C
!  Purpose: This module does the initial setting of all boundary       C
!           conditions                                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check whether the velocity components have the correct     C
!           sign at MASS_INFLOW cells                                  C
!  Author: M. Syamlal                                 Date: 08-JUL-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, BC_DT_0, TIME, BC_Jet_g0,
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        BC_I_e, BC_PLANE, BC_EP_g, BC_P_g, BC_T_g,    C
!                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g,     C
!                        MMAX, BC_ROP_s, BC_U_s, BC_V_s, BC_W_s        C
!  Variables modified: BC_TIME, BC_V_g, I, J, K, IJK, EP_g, P_g, T_g,  C
!                      T_s, U_g, V_g, W_g, ROP_s, U_s, V_s, W_s,C
!                      M                                               C
!                                                                      C
!  Local variables: L, IJK1, IJK2, IJK3                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_BC0 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE compar       !//
      USE mpi_utility  !//
      USE physprop
      USE bc
      USE fldvar
      USE indices
      USE run
      USE funits 
      USE scales 
      USE boundfunijk 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
! 
!                      Local index for boundary condition 
      INTEGER          L 
! 
!                      indices 
      INTEGER          I, J, K, IJK, M, N 
! 
!                      Local index for setting U velocity b.c. 
      INTEGER          IJK1 
! 
!                      Local index for setting V velocity b.c. 
      INTEGER          IJK2 
! 
!                      Local index for setting W velocity b.c. 
      INTEGER          IJK3 
! 
!                      Dummy variable for gas pressure 
      DOUBLE PRECISION PJ 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG, CALC_MW 
!-----------------------------------------------
!----------------------------------------------
!   Temporary variable //SP - To gather fluid_at
!----------------------------------------------
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG_G
!//SP  Logical function to identify a fluid cell in global coordiantes
      LOGICAL          FLUID_AT_G

      INCLUDE 'sc_p_g1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'sc_p_g2.inc'
! Global function to determine FLUID_AT_G
      FLUID_AT_G(IJK)    = FLAG_G(IJK) .EQ. 1
!
!  Setup for cyclic boundary conditions
!
!//? 1026 check the pointer location for ijk1 on all PEs and adjust loop counter if necessary
!//SP Modifications so that the present computation is performed only on the root processor
!//SP and broadcasted everywhere!!

    IF(myPE.eq.root) then

!//SP Allocate temp var
      ALLOCATE(FLAG_G(IJKMAX3))
      CALL GATHER(FLAG, FLAG_G, root)
      IJK1 = FUNIJK_GL(IMAX1/2 + 1,JMAX1,KMAX1/2 + 1)
!      write(*,"('(PE ',I2,'): IJK1 = ',I5)") myPE,IJK1 !//AIKEPARDBG

!//SP Exact implementation as in the serial code. In the serial version CYCLE has to be
!     replaced by EXIT to have the same meaning as in the original version
!
      DO IJK = IJK1, ijkmax3 
         IF (FLUID_AT_G(IJK)) EXIT  
      END DO 
      
      IF (CYCLIC) THEN 
         IJK_P_G = IJK 
      ELSE 
         PJ = UNDEFINED 
         DO L = 1, DIMENSION_BC 
            IF (BC_DEFINED(L) .AND. BC_TYPE(L)=='P_OUTFLOW') PJ = BC_P_G(L) 
         END DO 
         IF (PJ == UNDEFINED) THEN 
            IJK_P_G = IJK 
         ELSE 
            IJK_P_G = UNDEFINED_I 
         ENDIF 
      ENDIF 

    ELSE ! myPE.eq.root

      ALLOCATE(FLAG_G(1))
      CALL GATHER(FLAG, FLAG_G, root)

    ENDIF ! myPE.eq.root

!//SP Deallocate storage

      DEALLOCATE(FLAG_G)

!//SP Broadcast the values

      CALL BCAST(PJ)
      CALL BCAST(IJK_P_G)
!
!  Set the boundary conditions.
!
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L)=='FREE_SLIP_WALL' .OR. BC_TYPE(L)=='NO_SLIP_WALL'&
                .OR. BC_TYPE(L)=='PAR_SLIP_WALL') CYCLE  
!
!         Initialization for time dependent mass inflow
!
            BC_JET_G(L) = UNDEFINED 
            IF (BC_DT_0(L) /= UNDEFINED) THEN 
               BC_TIME(L) = TIME + BC_DT_0(L) 
               BC_OUT_N(L) = 0 
               BC_MOUT_G(L) = ZERO 
               BC_VOUT_G(L) = ZERO 
               M = 1 
               IF (MMAX > 0) THEN 
                  BC_MOUT_S(L,:MMAX) = ZERO 
                  BC_VOUT_S(L,:MMAX) = ZERO 
                  M = MMAX + 1 
               ENDIF 
               BC_JET_G(L) = BC_JET_G0(L) 
               IF (BC_JET_G(L) /= UNDEFINED) THEN 
                  SELECT CASE (BC_PLANE(L))  
                  CASE ('W')  
                     BC_U_G(L) = BC_JET_G(L) 
                  CASE ('E')  
                     BC_U_G(L) = BC_JET_G(L) 
                  CASE ('S')  
                     BC_V_G(L) = BC_JET_G(L) 
                  CASE ('N')  
                     BC_V_G(L) = BC_JET_G(L) 
                  CASE ('B')  
                     BC_W_G(L) = BC_JET_G(L) 
                  CASE ('T')  
                     BC_W_G(L) = BC_JET_G(L) 
                  END SELECT 
               ENDIF 
            ELSE 
               BC_TIME(L) = UNDEFINED 
            ENDIF 
!
!         set the field variables with bc variables
!
            DO K = BC_K_B(L), BC_K_T(L) 
               DO J = BC_J_S(L), BC_J_N(L) 
                  DO I = BC_I_W(L), BC_I_E(L) 	
!// 360 1025 Check if current i,j,k resides on this PE		  	  
!// SP - Changed to two layers
		  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                     IJK = BOUND_FUNIJK(I,J,K) 
!//AIKEPARDBG
!              write(UNIT_LOG,"('P1:set_bc0 for EP_G: IJK= ',I5,'  I= ',I4,&
!	      & ' J=',I4,' K=',I4,' BC_EP_G=',E12.4,' EP_G=',E12.4)") &
!	          IJK,I,J,K,BC_EP_G(L),EP_G(IJK) !//AIKEPARDBG
		     
                     IF (.NOT.WALL_AT(IJK)) THEN 
!
                        IF (P_OUTFLOW_AT(IJK) .OR. OUTFLOW_AT(IJK)) THEN 
                           P_STAR(IJK) = ZERO 
                           P_G(IJK) = SCALE(BC_P_G(L)) 
!
                           IF (BC_EP_G(L) /= UNDEFINED) EP_G(IJK) = BC_EP_G(L) 
                           IF (BC_T_G(L) /= UNDEFINED) T_G(IJK) = BC_T_G(L) 
                           N = 1 
                           IF (NMAX(0) > 0) THEN 
                              WHERE (BC_X_G(L,:NMAX(0)) /= UNDEFINED) X_G(IJK,:&
                                 NMAX(0)) = BC_X_G(L,:NMAX(0)) 
                              N = NMAX(0) + 1 
                           ENDIF 
                           DO M = 1, MMAX 
                              IF (BC_ROP_S(L,M) /= UNDEFINED) ROP_S(IJK,M) = &
                                 BC_ROP_S(L,M) 
                              IF(BC_T_S(L,M)/=UNDEFINED)T_S(IJK,M)=BC_T_S(L,M) 
                              IF (BC_THETA_M(L,M) /= UNDEFINED) THETA_M(IJK,M)&
                                  = BC_THETA_M(L,M) 
                              N = 1 
                              IF (NMAX(M) > 0) THEN 
                                 WHERE (BC_X_S(L,M,:NMAX(M)) /= UNDEFINED) X_S(&
                                    IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M)) 
                                 N = NMAX(M) + 1 
                              ENDIF 
                           END DO 
                        ELSE 
!
                           P_STAR(IJK) = ZERO 
!
                           EP_G(IJK) = BC_EP_G(L) 
                           P_G(IJK) = SCALE(BC_P_G(L)) 
                           T_G(IJK) = BC_T_G(L) 
                           N = 1 
                           IF (NMAX(0) > 0) THEN 
                              X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0)) 
                              N = NMAX(0) + 1 
                           ENDIF 
                           DO M = 1, MMAX 
                              ROP_S(IJK,M) = BC_ROP_S(L,M) 
                              T_S(IJK,M) = BC_T_S(L,M) 
                              THETA_M(IJK,M) = BC_THETA_M(L,M) 
                              N = 1 
                              IF (NMAX(M) > 0) THEN 
                                 X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M)) 
                                 N = NMAX(M) + 1 
                              ENDIF 
                           END DO 
                           IJK1 = IJK 
                           IJK2 = IJK 
                           IJK3 = IJK 
                           SELECT CASE (BC_PLANE(L))  
                           CASE ('W')  
                              IJK1 = BOUND_FUNIJK(I - 1,J,K) 
                           CASE ('S')  
                              IJK2 = BOUND_FUNIJK(I,J - 1,K) 
                           CASE ('B')  
                              IJK3 = BOUND_FUNIJK(I,J,K - 1) 
                           END SELECT 
!
!             When the boundary plane is W, S, or B the velocity components
!             need to be set for both sides of the boundary cell.
                           U_G(IJK) = BC_U_G(L) 
                           V_G(IJK) = BC_V_G(L) 
                           W_G(IJK) = BC_W_G(L) 
                           U_G(IJK1) = BC_U_G(L) 
                           V_G(IJK2) = BC_V_G(L) 
                           W_G(IJK3) = BC_W_G(L) 
!
                           M = 1 
                           IF (MMAX > 0) THEN 
                              U_S(IJK,:MMAX) = BC_U_S(L,:MMAX) 
                              V_S(IJK,:MMAX) = BC_V_S(L,:MMAX) 
                              W_S(IJK,:MMAX) = BC_W_S(L,:MMAX) 
                              U_S(IJK1,:MMAX) = BC_U_S(L,:MMAX) 
                              V_S(IJK2,:MMAX) = BC_V_S(L,:MMAX) 
                              W_S(IJK3,:MMAX) = BC_W_S(L,:MMAX) 
                              M = MMAX + 1 
                           ENDIF 
!                                                ! 'MASS_OUTFLOW'
                           IF (.NOT.MASS_OUTFLOW_AT(IJK)) THEN 
!
                              IF (MW_AVG == UNDEFINED) MW_MIX_G(IJK) = CALC_MW(&
                                 X_G,DIMENSION_3,IJK,NMAX(0),MW_G) 
                              IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(MW_MIX_G&
                                 (IJK),P_G(IJK),T_G(IJK)) 
                              ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
                           ENDIF 
                        ENDIF 
                     ENDIF 
!//AIKEPARDBG		     
!              write(UNIT_LOG,"('P2:set_bc0 for EP_G: IJK= ',I5,'  I= ',I4,&
!	      & ' J=',I4,' K=',I4,' BC_EP_G=',E12.4,' EP_G=',E12.4)") &
!	          IJK,I,J,K,BC_EP_G(L),EP_G(IJK)  !//AIKEPARDBG
		     
                  END DO 
               END DO 
            END DO 
         ENDIF 
      END DO 
!//? check for the variables that require communication due to the cyclic BC along k direction 
      RETURN  
      END SUBROUTINE SET_BC0 
