!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for U_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.  The center        C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE SOURCE_U_S(A_M, B_M, IER) 
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
      USE compar        !//d
      USE sendrecv    !// 400
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
      INTEGER          I, IJK, IJKE, IJKM, IPJK, IPJKM 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Internal surface number 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
! 
!                      average volume fraction 
      DOUBLE PRECISION EPSA 
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
! 
!                      Average density difference 
      DOUBLE PRECISION dro1, dro2, droa 
! 
!                      Average quantities 
      DOUBLE PRECISION wse, EPMUGA 
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
      DOUBLE PRECISION V0, Vmt, Vbf, Vcf, Vtza 
! 
!                      error message 
      CHARACTER*80     LINE(2) 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
!
!// 350 1213 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

!$omp  parallel do private( IJK, IJKE, ISV, Sdp, Sdps, V0, Vmt, Vbf, &
!$omp&  I,PGE,DRO1,DRO2,DROA, IJKM,IPJK,IPJKM,  WSE,VCF,EPMUGA,VTZA, &
!$omp&  EPSA, ROPSA, LINE) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3 
!
!           Wall or impermeable internal surface
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
               IF (IP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (SIP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_E(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 
!
!           dilute flow
               ELSE IF (EPSA <= DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
!
                  IF (EP_S(WEST_OF(IJK),M) > DIL_EP_S) THEN 
                     A_M(IJK,W,M) = ONE 
                  ELSE IF (EP_S(EAST_OF(IJK),M) > DIL_EP_S) THEN 
                     A_M(IJK,E,M) = ONE 
                  ELSE 
                     B_M(IJK,M) = -U_S(IJK,M) 
                  ENDIF 
!
!           Normal case
               ELSE 
!
!           Surface forces
!
!             Pressure terms
                  PGE = P_G(IJKE) 
                  IF (CYCLIC_X_PD) THEN 
                     IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X 
                  ENDIF 
                  IF (MODEL_B) THEN 
                     SDP = ZERO 
!
                  ELSE 
                     SDP = -P_SCALE*EPSA*(PGE - P_G(IJK))*AYZ(IJK) 
!
                  ENDIF 
!
                  IF (CLOSE_PACKED(M)) THEN 
                     SDPS = -((P_S(IJKE,M)-P_S(IJK,M))+(P_STAR(IJKE)-P_STAR(IJK&
                        )))*AYZ(IJK) 
                  ELSE 
                     SDPS = -(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK) 
                  ENDIF 
!
!             Shear stress terms
!
!           Volumetric forces
                  ROPSA = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I) 
!
!             Previous time step
                  V0 = AVG_X(ROP_SO(IJK,M),ROP_SO(IJKE,M),I)*ODT 
!
!             Interphase mass transfer
                  VMT = AVG_X(SUM_R_S(IJK,M),SUM_R_S(IJKE,M),I) 
!
!             Body force
                  IF (MODEL_B) THEN 
                     DRO1 = (RO_S(M)-RO_G(IJK))*EP_S(IJK,M) 
                     DRO2 = (RO_S(M)-RO_G(IJKE))*EP_S(IJKE,M) 
                     DROA = AVG_X(DRO1,DRO2,I) 
!
                     VBF = DROA*BFX_S(IJK,M) 
!
                  ELSE 
                     VBF = ROPSA*BFX_S(IJK,M) 
                  ENDIF 
!
!
!           Special terms for cylindrical coordinates
                  IF (CYLINDRICAL) THEN 
!
!             centrifugal force
                     IJKM = KM_OF(IJK) 
                     IPJK = IP_OF(IJK) 
                     IPJKM = IP_OF(IJKM) 
!//? make sure W_G(IJKM) for k-direction decomposition is up to date on PEs
!//I? make sure W_G(IPJK) for i-direction decomp. is up to date on PEs	       		     
                     WSE = AVG_X(HALF*(W_S(IJK,M)+W_S(IJKM,M)),HALF*(W_S(IPJK,M&
                        )+W_S(IPJKM,M)),I) 
                     VCF = ROPSA*WSE**2*OX_E(I) 
!
!             Tau_zz/X
                     EPMUGA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                     VTZA = 2.*EPMUGA*OX_E(I)*OX_E(I) 
                  ELSE 
                     VCF = ZERO 
                     VTZA = ZERO 
                  ENDIF 
!
!             Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(&
                     IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+ZMAX(VMT)+VTZA)*&
                     VOL_U(IJK)) 
                  B_M(IJK,M) = -(SDP + SDPS + TAU_U_S(IJK,M)+((V0+ZMAX((-VMT)))&
                     *U_SO(IJK,M)+VBF+VCF)*VOL_U(IJK))+B_M(IJK,M) 
               ENDIF 
            END DO 
            CALL SOURCE_U_S_BC (A_M, B_M, M, IER) 
         ENDIF 

!//? verify the location of the COMM
!//? verify whether in the do M=1,MMAX lopp for each M, or out of the loop
!!!      CALL SEND_RECV(A_M, 2)
!!!      CALL SEND_RECV(B_M, 2)

      END DO 

      
      RETURN  
      END SUBROUTINE SOURCE_U_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s_BC(A_m, B_m, M, IER)                        C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
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
      SUBROUTINE SOURCE_U_S_BC(A_M, B_M, M, IER) 
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
      USE compar        !//d
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
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
!
!  Set the default boundary conditions
!
!//? see if the questions in source_u_g for similar location required any
!//? changes. If yes then implement it here also.
      IF (DO_K) THEN 
         K1 = 1 
!// 350 1208 change do loop limits: 1,jmax2->jmin3,jmax3	 	 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
!// 360 1208 Check if current i,j,k resides on this PE
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
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
!// 350 1208 change do loop limits: 1,jmax2->jmin3,jmax3	 	 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
!// 360 1208 Check if current i,j,k resides on this PE
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = -ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
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
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft DO_K branch in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG

      J1 = 1 
!// 350 1208 change do loop limits: 1,jmax2->jmin3,jmax3	 	 	       
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
!// 360 1208 Check if current i,j,k resides on this PE
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = -ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
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
!// 350 1208 change do loop limits: 1,jmax2->jmin3,jmax3	 	 	       
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
!// 360 1208 Check if current i,j,k resides on this PE
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
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

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft J branch in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1208 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = -ONE 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = -ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = -ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = -ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 1 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ELSE                              !Johnson and Jackson partial slip 
!
!//? need to go over the subroutine to see any further modifications necessary?
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1208 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = ONE 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 2 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1208 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     			
                           IJK = FUNIJK(I,J,K) 
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
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,N,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J)) 
                                 A_M(IJK,N,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,S,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,S,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,T,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*&
                                    OX_E(I)) 
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,B,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX_E(I)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 3 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1208 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 4 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1208 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IM = IM1(I) 
                           IMJK = IM_OF(IJK) 
                           A_M(IMJK,E,M) = ZERO 
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM)) 
                           A_M(IMJK,N,M) = ZERO 
                           A_M(IMJK,S,M) = ZERO 
                           A_M(IMJK,T,M) = ZERO 
                           A_M(IMJK,B,M) = ZERO 
                           A_M(IMJK,0,M) = -ONE 
                           B_M(IMJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 5 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ELSE IF (BC_PLANE(L) == 'E') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
!// 360 1229 Check if current i,j,k resides on this PE		     
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     

                           IJK = FUNIJK(I,J,K) 
!
                           IP = IP1(I) 
                           IPJK = IP_OF(IJK) 
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I) 
                           A_M(IPJK,W,M) = ZERO 
                           A_M(IPJK,N,M) = ZERO 
                           A_M(IPJK,S,M) = ZERO 
                           A_M(IPJK,T,M) = ZERO 
                           A_M(IPJK,B,M) = ZERO 
                           A_M(IPJK,0,M) = -ONE 
                           B_M(IPJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 6 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
		  
               ENDIF 
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
!// 360 1229 Check if current i,j,k resides on this PE		     
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -U_S(IJK,M) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_S(IJKW,M) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): aft Part 7 in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
	       
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_U_S_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_U_s(I1, I2, J1, J2, K1, K2, L, M, A_m, b_m)      C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava,                Date: 14-APR-98  C
!          Princeton University                                        C
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
      SUBROUTINE JJ_BC_U_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                       JM, KM
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
!                      coefficients for granular bc
      DOUBLE PRECISION hw, gw, cw
!      
 !-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!
      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
!// 360 1229 Check if current i,j,k resides on this PE		     
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     

               IJK = FUNIJK(I,J,K) 
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
!
               IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,N,M) = -ONE 
                  ELSE 
! start anuj 4/20
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, NORTH_OF(IJK), 'N', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'U', M, &
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,N,M) = -(HALF*HW - ODY_N(J)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(J)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,S,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, SOUTH_OF(IJK), 'S', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'U', M, &
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,S,M) = -(HALF*HW - ODY_N(JM)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(JM)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,T,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, TOP_OF(IJK), 'T', 'U', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'U', M, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,B,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, BOTTOM_OF(IJK), 'B', 'U', L&
                           , M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'U', M, &
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 

!//AIKEPARDBG
!    write(*,"('(PE ',I2,'): eof JJ_BC_U_S in source_u_s')") myPE  !//AIKEPARDBG
!    write(*,"('(PE ',I2,'): eof JJ_BC_U_S in source_u_s')") myPE  !//AIKEPARDBG
!    call mfix_exit(myPE)     !//AIKEPARDBG
      
      RETURN  
      END SUBROUTINE JJ_BC_U_S 
