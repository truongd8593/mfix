!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_WALL_BC(IER)                                       C
!  Purpose: Set wall boundary conditions                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Add calculations for mass outflow boundary condition       C
!  Author: M. Syamlal                                 Date: 23-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose:Revised for MFIX 2.0.  This sub routine is different from   C
!          old set_wall_bc.
!  Author: M. Syamlal                                 Date: 18-JUL-96  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_I_w, BC_I_e, BC_J_s, BC_J_n,   C
!                        BC_K_b, BC_K_t, BC_TYPE, TIME, DT, BC_TIME,   C
!                        BC_V_g, BC_V_gh, BC_V_gl, BC_DT_l, BC_DT_h,   C
!                        BC_PLANE, IMAX2, JMAX2, KMAX2                 C
!  Variables modified: BC_V_g, BC_TIME, I, J, K, IJK, V_g              C
!                                                                      C
!  Local variables: L, IJK2, I1, I2, J1, J2, K1, K2                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_WALL_BC(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE bc
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE funits 
      USE compar    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      error index 
      INTEGER          IER 
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
!                      Index for setting V velocity b.c. 
      INTEGER          IJK2 
! 
!                      indices 
      INTEGER          IJK, IPJK, M 
! 
!                      Starting I index 
      INTEGER          I1 
! 
!                      Ending I index 
      INTEGER          I2 
! 
!                      Starting J index 
      INTEGER          J1 
! 
!                      Ending J index 
      INTEGER          J2 
! 
!                      Starting K index 
      INTEGER          K1 
! 
!                      Ending K index 
      INTEGER          K2 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Set the boundary conditions
!
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
!
!  The range of boundary cells
!
            I1 = BC_I_W(L) 
            I2 = BC_I_E(L) 
            J1 = BC_J_S(L) 
            J2 = BC_J_N(L) 
            K1 = BC_K_B(L) 
            K2 = BC_K_T(L) 
!
            SELECT CASE (BC_TYPE(L))  
            CASE ('FREE_SLIP_WALL')  
!
!  Set velocities for the range of boundary cells.  Use 1.0 as the sign
!  to make the velocity at the cell equal to that at the fluid cell.
!
               CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, BC_JJ_PS(L), ONE) 
!
            CASE ('NO_SLIP_WALL')  
!
!  Set velocities for the range of boundary cells.  Use -1.0 as the sign
!  to make the velocity at the cell equal to that at the fluid cell.
!
               CALL SET_WALL_BC1 (I1, I2, J1, J2, K1, K2, BC_JJ_PS(L), (-ONE)) 
!
            CASE ('PAR_SLIP_WALL')  
!             updating the boundary velocity may improve convergence
            END SELECT 
         ENDIF 
      END DO 
      
      K1 = 1 

      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
            IF(K1.NE.KSTART2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1&
               , 0, (-ONE)) 
         END DO 
      END DO 

      K1 = KMAX2 
      DO J1 = JSTART3, JEND3
         DO I1 = ISTART3, IEND3
!// 
            IF(K1.NE.KEND2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
	    
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1&
               , 0, (-ONE)) 
         END DO 
      END DO 
      J1 = 1 
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
!//
            IF(J1.NE.JSTART2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
	    
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1&
               , 0, (-ONE)) 
         END DO 
      END DO 
      J1 = JMAX2 
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
!//
            IF(J1.NE.JEND2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
	    
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1&
               , 0, (-ONE)) 
         END DO 
      END DO 
      I1 = 1 
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
!//
            IF(I1.NE.ISTART2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
            IJK = FUNIJK(I1,J1,K1) 

!
            IF (DEFAULT_WALL_AT(IJK)) THEN 
!
               IF (NS_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1, &
                  0, (-ONE)) 
               IF(FS_WALL_AT(IJK))CALL SET_WALL_BC1(I1,I1,J1,J1,K1,K1,0,ONE) 
            ENDIF 
!
!  For cylindrical coordinates the azimuthal component should be zero at center
!
            IF (CYLINDRICAL .AND. XMIN==ZERO) THEN 
               IPJK = IP_OF(IJK) 
               W_G(IJK) = -W_G(IPJK) 
               M = 1 
               IF (MMAX > 0) THEN 
                  W_S(IJK,:MMAX) = -W_S(IPJK,:MMAX) 
                  M = MMAX + 1 
               ENDIF 
            ENDIF 
         END DO 
      END DO 
      I1 = IMAX2 
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
!//
            IF(I1.NE.IEND2)   EXIT
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
	    
            IJK = FUNIJK(I1,J1,K1) 
            IF (DEFAULT_WALL_AT(IJK)) CALL SET_WALL_BC1 (I1, I1, J1, J1, K1, K1&
               , 0, (-ONE)) 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE SET_WALL_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_WALL_BC1(II1, II2, JJ1, JJ2, KK1, KK2, BC_JJ_PSL, &C
!                                                SIGN)                 C
!  Purpose: Set U, V, and W components for the specified cells by      C
!           copying the same or negative values from near by fluid cellC
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                 (name changed to set_wall_bc1)                       C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: V_s, W_s, U_s                                 C
!                                                                      C
!  Variables modified: I, J, K, V_g, W_g, U_g                          C
!                                                                      C
!  Local variables: SIGN, LWALL, LFLUID, I1, I2, J1, J2, K1, K2        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_WALL_BC1(II1, II2, JJ1, JJ2, KK1, KK2, &
                                           BC_JJ_PSL, SIGN) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE bc 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run 
      USE funits 
      USE compar  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Starting I index 
      INTEGER          II1 
! 
!                      Ending I index 
      INTEGER          II2 
! 
!                      Starting J index 
      INTEGER          JJ1 
! 
!                      Ending J index 
      INTEGER          JJ2 
! 
!                      Starting K index 
      INTEGER          KK1 
! 
!                      Ending K index 
      INTEGER          KK2 
  
!                      Johnson-Jackson boundary condition: 0= no, 1=yes 
      INTEGER          BC_JJ_PSL 
! 
!                      Sign with legal values +1 or -1 
      DOUBLE PRECISION SIGN 
! 
!                      Local indices near wall cell 
      INTEGER          I, J, K 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP 
      INTEGER          I1, I2, J1, J2, K1, K2
  
!                      Locall index for a fluid cell near the wall cell 
      INTEGER          LFLUID 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!// Limit I1, I2 and all to local processor first ghost layer

      I1 = II1
      I2 = II2
      J1 = JJ1
      J2 = JJ2
      K1 = KK1
      K2 = KK2

      IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
      IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
      IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
      IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
      IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
      IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
               IJK = FUNIJK(I,J,K) 
               IF (WALL_AT(IJK)) THEN 
                  IMJK = IM_OF(IJK) 
                  IJMK = JM_OF(IJK) 
                  IJKM = KM_OF(IJK) 
                  IPJK = IP_OF(IJK) 
                  IJPK = JP_OF(IJK) 
                  IJKP = KP_OF(IJK) 
!
!         Fluid cell at West
!
                  IF (.NOT.WALL_AT(IMJK)) THEN 
                     LFLUID = IMJK 
!                                                   Wall cell at North
                     IF (WALL_AT(IJPK)) THEN 
                        V_G(IJK) = SIGN*V_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(V_S,IJK,SIGN,V_S,LFLUID) 
                     ENDIF 
!                                                   Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN 
                        W_G(IJK) = SIGN*W_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(W_S,IJK,SIGN,W_S,LFLUID) 
                     ENDIF 
                  ENDIF 
!
!         Fluid cell at East
!
                  IF (.NOT.WALL_AT(IPJK)) THEN 
                     LFLUID = IPJK 
!                                                   Wall cell at North
                     IF (WALL_AT(IJPK)) THEN 
                        V_G(IJK) = SIGN*V_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(V_S,IJK,SIGN,V_S,LFLUID) 
                     ENDIF 
!                                                   Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN 
                        W_G(IJK) = SIGN*W_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(W_S,IJK,SIGN,W_S,LFLUID) 
                     ENDIF 
                  ENDIF 
!
!         Fluid cell at South
!
                  IF (.NOT.WALL_AT(IJMK)) THEN 
                     LFLUID = IJMK 
!                                                   Wall cell at East
                     IF (WALL_AT(IPJK)) THEN 
                        U_G(IJK) = SIGN*U_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(U_S,IJK,SIGN,U_S,LFLUID) 
                     ENDIF 
!                                                   Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN 
                        W_G(IJK) = SIGN*W_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(W_S,IJK,SIGN,W_S,LFLUID) 
                     ENDIF 
                  ENDIF 
!
!         Fluid cell at North
!
                  IF (.NOT.WALL_AT(IJPK)) THEN 
                     LFLUID = IJPK 
!                                                   Wall cell at East
                     IF (WALL_AT(IPJK)) THEN 
                        U_G(IJK) = SIGN*U_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(U_S,IJK,SIGN,U_S,LFLUID) 
                     ENDIF 
!                                                   Wall cell at Top
                     IF (WALL_AT(IJKP)) THEN 
                        W_G(IJK) = SIGN*W_G(LFLUID) 
                        IF(BC_JJ_PSL==0)CALL EQUAL(W_S,IJK,SIGN,W_S,LFLUID) 
                     ENDIF 
                  ENDIF 
                  IF (DO_K) THEN 
!
!           Fluid cell at Bottom
!
                     IF (.NOT.WALL_AT(IJKM)) THEN 
                        LFLUID = IJKM 
!                                                   Wall cell at East
                        IF (WALL_AT(IPJK)) THEN 
                           U_G(IJK) = SIGN*U_G(LFLUID) 
                           IF (BC_JJ_PSL == 0) CALL EQUAL (U_S, IJK, SIGN, U_S&
                              , LFLUID) 
                        ENDIF 
!                                                   Wall cell at North
                        IF (WALL_AT(IJPK)) THEN 
                           V_G(IJK) = SIGN*V_G(LFLUID) 
                           IF (BC_JJ_PSL == 0) CALL EQUAL (V_S, IJK, SIGN, V_S&
                              , LFLUID) 
                        ENDIF 
                     ENDIF 
!
!           Fluid cell at Top
!
                     IF (.NOT.WALL_AT(IJKP)) THEN 
                        LFLUID = IJKP 
!                                                   Wall cell at East
                        IF (WALL_AT(IPJK)) THEN 
                           U_G(IJK) = SIGN*U_G(LFLUID) 
                           IF (BC_JJ_PSL == 0) CALL EQUAL (U_S, IJK, SIGN, U_S&
                              , LFLUID) 
                        ENDIF 
!                                                   Wall cell at North
                        IF (WALL_AT(IJPK)) THEN 
                           V_G(IJK) = SIGN*V_G(LFLUID) 
                           IF (BC_JJ_PSL == 0) CALL EQUAL (V_S, IJK, SIGN, V_S&
                              , LFLUID) 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE SET_WALL_BC1 


!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kstart3,kend3      
!// 360 Check if i,j,k resides on current processor
!//     Limit I1, I2 and all to local processor first ghost layer
