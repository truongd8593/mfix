!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GEOMETRY                                           C
!  Purpose: Calculate X, X_E,  oX, oX_E                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Fix bugs                                                   C
!  Author: M. Syamlal                                 Date: 10-FEB-92  C
!  Revision Number: 2                                                  C
!  Purpose: Include logic for Variable Grid Spacing capability         C
!  Author: W. Rogers                                  Date: 06-APR-92  C
!  Revision Number: 3                                                  C
!  Purpose: Add oX, and oX_E calculations                              C
!  Author: M. Syamlal                                 Date: 8-MAY-92   C
!  Revision Number: 4                                                  C
!  Purpose: Add FX, FX_bar, FX_E, FX_E_bar, FY_N, FY_N_bar, FZ_T, and  C
!           FZ_T_bar calculations                                      C
!  Author: M. Syamlal                                 Date: 31-AUG-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: COORDINATES, IMAX2, DX, JMAX2, DY, KMAX2,     C
!                        DZ,                                           C
!                                                                      C
!  Variables modified: X, X_E, I,                                      C
!                      J, K,  oX,                                      C
!                      oX_E, FX, FX_bar, FX_E, FX_E_bar, FY_N,         C
!                      FY_N_bar, FZ_T, FZ_T_bar                        C
!                                                                      C
!  Local variables: DX_E, DY_N, DZ_T                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_GEOMETRY 
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
      USE run
      USE geometry
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
!                      Indices 
      INTEGER          I, J, K 
! 
!     X-direction dimension of U-momentum cell 
      DOUBLE PRECISION DX_E 
! 
!     Y-direction dimension of V-momentum cell 
      DOUBLE PRECISION DY_N 
! 
!     Z-direction dimension of W-momentum cell 
      DOUBLE PRECISION DZ_T 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
!
!
      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE. 
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE. 
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE. 
      CYCLIC = CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z 
      IF (CYLINDRICAL .AND. COMPARE(ZLENGTH,8.D0*ATAN(ONE)) .AND. DO_K) &
         CYCLIC_Z = .TRUE. 
!
      IF (CYCLIC_X) THEN 
         DX(1) = DX(IMAX1) 
         DX(IMAX2) = DX(IMIN1) 
      ENDIF 
      IF (CYCLIC_Y) THEN 
         DY(1) = DY(JMAX1) 
         DY(JMAX2) = DY(JMIN1) 
      ENDIF 
!//D 300 0912 For CYCLIC_Z, dz(1) of PE 0 = dz(kmax1) of PE 1
!//D          As dz() is the global array for all PEs no modification necessary
      IF (CYCLIC_Z) THEN 
         DZ(1) = DZ(KMAX1) 
         DZ(KMAX2) = DZ(KMIN1) 
      ENDIF 
!
      IF (COORDINATES == 'CARTESIAN') THEN 
         I = 1 
         IF (IMAX2 > 0) THEN 
            X(:IMAX2) = ONE 
            X_E(:IMAX2) = ONE 
            OX(:IMAX2) = ONE 
            OX_E(:IMAX2) = ONE 
            ODX(:IMAX2) = ONE/DX(:IMAX2) 
            I = IMAX2 + 1 
         ENDIF 
      ELSE IF (CYLINDRICAL) THEN 
         IF (XMIN == ZERO) THEN 
            ODX(1) = ONE/DX(1) 
            OX(1) = UNDEFINED 
            OX_E(1) = UNDEFINED 
            IF (DO_I) THEN 
               X(1) = -HALF*DX(1) 
               X_E(1) = 0.0 
            ELSE 
               X(1) = HALF*DX(1) 
               X_E(1) = DX(1) 
            ENDIF 
         ELSE 
            IF (DO_I) THEN 
               X_E(1) = XMIN 
               X(1) = XMIN - HALF*DX(1) 
            ELSE 
               X_E(1) = XMIN + DX(1) 
               X(1) = XMIN + HALF*DX(1) 
            ENDIF 
            OX(1) = ONE/X(1) 
            OX_E(1) = ONE/X_E(1) 
            ODX(1) = ONE/DX(1) 
         ENDIF 
!                                                #1 add the DO_I IF block
         IF (DO_I) THEN 
            DO I = IMIN1, IMAX2 
               X(I) = X(I-1) + (DX(I-1)+DX(I))/2. 
               X_E(I) = X_E(I-1) + DX(I) 
               OX(I) = ONE/X(I) 
               OX_E(I) = ONE/X_E(I) 
               ODX(I) = ONE/DX(I) 
            END DO 
         ENDIF 
      ENDIF 
!
      J = 1 
      IF (JMAX2 > 0) THEN 
         ODY(:JMAX2) = ONE/DY(:JMAX2) 
         J = JMAX2 + 1 
      ENDIF 

!// 200 0920 Changed the limit from KMAX2--> KMAX3
      DO K = 1, KMAX3 
!
         IF (K == 1) THEN 
            Z(K) = ZERO - HALF*DZ(K) 
            Z_T(K) = ZERO 
!// 200 0920 added initializations to take care of z(KMIN3)	    
	    Z(K-1) =Z(K)
            Z_T(K-1) = Z_T(K) 	    
!
!// 200 0920 added initializations to take care of z(KMAX3)	    
         ELSE IF (K == KMAX3) THEN
	    Z(K) =Z(K-1)
            Z_T(K) = Z_T(K-1) 	    
	    
         ELSE
            Z(K) = Z_T(K-1) + HALF*DZ(K) 
            Z_T(K) = Z_T(K-1) + DZ(K) 
!
         ENDIF 
!
         ODZ(K) = ONE/DZ(K) 
      END DO 
      DX_E = HALF*(DX(1)+DX(IMIN1)) 
      DY_N = HALF*(DY(1)+DY(JMIN1)) 
      DZ_T = HALF*(DZ(1)+DZ(KMIN1)) 
!
!
      ODX_E(1) = ONE/DX_E 
      ODY_N(1) = ONE/DY_N 
      ODZ_T(1) = ONE/DZ_T 
      FX(1) = HALF 
      FX_BAR(1) = HALF 
      FX_E(1) = HALF 
      FX_E_BAR(1) = HALF 
      FY_N(1) = HALF 
      FY_N_BAR(1) = HALF 
      FZ_T(1) = HALF 
      FZ_T_BAR(1) = HALF 
!       ..........................................
!       Look at 2 through IMAX1 U-momentum cells
      IF (DO_I) THEN 
         DO I = IMIN1, IMAX1 
            DX_E = HALF*(DX(I+1)+DX(I)) 
            ODX_E(I) = ONE/DX_E 
!            FX(I)       = DX_E * X_E(I) /
!     &             (DX_E * X_E(I) + HALF * (DX(I) + DX(I-1)) * X_E(I-1))
            FX(I) = HALF 
            FX_BAR(I) = ONE - FX(I) 
!            FX_E(I)     = X(I+1) * DX(I+1) /
!     &                  ( X(I+1) * DX(I+1) + X(I) * DX(I) )
            FX_E(I) = DX(I+1)/(DX(I+1)+DX(I)) 
            FX_E_BAR(I) = ONE - FX_E(I) 
         END DO 
      ENDIF 
!
!       Look at 2 through JMAX1 V-momentum cells
      IF (DO_J) THEN 
         DO J = JMIN1, JMAX1 
            DY_N = HALF*(DY(J+1)+DY(J)) 
            ODY_N(J) = ONE/DY_N 
            FY_N(J) = DY(J+1)/(DY(J+1)+DY(J)) 
            FY_N_BAR(J) = ONE - FY_N(J) 
         END DO 
      ENDIF 
!
!       Look at 2 through KMAX1 W-momentum cells
      IF (DO_K) THEN 
!//D 300 0912 no changes in the limits as they run over ACTIVE cells ONLY
         DO K = KMIN1, KMAX1 
            DZ_T = HALF*(DZ(K+1)+DZ(K)) 
            ODZ_T(K) = ONE/DZ_T 
            FZ_T(K) = DZ(K+1)/(DZ(K+1)+DZ(K)) 
            FZ_T_BAR(K) = ONE - FZ_T(K) 
         END DO 
      ENDIF 

!       ..........................................
!       Look at last U-, V-, and W-momentum cells
      DX_E = DX(IMAX2) 
      DY_N = DY(JMAX2) 
      DZ_T = DZ(KMAX2) 
      ODX_E(IMAX2) = ONE/DX_E 
      ODY_N(JMAX2) = ONE/DY_N 
      ODZ_T(KMAX2) = ONE/DZ_T 
!//S2D for 2D/3D decomp. do similar add ons for JMAX3, IMAX3 as done in KMAX3 below      
      FX(IMAX2) = HALF 
      FX_BAR(IMAX2) = HALF 
      FX_E(IMAX2) = HALF 
      FX_E_BAR(IMAX2) = HALF 
      
      FY_N(JMAX2) = HALF 
      FY_N_BAR(JMAX2) = HALF      

      FZ_T(KMAX2) = HALF 
      FZ_T_BAR(KMAX2) = HALF 
!// 200 0920 need to update values for KMAX3 also             
      FZ_T(KMAX3) = HALF 
      FZ_T_BAR(KMAX3) = HALF 
      
      IF (CYCLIC_X) THEN 
         FX_E(1) = FX_E(IMAX1) 
         FX_E_BAR(1) = FX_E_BAR(IMAX1) 
      ENDIF 
      IF (CYCLIC_Y) THEN 
         FY_N(1) = FY_N(JMAX1) 
         FY_N_BAR(1) = FY_N_BAR(JMAX1) 
      ENDIF 
!//? should we update the additional ghosts for CYCLIC_Z?      
      IF (CYCLIC_Z) THEN 
         FZ_T(1) = FZ_T(KMAX1) 
         FZ_T_BAR(1) = FZ_T_BAR(KMAX1) 
      ENDIF 
!=====================================================================
!
      RETURN  
      END SUBROUTINE SET_GEOMETRY 
