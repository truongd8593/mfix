!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DRAG_ss(L, M, F_ss, IER)                               C
!  Purpose: This module computes the coefficient of drag between       C
!           solids phase m and solids phase l                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: D_p, MMAX, IJK, EP_g, C_f, e, ROP_s, M, RO_s  C
!  Variables modified: DRAG_ss                                         C
!                                                                      C
!  Local variables: DPSUM                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE DRAG_SS(L, M, F_SS, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
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
!                      Relative velocity between solids phase m and l 
      DOUBLE PRECISION VREL 
! 
!                      Indices 
      INTEGER          I, IJK, IMJK, IJMK, IJKM 
! 
!                      Index of solids phases 
      INTEGER          L, M 
! 
!                      Index for storing solids-solids drag coefficients 
!                      in the upper triangle of the matrix 
      INTEGER          LM 
! 
!                      Cell center value of U_sm 
      DOUBLE PRECISION USCM 
! 
!                      Cell center value of U_sl 
      DOUBLE PRECISION USCL 
! 
!                      Cell center value of V_sm 
      DOUBLE PRECISION VSCM 
! 
!                      Cell center value of V_sl 
      DOUBLE PRECISION VSCL 
! 
!                      Cell center value of W_sm 
      DOUBLE PRECISION WSCM 
! 
!                      Cell center value of W_sl 
      DOUBLE PRECISION WSCL 
! 
!                      Particle diameters 
      DOUBLE PRECISION D_pm, D_pl 
! 
!                      Sum of particle diameters 
      DOUBLE PRECISION DPSUM 
! 
!                      Drag coefficient 
      DOUBLE PRECISION CONST 
! 
!                      Drag array 
      DOUBLE PRECISION F_ss(DIMENSION_3, DIMENSION_LM) 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
      D_PM = D_P(M) 
      D_PL = D_P(L) 
      DPSUM = D_PL + D_PM 
      LM = FUNLM(L,M) 
      CONST = 3.*(ONE + C_E)*(PI/2. + C_F*PI*PI/8.)*DPSUM**2/(2.*PI*(RO_S(L)*&
         D_PL**3+RO_S(M)*D_PM**3)) 
!
!!$omp  parallel do private( I, L, M,  IJK, IMJK, IJMK, IJKM, &
!!$omp&  USCM, VSCM, WSCM, &
!!$omp&  VREL, USCL, VSCL, WSCL) &
!!$omp&  schedule(static)
      DO IJK = 1, IJKMAX2 
         IF (.NOT.WALL_AT(IJK)) THEN 
            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!
!         Calculate velocity components at i, j, k
            USCL = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I) 
            VSCL = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L)) 
            WSCL = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L)) 
            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 
!
!         magnitude of solids-solids relative velocity
!
            VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + (WSCL - WSCM)**2) 
!
!
            F_SS(IJK,LM) = CONST*ROP_S(IJK,L)*ROP_S(IJK,M)*G_0(IJK,L,M)*VREL 
!
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE DRAG_SS 
