!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_S(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for scalar quantities               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
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
      SUBROUTINE PARTIAL_ELIM_S(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
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
!                      Indices 
      INTEGER          IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a0, b0, a1, b1, F10, Saxf0, Saxf1 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (MMAX == 1) THEN 
!
!$omp  parallel do private( IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!$omp&  a0, b0, a1, b1, F10, Saxf0, Saxf1) &
!$omp&  schedule(static)
         DO IJK = 1, IJKMAX2 
            IF (FLUID_AT(IJK)) THEN 
               IF (VXF(IJK,1) > ZERO) THEN 
!
                  IJKW = WEST_OF(IJK) 
                  IJKS = SOUTH_OF(IJK)
                  IJKE = EAST_OF(IJK) 
                  IJKN = NORTH_OF(IJK) 
!
                  A0 = A_M(IJK,0,0) 
                  B0 = B_M(IJK,0) 
                  A1 = A_M(IJK,0,1) 
                  B1 = B_M(IJK,1) 
                  F10 = -VXF(IJK,1) 
!
                  SAXF0 = -(A_M(IJK,E,0)*VAR_G(IJKE)+A_M(IJK,W,0)*VAR_G(IJKW)+&
                     A_M(IJK,N,0)*VAR_G(IJKN)+A_M(IJK,S,0)*VAR_G(IJKS)) 
                  SAXF1 = -(A_M(IJK,E,1)*VAR_S(IJKE,1)+A_M(IJK,W,1)*VAR_S(IJKW,&
                     1)+A_M(IJK,N,1)*VAR_S(IJKN,1)+A_M(IJK,S,1)*VAR_S(IJKS,1)) 
!
                  IF (DO_K) THEN 
                     IJKB = BOTTOM_OF(IJK) 
                     IJKT = TOP_OF(IJK) 
                     SAXF0 = SAXF0 - (A_M(IJK,T,0)*VAR_G(IJKT)+A_M(IJK,B,0)*&
                        VAR_G(IJKB)) 
                     SAXF1 = SAXF1 - (A_M(IJK,T,1)*VAR_S(IJKT,1)+A_M(IJK,B,1)*&
                        VAR_S(IJKB,1)) 
                  ENDIF 
!
!             gas phase
                  A_M(IJK,0,0) = A0 + (A1*F10)/(A1 + F10) 
                  B_M(IJK,0) = B0 + F10*(SAXF1 + B1)/(A1 + F10) 
!
!             solids phase
                  A_M(IJK,0,1) = A1 + (A0*F10)/(A0 + F10) 
                  B_M(IJK,1) = B1 + F10*(SAXF0 + B0)/(A0 + F10) 
!
               ENDIF 
            ENDIF 
         END DO 
      ELSE 
         IER = 1 
         WRITE (LINE, *) 'Error: Cannot do partial elimination for M > 1' 
         CALL WRITE_ERROR ('PARTIAL_ELIM_S', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_U(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for X vector quantities
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
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
      SUBROUTINE PARTIAL_ELIM_U(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
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
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a0, b0, a1, b1, F10, Saxf0, Saxf1 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (MMAX == 1) THEN 
!
         IF (MOMENTUM_X_EQ(0) .AND. MOMENTUM_X_EQ(1)) THEN 
!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!$omp&  a0, b0, a1, b1, F10, Saxf0, Saxf1) &
!$omp&  schedule(static)
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_E(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     IMJK = IM_OF(IJK) 
                     IJMK = JM_OF(IJK) 
                     IPJK = IP_OF(IJK) 
                     IJPK = JP_OF(IJK) 
!
                     A0 = A_M(IJK,0,0) 
                     B0 = B_M(IJK,0) 
                     A1 = A_M(IJK,0,1) 
                     B1 = B_M(IJK,1) 
                     F10 = -VXF(IJK,1) 
!
                     SAXF0 = -(A_M(IJK,E,0)*VAR_G(IPJK)+A_M(IJK,W,0)*VAR_G(IMJK&
                        )+A_M(IJK,N,0)*VAR_G(IJPK)+A_M(IJK,S,0)*VAR_G(IJMK)) 
                     SAXF1 = -(A_M(IJK,E,1)*VAR_S(IPJK,1)+A_M(IJK,W,1)*VAR_S(&
                        IMJK,1)+A_M(IJK,N,1)*VAR_S(IJPK,1)+A_M(IJK,S,1)*VAR_S(&
                        IJMK,1)) 
!
                     IF (DO_K) THEN 
                        IJKM = KM_OF(IJK) 
                        IJKP = KP_OF(IJK) 
                        SAXF0 = SAXF0 - (A_M(IJK,T,0)*VAR_G(IJKP)+A_M(IJK,B,0)*&
                           VAR_G(IJKM)) 
                        SAXF1 = SAXF1 - (A_M(IJK,T,1)*VAR_S(IJKP,1)+A_M(IJK,B,1&
                           )*VAR_S(IJKM,1)) 
                     ENDIF 
!
!               gas phase
                     A_M(IJK,0,0) = A0 + (A1*F10)/(A1 + F10) 
                     B_M(IJK,0) = B0 + F10*(SAXF1 + B1)/(A1 + F10) 
!
!               solids phase
                     A_M(IJK,0,1) = A1 + (A0*F10)/(A0 + F10) 
                     B_M(IJK,1) = B1 + F10*(SAXF0 + B0)/(A0 + F10) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_X_EQ(0)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_E(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               gas phase
                     A_M(IJK,0,0) = A_M(IJK,0,0) + F10 
                     B_M(IJK,0) = B_M(IJK,0) + F10*VAR_S(IJK,1) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_X_EQ(1)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_E(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               solids phase
                     A_M(IJK,0,1) = A_M(IJK,0,1) + F10 
                     B_M(IJK,1) = B_M(IJK,1) + F10*VAR_G(IJK) 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      ELSE 
         IER = 1 
         WRITE (LINE, *) 'Error: Cannot do partial elimination for M > 1' 
         CALL WRITE_ERROR ('PARTIAL_ELIM_U', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_U 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_V(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for Y vector quantities
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
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
      SUBROUTINE PARTIAL_ELIM_V(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
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
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a0, b0, a1, b1, F10, Saxf0, Saxf1 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (MMAX == 1) THEN 
!
         IF (MOMENTUM_Y_EQ(0) .AND. MOMENTUM_Y_EQ(1)) THEN 
!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!$omp&  a0, b0, a1, b1, F10, Saxf0, Saxf1) &
!$omp&  schedule(static)
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_N(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     IMJK = IM_OF(IJK) 
                     IJMK = JM_OF(IJK) 
                     IPJK = IP_OF(IJK) 
                     IJPK = JP_OF(IJK) 
!
                     A0 = A_M(IJK,0,0) 
                     B0 = B_M(IJK,0) 
                     A1 = A_M(IJK,0,1) 
                     B1 = B_M(IJK,1) 
                     F10 = -VXF(IJK,1) 
!
                     SAXF0 = -(A_M(IJK,E,0)*VAR_G(IPJK)+A_M(IJK,W,0)*VAR_G(IMJK&
                        )+A_M(IJK,N,0)*VAR_G(IJPK)+A_M(IJK,S,0)*VAR_G(IJMK)) 
                     SAXF1 = -(A_M(IJK,E,1)*VAR_S(IPJK,1)+A_M(IJK,W,1)*VAR_S(&
                        IMJK,1)+A_M(IJK,N,1)*VAR_S(IJPK,1)+A_M(IJK,S,1)*VAR_S(&
                        IJMK,1)) 
!
                     IF (DO_K) THEN 
                        IJKM = KM_OF(IJK) 
                        IJKP = KP_OF(IJK) 
                        SAXF0 = SAXF0 - (A_M(IJK,T,0)*VAR_G(IJKP)+A_M(IJK,B,0)*&
                           VAR_G(IJKM)) 
                        SAXF1 = SAXF1 - (A_M(IJK,T,1)*VAR_S(IJKP,1)+A_M(IJK,B,1&
                           )*VAR_S(IJKM,1)) 
                     ENDIF 
!
!               gas phase
                     A_M(IJK,0,0) = A0 + (A1*F10)/(A1 + F10) 
                     B_M(IJK,0) = B0 + F10*(SAXF1 + B1)/(A1 + F10) 
!
!               solids phase
                     A_M(IJK,0,1) = A1 + (A0*F10)/(A0 + F10) 
                     B_M(IJK,1) = B1 + F10*(SAXF0 + B0)/(A0 + F10) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_Y_EQ(0)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_N(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               gas phase
                     A_M(IJK,0,0) = A_M(IJK,0,0) + F10 
                     B_M(IJK,0) = B_M(IJK,0) + F10*VAR_S(IJK,1) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_Y_EQ(1)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_N(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               solids phase
                     A_M(IJK,0,1) = A_M(IJK,0,1) + F10 
                     B_M(IJK,1) = B_M(IJK,1) + F10*VAR_G(IJK) 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      ELSE 
         IER = 1 
         WRITE (LINE, *) 'Error: Cannot do partial elimination for M > 1' 
         CALL WRITE_ERROR ('PARTIAL_ELIM_V', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_V 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_W(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for Z vector quantities
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
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
      SUBROUTINE PARTIAL_ELIM_W(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
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
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a0, b0, a1, b1, F10, Saxf0, Saxf1 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
      IF (MMAX == 1) THEN 
!
         IF (MOMENTUM_Z_EQ(0) .AND. MOMENTUM_Z_EQ(1)) THEN 
!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, &
!$omp&  a0, b0, a1, b1, F10, Saxf0, Saxf1) &
!$omp&  schedule(static)
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_T(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     IMJK = IM_OF(IJK) 
                     IJMK = JM_OF(IJK) 
                     IPJK = IP_OF(IJK) 
                     IJPK = JP_OF(IJK) 
!
                     A0 = A_M(IJK,0,0) 
                     B0 = B_M(IJK,0) 
                     A1 = A_M(IJK,0,1) 
                     B1 = B_M(IJK,1) 
                     F10 = -VXF(IJK,1) 
!
                     SAXF0 = -(A_M(IJK,E,0)*VAR_G(IPJK)+A_M(IJK,W,0)*VAR_G(IMJK&
                        )+A_M(IJK,N,0)*VAR_G(IJPK)+A_M(IJK,S,0)*VAR_G(IJMK)) 
                     SAXF1 = -(A_M(IJK,E,1)*VAR_S(IPJK,1)+A_M(IJK,W,1)*VAR_S(&
                        IMJK,1)+A_M(IJK,N,1)*VAR_S(IJPK,1)+A_M(IJK,S,1)*VAR_S(&
                        IJMK,1)) 
!
                     IF (DO_K) THEN 
                        IJKM = KM_OF(IJK) 
                        IJKP = KP_OF(IJK) 
                        SAXF0 = SAXF0 - (A_M(IJK,T,0)*VAR_G(IJKP)+A_M(IJK,B,0)*&
                           VAR_G(IJKM)) 
                        SAXF1 = SAXF1 - (A_M(IJK,T,1)*VAR_S(IJKP,1)+A_M(IJK,B,1&
                           )*VAR_S(IJKM,1)) 
                     ENDIF 
!
!               gas phase
                     A_M(IJK,0,0) = A0 + (A1*F10)/(A1 + F10) 
                     B_M(IJK,0) = B0 + F10*(SAXF1 + B1)/(A1 + F10) 
!
!               solids phase
                     A_M(IJK,0,1) = A1 + (A0*F10)/(A0 + F10) 
                     B_M(IJK,1) = B1 + F10*(SAXF0 + B0)/(A0 + F10) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_Z_EQ(0)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_T(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               gas phase
                     A_M(IJK,0,0) = A_M(IJK,0,0) + F10 
                     B_M(IJK,0) = B_M(IJK,0) + F10*VAR_S(IJK,1) 
                  ENDIF 
               ENDIF 
            END DO 
         ELSE IF (MOMENTUM_Z_EQ(1)) THEN 
            DO IJK = 1, IJKMAX2 
               IF (FLOW_AT_T(IJK)) THEN 
                  IF (VXF(IJK,1) > ZERO) THEN 
!
                     F10 = -VXF(IJK,1) 
!
!               solids phase
                     A_M(IJK,0,1) = A_M(IJK,0,1) + F10 
                     B_M(IJK,1) = B_M(IJK,1) + F10*VAR_G(IJK) 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      ELSE 
         IER = 1 
         WRITE (LINE, *) 'Error: Cannot do partial elimination for M > 1' 
         CALL WRITE_ERROR ('PARTIAL_ELIM_W', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_W 
