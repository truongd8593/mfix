!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_XSI(DISCR, PHI, U, V, W, XSI_e, XSI_n, XSI_t)     C
!  Purpose: Determine convection weighting factors for higher order    C
!           discretization.                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-MAR-97   C
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
      SUBROUTINE CALC_XSI(DISCR, PHI, U, V, W, XSI_E, XSI_N, XSI_T) 
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
      USE indices
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
!                      discretization method
      INTEGER          DISCR
!
!                      convected quantity
      DOUBLE PRECISION PHI(DIMENSION_3)
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
!
!                      Convection weighting factors
      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
                       XSI_t(DIMENSION_3)
!
!                      Indices
      INTEGER          IJK, IJKC, IJKD, IJKU, I, J, K

!
!                      Error message
      CHARACTER*80     LINE
!
!                      
      DOUBLE PRECISION DEN, DEN1, PHI_C
!
!                      down wind factor
      DOUBLE PRECISION dwf
!
!                      Courant number
      DOUBLE PRECISION cf
!
!                      cell widths for QUICKEST
      DOUBLE PRECISION oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: PHI_C_OF, XSI, MINMOD, VANLEER, &
         ULTRA_QUICK, QUICKEST, SUPERBEE, SMART, MUSCL 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      SELECT CASE (DISCR)                        !first order upwinding 
      CASE (:1)  
!
!$omp    parallel do private(IJK)
         DO IJK = 1, IJKMAX2 
!
            XSI_E(IJK) = XSI(U(IJK),ZERO) 
            XSI_N(IJK) = XSI(V(IJK),ZERO) 
            IF (DO_K) XSI_T(IJK) = XSI(W(IJK),ZERO) 
         END DO 
      CASE (2)                                   !Superbee 

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
         DO IJK = 1, IJKMAX2 
!
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = SUPERBEE(PHI_C) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = SUPERBEE(PHI_C) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               DWF = SUPERBEE(PHI_C) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (3)                                   !SMART 
!
!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
         DO IJK = 1, IJKMAX2 
!
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = SMART(PHI_C) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = SMART(PHI_C) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               DWF = SMART(PHI_C) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (4)                                   !ULTRA-QUICK 

!$omp    parallel do private(IJK, I,J,K, IJKC,IJKD,IJKU, PHI_C,DWF,CF)
         DO IJK = 1, IJKMAX2 
!
            I = I_OF(IJK) 
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            CF = ABS(U(IJK))*DT*ODX_E(I) 
!
            DWF = ULTRA_QUICK(PHI_C,CF) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            J = J_OF(IJK) 
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            CF = ABS(V(IJK))*DT*ODY_N(J) 
!
            DWF = ULTRA_QUICK(PHI_C,CF) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               K = K_OF(IJK) 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K) 
!
               DWF = ULTRA_QUICK(PHI_C,CF) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (5)                                   !QUICKEST 

!$omp    parallel do &
!$omp&   private(IJK,I,J,K, IJKC,IJKD,IJKU, &
!$omp&           ODXC,ODXUC, PHI_C,CF,DWF, &
!$omp&           ODYC,ODYUC,  ODZC,ODZUC )
         DO IJK = 1, IJKMAX2 
!
            I = I_OF(IJK) 
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
               ODXC = ODX(I) 
               ODXUC = ODX_E(IM1(I)) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
               ODXC = ODX(IP1(I)) 
               ODXUC = ODX_E(IP1(I)) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            CF = ABS(U(IJK))*DT*ODX_E(I) 
            DWF = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX_E(I)) 
!
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            J = J_OF(IJK) 
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
               ODYC = ODY(J) 
               ODYUC = ODY_N(JM1(J)) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
               ODYC = ODY(JP1(J)) 
               ODYUC = ODY_N(JP1(J)) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            CF = ABS(V(IJK))*DT*ODY_N(J) 
!
            DWF = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY_N(J)) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               K = K_OF(IJK) 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
                  ODZC = ODZ(K) 
                  ODZUC = ODZ_T(KM1(K)) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
                  ODZC = ODZ(KP1(K)) 
                  ODZUC = ODZ_T(KP1(K)) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               CF = ABS(W(IJK))*DT*OX(I)*ODZ_T(K) 
!
               DWF = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ_T(K)) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (6)                                   !MUSCL 

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF )
         DO IJK = 1, IJKMAX2 
!
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = MUSCL(PHI_C) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = MUSCL(PHI_C) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               DWF = MUSCL(PHI_C) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (7)                                   !Van Leer 


!$omp    parallel do private( IJK, IJKC,IJKD,IJKU,  PHI_C,DWF )
         DO IJK = 1, IJKMAX2 
!
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = VANLEER(PHI_C) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = VANLEER(PHI_C) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               DWF = VANLEER(PHI_C) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
            ENDIF 
         END DO 
      CASE (8)                                   !Minmod 

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF )
         DO IJK = 1, IJKMAX2 
!
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = MINMOD(PHI_C) 
!
            XSI_E(IJK) = XSI(U(IJK),DWF) 
!
            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU = SOUTH_OF(IJKC) 
            ELSE 
               IJKC = NORTH_OF(IJK) 
               IJKD = IJK 
               IJKU = NORTH_OF(IJKC) 
            ENDIF 
!
            PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
            DWF = MINMOD(PHI_C) 
!
            XSI_N(IJK) = XSI(V(IJK),DWF) 
!
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  IJKC = IJK 
                  IJKD = TOP_OF(IJK) 
                  IJKU = BOTTOM_OF(IJKC) 
               ELSE 
                  IJKC = TOP_OF(IJK) 
                  IJKD = IJK 
                  IJKU = TOP_OF(IJKC) 
               ENDIF 
!
               PHI_C = PHI_C_OF(PHI(IJKU),PHI(IJKC),PHI(IJKD)) 
!
               DWF = MINMOD(PHI_C) 
!
               XSI_T(IJK) = XSI(W(IJK),DWF) 
!
            ENDIF 
         END DO 
      CASE DEFAULT                               !Error 
         WRITE (LINE, '(A,I2,A)') 'DISCRETIZE = ', DISCR, ' not supported.' 
         CALL WRITE_ERROR ('CALC_XSI', LINE, 1) 
         STOP  
      END SELECT 
      RETURN  
      END SUBROUTINE CALC_XSI 
