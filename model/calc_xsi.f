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
      SUBROUTINE CALC_XSI(DISCR, PHI, U, V, W, XSI_E, XSI_N, XSI_T,incr) 
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
      USE vshear
      USE compar    !//d
!//SP
      USE sendrecv
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
      CHARACTER*80     LINE(1)
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

      INTEGER incr
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: PHI_C_OF, XSI, MINMOD, VANLEER, &
         ULTRA_QUICK, QUICKEST, SUPERBEE, SMART, MUSCL 
!-----------------------------------------------
      INCLUDE 'function.inc'


	IF (SHEAR) THEN			
! calculate XSI_E,XSI_N,XSI_T when periodic shear BCs are used

	call CXS(incr,DISCR,U,V,W,PHI,XSI_E,XSI_N,XSI_T)


	ELSE
!
!
      SELECT CASE (DISCR)                        !first order upwinding 
      CASE (:1)  
!
!$omp    parallel do private(IJK)
!//SP
         DO IJK = IJKSTART3, IJKEND3 
!
            XSI_E(IJK) = XSI(U(IJK),ZERO) 
            XSI_N(IJK) = XSI(V(IJK),ZERO) 
            IF (DO_K) XSI_T(IJK) = XSI(W(IJK),ZERO) 
         END DO 
      CASE (2)                                   !Superbee 

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU, PHI_C,DWF)
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
!//SP
         DO IJK = IJKSTART3, IJKEND3
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
      
      ENDIF
!//SP
      call send_recv(XSI_E,2)
      call send_recv(XSI_N,2)
      call send_recv(XSI_T,2)
      RETURN  
      END SUBROUTINE CALC_XSI 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE CXS calculates XSI_E, XSI_N,XSI_T when using periodic shear
!BC's.  Uses true velocities.


	SUBROUTINE CXS(INCR,DISCR,U,V,W,PHI,XSI_E,XSI_N,XSI_T)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE geometry
      USE indices
      USE vshear
      USE compar    !//d
      USE sendrecv
      IMPLICIT NONE

	INTEGER incr,IJK,IJKC,IJKD,IJKU,I,DISCR
	DOUBLE PRECISION V(DIMENSION_3), U(DIMENSION_3), W(DIMENSION_3)	
	DOUBLE PRECISION PHI(DIMENSION_3),DWFE,DWFN,DWFT
	DOUBLE PRECISION PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV
	DOUBLE PRECISION PHICW,PHIDW,PHIUW
	DOUBLE PRECISION XSI_E(DIMENSION_3), XSI_N(DIMENSION_3)
	DOUBLE PRECISION XSI_T(DIMENSION_3),SRT

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: XSI
         
!-----------------------------------------------
        INCLUDE 'function.inc'


	IF (INCR .eq. 2) THEN			!V momentum balance
	SRT=(2d0*V_sh/XLENGTH)

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!$omp&	PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!$omp&  shared(DISCR)
!//SP
        DO IJK = IJKSTART3, IJKEND3

	I=I_OF(IJK)
	V(IJK)=V(IJK)+VSH(IJK)

            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
		PHICU=PHI(IJKC)
		PHIDU=PHI(IJKD)+SRT*1d0/oDX_E(I)
		PHIUU=PHI(IJKU)-SRT*1d0/oDX_E(IM1(I))

	

            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
		PHICU=PHI(IJKC)+SRT*1d0/oDX_E(I)
		PHIDU=PHI(IJKD)
		PHIUU=PHI(IJKU)+SRT*1d0/oDX_E(I)&
		+SRT*1d0/oDX_E(IP1(I))
            ENDIF 

            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU= SOUTH_OF(IJKC) 
 		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
           ELSE 
               IJKC= NORTH_OF(IJK) 
               IJKD= IJK 
               IJKU= NORTH_OF(IJKC) 
		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
            ENDIF 
	
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
		PHICW=PHI(IJKC)
		PHIDW=PHI(IJKD)
		PHIUW=PHI(IJKU)
		
	    ELSE 
		PHICW=0d0
		PHIDW=0d0
		PHIUW=0d0
		W(IJK)=0d0
		DWFT=0d0
	    END IF

	CALL DW(U(IJK),V(IJK),W(IJK),IJK,PHICU,PHIDU,PHIUU,PHICV,&
	PHIDV,PHIUV,PHICW,PHIDW,PHIUW,DWFE,DWFN,DWFT,DISCR)
            XSI_E(IJK) = XSI(U(IJK),DWFE) 
            XSI_N(IJK) = XSI(V(IJK),DWFN) 
            XSI_T(IJK) = XSI(W(IJK),DWFT) 

	V(IJK)=V(IJK)-VSH(IJK)
	END DO


	ELSE IF (INCR .eq. 1) THEN  			!u momentum balance

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!$omp&	PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!$omp&  shared(DISCR)
!//SP   
        DO IJK = IJKSTART3, IJKEND3

	V(IJK)=V(IJK)+VSHE(IJK)
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
		PHICU=PHI(IJKC)
		PHIDU=PHI(IJKD)
		PHIUU=PHI(IJKU)

	

            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
		PHICU=PHI(IJKC)
		PHIDU=PHI(IJKD)
		PHIUU=PHI(IJKU)
            ENDIF 

            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU= SOUTH_OF(IJKC) 
 		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
           ELSE 
               IJKC= NORTH_OF(IJK) 
               IJKD= IJK 
               IJKU= NORTH_OF(IJKC) 
		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
            ENDIF 

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
		PHICW=PHI(IJKC)
		PHIDW=PHI(IJKD)
		PHIUW=PHI(IJKU)
	     ELSE		
		PHICW=0d0
		PHIDW=0d0
		PHIUW=0d0
		W(IJK)=0d0
		DWFT=0d0
	     END IF
	
	CALL DW(U(IJK),V(IJK),W(IJK),IJK,PHICU,PHIDU,PHIUU,PHICV,&
	PHIDV,PHIUV,PHICW,PHIDW,PHIUW,DWFE,DWFN,DWFT,DISCR)

            XSI_E(IJK) = XSI(U(IJK),DWFE) 
            XSI_N(IJK) = XSI(V(IJK),DWFN) 
            XSI_T(IJK) = XSI(W(IJK),DWFT) 

	V(IJK)=V(IJK)-VSHE(IJK)
	END DO


	ELSE IF (INCR .eq. 0) THEN			!scalars and w momentum

!$omp    parallel do private(IJK, IJKC,IJKD,IJKU,DWFE,DWFN,DWFT,&
!$omp&	PHICU,PHIDU,PHIUU,PHICV,PHIDV,PHIUV,PHICW,PHIDW,PHIUW,I)&
!$omp&  shared(DISCR)
!//SP   
        DO IJK = IJKSTART3, IJKEND3

	V(IJK)=V(IJK)+VSH(IJK)
            IF (U(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = EAST_OF(IJK) 
               IJKU = WEST_OF(IJKC) 
		PHICU=PHI(IJKC)
		PHIDU=PHI(IJKD)
		PHIUU=PHI(IJKU)

            ELSE 
               IJKC = EAST_OF(IJK) 
               IJKD = IJK 
               IJKU = EAST_OF(IJKC) 
		PHICU=PHI(IJKC)
		PHIDU=PHI(IJKD)
		PHIUU=PHI(IJKU)
            ENDIF 

            IF (V(IJK) >= ZERO) THEN 
               IJKC = IJK 
               IJKD = NORTH_OF(IJK) 
               IJKU= SOUTH_OF(IJKC) 
 		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
           ELSE 
               IJKC= NORTH_OF(IJK) 
               IJKD= IJK 
               IJKU= NORTH_OF(IJKC) 
		PHICV=PHI(IJKC)
		PHIDV=PHI(IJKD)
		PHIUV=PHI(IJKU)
            ENDIF 

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
		PHICW=PHI(IJKC)
		PHIDW=PHI(IJKD)
		PHIUW=PHI(IJKU)
	   ELSE
		PHICW=0d0
		PHIDW=0d0
		PHIUW=0d0
	 	W(IJK)=0d0
		DWFT=0d0
	   END IF				


	
	CALL DW(U(IJK),V(IJK),W(IJK),IJK,PHICU,PHIDU,PHIUU,PHICV,&
	PHIDV,PHIUV,PHICW,PHIDW,PHIUW,DWFE,DWFN,DWFT,DISCR)


           XSI_E(IJK) = XSI(U(IJK),DWFE) 
           XSI_N(IJK) = XSI(V(IJK),DWFN) 
           XSI_T(IJK) = XSI(W(IJK),DWFT) 

	V(IJK)=V(IJK)-VSH(IJK)
	END DO
	
	ELSE
	write(*,*) 'INCR ERROR'
	END IF

!//SP
        call send_recv(XSI_E,2)
        call send_recv(XSI_N,2)
        call send_recv(XSI_T,2)
!        not needed, because VSH is first added and then subtracted from V
!        so that there is no net change in V
!        call send_recv(V,2)  

	RETURN
	END SUBROUTINE CXS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE DW calculates DWFs for various discretization methods
!when period shear BCs are used

	SUBROUTINE DW(UU,VV,WW,IJK,PHICU,PHIDU,PHIUU,PHICV,&
	PHIDV,PHIUV,PHICW,PHIDW,PHIUW,DWFE,DWFN,DWFT,DISCR)

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE geometry
      USE indices
      USE vshear
      
      USE compar    !//d
      USE sendrecv
      
      IMPLICIT NONE

!                      discretization method
      INTEGER          DISCR
!

!                      Indices
      INTEGER          IJK, IJKC, IJKD, IJKU, I, J, K,incr

!
!                      
      DOUBLE PRECISION  PHI_C


      DOUBLE PRECISION PHICU, PHIDU, PHIUU, DWFE, UU
      DOUBLE PRECISION PHICV, PHIDV, PHIUV, DWFN, VV
      DOUBLE PRECISION PHICW, PHIDW, PHIUW, DWFT, WW
!
!                      cell widths for QUICKEST
      DOUBLE PRECISION oDXc, oDXuc, oDYc, oDYuc, oDZc, oDZuc
      DOUBLE PRECISION CF,DZUC
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: PHI_C_OF, XSI, MINMOD, VANLEER, &
         ULTRA_QUICK, QUICKEST, SUPERBEE, SMART, MUSCL 
!-----------------------------------------------
      INCLUDE 'function.inc'

      SELECT CASE (DISCR)                        !first order upwinding 
      CASE (:1)  

	    DWFE=0d0
	    DWFN=0d0
	    DWFT=0d0
      CASE (2)                                   !Superbee 
            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            DWFE = SUPERBEE(PHI_C) 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            DWFN = SUPERBEE(PHI_C) 	 	

	IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            DWFT = SUPERBEE(PHI_C) 
	END IF

      CASE (3)                                   !SMART 

            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            DWFE = SMART(PHI_C) 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            DWFN = SMART(PHI_C) 	 	

	IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            DWFT = SMART(PHI_C) 
	END IF
      CASE (4)                                   !ULTRA-QUICK 


	    I=I_OF(IJK)
 	    J=J_OF(IJK)	


            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            CF = ABS(UU)*DT*ODX_E(I) 
            DWFE = ULTRA_QUICK(PHI_C,CF)  

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            CF = ABS(VV)*DT*ODY_N(J) 
            DWFN = ULTRA_QUICK(PHI_C,CF) 
           	 	

	IF (DO_K) THEN
	    K=K_OF(IJK)
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            CF = ABS(WW)*DT*OX(I)*ODZ_T(K) 
            DWFT = ULTRA_QUICK(PHI_C,CF) 
 	END IF    
      CASE(5)					!QUICKEST

	    I=I_OF(IJK)
 	    J=J_OF(IJK)	

            IF (UU >= ZERO) THEN 
               ODXC = ODX(I) 
               ODXUC = ODX_E(IM1(I))
            ELSE 
               ODXC = ODX(IP1(I)) 
               ODXUC = ODX_E(IP1(I)) 
            ENDIF 

            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 	    
            CF = ABS(UU)*DT*ODX_E(I) 
            DWFE = QUICKEST(PHI_C,CF,ODXC,ODXUC,ODX_E(I)) 

            IF (VV >= ZERO) THEN 

               ODYC = ODY(J) 
               ODYUC = ODY_N(JM1(J)) 
            ELSE  
               ODYC = ODY(JP1(J)) 
               ODYUC = ODY_N(JP1(J)) 
            ENDIF 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV)  
            CF = ABS(VV)*DT*ODY_N(J) 
            DWFN = QUICKEST(PHI_C,CF,ODYC,ODYUC,ODY_N(J)) 

	IF (DO_K) THEN
	    K=K_OF(IJK)
            IF (WW >= ZERO) THEN 
               ODZC = ODZ(K) 
               ODZUC = ODZ_T(KM1(K)) 
            ELSE 
               ODZC = ODZ(KP1(K)) 
               DZUC = ODZ_T(KP1(K)) 
            ENDIF 

            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW)
            CF = ABS(WW)*DT*OX(I)*ODZ_T(K) 
            DWFT = QUICKEST(PHI_C,CF,ODZC,ODZUC,ODZ_T(K)) 
	END IF

      CASE (6)                                   !MUSCL 
            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            DWFE = MUSCL(PHI_C) 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            DWFN = MUSCL(PHI_C) 	 	

	IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            DWFT = MUSCL(PHI_C)
	END IF	

      CASE (7)                                   !Van Leer 

            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            DWFE = VANLEER(PHI_C) 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            DWFN = VANLEER(PHI_C) 	 	

	IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            DWFT = VANLEER(PHI_C)
	END IF

      CASE (8)                                   !Minmod 
            PHI_C = PHI_C_OF(PHIUU,PHICU,PHIDU) 
            DWFE = MINMOD(PHI_C) 

            PHI_C = PHI_C_OF(PHIUV,PHICV,PHIDV) 
            DWFN = MINMOD(PHI_C) 	 	

	IF (DO_K) THEN
            PHI_C = PHI_C_OF(PHIUW,PHICW,PHIDW) 
            DWFT = MINMOD(PHI_C)
	END IF
      CASE DEFAULT                               !Error 
         WRITE (*,*) 'DISCRETIZE = ', DISCR, ' not supported.' 

         STOP  

      END SELECT 
	RETURN
	END SUBROUTINE DW
