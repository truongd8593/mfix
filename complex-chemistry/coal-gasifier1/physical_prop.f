!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PHYSICAL_PROP(DENSITY, SIZE, SP_HEAT,IER)              C
!  Purpose: Calculate physical properties that vary with time          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!  Revision Number: 2                                                 C
!  Purpose: To incorporate heat capacity changes resulting from        C
!              change in composition -- MGAS chemistry                 C
!  Author: S. Venkatesan                              Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Gunn, D.J., "Transfer of heat or mass to particles in fixed and   C
!      fluidised beds," Int. J Heat Mass Transfer, Vol. 21, pp 467-476,C
!      (1978).                                                         C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!    Bird, R.B., W.E. Stewart, E.N. Lightfoot, Transport Phenomena,    C
!      John Wiley & Sons, New York, 1960.                              C
!    Bauer, R., and E.U. Schlunder, "Effective radial thermal          C
!      conductivity of packings in gas flow: Part II: Thermal          C
!      conductivity of the packing fraction without gas flow," Int.    C
!      Chem. Eng., Vol. 18, 189-204 (1978).                            C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE PHYSICAL_PROP(DENSITY, SIZE, SP_HEAT, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      USE constant
      Use usr
      Use compar
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
!                      Indices
      INTEGER          IJK, M, N
!
!                      Average molecular weight
      DOUBLE PRECISION MW
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), SIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'usrnlst.inc'
!
      DOUBLE PRECISION TGX, TSX, VMLEFT, DIFF, EP_g2
      DOUBLE PRECISION Sc1o3, UGC, VGC, WGC, USCM, VSCM, WSCM, VREL, Re
      INTEGER IMJK, IJMK, IJKM, I
 
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'cp_fun2.inc'
!
!  Fluid phase
!
!$omp  parallel do private(IJK, MW, N, TGX)  
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
!----------------------------------------------------------------------	  
           TGX  = T_g(IJK)
!
! 1.1      Density
!
!              Average molecular weight
!              -----------------------------------------------------
                  MW = ZERO 
                  N = 1 
                  IF (NMAX(0) > 0) THEN 
                     MW = MW + SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0))) 
                     N = NMAX(0) + 1 
                  ENDIF 
                  MW = ONE/MAX(MW,OMW_MAX) 
!              -----------------------------------------------------
                  MW_MIX_G(IJK) = MW 
                  RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK)) 
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK) 
!
!           Constant pressure spcific heat of air in cal/g.K
!
            C_pg(IJK) =  X_g(IJK,1)*CPO2(TGX)&
                       + X_g(IJK,2)*CPCO(TGX)  + X_g(IJK,3)*CPCO2(TGX)&
                       + X_g(IJK,4)*CPCH4(TGX) + X_g(IJK,5)*CPH2(TGX)&
                       + X_g(IJK,6)*CPH2O(TGX) + X_g(IJK,7)*CPN2(TGX)&
                       + X_g(IJK,8)*CPTAR(TGX)
!----------------------------------------------------------------------	  
         ENDIF 
      END DO 
      
      
      DO M = 1, MMAX 
!
!
!$omp  parallel do private(IJK,VMLEFT,TSX,EP_g2,DIFF,Sc1o3,&
!$omp                UGC,VGC,WGC,IMJK,IJMK,IJKM,USCM,VSCM,WSCM,VREL,&
!$omp                   Re,I)  
         DO IJK =  IJKSTART3, IJKEND3 
            IF (.NOT.WALL_AT(IJK)) THEN 
!----------------------------------------------------------------------	  
              TSX  = T_s(IJK, M)
!
!             Specific heat of solids (Coal) in cal/g.K
!             Perry & Chilton(1973) -- Table 3-201 on page 3-136
!
              IF(ROP_s(IJK,M) .GT. ZERO)THEN
                C_ps(IJK, M) = X_s(IJK,M,1)*CPFC(TSX) &
	                     + X_s(IJK,M,2)*CPVM(TSX) &
			     + X_s(IJK,M,3)*CPM(TSX)  &
			     + X_s(IJK,M,4)*CPA(TSX) &
			     + (X_s(IJK,M,5)+X_s(IJK,M,6)+X_s(IJK,M,7)&
                                +X_s(IJK,M,8)) * CPA(TSX)
              ELSE
                C_ps(IJK, M) =  CPA(TSX)
              ENDIF 	  
!
!             VMLEFT is from Coal Conversion Systems Technical
!             Data Book (1978), p. 17
!
              IF(PAVM .NE. ZERO) THEN
                IF(TSX .LE. 1223.0)THEN
                  VMLEFT = ((867.2 / (TSX - 273.) )**3.914)/100.
                  VMSTAR(IJK)  = DAFC * VMLEFT
                ELSE
                  VMSTAR(IJK) = ZERO
                ENDIF
              ELSE
                VMSTAR(IJK) = ZERO
              ENDIF
!
!             Calculate Sherwood number for solids phases (Gunn 1978)
!
              EP_g2 = EP_g(IJK) * EP_g(IJK)
              DIFF = 4.26 * ((T_g(IJK)/1800.)**1.75) * 1013000. / P_g(IJK)
              Sc1o3 = (MU_g(IJK)/(RO_g(IJK) * DIFF))**(1./3.)
	      
	      IMJK = IM_OF(IJK)
	      IJMK = JM_OF(IJK)
	      IJKM = KM_OF(IJK)
	      I    = I_OF(IJK)
	      
              UGC = AVG_X_E(U_g(IMJK), U_g(IJK), I)
              VGC = AVG_Y_N(V_g(IJMK), V_g(IJK))
              WGC = AVG_Z_T(W_g(IJKM), W_g(IJK))
!
              USCM = AVG_X_E(U_s(IMJK,M), U_s(IJK,M), I)
              VSCM = AVG_Y_N(V_s(IJMK,M), V_s(IJK,M))
              WSCM = AVG_Z_T(W_s(IJKM,M), W_s(IJK,M))
!
              VREL = SQRT((UGC - USCM)**2 + (VGC-VSCM)**2&
                                        + (WGC-WSCM)**2 )
              IF (FLUIDorP_FLOW_AT(IJK)) THEN
	      IF (MU_g(IJK).gt.ZERO) then
              Re = EP_g(IJK) * D_p(M) * VREL * RO_g(IJK) / MU_g(IJK)
              ELSE
              Re = LARGE_NUMBER
              END IF
              N_sh(IJK, M) = ( (7. - 10. * EP_g(IJK) + 5. * EP_g2)&
                             *(ONE + 0.7 * Re**0.2 * Sc1o3)&
                            + (1.33 - 2.4*EP_g(IJK) + 1.2*EP_g2)&
                             * Re**0.7 * Sc1o3 )
              ELSE
              N_sh(IJK, M) = 0.0
              ENDIF
              
!----------------------------------------------------------------------	  
           ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE PHYSICAL_PROP 
