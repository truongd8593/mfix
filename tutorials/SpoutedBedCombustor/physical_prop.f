!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PHYSICAL_PROP(DENSITY, PSIZE, SP_HEAT,IER)              C
!  Purpose: Calculate physical properties that vary with time          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
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
!                                                                      C
!  Variables referenced: ANY_SPECIES_EQ                                C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: TGX, TSX                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
      SUBROUTINE PHYSICAL_PROP(DENSITY, PSIZE, SP_HEAT, IER) 
 
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
      USE compar 
      USE funits    
      IMPLICIT NONE
      INCLUDE 'usrnlst.inc'
!
!  Local Variables
!
!
!
!                      Error index
      INTEGER          IER
      INTEGER		IJK, M
      DOUBLE PRECISION MW
      DOUBLE PRECISION TGX, TSX
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!
!  Statement functions
!
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'cp_fun2.inc'
!
!
!
!$omp  parallel do private(IJK, MW, N)  
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
!
! 1.1      Density
!
            IF (DENSITY(0)) THEN 
               IF (MW_AVG == UNDEFINED) THEN 
!              Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
!              -----------------------------------------------------
                  MW = ONE/SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
!              -----------------------------------------------------
                  MW_MIX_G(IJK) = MW
                  RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ELSE
                  RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF
            ENDIF

           TGX  = T_g(IJK)
!
!     Constant pressure spcific heat of air in cal/g.K
!     Perry & Chilton (1973) -- Table 3-174
!
           IF (SP_HEAT(0) .AND. C_PG0==UNDEFINED) THEN
            C_pg(IJK) =  X_g(IJK,1)*CPO2(TGX)  + X_g(IJK,2)*CPCO(TGX) & 
                       + X_g(IJK,3)*CPCO2(TGX) + X_g(IJK,4)*CPN2(TGX)
           ENDIF
           IF (UNITS == 'SI') C_PG(IJK) = 4183.925*C_PG(IJK)    !in J/kg K
         ENDIF 
      END DO 
      
!
      DO M = 1, MMAX 
!
!$omp  parallel do private(IJK)  
         DO IJK = IJKSTART3, IJKEND3 
            IF (.NOT.WALL_AT(IJK)) THEN 
      
              TSX  = T_s(IJK, M)
!
!     Specific heat of solids (Coal) in cal/g.K
!     Perry & Chilton(1973) -- Table 3-201 on page 3-136
!
              IF (SP_HEAT(M) .AND. C_PS0==UNDEFINED) THEN
                C_ps(IJK, M) = X_s(IJK,M,1) * CPFC(TSX) &
                         + X_s(IJK,M,2) * CPA(TSX)
              ENDIF
              IF (UNITS == 'SI') C_PS(IJK,M) = 4183.925*C_PS(IJK,M) !J/kg K
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE PHYSICAL_PROP
