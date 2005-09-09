!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: FEX                                                    C
!     Purpose: Provide the source terms of ODEs in ODEPACK solver and     C 
!              reaction rates                                             C
!                                                                         C
!     provide ydot and rxn_source_g and rxn_source_s                      C
!     first-order form ydot(dy/dt) = f(t, y)                              C
!     rxn_source_g = source terms of gas species                          C
!     rxn_source_s = source terms of solids species                       C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!      
      SUBROUTINE FEX(NEQ, T, Y, YDOT)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE run
      USE physprop
      USE toleranc
      USe param1
      USE usr
      USE mchem
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!----------------------------------------------- 
!
!                      NEQ(1)Number of equations, NEQ(2) IJK of the cell      
      INTEGER          NEQ(2)
      INTEGER          IJK
      DOUBLE PRECISION T, Y(NEQ(1)), YDOT(NEQ(1))
!
!                      Tempory arrays for Mw_mix_g       
      DOUBLE PRECISION MW_MIX_g_ISAT
!
!                      Pressure
      DOUBLE PRECISION P_g_ISAT
      DOUBLE PRECISION P_H2, P_SiH2, P_SiH4
!
!                      Diffusion
      DOUBLE PRECISION DIFF_SiH2, DIFF_SiH4
      DOUBLE PRECISION R_SiH2, R_SiH4
!
!                      Cp
      DOUBLE PRECISION C_pg_ISAT, C_ps_ISAT
      DOUBLE PRECISION TGX, TSX, XXX
!
!                      Reaction rates
      DOUBLE PRECISION RxnaF, RxnaB, RxnbF, RxnbB, RxncF, RxndF
!
!                      Heat of reaction
      DOUBLE PRECISION HORa, HORb, HORc, HORd
!
!                      Loop indices
      INTEGER          NL

      INCLUDE 'cp_fun1.inc'
      INCLUDE 'cp_fun2.inc'

!
!      The number of the cell
!
       IJK = NEQ(2)
!
!      Temperature
!
       TGX = Y(2)
       TSX = Y(9)
!
!      Calculate the pressure 
!
       MW_MIX_g_ISAT = ZERO
       DO NL = 1, NMAX(0)
          MW_MIX_g_ISAT = MW_MIX_g_ISAT + Y(NL+2)/MW_g(NL)
       END DO
       MW_MIX_g_ISAT = ONE/MAX(MW_MIX_g_ISAT,OMW_MAX)

       P_g_ISAT = Y(1)*(8314.56D4*Y(2))/MW_MIX_g_ISAT
!
!      partial pressure
!
       P_SiH4 = P_g_ISAT * 0.1D-3 * Y(3) * MW_MIX_g_ISAT/MW_g(1)
       P_SiH2 = P_g_ISAT * 0.1D-03 * Y(4) *MW_MIX_g_ISAT/MW_g(2)
       P_H2 = P_g_ISAT * 0.1D-3 * Y(5) * MW_MIX_g_ISAT/MW_g(3)
!
!      Diffusion
!
       R_SiH2 = (8314.7295d0/30.11d0)
       R_SiH4 = (8314.7295d0/32.12d0)
       DIFF_SiH2 = 3.45d-5*Y(2)**1.5
       DIFF_SiH4 = 3.45d-5*Y(2)**1.5
!
!      Heat capacibility
!
       C_pg_ISAT = Y(3)*CPSiH4(TGX) + Y(4)*CPSiH2(TGX)   &
        +Y(5)*CPH2(TGX) + Y(6)*CPSi2H6(TGX)              &
        +Y(7)*CPN2(TGX)
!
       IF( Y(8) .GT. ZERO)THEN
                C_ps_ISAT = Y(10)*CPSi(TSX) + Y(11)*CPAL2O3(TSX)
       ELSE
                C_ps_ISAT = CPAL2O3(TSX) 
       ENDIF

!  User input is required in sections 1 through 4.
!
!111111111111111111111111111111111111111111111111111111111111111111111111111111
!
! 1. Write the rates of various reactions(same as 1 in rrates.f):
!    Write the reaction rates for each of the reactions as RXNxF and RXNxB (both
!    quantities >= 0), where x identifies the reaction, F stands for forward
!    rate, and B stands for the backward rate.  The rates can be in
!    g-mole/(cm^3.s) or g/(cm^3.s).  For the sake of clarity, give the reaction
!    scheme and the units in a comment statement above the rate expression.
!    The volume (cm^3) is that of the computational cell.  Therefore, for
!    example, the rate term of a gas phase reaction will have a multiplicative
!    factor of epsilon. Note that X_g and X_s are mass fractions
!
!  a)   SiH4 <-> SiH2 + H2                  (g-mole/cm^3.s)

        RxnaF = (ONE-Y(8)) * 1.08D13 * EXP(-25600./Y(2)) &	
     &        * (Y(1)*Y(3)/MW_g(1)) 
        RxnaB = (ONE-Y(8)) * 1.26D12 * EXP(-2520./Y(2))  &
     &        * (Y(1)*Y(4)/MW_g(2))	&	
     &	        * (Y(1)*Y(5)/MW_g(3))

!  b)   SiH4 + SiH2 <-> Si2H6                  (g-mole/cm^3.s)

        RxnbF = (ONE-Y(8)) * 1.0D14 * EXP(-2300./Y(2))	 &
     & 	        * (Y(1)*y(3)/MW_g(1))		 &
     &	        * (y(1)*y(4)/MW_g(2))
        RxnbB = (ONE-Y(8)) * 5.62D15 * exp(-26200./Y(2)) &	 
     &	        * (Y(1)*Y(6)/MW_g(4)) 

!  c)   SiH4  -> Si(s) + 2H2                  (g-mole/cm^3.s)      
    
        RxncF  = (6.D0 * Y(8) / Y(12) )	&	
     &	        * 2.15D10 * EXP(-23016.0/Y(2))	&	
     &	        * (Y(1)*Y(3)/MW_g(1))		&
     &		/ (ONE + 0.034d0 * P_H2 + 7.6D-3&	
     &		* EXP(3954.0/y(2)) * P_SiH4)

!  d)   SiH2  -> Si(s) + H2                  (g-mole/cm^3.s) 

        RxndF  = DIFF_SiH2 * N_sh(IJK, 1) * P_SiH2 * 6.0d0 * Y(8) &	
     &               / ( ( Y(12) * MW_g(2) ) * (Y(12) * R_SiH2 * Y(2)))
!
!22222222222222222222222222222222222222222222222222222222222222222222222222222
! 2. Write the reaction rates of various species:
!    Obtain the rates of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.
!    the reactions rates for gas species n are added to get RXN_source_g(n) and
!    solids species n of m phases are added to get RXN_source_s(m, n). The source
!    terms of inert species are set zero
!
!  GAS SPECIES
!
!  (1) SiH4    
      RXN_source_g(1) = (-RxnaF+RxnaB-RxnbF+RxnbB-RxncF)*MW_g(1)
!  (2) SiH2
      RXN_source_g(2) = (RxnaF -RxnaB-RxnbF+RxnbB-RxndF)*MW_g(2)
!  (3) H2
      RXN_source_g(3) = (RxnaF-RxnaB+2.0d0*RxncF+RxndF)*MW_g(3)
!  (4) Si2H6
      RXN_source_g(4) = (RxnbF-RxnbB)*MW_g(4)
!  (5) N2
      RXN_source_g(5) = ZERO
!
!  SOLIDS SPECIES
!
!  (1) Si
      RXN_source_s(1,1) = (RxncF+RxndF)*MW_s(1,1)
!  (2) Al2O3
      RXN_source_s(1,2) = ZERO
!
!33333333333333333333333333333333333333333333333333333333333333333333333333333
! 3.  Determine the heat of reactions in cal/gmol at the
!     temperature T_g or T_s for RXNx (RXNxF and RXNxB). Note that 
!     for exothermic reactions
!     HORx will be negative.
!
      HORa = 57060.0D0
      HORb = -54350.0D0
      HORc = -8210.0D0
      HORd = -65360.0D0
    
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
! 4. Wring the ODEs need to solve in ODEPACK format  
!    (referring to the manual)
!
!     RO_g
      YDOT(1) = ONE/(ONE-Y(8))*(SUM(RXN_source_g(1:5)) + Y(1)/RO_s(1)*SUM(RXN_source_s(1,1:2)) )
!     T_g
      YDOT(2) = -(HORa * (RxnaF-RxnaB) + HORb*(RxnbF-RxnbB) )   &
                /Y(1)/(ONE-Y(8))/C_pg_ISAT
!     SiH4
      YDOT(3) = RXN_source_g(1)/((ONE-Y(8))*Y(1)) - Y(3)/((ONE-Y(8))*Y(1) ) &
                *(SUM(RXN_source_g(1:5)))
!     SiH2
      YDOT(4) = RXN_source_g(2)/((ONE-Y(8))*Y(1)) - Y(4)/((ONE-Y(8))*Y(1) ) &
                *(SUM(RXN_source_g(1:5)))
!     H2
      YDOT(5) = RXN_source_g(3)/((ONE-Y(8))*Y(1)) - Y(5)/((ONE-Y(8))*Y(1) ) &
                *(SUM(RXN_source_g(1:5)))
!     Si2H6
      YDOT(6) = RXN_source_g(4)/((ONE-Y(8))*Y(1)) - Y(6)/((ONE-Y(8))*Y(1) ) &
                *(SUM(RXN_source_g(1:5)))
!     N2
!      YDOT(7) = ZERO
      YDOT(7) = RXN_source_g(5)/((ONE-Y(8))*Y(1)) - Y(7)/((ONE-Y(8))*Y(1) ) &
                *(SUM(RXN_source_g(1:5)))
!     EP_s(1)
      YDOT(8) = SUM(RXN_source_s(1,1:2))/RO_S(1)
!     T_s(1)
      YDOT(9) = - (HORc * RxncF + HORd * RxndF) /RO_S(1)/C_ps_ISAT/Y(8)
!     Si
      YDOT(10) = RXN_source_s(1,1)/(Y(8)*RO_S(1)) -   &
                 Y(10)/(Y(8)*RO_S(1) )*SUM(RXN_source_s(1,1:2))
!     Al2O3
!      YDOT(11) = ZERO
      YDOT(11) = RXN_source_s(1,2)/(Y(8)*RO_S(1)) -   &
                 Y(11)/(Y(8)*RO_S(1) )*SUM(RXN_source_s(1,1:2))
!     D_P(1)      
      IF (CALL_GROW) THEN
         YDOT(12) = ONE/(Y(8)*RO_S(1))*(4.0d0/3.0d0*Y(12)*SUM(RXN_source_s(1,1:2)) &
                    - Y(12) * SUM(RXN_source_s(1,1:2)) )
      ELSE
         YDOT(12) = ZERO
      END IF

      RETURN
      END
