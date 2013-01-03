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
!
      USE mflux
      USE fldvar
      USE geometry
      USE sendrecv  
      USE compar 
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
!                      Time for integration
      DOUBLE PRECISION chem_Dt
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
      DOUBLE PRECISION TGX, TSX
!
!                      Reaction rates
      DOUBLE PRECISION RxnaF, RxnaB, RxnbF, RxnbB, RxncF, RxncB
      DOUBLE PRECISION gas_const
      DOUBLE PRECISION HORbg,HORbc,HORbt,HORtar
      DOUBLE PRECISION mf_bg_co, mf_bg_co2, mf_bg_ch4, mf_bg_h2, mf_bg_h2o
      DOUBLE PRECISION A_bg,A_bc,A_bt
      DOUBLE PRECISION E_bg,E_bc,E_bt
      DOUBLE PRECISION ZERO_gs_mfrac
      
!                      Loop indices
      INTEGER          NL
!----------------------------------------------- 
!
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'cp_fun2.inc'
!
!      The number of the cell
!
       IJK = NEQ(2)
!
!QX 
!
!      Temperature
       TGX = Y(2)
       TSX = Y(19)
       ZERO_gs_mfrac = 1.d-32

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
       C_pg_ISAT = Y(3)*CPO2(TGX) + Y(4)*CPN2(TGX) + Y(5)*CPH2O(TGX) +&
            Y(6)*CPCO(TGX) + Y(7)*CPCO2(TGX) + Y(8)*CPH2(TGX) + &
            Y(9)*CPCH4(TGX) + Y(10)*CPTAR(TGX) +Y(11)*CPTAR(TGX)

       IF( Y(18) .GT. 1.D-15 .and. &
            (Y(21)+Y(22)+Y(23)+Y(24)+Y(25)) .gt. 0.D0 )THEN
          C_ps_ISAT = Y(21)*CPBIO(TSX) + Y(22)*CPH2O(TSX) + &
               (Y(23)+Y(25))*CPBIOCHAR(TSX) + Y(24)*CPGAS(TSX)
       ELSE
                C_ps_ISAT = CPBIO(TSX)
       ENDIF
!end
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
!QX
! biomass devolatilization
       A_bg = 1.3D+8    
       E_bg = 1.403D+5/4.1868D0 
       mf_bg_co  = 0.270D0
       mf_bg_co2 = 0.386D0
       mf_bg_ch4 = 0.056D0
       mf_bg_h2  = 0.032D0
       mf_bg_h2o = 0.256D0
       A_bc = 1.08D+7        
       E_bc = 1.213D+5/4.1868D0
       A_bt = 2.D+8          
       E_bt = 1.331D+5/4.1868D0
       gas_const = 1.986d0
!  bio ---> CO + CO2 + CH4 + H2 + H2O + char + tar
       if(Y(21) .gt. ZERO_gs_mfrac) then
          RxnaF = Y(17)*Y(18)*Y(21)* A_bg*EXP(-E_bg/(gas_const*Y(19)))
       else
          RxnaF = ZERO
       endif

       if(Y(21) .gt. ZERO_gs_mfrac) then
          RxnbF = Y(17)*Y(18)*Y(21)* A_bc*EXP(-E_bc/(gas_const*Y(19)))
       else
          RxnbF = ZERO
       endif

       if(Y(21) .gt. ZERO_gs_mfrac) then
          RxncF = Y(17)*Y(18)*Y(21)* A_bt*EXP(-E_bt/(gas_const*Y(19)))
       else
          RxncF = ZERO
       endif
!end

!22222222222222222222222222222222222222222222222222222222222222222222222222222
! 2. Write the reaction rates of various species:
!    Obtain the rates of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.
!    the reactions rates for gas species n are added to get RXN_source_g(n) and
!    solids species n of m phases are added to get RXN_source_s(m, n). The source
!    terms of inert species are set zero
!
       RXN_source_g(:) = ZERO
       RXN_source_s(:,:) = ZERO

!  GAS SPECIES, O2, N2, H2O, CO, CO2, CH4, tar
       RXN_source_g(1) = ZERO

       RXN_source_g(2) = ( (-RxnaF-RxnbF-RxncF) * RO_ss(2,4)/RO_ss(2,1)  &
            + RxnbF*RO_ss(2,4)/RO_ss(2,3) )

       RXN_source_g(3) = mf_bg_h2o * RxnaF

       RXN_source_g(4) = mf_bg_co*RxnaF

       RXN_source_g(5) = mf_bg_co2 * RxnaF

       RXN_source_g(6) = mf_bg_h2*RxnaF

       RXN_source_g(7) = mf_bg_ch4*RxnaF 

       RXN_source_g(8) = RxncF

!SOLIDS SPECIES, bio, char, pores, ash 
       RXN_source_s(1,1) = ZERO

       RXN_source_s(2,1) = -RxnaF - RxnbF - RxncF

       RXN_source_s(2,3) = 0.9d0*RxnbF

       RXN_source_s(2,4) = -( (-RxnaF-RxnbF-RxncF) * RO_ss(2,4)/RO_ss(2,1)  &
            + (RxnbF)*RO_ss(2,4)/RO_ss(2,3) )

       RXN_source_s(2,5) = 0.1d0*RxnbF
!end
!33333333333333333333333333333333333333333333333333333333333333333333333333333
! 3.  Determine the heat of reactions in cal/gmol at the
!     temperature T_g or T_s for RXNx (RXNxF and RXNxB). Note that 
!     for exothermic reactions
!     HORx will be negative.
!
!QX
!-----------------------------------------------------------------------------
! Reaction heat
       HORbg = ZERO
       HORbc = ZERO
       HORbt = ZERO
       HORtar = ZERO
       HORbg = 150.d0/4.1816d0 * RxnaF
       HORbc = 150.d0/4.1868d0 * RxnbF
       HORbt = 150.d0/4.1868d0 * RxncF
!end
    
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
! 4. Wring the ODEs need to solve in ODEPACK format  
!    (referring to the manual)
!QX
!-----------------------------------------------------------------------------
! gas, RO_g, T_g, O2, N2, CO, CO2, H2, CH4, tar
      YDOT(1) = ONE/(ONE-Y(18)-Y(13)) * (SUM(RXN_source_g(1:9)))

      YDOT(2) = -(HORtar)/(Y(1)*(ONE-Y(18)-Y(13))*C_pg_ISAT)

      YDOT(3) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(1) - Y(3) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(4) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(1) - Y(4) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(6) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(4) - Y(6) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(7) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(5) - Y(7) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(8) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(6) - Y(8) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(9) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(7) - Y(9) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(10) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(8) - Y(10) * &
           SUM(RXN_source_g(1:9)) )

      YDOT(11) = ONE/(Y(1)*(ONE-Y(18)-Y(13))) * (RXN_source_g(9) - Y(11) * &
           SUM(RXN_source_g(1:9)) )
!     RO_S, EP_S, T_S
      YDOT(12) = ZERO
      YDOT(13) = ZERO
      YDOT(14) = ZERO
      YDOT(15) = ZERO
      YDOT(16) = ZERO
!     RO_S, EP_S, T_S, dp, bio, moisture, char, pore, ash
      if(Y(18) .eq. ZERO) then
         YDOT(17) = ZERO
      else
         YDOT(17) = ONE/Y(18) * (SUM(RXN_source_s(2,1:5)))
      endif

      YDOT(18) = ZERO

      if(Y(18) .eq. ZERO) then
         YDOT(19) = ZERO
      elseif(Y(18) .gt. 1.D-6 .and. C_ps_ISAT .gt. 0.D0) then
         YDOT(19) = -(HORbg + HORbc + HORbt) / ( Y(17)*Y(18)*C_ps_ISAT )
      endif

      YDOT(20) = ZERO

      if(Y(18) .eq. ZERO) then
         YDOT(21) = ZERO
      else
         YDOT(21) = ONE/(Y(17) * Y(18)) * ( RXN_source_s(2,1) - Y(21) * &
              (SUM(RXN_source_s(2,1:5))) )
      endif

      if(Y(18) .eq. ZERO) then
         YDOT(22) = ZERO
      else
         YDOT(22) = ONE/(Y(17) * Y(18)) * ( RXN_source_s(2,2) - Y(22) * &
              (SUM(RXN_source_s(2,1:5))) )
      endif

      if(Y(18) .eq. ZERO) then
         YDOT(23) = ZERO
      else
         YDOT(23) = ONE/(Y(17) * Y(18)) * ( RXN_source_s(2,3) - Y(23) * &
              (SUM(RXN_source_s(2,1:5))) )
      endif

      if(Y(18) .eq. ZERO) then
         YDOT(24) = ZERO
      else
         YDOT(24) = ONE/(Y(17) * Y(18)) * ( RXN_source_s(2,4) - Y(24) * &
              (SUM(RXN_source_s(2,1:5))) )
      endif

      if(Y(18) .eq. ZERO) then
         YDOT(25) = ZERO
      else
         YDOT(25) = ONE/(Y(17) * Y(18)) * ( RXN_source_s(2,5) - Y(25) * &
              (SUM(RXN_source_s(2,1:5))) )
      endif

      RETURN
      END
