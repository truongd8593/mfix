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
!
!                      Loop indices
      INTEGER          NL

!
!      The number of the cell
!
       IJK = NEQ(2)
!
!
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
!
!  SOLIDS SPECIES
!
!
!33333333333333333333333333333333333333333333333333333333333333333333333333333
! 3.  Determine the heat of reactions in cal/gmol at the
!     temperature T_g or T_s for RXNx (RXNxF and RXNxB). Note that 
!     for exothermic reactions
!     HORx will be negative.
!    
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
! 4. Wring the ODEs need to solve in ODEPACK format  
!    (referring to the manual)
!


      RETURN
      END
