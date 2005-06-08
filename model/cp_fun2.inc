!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: cp_fun2.inc                                            C
!  Purpose: This is an include file containing all CP functions        C
!                                                                      C
!  Author: S. Venkatesan                              Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: Noted at the appropriate location   C
!                                                                      C
!  Variables referenced: None.                                         C
!  Variables modified: None.                                           C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!                      Statement Function Definitions for the various
!                         constant-pressure heat capacities (Cp's)!
!
      CPO2(XXX)= (8.27 + 0.000258*XXX - 187700.0/XXX**2 )/32.
      CPCO(XXX)= ( 6.6 + 0.0012*XXX )/28.
      CPCO2(XXX)= ( 10.34 + 0.00274*XXX - 195500.0/XXX**2 )/44.
      CPCH4(XXX)= ( 5.34 + 0.0115*XXX )/16.
!      CPH2(XXX)= ( 6.62 + 0.00081*XXX )/2.
      CPH2O(XXX)= ( 8.22 + 0.00015*XXX + 0.00000134*XXX**2 )/18.
      CPH2S(XXX)= (7.2 + 0.0036 * XXX )/34.
!      CPN2(XXX)= ( 6.5 + 0.001*XXX )/28.
      CPNH3(XXX)= ( 6.7 + 0.0063*XXX )/17.
!     Perry & Chilton,  page 3-136 (tar oils - 200 C)
      CPTAR(XXX)= 0.45
!     Smith & Van Ness (1959), p. 122 (Table 5-1)
      CPC2H4(XXX)= (2.83 + 28.6E-3*XXX - 8.726E-6*XXX**2)/28.
      CPC2H6(XXX)= (2.247 + 38.2E-3*XXX - 11.05E-6*XXX**2)/30.
      CPC3H8(XXX)= (2.41 + 57.2E-3*XXX - 17.53E-6*XXX**2)/44.
      CPC6H6(XXX)= (-0.409 + 77.62E-3*XXX - 26.4E-6*XXX**2)/78.
!
!     SPECIFIC HEAT OF COAL -- IGT Coal COnversion Systems Technical Data 
!     Book (Sec. IA30.6, pp 3-4)
!
      CPFC(XXX)= -0.1315 + 1.341E-3*XXX - 1.087E-6*XXX*XXX &
                + 3.06E-10*XXX*XXX*XXX
      CPVM(XXX)= 0.1743 + 8.1E-4*XXX
      CPM(XXX)= 1.00763
      CPA(XXX)= 0.1442 + 1.4E-4*XXX
!
!
      CPSiH4(XXX)   = 14.13/32.12
      CPSiH2(XXX)   = 9.413/30.11
      CPH2(XXX)     = ( 6.62 + 0.00081*XXX )/2.02
      CPSi2H6(XXX)  = 43.68/62.22
      CPN2(XXX)     = ( 6.5 + 0.001*XXX )/28.01
      
      CPSi(XXX)     = 0.1442  !?
      CPAl2O3(XXX)  = 0.1442  !?
