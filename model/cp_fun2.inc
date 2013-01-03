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
      CPO2(XXX)= (8.27d0 + 0.000258d0*XXX - 187700.0d0/XXX**2 )/32.d0
      CPCO(XXX)= ( 6.6d0 + 0.0012d0*XXX )/28.d0
      CPCO2(XXX)= ( 10.34d0 + 0.00274d0*XXX - 195500.0d0/XXX**2 )/44.d0
      CPCH4(XXX)= ( 5.34d0 + 0.0115d0*XXX )/16.d0
!      CPH2(XXX)= ( 6.62d0 + 0.00081d0*XXX )/2.d0
      CPH2O(XXX)= ( 8.22d0 + 0.00015d0*XXX + 0.00000134d0*XXX**2 )/18.d0
      CPH2S(XXX)= (7.2d0 + 0.0036d0 * XXX )/34.d0
!      CPN2(XXX)= ( 6.5d0 + 0.001d0*XXX )/28.d0
      CPNH3(XXX)= ( 6.7d0 + 0.0063d0*XXX )/17.d0
!     Perry & Chilton,  page 3-136 (tar oils - 200 C)
!     Smith & Van Ness (1959), p. 122 (Table 5-1)
      CPC2H4(XXX)= (2.83d0 + 28.6D-3*XXX - 8.726D-6*XXX**2)/28.d0
      CPC2H6(XXX)= (2.247d0 + 38.2D-3*XXX - 11.05D-6*XXX**2)/30.d0
      CPC3H8(XXX)= (2.41d0 + 57.2D-3*XXX - 17.53D-6*XXX**2)/44.d0
      CPC6H6(XXX)= (-0.409d0 + 77.62D-3*XXX - 26.4D-6*XXX**2)/78.d0
!
!     SPECIFIC HEAT OF COAL -- IGT Coal COnversion Systems Technical Data 
!     Book (Sec. IA30.6, pp 3-4)
!
      CPFC(XXX)= -0.1315d0 + 1.341D-3*XXX - 1.087D-6*XXX*XXX &
                + 3.06D-10*XXX*XXX*XXX
      CPVM(XXX)= 0.1743d0 + 8.1D-4*XXX
      CPM(XXX)= 1.00763d0
      CPA(XXX)= 0.1442d0 + 1.4D-4*XXX
!
!
      CPSiH4(XXX)   = 14.13d0/32.12d0
      CPSiH2(XXX)   = 9.413d0/30.11d0
      CPH2(XXX)     = ( 6.62d0 + 0.00081d0*XXX )/2.02d0
      CPSi2H6(XXX)  = 43.68d0/62.22d0
      CPN2(XXX)     = ( 6.5d0 + 0.001d0*XXX )/28.01d0
      
      CPSi(XXX)     = 0.1442d0  !?
      CPAl2O3(XXX)  = 0.1442d0  !?
!QX

      CPGAS(XXX)    = 0.263d0 !(cal/g K)
      CPTAR(XXX)    = 0.598d0 !(cal/g K)
      CPBIO(XXX)    = 0.55d0  !(cal/g K)
      CPBIOCHAR(XXX)   = 1.6d0/4.1868d0 !(cal/g K)
