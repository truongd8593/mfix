!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DRAG_gs (M, IER)                                       C
!  Purpose: Calculate the gas-solids drag coefficient                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To use Tsuji drag correlation also if needed in DES        C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: Added switch function to Gidaspow - GIDASPOW_BLEND         C
!  Author: Charles E.A. Finney                        Date: 23-Mar-06  C
!  Reviewer: Sreekanth Pannala                        Date: 24-Mar-06  C
!                                                                      C
!  Literature/Document References:                                     C
!    Syamlal-O'Brien:                                                  C
!      Syamlal M, O'Brien TJ (1988). International Journal of          C
!        Multiphase Flow 14: 473-481.                                  C
!    Gidaspow:                                                         C
!      Ding J, Gidaspow D (1990). AIChE Journal 36: 523-538.           C
!    Gidaspow,blended (original source unknown):                       C
!      Lathouwers D, Bellan J (2000). Proceedings of the 2000 U.S. DOE C
!        Hydrogen Program Review NREL/CP-570-28890. Available from     C
!      http://www.eere.energy.gov/hydrogenandfuelcells/pdfs/28890k.pdf.C
!    Wen-Yu:                                                           C
!      Wen CY, Yu YH (1966). Chemical Engineering Progress Symposium   C
!        Series 62: 100-111.                                           C
!    Koch-Hill (modified):                                             C
!      Benyahia S, Syamlal M, O'Brien TJ (2006). Powder Technology     C
!        162: 166-174.                                                 C
!      Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics    C
!        448: 213-241.                                                 C
!      Hill RJ, Koch DL, Ladd JC (2001). Journal of Fluid Mechanics    C
!        448: 243-278.                                                 C
!                                                                      C
!  Variables referenced: EP_g, RO_g, MU_g, D_p                         C
!  Variables modified: DRAG_gs                                         C
!                                                                      C
!  Local variables: A, B, V_rm, Re                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE DRAG_GS(M, IER) 
!...  Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...  Switches: -xf
!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE usr
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar  
      USE drag  
      USE sendrecv 
      USE discretelement
      USE ur_facs 
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     Solids phase 
      INTEGER          M 
!     
!     Error index 
      INTEGER          IER 
!     
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1) 
!     defined in the data file.
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0 
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0 
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0 
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0 
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     
!     Indices 
      INTEGER          I,  IJK, IMJK, IJMK, IJKM, Im
!     
!     Cell center value of U_sm 
      DOUBLE PRECISION USCM 
!     
!     Cell center value of U_g 
      DOUBLE PRECISION UGC 
!     
!     Cell center value of V_sm 
      DOUBLE PRECISION VSCM 
!     
!     Cell center value of V_g 
      DOUBLE PRECISION VGC 
!     
!     Cell center value of W_sm 
      DOUBLE PRECISION WSCM 
!     
!     Cell center value of W_g 
      DOUBLE PRECISION WGC 
!     
!     Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION VREL 
!     
!     Reynolds number 
      DOUBLE PRECISION Re 
!     
!     Ratio of settling velocity of a multiparticle 
!     system to that of a single particle 
      DOUBLE PRECISION V_rm 
!     
!     Function of EP_g 
      DOUBLE PRECISION A 
!     
!     Function of EP_g 
      DOUBLE PRECISION B 
!     
!     Single sphere drag coefficient x Re 
      DOUBLE PRECISION C_DsxRe, C_DsxReT 
!     
!     single sphere drag coefficient 
      DOUBLE PRECISION C_d 
!     
!     drag coefficient 
      DOUBLE PRECISION DgA  

! --- Gidaspow switch function variables [ceaf 2006-03-23]
      DOUBLE PRECISION Ergun
      DOUBLE PRECISION WenYu
      DOUBLE PRECISION PHI_gs
! --- end Gidaspow switch function variables

!     
!     Gas Laminar viscosity redefined here to set
!     viscosity at pressure boundaries
      DOUBLE PRECISION Mu
!     
!     Gidaspow Reynolds number
      DOUBLE PRECISION Re_g
!
!***********************************************************
!     Declaration of variables relevant to the Koch and Hill
!     drag correlation, sof
!***********************************************************     
!     Stokes Drag Force
      DOUBLE PRECISION F_STOKES
!
!      zero Re function for low Reynolds number
       DOUBLE PRECISION F_0
!
!      inertial function for low Reynolds number
       DOUBLE PRECISION F_1
!
!      zero Re function for high Reynolds number
       DOUBLE PRECISION F_2
!
!      inertial function for high Reynolds number
       DOUBLE PRECISION F_3
!      
!      dimensionless drag force F
       DOUBLE PRECISION F
!      
!      transition Reynolds numbers
       DOUBLE PRECISION Re_Trans_1, Re_Trans_2
!      
!      solids volume fraction
       DOUBLE PRECISION phis
!      
!      weighting factor to compute F_0 and F_2
       DOUBLE PRECISION w, D_p_av, Y_i
!     
!     Hill and Koch Reynolds number
      DOUBLE PRECISION Re_kh
!
!     End of Koch and Hill variables declaration, sof
!***********************************************************
       DOUBLE PRECISION dxf, dxxf, dyf, dyyf, dxyf, dxxxf, dyyyf, dxxyf, dyyxf, &
                        dxxxxf, dyyyyf, dxxxyf, dyyyxf, dxxyyf, d5xf, d5yf, d4xyf, d4yxf, d3x2yf, d3y2xf, &
			d6xf, d6yf, d5xyf, d5yxf, d4x2yf, d4y2xf, d3x3yf, &
			d7xf, d7yf, d6xyf, d6yxf, d5x2yf, d5y2xf, d4x3yf, d4y3xf
       DOUBLE PRECISION term1(IMAX2), term2(IMAX2), term3(IMAX2), term4(IMAX2), &
                        term5(IMAX2), term6(IMAX2), term7(IMAX2), term8(IMAX2), &
                        term9(IMAX2), term10(IMAX2), term11(IMAX2), term12(IMAX2), &
                        term13(IMAX2), term14(IMAX2), term15(IMAX2), term16(IMAX2), &
                        term17(IMAX2), term18(IMAX2), term19(IMAX2), term20(IMAX2), &
			term21(IMAX2), term22(IMAX2), term23(IMAX2)
       DOUBLE PRECISION term24(IMAX2), term25(IMAX2), term26(IMAX2), term27(IMAX2), term28(IMAX2), term29(IMAX2)
       DOUBLE PRECISION term30(IMAX2), term31(IMAX2), term32(IMAX2), term33(IMAX2),term34(IMAX2)
       DOUBLE PRECISION term35(IMAX2), term36(IMAX2), term37(IMAX2), term38(IMAX2),term39(IMAX2)
       DOUBLE PRECISION term40(IMAX2), term41(IMAX2), term42(IMAX2), term43(IMAX2), term44(IMAX2), &
                        term45(IMAX2), term46(IMAX2), term47(IMAX2)
       Double Precision term48(IMAX2),term49(IMAX2),term50(IMAX2),term51(IMAX2),term52(IMAX2),term53(IMAX2), &
                        term54(IMAX2),term55(IMAX2),term56(IMAX2),term57(IMAX2),term58(IMAX2),term59(IMAX2), &
			term60(IMAX2),term61(IMAX2),term62(IMAX2)
       Double Precision epg, dp, Vrsm, rhoG, Drag_corr
!     
!     
!     Current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/RE ! Tsuji drag
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2 ! Dalla Valle (1948) 
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0 * Re**0.657D0)      ! Turton and
!     &          + 0.413D0 * Re**2.09D0 / (Re**1.09D0 + 16300.D0) ! Levenspiel (1986)
!     
!     
!     
!     !$omp  parallel do private( I,  IJK, IMJK, IJMK, IJKM, &
!     !$omp&  USCM, VSCM, WSCM, &
!     !$omp&  VREL, UGC, VGC, WGC, Re, V_rm, A, B) &
!     !$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUIDorP_FLOW_AT(IJK)) THEN 
!     
            I = I_OF(IJK) 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
!!!!!!!!!!!!!!!!!!!!

!     Computing drag force: BETA * (Vg-Vs) 
         
       
          Mu = MU_G(IJK) 
	  epg = Ep_g_mean(I)
	  dp = D_P(IJK,M)
	  rhoG = RO_g(IJK)
	  Vrsm = ABS(V_g_mean(I) - V_S_mean(I))
	  
f = 18.0*(1-epg)*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/(epg**2.65*dp**2)

dxf = -18.0*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**2.65*dp**2)-47.7*(1-epg)*Mu*Vrsm*(0.15*abs(Vrsm)**0.687* &
     (epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**3.65*dp**2)+1.8549*(1-epg)*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**2.65*dp*(epg*dp*rhoG/Mu)**0.313)

dxxf = 95.39999999999999*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**3.65*dp**2)+174.105*(1-epg)*Mu*Vrsm*(0.15*abs(Vrsm)** &
     0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**4.65*dp**2)-3.7098*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**2.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-9.830969999999999*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG/(epg**3.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-0.5805837*(1-epg)*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**2.65*Mu*(epg*dp*rhoG/Mu)**1.313)

dyf = 18.0*(1-epg)*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**2.65*dp**2)+1.8549*(1-epg)*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/(epg**2.65*abs(Vrsm)**1.313*dp**2)

dyyf = 5.5647*(1-epg)*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**1.313*dp**2)-2.4354837*(1-epg)*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)

dxyf = -18.0*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**2.65*dp**2)-47.7*(1-epg)*Mu*(0.15*abs(Vrsm)**0.687* &
     (epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**3.65*dp**2)-1.8549*Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**1.313*dp**2)-4.915484999999999*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+1.8549*(1-epg)*abs(Vrsm)** &
     0.687*rhoG/(epg**2.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1.2743163*(1-epg)*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/Mu)**0.313)

dxxxf = -522.3149999999999*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**4.65*dp**2)-809.58825*(1-epg)*Mu*Vrsm*(0.15*abs(Vrsm)** &
     0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**5.65*dp**2)+29.49291*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**3.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+53.82456074999999*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG/(epg**4.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1.7417511*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**2.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+4.615640414999999*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG**2/(epg**3.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+0.7623063981*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp*rhoG**3 &
     /(epg**2.65*Mu**2*(epg*dp*rhoG/Mu)**2.313)

dyyyf = 5.5647*(1-epg)*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**1.313*dp**2)-14.6129022*(1-epg)*Mu*Vrsm** &
     2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+8.068757498099998*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)

dxxyf = 95.39999999999999*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**3.65*dp**2)+174.105*(1-epg)*Mu*(0.15*abs(Vrsm)**0.687* &
     (epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**4.65*dp**2)+9.830969999999999*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+17.94152025*(1-epg)*Mu*Vrsm**2 &
     *(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-3.7098*abs(Vrsm)**0.687* &
     rhoG/(epg**2.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-9.830969999999999*(1-epg)*abs(Vrsm)**0.687*rhoG/ &
     (epg**3.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-2.5486326*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-6.753876389999999*(1-epg)*Vrsm**2*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-0.5805837*(1-epg)*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**2.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-0.3988610019*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/Mu)**1.313)

dyyxf = -5.5647*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**1.313*dp**2)-14.746455*(1-epg)*Mu*Vrsm* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+2.4354837*Mu*Vrsm**3*(epg* &
     dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+6.454031805*(1-epg)*Mu*Vrsm**3 &
     *(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)+3.822948900000001*(1-epg)* &
     Vrsm*rhoG/(epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1.6731773019*(1-epg)*Vrsm**3*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/Mu)**0.313)

dxxxxf = 3238.353*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**5.65*dp**2)+4574.173612500001*(1-epg)*Mu*Vrsm*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**6.65*dp**2)-215.298243*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**4.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-333.71227665*(1-epg)*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**5.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-18.46256166*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**3.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-33.69417502949999*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG**2/(epg**4.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-3.049225592399999*Vrsm*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**2.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-8.080447819859998*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**3.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-1.7632146988053*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp** &
     2*rhoG**4/(epg**2.65*Mu**3*(epg*dp*rhoG/Mu)**3.313)

dyyyyf = -36.53225550000001*(1-epg)*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+80.68757498099998*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-42.86930858740529*(1-epg)* &
     Mu*Vrsm**5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)

dxxxyf = -522.3149999999999*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**4.65*dp**2)-809.58825*(1-epg)*Mu*(0.15*abs(Vrsm)** &
     0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**5.65*dp**2)-53.82456075*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-83.42806916250001*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+29.49291*abs(Vrsm)**0.687* &
     rhoG/(epg**3.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+53.82456074999999*(1-epg)*abs(Vrsm)**0.687*rhoG/ &
     (epg**4.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+20.26162917*Vrsm**2*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+36.97747323525*(1-epg)*Vrsm**2*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1.7417511*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**2.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+4.615640414999999*(1-epg)*abs(Vrsm)**0.687*rhoG** &
     2/(epg**3.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+1.1965830057*Vrsm**2*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+3.170944965105*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+0.7623063981*(1-epg)*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**2.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+0.5237044954947*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/Mu)**2.313)

dyyyxf = -5.5647*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**1.313*dp**2)-14.746455*(1-epg)*Mu*(epg* &
     dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+14.6129022*Mu*Vrsm**2*(epg* &
     dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+38.72419083*(1-epg)*Mu*Vrsm**2 &
     *(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-8.068757498099998*Mu*Vrsm** &
     4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-21.382207369965*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+3.822948900000001*(1-epg)* &
     rhoG/(epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-10.0390638114*(1-epg)*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+5.543236401194699*(1-epg)*Vrsm**4*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/Mu)**0.313)

dxxyyf = 29.49291*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+53.82456075*(1-epg)*Mu*Vrsm* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-12.90806361*Mu*Vrsm**3*(epg* &
     dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-23.55721608825*(1-epg)*Mu* &
     Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)-7.645897800000001*Vrsm*rhoG/ &
     (epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-20.26162917*(1-epg)*Vrsm*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+3.346354603800001*Vrsm**3*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+8.867839700069998*(1-epg)*Vrsm**3*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1.1965830057*(1-epg)*Vrsm*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+0.5237044954947*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/Mu)**1.313)

d5xf = -22870.86806250001*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**6.65*dp**2)-30418.25452312501*(1-epg)*Mu*Vrsm*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**7.65*dp**2)+1668.56138325*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**5.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+2356.842953840625*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG/(epg**6.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+168.4708751475*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**4.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+261.1298564786249*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG**2/(epg**5.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+40.40223909929999*Vrsm*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**3.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+73.73408635622246*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**4.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+8.816073494026497*Vrsm*abs(Vrsm)**0.687*dp**2*rhoG**4 &
     /(epg**2.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+23.36259475917022*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**2 &
     *rhoG**4/(epg**3.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+5.841530297141957*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**3 &
     *rhoG**5/(epg**2.65*Mu**4*(epg*dp*rhoG/Mu)**4.313)

d5yf = -36.53225550000001*(1-epg)*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+363.0940874144999*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-643.0396288110793*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+313.5032536996948*(1-epg)* &
     Mu*Vrsm**6*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp**2)

d4xyf = 3238.353*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**5.65*dp**2)+4574.173612500001*(1-epg)*Mu*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**6.65*dp**2)+333.71227665*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+471.3685907681251*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**1.313*dp**2)-215.298243*abs(Vrsm)**0.687* &
     rhoG/(epg**4.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-333.71227665*(1-epg)*abs(Vrsm)**0.687*rhoG/ &
     (epg**5.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-147.909892941*Vrsm**2*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-229.26033405855*(1-epg)*Vrsm**2*rhoG/ &
     (epg**5.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-18.46256166*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**3.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-33.69417502949999*(1-epg)*abs(Vrsm)**0.687*rhoG** &
     2/(epg**4.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-12.68377986042*Vrsm**2*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-23.14789824526649*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-3.049225592399999*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**2.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-8.080447819859998*(1-epg)*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**3.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-2.0948179819788*Vrsm**2*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-5.551267652243819*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-1.7632146988053*(1-epg)*abs(Vrsm)**0.687*dp**2* &
     rhoG**4/(epg**2.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)-1.211328498079241*(1-epg)*Vrsm**2*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/Mu)**3.313)

d4yxf = 36.53225550000001*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+96.81047707500001*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-80.68757498099998*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-213.8220736996499*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+42.86930858740529*Mu*Vrsm** &
     5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+113.603667756624*(1-epg)* &
     Mu*Vrsm**5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**7.313*dp**2)-25.09765952850001*(1-epg)* &
     Vrsm*rhoG/(epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+55.43236401194699*(1-epg)*Vrsm**3*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-29.45121499954744*(1-epg)*Vrsm**5*rhoG/ &
     (epg**2.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/Mu)**0.313)

d3x2yf = -161.47368225*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-250.2842074875*(1-epg)*Mu* &
     Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+70.67164826474999*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+109.5410548103625*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**3.313*dp**2)+60.78488751*Vrsm*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+110.93241970575*(1-epg)*Vrsm*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-26.60351910021*Vrsm**3*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-48.55142235788324*(1-epg)*Vrsm**3*rhoG/ &
     (epg**4.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+3.589749017099999*Vrsm*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+9.512834895314999*(1-epg)*Vrsm*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-1.5711134864841*Vrsm**3*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-4.163450739182864*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+1.5711134864841*(1-epg)*Vrsm*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-0.68762400258454*(1-epg)*Vrsm**3*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/Mu)**2.313)

d3y2xf = 29.49291*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**1.313*dp**2)+53.82456075*(1-epg)*Mu* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-77.44838166*Mu*Vrsm**2*(epg* &
     dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-141.3432965295*(1-epg)*Mu* &
     Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+42.76441473992999*Mu*Vrsm** &
     4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+78.04505690037223*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**5.313*dp**2)-7.645897800000001*rhoG/ &
     (epg**2.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-20.26162917*(1-epg)*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+20.0781276228*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+53.20703820042*(1-epg)*Vrsm**2*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-11.0864728023894*Vrsm**4*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-29.37915292633191*(1-epg)*Vrsm**4*rhoG/ &
     (epg**3.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1.1965830057*(1-epg)*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+3.1422269729682*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-1.735032993573941*(1-epg)*Vrsm**4*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/Mu)**1.313)

d6xf = 182509.52713875*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**7.65*dp**2)+232699.6471019063*(1-epg)*Mu*Vrsm*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**8.65*dp**2)-14141.05772304375*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**6.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-18807.60677164819*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG/(epg**7.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1566.779138871749*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**5.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-2213.075533656346*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG**2/(epg**6.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-442.4045181373348*Vrsm*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**4.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-685.7270031128689*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**5.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-140.1755685550213*Vrsm*abs(Vrsm)**0.687*dp**2*rhoG**4 &
     /(epg**3.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)-255.8204126129138*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**2 &
     *rhoG**4/(epg**4.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)-35.04918178285174*Vrsm*abs(Vrsm)**0.687*dp**3*rhoG**5 &
     /(epg**2.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)-92.88033172455711*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**3 &
     *rhoG**5/(epg**3.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)-25.19452017157326*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**4 &
     *rhoG**6/(epg**2.65*Mu**5*(epg*dp*rhoG/Mu)**5.313)

d6yf = 847.2195373004998*(1-epg)*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-4501.277401677555*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+6583.568327693592*(1-epg)* &
     Mu*Vrsm**5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp** &
     2)-2919.655801705258*(1-epg)*Mu*Vrsm**7*(epg*dp*rhoG/ &
     Mu)**0.687/(epg**2.65*abs(Vrsm)**11.313*dp**2)

d5xyf = -22870.86806250001*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**6.65*dp**2)-30418.25452312501*(1-epg)*Mu*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**7.65*dp**2)-2356.842953840626*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**1.313*dp**2)-3134.601128608032*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**7.65*abs(Vrsm)**1.313*dp**2)+1668.56138325*abs(Vrsm)** &
     0.687*rhoG/(epg**5.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+2356.842953840625*(1-epg)*abs(Vrsm)**0.687*rhoG/ &
     (epg**6.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1146.30167029275*Vrsm**2*rhoG/ &
     (epg**5.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1619.151109288509*(1-epg)*Vrsm**2*rhoG/ &
     (epg**6.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+168.4708751475*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**4.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+261.1298564786249*(1-epg)*abs(Vrsm)**0.687*rhoG** &
     2/(epg**5.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+115.7394912263325*Vrsm**2*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+179.3962114008153*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**5.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+40.40223909929999*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**3.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+73.73408635622246*(1-epg)*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**4.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+27.7563382612191*Vrsm**2*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)+50.65531732672483*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**4.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)+8.816073494026497*abs(Vrsm)**0.687*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+23.36259475917022*(1-epg)*abs(Vrsm)**0.687*dp**2* &
     rhoG**4/(epg**3.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+6.056642490396204*Vrsm**2*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)+16.05010259954994*(1-epg)*Vrsm**2*dp**2*rhoG**4/ &
     (epg**3.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)+5.841530297141957*(1-epg)*abs(Vrsm)**0.687*dp**3* &
     rhoG**5/(epg**2.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)+4.013131314136524*(1-epg)*Vrsm**2*dp**3*rhoG**5/ &
     (epg**2.65*Mu**4*abs(Vrsm)**1.313*(epg*dp*rhoG/Mu)**4.313)

d5yxf = 36.53225550000001*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**3.313*dp**2)+96.81047707500001*(1-epg)* &
     Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-363.0940874144999*Mu*Vrsm** &
     2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-962.1993316484249*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+643.0396288110793*Mu*Vrsm** &
     4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+1704.05501634936*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**7.313*dp**2)-313.5032536996948*Mu*Vrsm** &
     6*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp** &
     2)-830.7836223041912*(1-epg)*Mu*Vrsm**6*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**9.312999999999999*dp** &
     2)-25.09765952850001*(1-epg)*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+249.4456380537615*(1-epg)*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-441.7682249932115*(1-epg)*Vrsm**4*rhoG/ &
     (epg**2.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+215.3767352916904*(1-epg)*Vrsm**6*rhoG/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp*(epg*dp*rhoG/ &
     Mu)**0.313)

d4x2yf = 1001.13682995*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+1414.105772304375*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**1.313*dp**2)-438.16421924145*Mu*Vrsm**3* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**3.313*dp**2)-618.9069596785482*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**3.313*dp**2)-443.729678823*Vrsm*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-687.78100217565*(1-epg)*Vrsm*rhoG/ &
     (epg**5.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+194.205689431533*Vrsm**3*rhoG/ &
     (epg**4.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+301.0188186188761*(1-epg)*Vrsm**3*rhoG/ &
     (epg**5.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-38.05133958126*Vrsm*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-69.44369473579948*(1-epg)*Vrsm*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+16.65380295673146*Vrsm**3*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+30.3931903960349*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-6.284453945936399*Vrsm*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-16.65380295673146*(1-epg)*Vrsm*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)+2.750496010338164*Vrsm**3*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)+7.288814427396134*(1-epg)*Vrsm**3*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)-3.633985494237723*(1-epg)*Vrsm*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)+1.590474317978043*(1-epg)*Vrsm**3*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**3.313*(epg*dp*rhoG/Mu)**3.313)

d4y2xf = -193.62095415*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-353.35824132375*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+427.6441473992999*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+780.4505690037222*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**5.313*dp**2)-227.207335513248*Mu*Vrsm**5* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**7.313*dp**2)-414.6533873116776*(1-epg)* &
     Mu*Vrsm**5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**7.313*dp**2)+50.19531905700001*Vrsm*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+133.01759550105*(1-epg)*Vrsm*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-110.864728023894*Vrsm**3*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-293.791529263319*(1-epg)*Vrsm**3*rhoG/ &
     (epg**3.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+58.90242999909487*Vrsm**5*rhoG/ &
     (epg**2.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+156.0914394976014*(1-epg)*Vrsm**5*rhoG/ &
     (epg**3.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+7.855567432420501*(1-epg)*Vrsm*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-17.35032993573941*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)+9.218230294858346*(1-epg)*Vrsm**5*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**7.313*(epg*dp*rhoG/Mu)**1.313)

d3x3yf = -161.47368225*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**1.313*dp**2)-250.2842074875*(1-epg)*Mu* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+424.0298895884999*Mu*Vrsm** &
     2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+657.246328862175*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**3.313*dp**2)-234.1351707011167*Mu*Vrsm** &
     4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**5.313*dp**2)-362.909514586731*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**5.313*dp**2)+60.78488751*rhoG/ &
     (epg**3.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+110.93241970575*(1-epg)*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-159.62111460126*Vrsm**2*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-291.3085341472995*(1-epg)*Vrsm**2*rhoG/ &
     (epg**4.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+88.13745877899571*Vrsm**4*rhoG/ &
     (epg**3.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+160.8508622716672*(1-epg)*Vrsm**4*rhoG/ &
     (epg**4.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+3.589749017099999*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+9.512834895314999*(1-epg)*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-9.426680918904598*Vrsm**2*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-24.98070443509719*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+5.205098980721822*Vrsm**4*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)+13.79351229891283*(1-epg)*Vrsm**4*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)+1.5711134864841*(1-epg)*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-4.125744015507246*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)+2.278098320562584*(1-epg)*Vrsm**4*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**5.313*(epg*dp*rhoG/Mu)**2.313)

d7xf = -1628897.529713344*Mu*Vrsm*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**8.65*dp**2)-2012851.94743149*(1-epg)*Mu*Vrsm*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**9.65*dp**2)+131653.2474015373*Vrsm*abs(Vrsm)**0.687*rhoG/ &
     (epg**7.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+167857.8904369601*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG/(epg**8.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+15491.52873559442*Vrsm*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**6.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+20603.73321834059*(1-epg)*Vrsm*abs(Vrsm)**0.687* &
     rhoG**2/(epg**7.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)+4800.089021790083*Vrsm*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**5.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+6780.125743278491*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**6.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)+1790.742888290397*Vrsm*abs(Vrsm)**0.687*dp**2*rhoG**4 &
     /(epg**4.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+2775.651476850115*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**2 &
     *rhoG**4/(epg**5.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)+650.1623220718998*Vrsm*abs(Vrsm)**0.687*dp**3*rhoG**5 &
     /(epg**3.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)+1186.546237781217*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**3 &
     *rhoG**5/(epg**4.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)+176.3616412010128*Vrsm*abs(Vrsm)**0.687*dp**4*rhoG**6 &
     /(epg**2.65*Mu**5*(epg*dp*rhoG/ &
     Mu)**5.313)+467.3583491826839*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**4 &
     *rhoG**6/(epg**3.65*Mu**5*(epg*dp*rhoG/ &
     Mu)**5.313)+133.8584856715687*(1-epg)*Vrsm*abs(Vrsm)**0.687*dp**5 &
     *rhoG**7/(epg**2.65*Mu**6*(epg*dp*rhoG/Mu)**6.313)

d7yf = 847.2195373004998*(1-epg)*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-18005.10960671022*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+65835.68327693592*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp** &
     2)-81750.36244774722*(1-epg)*Mu*Vrsm**6*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**11.313*dp**2)+33030.06608469158* &
     (1-epg)*Mu*Vrsm**8*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**13.313*dp**2)

d6xyf = 182509.52713875*Mu*(0.15*abs(Vrsm)**0.687*(epg*dp*rhoG/ &
     Mu)**0.687+1)/ &
     (epg**7.65*dp**2)+232699.6471019063*(1-epg)*Mu*(0.15* &
     abs(Vrsm)**0.687*(epg*dp*rhoG/Mu)**0.687+1)/ &
     (epg**8.65*dp**2)+18807.60677164819*Mu*Vrsm**2*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**7.65*abs(Vrsm)**1.313*dp**2)+23979.69863385145*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**8.65*abs(Vrsm)**1.313*dp**2)-14141.05772304375*abs(Vrsm)**0.687 &
     *rhoG/(epg**6.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-18807.60677164819*(1-epg)*abs(Vrsm)**0.687*rhoG/ &
     (epg**7.65*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-9714.906655731058*Vrsm**2*rhoG/ &
     (epg**6.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-12920.82585212231*(1-epg)*Vrsm**2*rhoG/ &
     (epg**7.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1566.779138871749*abs(Vrsm)**0.687*rhoG**2/ &
     (epg**5.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-2213.075533656346*(1-epg)*abs(Vrsm)**0.687*rhoG** &
     2/(epg**6.65*Mu*(epg*dp*rhoG/ &
     Mu)**1.313)-1076.377268404892*Vrsm**2*rhoG**2/ &
     (epg**5.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-1520.38289162191*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**6.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-442.4045181373348*abs(Vrsm)**0.687*dp*rhoG**3/ &
     (epg**4.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-685.7270031128689*(1-epg)*abs(Vrsm)**0.687*dp* &
     rhoG**3/(epg**5.65*Mu**2*(epg*dp*rhoG/ &
     Mu)**2.313)-303.9319039603491*Vrsm**2*dp*rhoG**3/ &
     (epg**4.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-471.094451138541*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**5.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-140.1755685550213*abs(Vrsm)**0.687*dp**2*rhoG**4/ &
     (epg**3.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)-255.8204126129138*(1-epg)*abs(Vrsm)**0.687*dp**2* &
     rhoG**4/(epg**4.65*Mu**3*(epg*dp*rhoG/ &
     Mu)**3.313)-96.30061559729964*Vrsm**2*dp**2*rhoG**4/ &
     (epg**3.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)-175.7486234650718*(1-epg)*Vrsm**2*dp**2*rhoG**4/ &
     (epg**4.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)-35.04918178285174*abs(Vrsm)**0.687*dp**3*rhoG**5/ &
     (epg**2.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)-92.88033172455711*(1-epg)*abs(Vrsm)**0.687*dp**3* &
     rhoG**5/(epg**3.65*Mu**4*(epg*dp*rhoG/ &
     Mu)**4.313)-24.07878788481914*Vrsm**2*dp**3*rhoG**5/ &
     (epg**2.65*Mu**4*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**4.313)-63.80878789477074*(1-epg)*Vrsm**2*dp**3*rhoG**5/ &
     (epg**3.65*Mu**4*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**4.313)-25.19452017157326*(1-epg)*abs(Vrsm)**0.687*dp**4* &
     rhoG**6/(epg**2.65*Mu**5*(epg*dp*rhoG/ &
     Mu)**5.313)-17.30863535787083*(1-epg)*Vrsm**2*dp**4*rhoG**6/ &
     (epg**2.65*Mu**5*abs(Vrsm)**1.313*(epg*dp*rhoG/Mu)**5.313)

d6yxf = -847.2195373004998*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**5.313*dp**2)-2245.131773846325*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+4501.277401677555*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**7.313*dp**2)+11928.38511444552*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**7.313*dp**2)-6583.568327693592*Mu*Vrsm** &
     5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp** &
     2)-17446.45606838802*(1-epg)*Mu*Vrsm**5*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**9.312999999999999*dp** &
     2)+2919.655801705258*Mu*Vrsm**7*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**2.65*abs(Vrsm)**11.313*dp**2)+7737.087874518932* &
     (1-epg)*Mu*Vrsm**7*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**11.313*dp**2)+582.0398221254435* &
     (1-epg)*Vrsm*rhoG/(epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-3092.377574952481*(1-epg)*Vrsm**3*rhoG/ &
     (epg**2.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+4522.911441125498*(1-epg)*Vrsm**5*rhoG/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-2005.803535771512*(1-epg)*Vrsm**7*rhoG/ &
     (epg**2.65*abs(Vrsm)**11.313*dp*(epg*dp*rhoG/Mu)**0.313)

d5x2yf = -7070.528861521876*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**1.313*dp**2)-9403.803385824096*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**7.65*abs(Vrsm)**1.313*dp**2)+3094.534798392741*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**3.313*dp**2)+4115.731281862346*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**7.65*abs(Vrsm)**3.313*dp**2)+3438.90501087825*Vrsm*rhoG/ &
     (epg**5.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+4857.453327865528*(1-epg)*Vrsm*rhoG/ &
     (epg**6.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1505.094093094381*Vrsm**3*rhoG/ &
     (epg**5.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-2125.945406495813*(1-epg)*Vrsm**3*rhoG/ &
     (epg**6.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+347.2184736789974*Vrsm*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+538.188634202446*(1-epg)*Vrsm*rhoG**2/ &
     (epg**5.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-151.9659519801745*Vrsm**3*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-235.5472255692705*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**5.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+83.26901478365728*Vrsm*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)+151.9659519801745*(1-epg)*Vrsm*dp*rhoG**3/ &
     (epg**4.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-36.44407213698067*Vrsm**3*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)-66.51043164998971*(1-epg)*Vrsm**3*dp*rhoG**3/ &
     (epg**4.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)+18.16992747118861*Vrsm*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)+48.15030779864982*(1-epg)*Vrsm*dp**2*rhoG**4/ &
     (epg**3.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)-7.952371589890215*Vrsm**3*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**3.313)-21.07378471320907*(1-epg)*Vrsm**3*dp**2*rhoG**4/ &
     (epg**3.65*Mu**3*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**3.313)+12.03939394240958*(1-epg)*Vrsm*dp**3*rhoG**5/ &
     (epg**2.65*Mu**4*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**4.313)-5.269241415461257*(1-epg)*Vrsm**3*dp**3*rhoG**5/ &
     (epg**2.65*Mu**4*abs(Vrsm)**3.313*(epg*dp*rhoG/Mu)**4.313)

d5y2xf = -193.62095415*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**3.313*dp**2)-353.35824132375*(1-epg)* &
     Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+1924.39866329685*Mu*Vrsm**2* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**5.313*dp**2)+3512.02756051675*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**5.313*dp**2)-3408.11003269872*Mu*Vrsm**4* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**7.313*dp**2)-6219.800809675164*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**7.313*dp**2)+1661.567244608383*Mu*Vrsm** &
     6*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**3.65*abs(Vrsm)**9.312999999999999*dp** &
     2)+3032.360221410298*(1-epg)*Mu*Vrsm**6*(epg*dp*rhoG/ &
     Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**9.312999999999999*dp** &
     2)+50.19531905700001*rhoG/ &
     (epg**2.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+133.01759550105*(1-epg)*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-498.891276107523*Vrsm**2*rhoG/ &
     (epg**2.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1322.061881684936*(1-epg)*Vrsm**2*rhoG/ &
     (epg**3.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+883.5364499864231*Vrsm**4*rhoG/ &
     (epg**2.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+2341.371592464021*(1-epg)*Vrsm**4*rhoG/ &
     (epg**3.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-430.7534705833808*Vrsm**6*rhoG/ &
     (epg**2.65*abs(Vrsm)**9.312999999999999*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-1141.496697045959*(1-epg)*Vrsm**6*rhoG/ &
     (epg**3.65*abs(Vrsm)**9.312999999999999*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+7.855567432420501*(1-epg)*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-78.07648471082733*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)+138.2734544228752*(1-epg)*Vrsm**4*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**7.313*(epg*dp*rhoG/ &
     Mu)**1.313)-67.41291814629908*(1-epg)*Vrsm**6*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**9.312999999999999*(epg*dp*rhoG/ &
     Mu)**1.313)

d4x3yf = 1001.13682995*Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**1.313*dp**2)+1414.105772304375*(1-epg)* &
     Mu*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**1.313*dp**2)-2628.9853154487*Mu*Vrsm**2* &
     (epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**3.313*dp**2)-3713.441758071289*(1-epg)* &
     Mu*Vrsm**2*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**3.313*dp**2)+1451.638058346924*Mu*Vrsm** &
     4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**5.313*dp**2)+2050.43875741503*(1-epg)* &
     Mu*Vrsm**4*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**6.65*abs(Vrsm)**5.313*dp**2)-443.729678823*rhoG/ &
     (epg**4.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-687.78100217565*(1-epg)*rhoG/ &
     (epg**5.65*abs(Vrsm)**1.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1165.234136589198*Vrsm**2*rhoG/ &
     (epg**4.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1806.112911713257*(1-epg)*Vrsm**2*rhoG/ &
     (epg**5.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-643.4034490866686*Vrsm**4*rhoG/ &
     (epg**4.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-997.2753460843366*(1-epg)*Vrsm**4*rhoG/ &
     (epg**5.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-38.05133958126*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)-69.44369473579948*(1-epg)*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**1.313)+99.92281774038875*Vrsm**2*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+182.3591423762095*(1-epg)*Vrsm**2*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-55.17404919565131*Vrsm**4*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)-100.6926397820636*(1-epg)*Vrsm**4*rhoG**2/ &
     (epg**4.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)-6.284453945936399*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)-16.65380295673146*(1-epg)*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**2.313)+16.50297606202899*Vrsm**2*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)+43.73288656437681*(1-epg)*Vrsm**2*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)-9.112393282250336*Vrsm**4*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**2.313)-24.14784219796339*(1-epg)*Vrsm**4*dp*rhoG**3/ &
     (epg**3.65*Mu**2*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**2.313)-3.633985494237723*(1-epg)*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**1.313*(epg*dp*rhoG/ &
     Mu)**3.313)+9.542845907868259*(1-epg)*Vrsm**2*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**3.313)-5.269241415461258*(1-epg)*Vrsm**4*dp**2*rhoG**4/ &
     (epg**2.65*Mu**3*abs(Vrsm)**5.313*(epg*dp*rhoG/Mu)**3.313)

d4y3xf = 1060.07472397125*Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**3.313*dp**2)+1643.115822155437*(1-epg)* &
     Mu*Vrsm*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**3.313*dp**2)-2341.351707011167*Mu*Vrsm** &
     3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**5.313*dp**2)-3629.09514586731*(1-epg)* &
     Mu*Vrsm**3*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**5.313*dp**2)+1243.960161935033*Mu*Vrsm** &
     5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**4.65*abs(Vrsm)**7.313*dp**2)+1928.138250999301*(1-epg)* &
     Mu*Vrsm**5*(epg*dp*rhoG/Mu)**0.687/ &
     (epg**5.65*abs(Vrsm)**7.313*dp**2)-399.05278650315*Vrsm*rhoG/ &
     (epg**3.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-728.2713353682487*(1-epg)*Vrsm*rhoG/ &
     (epg**4.65*abs(Vrsm)**3.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+881.3745877899571*Vrsm**3*rhoG/ &
     (epg**3.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)+1608.508622716672*(1-epg)*Vrsm**3*rhoG/ &
     (epg**4.65*abs(Vrsm)**5.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-468.2743184928042*Vrsm**5*rhoG/ &
     (epg**3.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-854.6006312493676*(1-epg)*Vrsm**5*rhoG/ &
     (epg**4.65*abs(Vrsm)**7.313*dp*(epg*dp*rhoG/ &
     Mu)**0.313)-23.56670229726149*Vrsm*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)-62.45176108774297*(1-epg)*Vrsm*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**1.313)+52.05098980721822*Vrsm**3*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)+137.9351229891283*(1-epg)*Vrsm**3*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**1.313)-27.65469088457504*Vrsm**5*rhoG**2/ &
     (epg**2.65*Mu*abs(Vrsm)**7.313*(epg*dp*rhoG/ &
     Mu)**1.313)-73.28493084412385*(1-epg)*Vrsm**5*rhoG**2/ &
     (epg**3.65*Mu*abs(Vrsm)**7.313*(epg*dp*rhoG/ &
     Mu)**1.313)-10.31436003876812*(1-epg)*Vrsm*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**3.313*(epg*dp*rhoG/ &
     Mu)**2.313)+22.78098320562584*(1-epg)*Vrsm**3*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**5.313*(epg*dp*rhoG/ &
     Mu)**2.313)-12.10353637714901*(1-epg)*Vrsm**5*dp*rhoG**3/ &
     (epg**2.65*Mu**2*abs(Vrsm)**7.313*(epg*dp*rhoG/Mu)**2.313)
	 
	 term1(I) = 0.5d0 * dxxf
	 term2(I) = 0.5d0 * dyyf
	 term3(I) = dxyf
	 term4(I) = dxf
	 term5(I) = dyf
	 term6(I) = ONE/6d0 * dxxxf
	 term7(I) = ONE/6d0 * dyyyf
	 term8(I) = 0.5d0 * dxxyf
	 term9(I) = 0.5d0 * dyyxf
	 term10(I) = ONE/12d0 * dxxxxf
	 term11(I) = ONE/12d0 * dyyyyf
	 term12(I) = ONE/6d0 * dxxxyf
	 term13(I) = ONE/6d0 * dyyyxf
	 term14(I) = 0.25d0 * dxxyyf
	 
	 term15(I) = ONE/120d0 * d5xf
	 term16(I) = ONE/120d0 * d5yf
	 term17(I) = ONE/24d0 * d4xyf
	 term18(I) = ONE/24d0 * d4yxf
	 term19(I) = ONE/12d0 * d3x2yf
	 term20(I) = ONE/12d0 * d3y2xf
	 
	 term21(I) = ONE/720d0 * d6xf
	 term22(I) = ONE/720d0 * d6yf
	 term23(I) = ONE/120d0 * d5xyf
	 term24(I) = ONE/120d0 * d5yxf
	 term25(I) = ONE/48d0 * d4x2yf
	 term26(I) = ONE/48d0 * d4y2xf
	 term27(I) = ONE/36d0 * d3x3yf
	 
	 term28(I) = ONE/5040d0 * d7xf
	 term29(I) = ONE/5040d0 * d7yf
	 term30(I) = ONE/720d0 * d6xyf
	 term31(I) = ONE/720d0 * d6yxf
	 term32(I) = ONE/240d0 * d5x2yf
	 term33(I) = ONE/240d0 * d5y2xf
	 term34(I) = ONE/144d0 * d4x3yf
	 term35(I) = ONE/144d0 * d4y3xf
	 
	 term1(I) = term1(I) * (epg2(I))
	 
	 term2(I) = term2(I) * (Vr2(I))
	 
	 term3(I) = term3(I) * (epgVr(I))
	 
	 term4(I) = term4(I) * ZERO !epgVrs(I)
	 term5(I) = term5(I) * ZERO !VrVrs(I)
	 
	 term6(I) = term6(I) * (epg3(I))
	 term7(I) = term7(I) * (Vr3(I))
	 term8(I) = term8(I) * (epg2Vr(I))
	 term9(I) = term9(I) * (epgVr2(I))
	 
	 term10(I) = term10(I) * (epg4(I))
	 term11(I) = term11(I) * (Vr4(I))
	 term12(I) = term12(I) * (epg3Vr(I))
	 term13(I) = term13(I) * (epgVr3(I))
	 term14(I) = term14(I) * (epg2Vr2(I))
	 
	 term15(I) = term15(I) * (epg5(I))
	 term16(I) = term16(I) * (Vr5(I))
	 term17(I) = term17(I) * (epg4Vr(I))
	 term18(I) = term18(I) * (epgVr4(I))
	 term19(I) = term19(I) * (epg3Vr2(I))
	 term20(I) = term20(I) * (epg2Vr3(I))
	 
	 term21(I) = term21(I) * (epg6(I))
	 term22(I) = term22(I) * (Vr6(I))
	 term23(I) = term23(I) * (epg5Vr(I))
	 term24(I) = term24(I) * (epgVr5(I))
	 term25(I) = term25(I) * (epg4Vr2(I))
	 term26(I) = term26(I) * (epg2Vr4(I))
	 term27(I) = term27(I) * (epg3Vr3(I))
	 
	 term28(I) = term28(I) * (epg7(I))
	 term29(I) = term29(I) * (Vr7(I))
	 term30(I) = term30(I) * (epg6Vr(I))
	 term31(I) = term31(I) * (epgVr6(I))
	 term32(I) = term32(I) * (epg5Vr2(I))
	 term33(I) = term33(I) * (epg2Vr5(I))
	 term34(I) = term34(I) * (epg4Vr3(I))
	 term35(I) = term35(I) * (epg3Vr4(I))

 Drag_corr = f+term1(I)+term2(I)+term3(I)+term6(I)+term7(I)+term8(I)+term9(I)+&
	 term10(I)+term11(I)+term12(I)+term13(I)+term14(I)+&
	 term15(I)+term16(I)+term17(I)+term18(I)+term19(I)+term20(I)+&
	 term21(I)+term22(I)+term23(I)+term24(I)+term25(I)+term26(I)+term27(I)+&
	 term28(I)+term29(I)+term30(I)+term31(I)+term32(I)+term33(I)+term34(I)+term35(I)

! 	write(*,'(I3, 8(1X,E13.6))') I, Drag_corr
            
!     Calculate velocity components at i, j, k
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I) 
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK)) 
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK)) 
            IF(.NOT.DES_CONTINUUM_COUPLED) THEN
               USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
               VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M)) 
               WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M)) 
            ELSE
               USCM = DES_U_S(IJK,M)
               VSCM = DES_V_S(IJK,M)
               WSCM = DES_W_S(IJK,M)
            ENDIF
!     
!     magnitude of gas-solids relative velocity
!     
            VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2) 
!     
!     ! Laminar viscosity at a pressure boundary is given the value of the fluid cell next to it.
!     ! This applies just to the calculation of the drag, in other routines the value of viscosity
!     ! at a pressure boundary has always a zero value.
!     
            IF (P_OUTFLOW_AT(IJK)) THEN
               IF( FLUID_AT(EAST_OF(IJK) )) THEN
                  Mu = MU_G(EAST_OF(IJK))
               ELSE IF ( FLUID_AT(WEST_OF(IJK)) ) THEN
                  Mu = MU_G(WEST_OF(IJK))
               ELSE IF ( FLUID_AT(NORTH_OF(IJK)) ) THEN
                  Mu = MU_G(NORTH_OF(IJK))
               ELSE IF ( FLUID_AT(SOUTH_OF(IJK)) ) THEN
                  Mu = MU_G(SOUTH_OF(IJK))
               ELSE IF ( FLUID_AT(TOP_OF(IJK)) ) THEN
                  Mu = MU_G(TOP_OF(IJK))
               ELSE IF ( FLUID_AT(BOTTOM_OF(IJK)) ) THEN
                  Mu = MU_G(BOTTOM_OF(IJK))
               ENDIF
            ELSE
               Mu = MU_G(IJK)
            ENDIF
            
            

!     Reynolds number
            if(Mu > ZERO)then
               RE = D_P(IJK,M)*VREL*RO_G(IJK)/Mu

!     Note the presence of gas volume fraction in ROP_G
               RE_G = D_P(IJK,M)*VREL*ROP_G(IJK)/Mu

!     Note Reynolds' number for Hill and Koch has an additional factor of 1/2 & ep_g
               RE_kh = 0.5D0*D_P(IJK,M)*VREL*ROP_G(IJK)/Mu
            else 
               RE = LARGE_NUMBER 
               RE_G = LARGE_NUMBER
               RE_kh = LARGE_NUMBER
            endif
!     
!     To select one of the following models uncomment (delete) lower
!     case c's.
!     
!---------------Begin Syamlal and O'Brien ---------------------------
!     
!     Calculate V_rm
!     
            IF(TRIM(DRAG_TYPE).EQ.'SYAM_OBRIEN') then

               IF (EP_s(IJK,M) <= ZERO) THEN 
                  F_gstmp = ZERO 
               ELSE IF (EP_G(IJK) == ZERO) THEN 
                  F_gstmp = ZERO 
               ELSE 
                  A = EP_G(IJK)**4.14D0 
                  IF (EP_G(IJK) <= 0.85D0) THEN 
                     B = drag_c1*EP_G(IJK)**1.28D0 
                  ELSE 
                     B = EP_G(IJK)**drag_d1
                  ENDIF 
                  V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) 
!------------------Begin cluster correction --------------------------
!     uncomment the following four lines and comment the above line V_RM=... 
!     V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) & 
!     * ( ONE + C(1) * exp( -a2*(Re - Re_c)**2 &
!     - a3*(EP_g(IJK)-ep_c)**2 &
!     )       * Re * (1. - EP_g(IJK))                )
!------------------End cluster correction ----------------------------
!     
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
!     
                  IF(TSUJI_DRAG) THEN
                     IF(EP_G(IJK).LE.0.8D0) THEN
                        F_gstmp = (Mu*EP_S(IJK,M)/(D_P(IJK,M)**2))*&
                        (150D0*(EP_S(IJK,M)/EP_G(IJK)) + 1.75D0*RE)
                     ELSE IF(EP_G(IJK).GT.0.8D0) THEN
                        IF(RE*EP_G(IJK).GT.1000D0) THEN
                           F_gstmp = 0.75D0*0.43D0*Mu*EP_S(IJK,M)*RE/(D_P(IJK,M)**2 *&
                           EP_G(IJK)**1.7D0)
                        ELSE IF(RE*EP_G(IJK).LE.1000D0) THEN
                           F_gstmp = 0.75D0*C_DSXRET(RE*EP_G(IJK))*Mu*EP_S(IJK,M)*&
                           RE/(D_P(IJK,M)**2 *EP_G(IJK)**1.7D0)
                        END IF
                     END IF 
                  ELSE IF(MODEL_B) THEN 
                     F_gstmp = 0.75D0*Mu*EP_S(IJK,M)*C_DSXRE(RE/V_RM)/(&
                     V_RM*D_P(IJK,M)*D_P(IJK,M)) 
                  ELSE
                     F_gstmp = 0.75D0*Mu*EP_S(IJK,M)*EP_G(IJK)*C_DSXRE(RE&
                     /V_RM)/(V_RM*D_P(IJK,M)*D_P(IJK,M)) 
                  ENDIF 
               ENDIF 
!---------------End Syamlal and O'Brien ---------------------------
!     
!--------------------------Begin Gidaspow --------------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') then
               IF(EP_g(IJK) .LE. 0.8D0) THEN
                  DgA = 150D0 * (ONE - EP_g(IJK)) * Mu &
                  / ( EP_g(IJK) * D_p(IJK,M)**2 ) &
                  + 1.75D0 * RO_g(IJK) * VREL / D_p(IJK,M)
               ELSE
                  IF(Re_G .LE. 1000D0)THEN
                     C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
                  ELSE
                     C_d = 0.44D0
                  ENDIF
                  DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                  /D_p(IJK,M)
               ENDIF
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!--------------------------End Gidaspow --------------------------
!     
!-----------------------Begin Gidaspow_blend ---------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND') then
!              Dense phase - EP_g < 0.8
               Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu &
               / ( EP_g(IJK) * D_p(IJK,M)**2 ) &
               + 1.75D0 * RO_g(IJK) * VREL / D_p(IJK,M)
!
!              Dilute phase - EP_g >= 0.8
               IF(Re_G .LE. 1000D0)THEN
                  C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
               ELSE
                  C_d = 0.44D0
               ENDIF
               WenYu = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
               /D_p(IJK,M)
!
!              Switch function
               PHI_gs = ATAN(150D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
!              Blend the models
               DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!-----------------------End Gidaspow_blend -----------------------
!     
!--------------------------Begin WEN_YU --------------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU') then
                IF(Re_G .LE. 1000D0)THEN
                   C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
                ELSE
                   C_d = 0.44D0
                ENDIF
                DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                  /D_p(IJK,M)
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
		  F_gstmp = Drag_corr/Vrsm
		!  write(*,'(I5, 8(1X,G13.6))') I, (VGC - VSCM)
               ENDIF
               
!--------------------------End WEN_YU ----------------------------
!     
!--------------------Begin Koch & Hill (2001) --------------------
!     
!!!   Added by Clay Sutton (Lehigh University) 7-14-04
!!!   
!!!   MODIFICATIONS:
!!!   
!!!   1) Declared new variables F_STOKES, F_0, F_1, F_3
!!!   
!!!   2) Added new drag closure lines
!!!
!!!  Clay's implementation was modified by Sof (01-21-2005)
!!!  for a report explaining these changes contact sof@fluent.com
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') then
!     
              F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)*EP_g(IJK)/D_p(IJK,M)**2
	       
	      phis = ONE-EP_G(IJK) ! EP_s(IJK,M) for polydisperse systems (sof --> 03-27-2007)
	      w = EXP(-10.0D0*(0.4D0-phis)/phis)
	   
	      IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
	        F_0 = (1.0D0-w) *                                           &
	              (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
	                *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*      &
	  	        phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *   &
			10.0D0*phis/(1.0D0-phis)**3
	               
	      ELSE IF(phis >= 0.4D0) THEN
	        F_0 = 10.0D0*phis/(1.0D0-phis)**3
	      ENDIF
	   
	      IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
	        F_1 = dsqrt(2.0D0/phis) / 40.0D0
	      ELSE IF(phis > 0.1D0) THEN
	        F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
	      ENDIF
	   
	      IF(phis < 0.4D0) THEN
	        F_2 = (1.0D0-w) *                                           &
	              (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
	                *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*      &
	  	        phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *  &
			10.0D0*phis/(1.0D0-phis)**3
	   
	      ELSE
	        F_2 = 10.0D0*phis/(1.0D0-phis)**3
	      ENDIF
	   
	      IF(phis < 0.0953D0) THEN
	        F_3 = 0.9351D0*phis + 0.03667D0
	      ELSE
	        F_3 = 0.0673D0 + 0.212D0*phis +0.0232D0/(1.0-phis)**5
	      ENDIF
	   
	      Re_Trans_1 = (F_2 - 1.0D0)/(3.0D0/8.0D0 - F_3)
	      Re_Trans_2 = (F_3 + dsqrt(F_3*F_3 - 4.0D0*F_1 &
	                 *(F_0-F_2))) / (2.0D0*F_1)
	   
	      IF(phis <= 0.01D0 .AND. Re_kh <= Re_Trans_1) THEN
	        F = 1.0D0 + 3.0D0/8.0D0*Re_kh
	   
	      ELSE IF(phis > 0.01D0 .AND. Re_kh <= Re_Trans_2) THEN
	        F = F_0 + F_1*Re_kh*Re_kh
	   
	  
	      ELSE IF(phis <= 0.01D0 .AND. Re_kh > Re_Trans_1 .OR.         &
	           phis >  0.01D0 .AND. Re_kh > Re_Trans_2) THEN
	        F = F_2 + F_3*Re_kh
	  
	      ELSE
	        F = zero
	      ENDIF
	   
!  This is a check for phis (or eps_(ijk,m)) to be within physical range
	      IF(phis <= ZERO .OR. phis > ONE) F = zero
	   
	      DgA = F * F_STOKES
!!!   
!!!   Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
              IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
              ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
              ENDIF
!     
!---------------------End Koch & Hill (2001) ----------------------
! 
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'BVK') then
!     
              F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)**2/D_p(IJK,M)**2 ! eq(9) BVK J. fluid. Mech. 528, 2005
	       
	      phis = ONE-EP_g(IJK)
	      D_p_av = ZERO
	      DO Im = 1, MMAX
	         D_p_av = D_p_av + EP_S(IJK,Im)/ (D_p(IJK,Im)*phis)
	      ENDDO 
	      IF(D_p_av > ZERO) D_p_av = ONE / D_p_av

	      Y_i = D_p(IJK,M)/D_p_av
	      
	      RE = D_p_av*VREL*ROP_G(IJK)/Mu
	      
	      F = 10d0 * phis / EP_g(IJK)**2 + EP_g(IJK)**2 * (ONE+1.5d0*DSQRT(phis))
	      F = F + 0.413d0*Re/(24d0*EP_g(IJK)**2) * (ONE/EP_G(IJK) + 3d0*EP_G(IJK) &
	          *phis + 8.4d0/Re**0.343) / (ONE+10**(3d0*phis)/Re**(0.5+2*phis))
	      
	      F = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) * F
	      
	      IF(Re == ZERO) F = ZERO
	      DgA = F * F_STOKES
!!!   
!!!   Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
              IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
              ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
              ENDIF
!     
!---- End Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62 (Jan 2007) -----
! 
            ELSE
              CALL START_LOG 
              IF(.not.DMP_LOG)call open_pe_log(ier)
	      if(mype == pe_io) WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              CALL END_LOG 
              call mfix_exit(myPE)  
            ENDIF
      
            F_gs(IJK, M) = (ONE - UR_F_gs) * F_gs(IJK, M) + UR_F_gs * F_gstmp
         
	 ELSE 
            F_gs(IJK, M) = ZERO 
         ENDIF 

      END DO
      
      RETURN  
      END SUBROUTINE DRAG_GS 

!//   Comments on the modifications for DMP version implementation      
!//   001 Include header file and common declarations for parallelization
!//   350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!     
