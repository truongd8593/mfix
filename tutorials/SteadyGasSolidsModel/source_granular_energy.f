!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_GRANULAR_ENERGY(sourcelhs,sourcerhs,IJK,M,IER)  C
!  Purpose: Calculate the source terms in the granular energy equation C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 04-FEB-98  C
!  Reviewer: M. Syamlal                               Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!     Local variables: sourcelhs, sourcerhs                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE run
      USE usr
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE compar        !//d
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
      INTEGER          IJK, I, J, K, M, MM
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs, shearProd
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs, collDiss
!
!                      Particle relaxation time
      DOUBLE PRECISION Tau_12_st
!
!                      Sum of eps*G_0
      DOUBLE PRECISION SUM_EpsGo
!
!                      Slip velocity
      DOUBLE PRECISION VSLIP
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
       Double Precision term63(IMAX2),term64(IMAX2),term65(IMAX2),term66(IMAX2),term67(IMAX2),term68(IMAX2), &
                        term69(IMAX2),term70(IMAX2),term71(IMAX2),term72(IMAX2),term73(IMAX2),term74(IMAX2), &
			term75(IMAX2),term76(IMAX2),term77(IMAX2),term78(IMAX2),term79(IMAX2),term80(IMAX2), &
			term81(IMAX2),term82(IMAX2),term83(IMAX2),term84(IMAX2),term85(IMAX2),term86(IMAX2), &
			term87(IMAX2),term88(IMAX2),term89(IMAX2),term90(IMAX2),term91(IMAX2),term92(IMAX2), &
			term93(IMAX2),term94(IMAX2),term95(IMAX2),term96(IMAX2),term97(IMAX2),term98(IMAX2), &
			grad_s_avg(IMAX2)
!                       Functions to compute dissipations, eq 42 of Hrenya-Sinclair paper
       DOUBLE PRECISION f1, f2, df1, df2, d2f1, d2f2, d3f1, d3f2, d4f1, x1, &
                        d4f2, d5f1, d5f2, d6f1, d6f2, d7f1, d7f2, d8f1, d8f2, d9f1, d9f2
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
      INCLUDE 's_pr1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 's_pr2.inc'
!
      I = I_OF(IJK)
      J = J_OF(IJK) 
      K = K_OF(IJK) 
!
	 x1 = ONE - Ep_g_mean(I) 
	 f1  = 5.0d0*SQRT_Pi/(96.0d0*Eta) * (ONE + 1.6d0*Eta*G_0(IJK,M,M)*x1)**2/G_0(IJK,M,M) + &
	                                1.6d0/SQRT_Pi*Eta*G_0(IJK,M,M)*x1**2
         
	 grad_s_avg(I) = ( V_s_mean(I+1) - V_s_mean(I-1) )/2d0
	 
	 df1 = (-5.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+16.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/(5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/(5.0d0*SQRT_Pi)
	 
	 d2f1 = ((-5.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+5.0d0*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1)))+16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/(5.0d0*SQRT_Pi)+32.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     (5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/(5.0d0*SQRT_Pi))/2d0
	 
	 d3f1 = ((-5.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-5.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-5.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1)))+48.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/(5.0d0*SQRT_Pi)+48.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     (5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/(5.0d0*SQRT_Pi))/6d0
	 
	 d4f1 = ((-5.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+5.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-15.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+5.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+(-5.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**3*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-5.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-5.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1)))+5.0d0*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+96.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/(5.0d0*SQRT_Pi)+64.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     (5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)/(5.0d0*SQRT_Pi))/24d0
	 
	 d5f1 = ((-5.0d0)*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/ &
     (one-x1)**6)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+25.0d0*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-25.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-75.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+25.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**5*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**6)+(-25.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(6.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+25.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-75.0d0)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+25.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**5)+(-25.0d0)*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+25.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-25.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-25.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**3*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-25.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1)))+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1)))+32*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     SQRT_Pi+16*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/SQRT_Pi+8.0d0*Eta*x1**2*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/(5.0d0*SQRT_Pi))/120d0
	 
	 d6f1 = ((-5.0d0)*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/(one-x1)**8+12600.0d0/ &
     (one-x1)**7)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+25.0d0*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+25.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-75.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-75.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-75.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+25*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+225.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+(-375.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**6)+75.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**6*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**7)+(-5.0d0)*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+25.0d0*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-75.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+(-225.0d0)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+150*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**5)-75*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**5*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**6)+(-25.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+75.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-225.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+75.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+(-25.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-25.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+25.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-5.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(48.0d0*Eta*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)/5.0d0+8.0d0*Eta*x1*(10080.0d0*x1**2/ &
     (one-x1)**9+22680.0d0*x1/(one-x1)**8+12600.0d0/(one-x1)**7)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-25.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+75.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-225.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+75.0d0*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**4*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+(-25.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+75.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-75.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-25.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+5.0d0*(8*Eta*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-75.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+75.0d0*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)**2*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-25.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+25.0d0*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+25.0d0*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)**2*SQRT_Pi/ &
     (24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+48*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/SQRT_Pi+96.0d0*Eta*x1*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)/ &
     (5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)/(5.0d0*SQRT_Pi))/720d0
	 
	 d7f1 = ((-5.0d0)*(90720.0d0*x1**2/(one-x1)**10+201600.0d0*x1/(one-x1)**9+110880.0d0/ &
     (one-x1)**8)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(96.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+35.0d0*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+35.0d0*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+175.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-105.0d0)*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-525.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+(-175.0d0)*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)**2*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-525.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+175.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+525.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+525.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**3*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**5)+(-875.0d0)*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**6)+(-2625.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**3*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**6)+1575.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**5*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+1)**2*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**7)+(-525.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**7*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)**2*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**8)+(-35.0d0)*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+35.0d0*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+175.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+175.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-525.0d0)*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-525.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-525.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**3*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+350*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**3*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+1575.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**5)+(-2625.0d0)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)**4*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**6)+525*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**6*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**7)+(-35.0d0)*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/ &
     (one-x1)**6)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+175.0d0*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+175.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/ &
     (4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-525.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+(-1575.0d0)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+525*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+(-525.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**5*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**6)+(-175.0d0)*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+175.0d0*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(6.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+175.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-525.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+175.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**4*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**5)+(-175.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+175.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-175.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+(-35.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+35.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+(-35.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(48.0d0*Eta*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)/5.0d0+8.0d0*Eta*x1*(10080.0d0*x1**2/ &
     (one-x1)**9+22680.0d0*x1/(one-x1)**8+12600.0d0/(one-x1)**7)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+5.0d0*(56.0d0*Eta*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)/5.0d0+8.0d0*Eta*x1*(90720.0d0*x1**2/ &
     (one-x1)**10+201600.0d0*x1/(one-x1)**9+110880.0d0/(one-x1)**8)/ &
     5.0d0)*(8.0d0*Eta*x1*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+1)*SQRT_Pi/(48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-35.0d0)*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+175.0d0*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+175.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-525.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**4)+(-1575.0d0)*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)**2*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/ &
     (8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**4)+525*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**3*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)**2*SQRT_Pi/(Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+(-525.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**5*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)**2*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**6)+(-175.0d0)*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+175.0d0*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**3)+525.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)**2*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-1575.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+525.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**4*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**5)+(-175.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(12.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+175.0d0*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-175.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**3*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(2.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-175.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(32.0d0*Eta*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/ &
     (one-x1)**5+42.0d0/(one-x1)**4)/5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/ &
     (one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/(one-x1)**5)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+175.0d0*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**2*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-35.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/ &
     5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+35.0d0*(48.0d0*Eta*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0+8.0d0*Eta*x1*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)/5.0d0)*(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))/5.0d0+8.0d0*Eta*x1*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/5.0d0)*SQRT_Pi/ &
     (48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+(-175.0d0)*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)**2*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+525.0d0*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**3)+(-525.0d0)*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)**3*(16.0d0*Eta*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)/ &
     5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/ &
     (one-x1)**3)/5.0d0)**2*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**4)+(-175.0d0)*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/(8.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**2)+175.0d0*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)**2*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/(4.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1))**3)+(-175.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/ &
     (one-x1)**3+2.5d0/(one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/ &
     (16.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1))**2)+35.0d0*(8*Eta*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/ &
     (one-x1)**6+240.0d0/(one-x1)**5)+8.0d0*Eta*x1*(1260.0d0*x1**2/ &
     (one-x1)**8+2880.0d0*x1/(one-x1)**7+1620.0d0/(one-x1)**6)/ &
     5.0d0)*(16.0d0*Eta*(1.5d0*x1**2/(one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/ &
     (one-x1)**2)/5.0d0+8.0d0*Eta*x1*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0)*SQRT_Pi/(16.0d0*Eta*(0.5d0*x1**2/ &
     (one-x1)**3+1.5d0*x1/(one-x1)**2+one/(one-x1)))+(-175.0d0)*(1.5d0*x1**2/ &
     (one-x1)**4+4.0d0*x1/(one-x1)**3+2.5d0/(one-x1)**2)*(24.0d0*Eta*(6.0d0*x1**2/ &
     (one-x1)**5+15.0d0*x1/(one-x1)**4+9.0d0/(one-x1)**3)/ &
     5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/(one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/ &
     (one-x1)**4)/5.0d0)**2*SQRT_Pi/(24.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/ &
     (one-x1)**2+one/(one-x1))**2)+175.0d0*(32.0d0*Eta*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/ &
     5.0d0+8.0d0*Eta*x1*(180.0d0*x1**2/(one-x1)**7+420.0d0*x1/(one-x1)**6+240.0d0/ &
     (one-x1)**5)/5.0d0)*(24.0d0*Eta*(6.0d0*x1**2/(one-x1)**5+15.0d0*x1/ &
     (one-x1)**4+9.0d0/(one-x1)**3)/5.0d0+8.0d0*Eta*x1*(30.0d0*x1**2/ &
     (one-x1)**6+72.0d0*x1/(one-x1)**5+42.0d0/(one-x1)**4)/5.0d0)*SQRT_Pi/ &
     (48.0d0*Eta*(0.5d0*x1**2/(one-x1)**3+1.5d0*x1/(one-x1)**2+one/ &
     (one-x1)))+336.0d0*Eta*(1260.0d0*x1**2/(one-x1)**8+2880.0d0*x1/ &
     (one-x1)**7+1620.0d0/(one-x1)**6)/ &
     (5.0d0*SQRT_Pi)+112.0d0*Eta*x1*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/ &
     (one-x1)**8+12600.0d0/(one-x1)**7)/ &
     (5.0d0*SQRT_Pi)+8.0d0*Eta*x1**2*(90720.0d0*x1**2/ &
     (one-x1)**10+201600.0d0*x1/(one-x1)**9+110880.0d0/(one-x1)**8)/ &
     (5.0d0*SQRT_Pi))/5040d0
	 
	 f2  = DSQRT(Th_mean(I))
	 df2 = 0.5d0/DSQRT(Th_mean(I))
	 d2f2 = -0.25d0/Th_mean(I)**1.5d0/2d0
	 d3f2 = 1.9375d0*0.375d0/Th_mean(I)**2.5d0/6d0  !modified by Euler Transormation to conv. div. series
	 d4f2 = -1.625d0*0.9375d0/Th_mean(I)**3.5d0/24d0
	 d5f2 = 0.875d0*3.28125d0/Th_mean(I)**4.5d0/120d0
	 d6f2 = -0.375d0*14.765625d0/Th_mean(I)**5.5d0/720d0
	 d7f2 = 0.0625d0*81.210938d0/Th_mean(I)**6.5d0/5040d0
	 
	 term1(I) = f1 * f2
	 term2(I) = df1 * df2
	 term3(I) = f1 * d2f2
	 term4(I) = f2 * d2f1
	 term5(I) = df1 * d2f2
	 term6(I) = df2 * d2f1
	 term7(I) = d2f1 * d2f2
	 term8(I) = f1 * d3f2
	 term9(I) = df1 * d3f2
	 term10(I) = d2f1 * d3f2
	 term11(I) = f2 * d3f1
	 term12(I) = df2 * d3f1
	 term13(I) = d2f2 * d3f1
	 term14(I) = d3f1 * d3f2
	 term15(I) = f2 * d4f1
	 term16(I) = df2 * d4f1
	 term17(I) = d2f2 * d4f1
	 term18(I) = d3f2 * d4f1
	 term19(I) = d4f2 * d4f1
	 term20(I) = f1 * d4f2
	 term21(I) = df1 * d4f2
	 term22(I) = d2f1 * d4f2
	 term23(I) = d3f1 * d4f2
	 term24(I) = f2 * d5f1
	 term25(I) = df2 * d5f1
	 term26(I) = d2f2 * d5f1
	 term27(I) = d3f2 * d5f1
	 term28(I) = d4f2 * d5f1
	 term29(I) = d5f2 * d5f1
	 term30(I) = f1 * d5f2
	 term31(I) = df1 * d5f2
	 term32(I) = d2f1 * d5f2
	 term33(I) = d3f1 * d5f2
	 term34(I) = d4f1 * d5f2
	 
	 term35(I) = f2 * d6f1
	 term36(I) = df2 * d6f1
	 term37(I) = d2f2 * d6f1
	 term38(I) = d3f2 * d6f1
	 term39(I) = d4f2 * d6f1
	 term40(I) = d5f2 * d6f1
	 term41(I) = d6f2 * d6f1
	 term42(I) = f1 * d6f2
	 term43(I) = df1 * d6f2
	 term44(I) = d2f1 * d6f2
	 term45(I) = d3f1 * d6f2
	 term46(I) = d4f1 * d6f2
	 term47(I) = d5f1 * d6f2
	 
	 term48(I) = f2   * d7f1
	 term49(I) = df2  * d7f1
	 term50(I) = d2f2 * d7f1
	 term51(I) = d3f2 * d7f1
	 term52(I) = d4f2 * d7f1
	 term53(I) = d5f2 * d7f1
	 term54(I) = d6f2 * d7f1
	 term55(I) = d7f2 * d7f1
	 term56(I) = f1   * d7f2
	 term57(I) = df1  * d7f2
	 term58(I) = d2f1 * d7f2
	 term59(I) = d3f1 * d7f2
	 term60(I) = d4f1 * d7f2
	 term61(I) = d5f1 * d7f2
	 term62(I) = d6f1 * d7f2
	 
	 term63(I) = f1 * df2
	 term64(I) = f2 * df1
!!
!!!!!! see ../turbulence_Prod.f  for details !!!!	 
	 
	 term1(I) = term1(I)*grad_s_avg(I)**2*64d0 + term1(I)*gsSqr(I)*64d0
	 
	 term2(I) = term2(I) * eps_Th(I)*grad_s_avg(I)**2*64d0 + term2(I) * eps_Thg(I)*grad_s_avg(I)*128d0 + term2(I) * eps_Thg2(I)*64d0
	 
	 term3(I) = term3(I) * Th2(I)*grad_s_avg(I)**2*64d0 + term3(I) * Th2g(I)*grad_s_avg(I)*128d0 + term3(I) * Th2g2(I)*64d0
	 term4(I) = term4(I) * eps2(I)*grad_s_avg(I)**2*64d0 + term4(I) * eps2g(I)*grad_s_avg(I)*128d0 + term4(I) * eps2g2(I)*64d0
	 term5(I) = term5(I) * eps_Th2(I)*grad_s_avg(I)**2*64d0 + term5(I) * eps_Th2g(I)*grad_s_avg(I)*128d0 + term5(I) * eps_Th2g2(I)*64d0
	 term6(I) = term6(I) * Th_eps2(I)*grad_s_avg(I)**2*64d0 + term6(I) * Th_eps2g(I)*grad_s_avg(I)*128d0 + term6(I) * Th_eps2g2(I)*64d0
	 term7(I) = term7(I) * eps2_Th2(I)*grad_s_avg(I)**2*64d0 + term7(I) * eps2_Th2g(I)*grad_s_avg(I)*128d0 + term7(I) * eps2_Th2g2(I)*64d0
	 
	 term8(I) = term8(I) * Th3(I)*grad_s_avg(I)**2*64d0 + term8(I) * Th3g(I)*grad_s_avg(I)*128d0 + term8(I) * Th3g2(I)*64d0
	 term9(I) = term9(I) * eps_Th3(I)*grad_s_avg(I)**2*64d0 + term9(I) * eps_Th3g(I)*grad_s_avg(I)*128d0 + term9(I) * eps_Th3g2(I)*64d0
	 term10(I) = term10(I) * eps2_Th3(I)*grad_s_avg(I)**2*64d0 + term10(I) * eps2_Th3g(I)*grad_s_avg(I)*128d0 + term10(I) * eps2_Th3g2(I)*64d0
	 term11(I) = term11(I) * eps3(I)*grad_s_avg(I)**2*64d0 + term11(I) * eps3g(I)*grad_s_avg(I)*128d0 + term11(I) * eps3g2(I)*64d0
	 term12(I) = term12(I) * Th_eps3(I)*grad_s_avg(I)**2*64d0 + term12(I) * Th_eps3g(I)*grad_s_avg(I)*128d0 + term12(I) * Th_eps3g2(I)*64d0
	 term13(I) = term13(I) * Th2_eps3(I)*grad_s_avg(I)**2*64d0 + term13(I) * Th2_eps3g(I)*grad_s_avg(I)*128d0 + term13(I) * Th2_eps3g2(I)*64d0
	 term14(I) = term14(I) * Th3_eps3(I)*grad_s_avg(I)**2*64d0 + term14(I) * Th3_eps3g(I)*grad_s_avg(I)*128d0 + term14(I) * Th3_eps3g2(I)*64d0
	 term15(I) = term15(I) * eps4(I)*grad_s_avg(I)**2*64d0 + term15(I) * eps4g(I)*grad_s_avg(I)*128d0 + term15(I) * eps4g2(I)*64d0
	 term16(I) = term16(I) * Th1_eps4(I)*grad_s_avg(I)**2*64d0 + term16(I) * Th1_eps4g(I)*grad_s_avg(I)*128d0 + term16(I) * Th1_eps4g2(I)*64d0
	 term17(I) = term17(I) * Th2_eps4(I)*grad_s_avg(I)**2*64d0 + term17(I) * Th2_eps4g(I)*grad_s_avg(I)*128d0 + term17(I) * Th2_eps4g2(I)*64d0
	 term18(I) = term18(I) * Th3_eps4(I)*grad_s_avg(I)**2*64d0 + term18(I) * Th3_eps4g(I)*grad_s_avg(I)*128d0 + term18(I) * Th3_eps4g2(I)*64d0
	 term19(I) = term19(I) * Th4_eps4(I)*grad_s_avg(I)**2*64d0 + term19(I) * Th4_eps4g(I)*grad_s_avg(I)*128d0 + term19(I) * Th4_eps4g2(I)*64d0
	 term20(I) = term20(I) * Th4(I)*grad_s_avg(I)**2*64d0 + term20(I) * Th4g(I)*grad_s_avg(I)*128d0 + term20(I) * Th4g2(I)*64d0
	 term21(I) = term21(I) * eps1_Th4(I)*grad_s_avg(I)**2*64d0 + term21(I) * eps1_Th4g(I)*grad_s_avg(I)*128d0 + term21(I) * eps1_Th4g2(I)*64d0
	 term22(I) = term22(I) * eps2_Th4(I)*grad_s_avg(I)**2*64d0 + term22(I) * eps2_Th4g(I)*grad_s_avg(I)*128d0 + term22(I) * eps2_Th4g2(I)*64d0
	 term23(I) = term23(I) * eps3_Th4(I)*grad_s_avg(I)**2*64d0 + term23(I) * eps3_Th4g(I)*grad_s_avg(I)*128d0 + term23(I) * eps3_Th4g2(I)*64d0
	 
	 term24(I) = term24(I) * eps5(I)*grad_s_avg(I)**2*64d0 + term24(I) * eps5g(I)*grad_s_avg(I)*128d0 + term24(I) * eps5g2(I)*64d0
	 term25(I) = term25(I) * Th1_eps5(I)*grad_s_avg(I)**2*64d0 + term25(I) * Th1_eps5g(I)*grad_s_avg(I)*128d0 + term25(I) * Th1_eps5g2(I)*64d0
	 term26(I) = term26(I) * Th2_eps5(I)*grad_s_avg(I)**2*64d0 + term26(I) * Th2_eps5g(I)*grad_s_avg(I)*128d0 + term26(I) * Th2_eps5g2(I)*64d0
	 term27(I) = term27(I) * Th3_eps5(I)*grad_s_avg(I)**2*64d0 + term27(I) * Th3_eps5g(I)*grad_s_avg(I)*128d0 + term27(I) * Th3_eps5g2(I)*64d0
	 term28(I) = term28(I) * Th4_eps5(I)*grad_s_avg(I)**2*64d0 + term28(I) * Th4_eps5g(I)*grad_s_avg(I)*128d0 + term28(I) * Th4_eps5g2(I)*64d0
	 term29(I) = term29(I) * Th5_eps5(I)*grad_s_avg(I)**2*64d0 + term29(I) * Th5_eps5g(I)*grad_s_avg(I)*128d0 + term29(I) * Th5_eps5g2(I)*64d0
	 term30(I) = term30(I) * Th5(I)*grad_s_avg(I)**2*64d0 + term30(I) * Th5g(I)*grad_s_avg(I)*128d0 + term30(I) * Th5g2(I)*64d0
	 term31(I) = term31(I) * eps1_Th5(I)*grad_s_avg(I)**2*64d0 + term31(I) * eps1_Th5g(I)*grad_s_avg(I)*128d0 + term31(I) * eps1_Th5g2(I)*64d0
	 term32(I) = term32(I) * eps2_Th5(I)*grad_s_avg(I)**2*64d0 + term32(I) * eps2_Th5g(I)*grad_s_avg(I)*128d0 + term32(I) * eps2_Th5g2(I)*64d0
	 term33(I) = term33(I) * eps3_Th5(I)*grad_s_avg(I)**2*64d0 + term33(I) * eps3_Th5g(I)*grad_s_avg(I)*128d0 + term33(I) * eps3_Th5g2(I)*64d0
	 term34(I) = term34(I) * eps4_Th5(I)*grad_s_avg(I)**2*64d0 + term34(I) * eps4_Th5g(I)*grad_s_avg(I)*128d0 + term34(I) * eps4_Th5g2(I)*64d0
	 
	 term35(I) = term35(I) * eps6(I)*grad_s_avg(I)**2*64d0 + term35(I) * eps6g(I)*grad_s_avg(I)*128d0 + term35(I) * eps6g2(I)*64d0
	 term36(I) = term36(I) * Th1_eps6(I)*grad_s_avg(I)**2*64d0 + term36(I) * Th1_eps6g(I)*grad_s_avg(I)*128d0 + term36(I) * Th1_eps6g2(I)*64d0
	 term37(I) = term37(I) * Th2_eps6(I)*grad_s_avg(I)**2*64d0 + term37(I) * Th2_eps6g(I)*grad_s_avg(I)*128d0 + term37(I) * Th2_eps6g2(I)*64d0
	 term38(I) = term38(I) * Th3_eps6(I)*grad_s_avg(I)**2*64d0 + term38(I) * Th3_eps6g(I)*grad_s_avg(I)*128d0 + term38(I) * Th3_eps6g2(I)*64d0
	 term39(I) = term39(I) * Th4_eps6(I)*grad_s_avg(I)**2*64d0 + term39(I) * Th4_eps6g(I)*grad_s_avg(I)*128d0 + term39(I) * Th4_eps6g2(I)*64d0
	 term40(I) = term40(I) * Th5_eps6(I)*grad_s_avg(I)**2*64d0 + term40(I) * Th5_eps6g(I)*grad_s_avg(I)*128d0 + term40(I) * Th5_eps6g2(I)*64d0
	 term41(I) = term41(I) * Th6_eps6(I)*grad_s_avg(I)**2*64d0 + term41(I) * Th6_eps6g(I)*grad_s_avg(I)*128d0 + term41(I) * Th6_eps6g2(I)*64d0
	 term42(I) = term42(I) * Th6(I)*grad_s_avg(I)**2*64d0 + term42(I) * Th6g(I)*grad_s_avg(I)*128d0+ term42(I) * Th6g2(I)*64d0
	 term43(I) = term43(I) * eps1_Th6(I)*grad_s_avg(I)**2*64d0 + term43(I) * eps1_Th6g(I)*grad_s_avg(I)*128d0 + term43(I) * eps1_Th6g2(I)*64d0
	 term44(I) = term44(I) * eps2_Th6(I)*grad_s_avg(I)**2*64d0 + term44(I) * eps2_Th6g(I)*grad_s_avg(I)*128d0 + term44(I) * eps2_Th6g2(I)*64d0
	 term45(I) = term45(I) * eps3_Th6(I)*grad_s_avg(I)**2*64d0 + term45(I) * eps3_Th6g(I)*grad_s_avg(I)*128d0 + term45(I) * eps3_Th6g2(I)*64d0
	 term46(I) = term46(I) * eps4_Th6(I)*grad_s_avg(I)**2*64d0 + term46(I) * eps4_Th6g(I)*grad_s_avg(I)*128d0 + term46(I) * eps4_Th6g2(I)*64d0
	 term47(I) = term47(I) * eps5_Th6(I)*grad_s_avg(I)**2*64d0 + term47(I) * eps5_Th6g(I)*grad_s_avg(I)*128d0 + term47(I) * eps5_Th6g2(I)*64d0
	 
	 term48(I) = term48(I) * eps7(I)*grad_s_avg(I)**2*64d0 + term48(I) * eps7g(I)*grad_s_avg(I)*128d0
	 term49(I) = term49(I) * Th1_eps7(I)*grad_s_avg(I)**2*64d0 + term49(I) * Th1_eps7g(I)*grad_s_avg(I)*128d0 + term49(I) * Th1_eps7g2(I)*64d0
	 term50(I) = term50(I) * Th2_eps7(I)*grad_s_avg(I)**2*64d0 + term50(I) * Th2_eps7g(I)*grad_s_avg(I)*128d0 + term50(I) * Th2_eps7g2(I)*64d0
	 term51(I) = term51(I) * Th3_eps7(I)*grad_s_avg(I)**2*64d0 + term51(I) * Th3_eps7g(I)*grad_s_avg(I)*128d0 + term51(I) * Th3_eps7g2(I)*64d0
	 term52(I) = term52(I) * Th4_eps7(I)*grad_s_avg(I)**2*64d0 + term52(I) * Th4_eps7g(I)*grad_s_avg(I)*128d0 + term52(I) * Th4_eps7g2(I)*64d0
	 term53(I) = term53(I) * Th5_eps7(I)*grad_s_avg(I)**2*64d0 + term53(I) * Th5_eps7g(I)*grad_s_avg(I)*128d0 + term53(I) * Th5_eps7g2(I)*64d0
	 term54(I) = term54(I) * Th6_eps7(I)*grad_s_avg(I)**2*64d0 + term54(I) * Th6_eps7g(I)*grad_s_avg(I)*128d0 + term54(I) * Th6_eps7g2(I)*64d0
	 term55(I) = term55(I) * Th7_eps7(I)*grad_s_avg(I)**2*64d0 + term55(I) * Th7_eps7g(I)*grad_s_avg(I)*128d0 + term55(I) * Th7_eps7g2(I)*64d0
	 term56(I) = term56(I) * Th7(I)*grad_s_avg(I)**2*64d0 + term56(I) * Th7g(I)*grad_s_avg(I)*128d0 + term56(I) * Th7g2(I)*64d0
	 term57(I) = term57(I) * eps1_Th7(I)*grad_s_avg(I)**2*64d0 + term57(I) * eps1_Th7g(I)*grad_s_avg(I)*128d0 + term57(I) * eps1_Th7g2(I)*64d0
	 term58(I) = term58(I) * eps2_Th7(I)*grad_s_avg(I)**2*64d0 + term58(I) * eps2_Th7g(I)*grad_s_avg(I)*128d0 + term58(I) * eps2_Th7g2(I)*64d0
	 term59(I) = term59(I) * eps3_Th7(I)*grad_s_avg(I)**2*64d0 + term59(I) * eps3_Th7g(I)*grad_s_avg(I)*128d0 + term59(I) * eps3_Th7g2(I)*64d0
	 term60(I) = term60(I) * eps4_Th7(I)*grad_s_avg(I)**2*64d0 + term60(I) * eps4_Th7g(I)*grad_s_avg(I)*128d0 + term60(I) * eps4_Th7g2(I)*64d0
	 term61(I) = term61(I) * eps5_Th7(I)*grad_s_avg(I)**2*64d0 + term61(I) * eps5_Th7g(I)*grad_s_avg(I)*128d0 + term61(I) * eps5_Th7g2(I)*64d0
	 term62(I) = term62(I) * eps6_Th7(I)*grad_s_avg(I)**2*64d0 + term62(I) * eps6_Th7g(I)*grad_s_avg(I)*128d0 + term62(I) * eps6_Th7g2(I)*64d0
	 
	 term63(I) = term63(I) * Thg(I)*grad_s_avg(I)*128d0 + term63(I) * Thg2(I) * 64d0
	 term64(I) = term64(I) * epsg(I)*grad_s_avg(I)*128d0 + term64(I) * epsg(I) * 64d0
	 
	 
	 shearProd = (  term1(I)+ term2(I)+ term3(I)+term4(I)+term5(I)+term6(I)+term7(I)+ &
	                term8(I)+term9(I)+term10(I)+term11(I)+term12(I)+term13(I)+term14(I)+ &
			term15(I)+term16(I)+term17(I)+term18(I)+term19(I)+term20(I)+term21(I)+term22(I)+term23(I)+ &
			term24(I)+term25(I)+term26(I)+term27(I)+term28(I)+term29(I)+term30(I)+term31(I)+term32(I)+ &
			term33(I)+term34(I)+ &
			term35(I)+term36(I)+term37(I)+term38(I)+term39(I)+term40(I)+term41(I)+term42(I)+term43(I)+ &
			term44(I)+term45(I)+term46(I)+term47(I)+ &
			term48(I)+term49(I)+term50(I)+term51(I)+term52(I)+term53(I)+term54(I)+term55(I)+term56(I)+ &
			term57(I)+term58(I)+term59(I)+term60(I)+term61(I)+term62(I)+term63(I)+term64(I)  )
	 
	 shearProd = shearProd * RO_s(M)*D_P(IJK,M)
	 
! this is what was in the transient code	 
!	 shearProd =  (5.0d0*SQRT_Pi/(96.0d0*Eta) * (ONE + 1.6d0*Eta*G_0(IJK,M,M)*ep_s(ijk,m))**2/G_0(IJK,M,M) + &
!	                                1.6d0/SQRT_Pi*Eta*G_0(IJK,M,M)*ep_s(ijk,m)**2) * DSQRT(theta_m(ijk,m)) * TRD_S2(IJK,M)
	 
	 !2d0*MU_S_C(IJK,M)* TRD_S2(IJK,M)
      
      
      
      
      
      SOURCERHS = shearProd * VOL(IJK)  
!
!
!---------------------------------------------------------------------
!  The following lines are commented out since Kphi_s has been set to zero
!  in subroutine calc_mu_s.  To activate the feature uncomment the following
!  lines and the lines in calc_mu_s.
!
!
!      sourcerhs=sourcerhs+
!     &    AVG_X_h(Kphi_s(IJK,M),Kphi_s(EAST_OF(IJK),M),I)
!     &    *(EP_s(EAST_OF(IJK),M)-EP_s(IJK,M))*oDX_E(I)
!     &    *AYZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(WEST_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,I_OF(WEST_OF(IJK)))*(EP_s(IJK,M)-EP_s(WEST_OF(IJK),M))
!     &    *oDX_E(I_OF(WEST_OF(IJK)))*AYZ(WEST_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Y_h(Kphi_s(IJK,M),Kphi_s(NORTH_OF(IJK),M),J)
!     &    *(EP_s(NORTH_OF(IJK),M)-EP_s(IJK,M))*oDY_N(J)
!     &    *AXZ(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_X_h(Kphi_s(SOUTH_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,J_OF(SOUTH_OF(IJK)))*(EP_s(IJK,M)-EP_s(SOUTH_OF(IJK),M))
!     &    *oDY_N(J_OF(SOUTH_OF(IJK)))
!     &    *AXZ(JM_OF(IJK))
!
!      sourcerhs=sourcerhs+
!     &    AVG_Z_h(Kphi_s(IJK,M),Kphi_s(TOP_OF(IJK),M),K)
!     &    *(EP_s(TOP_OF(IJK),M)-EP_s(IJK,M))
!     &    *oX(I)*oDZ_T(K)
!     &    *AXY(IJK)
!
!      sourcerhs=sourcerhs-
!     &    AVG_Z_h(Kphi_s(BOTTOM_OF(IJK),M),Kphi_s(IJK,M)
!     &    ,K_OF(BOTTOM_OF(IJK)))*(EP_s(IJK,M)-
!     &    EP_s(BOTTOM_OF(IJK),M))
!     &    *oX(I)*oDZ_T(K_OF(BOTTOM_OF(IJK)))
!     &    *AXY(KM_OF(IJK))
!---------------------------------------------------------------------
!
! Changes needed for multitype particles, sof June 16 2005
! Sum of eps*G_0 is used instead of Eps*G_0
!
      SUM_EpsGo = ZERO
      DO MM = 1, MMAX
        SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,MM,MM)
      ENDDO
      
	 x1 = ONE - Ep_g_mean(I)
         f1  = x1**2 * G_0(IJK,M,M)
	 
	 df1 = 2d0*x1*(0.5d0*x1**2/(1d0-x1)**3+1.5d0*x1/(1d0-x1)**2+1d0/(1d0-x1)) + &
	       x1**2*(1.5d0*x1**2/(1d0-x1)**4+4.0d0*x1/(1d0-x1)**3+2.5d0/(1d0-x1)**2)
	 d2f1 = (2d0*(0.5d0*x1**2/(1d0-x1)**3+1.5d0*x1/(1d0-x1)**2+1d0/(1d0-x1))+4d0*x1*(1.5d0*x1**2/(1d0-x1) &
	        **4+4.0d0*x1/(1d0-x1)**3+2.5d0/(1d0-x1)**2)+x1**2*(6.0d0*x1**2/(1d0-x1)**5+15.0d0*x1/(    &
		1d0-x1)**4+9.0d0/(1d0-x1)**3))/2d0
	 d3f1 = (6d0*(1.5d0*x1**2/(1d0-x1)**4+4.0d0*x1/(1d0-x1)**3+2.5d0/(1d0-x1)**2)+6d0*x1*(6.0d0*x1**2/(1d0 &
	        -x1)**5+15.0*x1/(1d0-x1)**4+9.0d0/(1d0-x1)**3)+x1**2*(30.0d0*x1**2/(1d0-x1)**6+72.0d0  &
		*x1/(1d0-x1)**5+42.0d0/(1d0-x1)**4))/6d0
	 d4f1 = (12d0*(6.0d0*x1**2/(1d0-x1)**5+15.0d0*x1/(1d0-x1)**4+9.0d0/(1d0-x1)**3)+8d0*x1*(30.0d0*x1**2 &
	        /(1d0-x1)**6+72.0d0*x1/(1d0-x1)**5+42.0d0/(1d0-x1)**4)+x1**2*(180.0d0*x1**2/(1d0-x1)**7 &
		+420.0d0*x1/(1d0-x1)**6+240.0d0/(1d0-x1)**5))/24d0
	 d5f1 = (20d0*(30.0d0*x1**2/(1d0-x1)**6+72.0d0*x1/(1d0-x1)**5+42.0d0/(1d0-x1)**4)+10d0*x1*(180.0d0* &
	        x1**2/(1d0-x1)**7+420.0d0*x1/(1d0-x1)**6+240.0d0/(1d0-x1)**5)+x1**2*(1260.0d0*x1**2 &
		/(1d0-x1)**8+2880.0d0*x1/(1d0-x1)**7+1620.0d0/(1d0-x1)**6))/120d0
	 d6f1 = (30d0*(180.0d0*x1**2/(1d0-x1)**7+420.0d0*x1/(1d0-x1)**6+240.0d0/(1d0-x1)**5)+12d0*x1*(1260.0d0 &
               *x1**2/(1d0-x1)**8+2880.0d0*x1/(1d0-x1)**7+1620.0d0/(1d0-x1)**6)+x1**2*(10080.0d0 &
               *x1**2/(1d0-x1)**9+22680.0d0*x1/(1d0-x1)**8+12600.0d0/(1d0-x1)**7))/720d0
	 d7f1 = (42d0*(1260.0d0*x1**2/(1d0-x1)**8+2880.0d0*x1/(1d0-x1)**7+1620.0d0/(1d0-x1)**6)+14d0*x1*( &
                10080.0d0*x1**2/(1d0-x1)**9+22680.0d0*x1/(1d0-x1)**8+12600.0d0/(1d0-x1)**7)+x1**2 &
               *(90720.0d0*x1**2/(1d0-x1)**10+201600.0d0*x1/(1d0-x1)**9+110880.0d0/(1d0-x1)**8))/5040d0
	 d8f1 = (56d0*(10080.0d0*x1**2/(one-x1)**9+22680.0d0*x1/(one-x1)**8+12600.0d0/ &
     (one-x1)**7)+16d0*x1*(90720.0d0*x1**2/(one-x1)**10+201600.0d0*x1/ &
     (one-x1)**9+110880.0d0/(one-x1)**8)+x1**2*(907200.0d0*x1**2/ &
     (one-x1)**11+1995840.0d0*x1/(one-x1)**10+1088640.0d0/(one-x1)**9))/40320d0
	 d9f1 = (72d0*(90720.0d0*x1**2/(one-x1)**10+201600.0d0*x1/(one-x1)**9+110880.0d0/ &
     (one-x1)**8)+18d0*x1*(907200.0d0*x1**2/(one-x1)**11+1995840.0d0*x1/ &
     (one-x1)**10+1088640.0d0/(one-x1)**9)+x1**2*(9979200.0d0*x1**2/ &
     (one-x1)**12+2.17728D+7*x1/(one-x1)**11+1.17936D+7/(one-x1)**10))/362880d0
	 
	 f2  = Th_mean(I)**1.5
	 df2 = 1.5d0*DSQRT(Th_mean(I))
	 d2f2 = 0.75d0 / DSQRT(Th_mean(I))/2d0 ! this !2 (factorial of 2)
	 d3f2 = -0.375d0/Th_mean(I)**1.5/6d0 !   this !3 (1.+1./2.+1./4.+1./8.+1./16.)*
	 d4f2 = 0.5625d0/Th_mean(I)**2.5/24d0 !  this !4 (1./2.+2./4.+3./8.+4./16.)*
	 d5f2 = -(1.+1./2.+1./4.+1./8.+1./16.)*1.40625d0/Th_mean(I)**3.5/120d0 !  this !5 (1./4.+3./8.+4./16.)*
	 d6f2 = (1./2.+2./4.+3./8.+4./16.)*4.921875d0/Th_mean(I)**4.5/720d0 !  this !6 (1./8.+4./16.)*
	 d7f2 = -(1./4.+3./8.+4./16.)*22.1484375d0/Th_mean(I)**5.5/5040d0 !  this !7 (1./16.)*
	 d8f2 = (1./8.+4./16.)*121.81641d0/Th_mean(I)**6.5/40320d0
	 d9f2 = -(1./16.)*791.80664d0/Th_mean(I)**7.5/362880d0
	 
	 term1(I) = f1 * f2
	 
	 term2(I) = df1 * df2
	 
	 term3(I) = f1 * d2f2
	 term4(I) = f2 * d2f1
	 term5(I) = df1 * d2f2
	 term6(I) = df2 * d2f1
	 term7(I) = d2f1 * d2f2
	 
	 term8(I) = f1 * d3f2
	 term9(I) = df1 * d3f2
	 term10(I) = d2f1 * d3f2
	 term11(I) = f2 * d3f1
	 term12(I) = df2 * d3f1
	 term13(I) = d2f2 * d3f1
	 term14(I) = d3f1 * d3f2
	 
	 term15(I) = f2 * d4f1
	 term16(I) = df2 * d4f1
	 term17(I) = d2f2 * d4f1
	 term18(I) = d3f2 * d4f1
	 term19(I) = d4f2 * d4f1
	 term20(I) = f1 * d4f2
	 term21(I) = df1 * d4f2
	 term22(I) = d2f1 * d4f2
	 term23(I) = d3f1 * d4f2
	 
	 term24(I) = f2 * d5f1
	 term25(I) = df2 * d5f1
	 term26(I) = d2f2 * d5f1
	 term27(I) = d3f2 * d5f1
	 term28(I) = d4f2 * d5f1
	 term29(I) = d5f2 * d5f1
	 term30(I) = f1 * d5f2
	 term31(I) = df1 * d5f2
	 term32(I) = d2f1 * d5f2
	 term33(I) = d3f1 * d5f2
	 term34(I) = d4f1 * d5f2
	 
	 term35(I) = f2 * d6f1
	 term36(I) = df2 * d6f1
	 term37(I) = d2f2 * d6f1
	 term38(I) = d3f2 * d6f1
	 term39(I) = d4f2 * d6f1
	 term40(I) = d5f2 * d6f1
	 term41(I) = d6f2 * d6f1
	 term42(I) = f1 * d6f2
	 term43(I) = df1 * d6f2
	 term44(I) = d2f1 * d6f2
	 term45(I) = d3f1 * d6f2
	 term46(I) = d4f1 * d6f2
	 term47(I) = d5f1 * d6f2
	 
	 term48(I) = f2   * d7f1
	 term49(I) = df2  * d7f1
	 term50(I) = d2f2 * d7f1
	 term51(I) = d3f2 * d7f1
	 term52(I) = d4f2 * d7f1
	 term53(I) = d5f2 * d7f1
	 term54(I) = d6f2 * d7f1
	 term55(I) = d7f2 * d7f1
	 term56(I) = f1   * d7f2
	 term57(I) = df1  * d7f2
	 term58(I) = d2f1 * d7f2
	 term59(I) = d3f1 * d7f2
	 term60(I) = d4f1 * d7f2
	 term61(I) = d5f1 * d7f2
	 term62(I) = d6f1 * d7f2
	 
	 term63(I) = f2   * d8f1
	 term64(I) = df2  * d8f1
	 term65(I) = d2f2 * d8f1
	 term66(I) = d3f2 * d8f1
	 term67(I) = d4f2 * d8f1
	 term68(I) = d5f2 * d8f1
	 term69(I) = d6f2 * d8f1
	 term70(I) = d7f2 * d8f1
	 term71(I) = d8f2 * d8f1
	 term72(I) = f1   * d8f2
	 term73(I) = df1  * d8f2
	 term74(I) = d2f1 * d8f2
	 term75(I) = d3f1 * d8f2
	 term76(I) = d4f1 * d8f2
	 term77(I) = d5f1 * d8f2
	 term78(I) = d6f1 * d8f2
	 term79(I) = d7f1 * d8f2
	 
	 term80(I) = f2   * d9f1
	 term81(I) = df2  * d9f1
	 term82(I) = d2f2 * d9f1
	 term83(I) = d3f2 * d9f1
	 term84(I) = d4f2 * d9f1
	 term85(I) = d5f2 * d9f1
	 term86(I) = d6f2 * d9f1
	 term87(I) = d7f2 * d9f1
	 term88(I) = d8f2 * d9f1
	 term89(I) = d9f2 * d9f1
	 term90(I) = f1   * d9f2
	 term91(I) = df1  * d9f2
	 term92(I) = d2f1 * d9f2
	 term93(I) = d3f1 * d9f2
	 term94(I) = d4f1 * d9f2
	 term95(I) = d5f1 * d9f2
	 term96(I) = d6f1 * d9f2
	 term97(I) = d7f1 * d9f2
	 term98(I) = d8f1 * d9f2
	 
	 term2(I) = term2(I) * eps_Th(I)
	 
	 term3(I) = term3(I) * Th2(I)
	 term4(I) = term4(I) * eps2(I)
	 term5(I) = term5(I) * eps_Th2(I)
	 term6(I) = term6(I) * Th_eps2(I)
	 term7(I) = term7(I) * eps2_Th2(I)
	 
	 term8(I) = term8(I) * Th3(I)
	 term9(I) = term9(I) * eps_Th3(I)
	 term10(I) = term10(I) * eps2_Th3(I)
	 term11(I) = term11(I) * eps3(I)
	 term12(I) = term12(I) * Th_eps3(I)
	 term13(I) = term13(I) * Th2_eps3(I)
	 term14(I) = term14(I) * Th3_eps3(I)
	 term15(I) = term15(I) * eps4(I)
	 term16(I) = term16(I) * Th1_eps4(I)
	 term17(I) = term17(I) * Th2_eps4(I)
	 term18(I) = term18(I) * Th3_eps4(I)
	 term19(I) = term19(I) * Th4_eps4(I)
	 term20(I) = term20(I) * Th4(I)
	 term21(I) = term21(I) * eps1_Th4(I)
	 term22(I) = term22(I) * eps2_Th4(I)
	 term23(I) = term23(I) * eps3_Th4(I)
	 
	 term24(I) = term24(I) * eps5(I)
	 term25(I) = term25(I) * Th1_eps5(I)
	 term26(I) = term26(I) * Th2_eps5(I)
	 term27(I) = term27(I) * Th3_eps5(I)
	 term28(I) = term28(I) * Th4_eps5(I)
	 term29(I) = term29(I) * Th5_eps5(I)
	 term30(I) = term30(I) * Th5(I)
	 term31(I) = term31(I) * eps1_Th5(I)
	 term32(I) = term32(I) * eps2_Th5(I)
	 term33(I) = term33(I) * eps3_Th5(I)
	 term34(I) = term34(I) * eps4_Th5(I)
	 
	 term35(I) = term35(I) * eps6(I)
	 term36(I) = term36(I) * Th1_eps6(I)
	 term37(I) = term37(I) * Th2_eps6(I)
	 term38(I) = term38(I) * Th3_eps6(I)
	 term39(I) = term39(I) * Th4_eps6(I)
	 term40(I) = term40(I) * Th5_eps6(I)
	 term41(I) = term41(I) * Th6_eps6(I)
	 term42(I) = term42(I) * Th6(I)
	 term43(I) = term43(I) * eps1_Th6(I)
	 term44(I) = term44(I) * eps2_Th6(I)
	 term45(I) = term45(I) * eps3_Th6(I)
	 term46(I) = term46(I) * eps4_Th6(I)
	 term47(I) = term47(I) * eps5_Th6(I)
	 
	 term48(I) = term48(I) * eps7(I)
	 term49(I) = term49(I) * Th1_eps7(I)
	 term50(I) = term50(I) * Th2_eps7(I)
	 term51(I) = term51(I) * Th3_eps7(I)
	 term52(I) = term52(I) * Th4_eps7(I)
	 term53(I) = term53(I) * Th5_eps7(I)
	 term54(I) = term54(I) * Th6_eps7(I)
	 term55(I) = term55(I) * Th7_eps7(I)
	 term56(I) = term56(I) * Th7(I)
	 term57(I) = term57(I) * eps1_Th7(I)
	 term58(I) = term58(I) * eps2_Th7(I)
	 term59(I) = term59(I) * eps3_Th7(I)
	 term60(I) = term60(I) * eps4_Th7(I)
	 term61(I) = term61(I) * eps5_Th7(I)
	 term62(I) = term62(I) * eps6_Th7(I)
	 
	 term63(I) = term63(I) * eps8(I)
	 term64(I) = term64(I) * Th1_eps8(I)
	 term65(I) = term65(I) * Th2_eps8(I)
	 term66(I) = term66(I) * Th3_eps8(I)
	 term67(I) = term67(I) * Th4_eps8(I)
	 term68(I) = term68(I) * Th5_eps8(I)
	 term69(I) = term69(I) * Th6_eps8(I)
	 term70(I) = term70(I) * Th7_eps8(I)
	 term71(I) = term71(I) * Th8_eps8(I)
	 term72(I) = term72(I) * Th8(I)
	 term73(I) = term73(I) * eps1_Th8(I)
	 term74(I) = term74(I) * eps2_Th8(I)
	 term75(I) = term75(I) * eps3_Th8(I)
	 term76(I) = term76(I) * eps4_Th8(I)
	 term77(I) = term77(I) * eps5_Th8(I)
	 term78(I) = term78(I) * eps6_Th8(I)
	 term79(I) = term79(I) * eps7_Th8(I)
	 
	 term80(I) = term80(I) * eps9(I)
	 term81(I) = term81(I) * Th1_eps9(I)
	 term82(I) = term82(I) * Th2_eps9(I)
	 term83(I) = term83(I) * Th3_eps9(I)
	 term84(I) = term84(I) * Th4_eps9(I)
	 term85(I) = term85(I) * Th5_eps9(I)
	 term86(I) = term86(I) * Th6_eps9(I)
	 term87(I) = term87(I) * Th7_eps9(I)
	 term88(I) = term88(I) * Th8_eps9(I)
	 term89(I) = term89(I) * Th9_eps9(I)
	 term90(I) = term90(I) * Th9(I)
	 term91(I) = term91(I) * eps1_Th9(I)
	 term92(I) = term92(I) * eps2_Th9(I)
	 term93(I) = term93(I) * eps3_Th9(I)
	 term94(I) = term94(I) * eps4_Th9(I)
	 term95(I) = term95(I) * eps5_Th9(I)
	 term96(I) = term96(I) * eps6_Th9(I)
	 term97(I) = term97(I) * eps7_Th9(I)
	 term98(I) = term98(I) * eps8_Th9(I)
      
      collDiss = term1(I)+term2(I)+term3(I)+term4(I)+term5(I)+term6(I)+term7(I)+ &
	         term8(I)+term9(I)+term10(I)+term11(I)+term12(I)+term13(I)+term14(I)+ &
		 term15(I)+term16(I)+term17(I)+term18(I)+term19(I)+term20(I)+term21(I)+term22(I)+term23(I)+ &
		 term24(I)+term25(I)+term26(I)+term27(I)+term28(I)+term29(I)+term30(I)+term31(I)+term32(I)+ &
		 term33(I)+term34(I)+ &
		 term35(I)+term36(I)+term37(I)+term38(I)+term39(I)+term40(I)+term41(I)+term42(I)+term43(I)+ &
		 term44(I)+term45(I)+term46(I)+term47(I)+ &
		 term48(I)+term49(I)+term50(I)+term51(I)+term52(I)+term53(I)+term54(I)+term55(I)+term56(I)+ &
		 term57(I)+term58(I)+term59(I)+term60(I)+term61(I)+term62(I)+ &
		 term63(I)+term64(I)+term65(I)+term66(I)+term67(I)+term68(I)+term69(I)+term70(I)+term71(I)+ &
	         term72(I)+term73(I)+term74(I)+term75(I)+term76(I)+term77(I)+term78(I)+term79(I)+ &
		 term80(I)+term81(I)+term82(I)+term83(I)+term84(I)+term85(I)+term86(I)+term87(I)+term88(I)+ &
		 term89(I)+term90(I)+term91(I)+term92(I)+term93(I)+term94(I)+term95(I)+term96(I)+term97(I)+ &
		 term98(I)
      collDiss = collDiss*(48d0/DSQRT(PI))*ETA*(ONE-ETA)*RO_s(M)/D_P(IJK,M)/Th_mean(I)
   ! steady state expression was the following: 
!      collDiss = (48d0/DSQRT(PI))*ETA*(ONE-ETA)*ROP_S(IJK,M)*SUM_EpsGo*DSQRT(THETA_M(IJK,M))/D_P(IJK,M)
     
      SOURCELHS = collDiss * VOL(IJK)
!      write(*,'(I5, 8(1X,G13.6))') I, collDiss*THETA_M(IJK,M), shearProd
!
      RETURN  
      END SUBROUTINE SOURCE_GRANULAR_ENERGY 
!-----------------------------------------------  
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: SOURCE_IA_NONEP_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER)
!  Purpose: This module computes the source term of the species
!           granular energy equation
!
!  Literature/Document References:    
!	Iddir, Y.H., "Modeling of the multiphase mixture of particles 
!	  using the kinetic theory approach," PhD Thesis, Illinois
!	  Institute of Technology, Chicago, Illinois, 2004
!    Iddir, Y.H., & H. Arastoopour, "Modeling of multitype particle
!      flow using the kinetic theory approach," AIChE J., Vol 51,
!      No 6, June 2005
!                                                         
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      SUBROUTINE SOURCE_IA_NONEP_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK, M, IER) 
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE residual
      use kintheory
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!                      Error index
      INTEGER          IER
!     
!                      Indices
      INTEGER          IJK, I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                       IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
!
!                      phase index 
      INTEGER          M, L, LM
!
!                      velocities
      DOUBLE PRECISION UsM_e, UsM_w, VsM_n, VsM_s, WsM_t, WsM_b,&
                       UsL_e, UsL_w, VsL_n, VsL_s, WsL_t, WsL_b,&
                       UsM_p, VsM_p, WsM_P, UsL_p, VsL_p, WsL_p
!
!                      number densities
      DOUBLE PRECISION NU_PM_E, NU_PM_W, NU_PM_N, NU_PM_S, NU_PM_T,&
                       NU_PM_B, NU_PM_p,&
                       NU_PL_E, NU_PL_W, NU_PL_N, NU_PL_S, NU_PL_T,&
                       NU_PL_B, NU_PL_p
!
!                      temperature of species L
      DOUBLE PRECISION T_PL_E, T_PL_W, T_PL_N, T_PL_S, T_PL_T,&
                       T_PL_B, T_PL_p
!
!                      particle characteristics
      DOUBLE PRECISION M_PM, M_PL, D_PM, D_PL
!
!                      coefficients in heat flux term
      DOUBLE PRECISION Knu_sL_e, Knu_sL_w, Knu_sL_n, Knu_sL_s, Knu_sL_t,&
                       Knu_sL_b, Knu_sM_e, Knu_sM_w, Knu_sM_n, Knu_sM_s,&
                       Knu_sM_t, Knu_sM_b,&
                       Kvel_s_e, Kvel_s_w, Kvel_s_n, Kvel_s_s, Kvel_s_t,&
                       Kvel_s_b,&
                       Kth_sL_e, Kth_sL_w, Kth_sL_n, Kth_sL_s, Kth_sL_t,&
                       Kth_sL_b
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs, S10_rhs, S15_rhs, S16_rhs,&
                       S11a_sum_rhs, S11b_sum_rhs, S11c_sum_rhs, S12a_sum_rhs,&
                       S12b_sum_rhs, S12c_sum_rhs, &
                       S13_sum_lhs, S13_sum_rhs, &
                       S14a_sum_rhs, S14b_sum_rhs, S14c_sum_rhs,&
                       S17_sum_rhs, S18_sum_rhs, &
                       S9_sum_rhs, s21a_sum_rhs, S21b_sum_rhs, S21c_sum_rhs
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs, S10_lhs, S16_lhs,&
                       S11a_sum_lhs, S11b_sum_lhs, S11c_sum_lhs, S12a_sum_lhs,&
                       S12b_sum_lhs, S12c_sum_lhs, &
                       S17_sum_lhs, S18_sum_lhs, S20_sum_lhs

!
!----------------------------------------------- 
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------
!
      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK) 
      IJKE = EAST_OF(IJK) 
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK) 
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK) 
      IJKB = BOTTOM_OF(IJK) 
!
!     initialize summation variables
      S9_sum_rhs = ZERO
      S11a_sum_rhs = ZERO
      S12a_sum_rhs = ZERO 
      S13_sum_rhs = ZERO  
      S14a_sum_rhs = ZERO
      S14b_sum_rhs = ZERO
      S14c_sum_rhs = ZERO   
      S17_sum_rhs = ZERO
      S18_sum_rhs = ZERO
      S21a_sum_rhs = ZERO
      S21b_sum_rhs = ZERO
      S21c_sum_rhs = ZERO

      S11a_sum_lhs = ZERO
      S11b_sum_lhs = ZERO
      S11c_sum_lhs = ZERO
      S12a_sum_lhs = ZERO
      S12b_sum_lhs = ZERO
      S12c_sum_lhs = ZERO
      S13_sum_lhs = ZERO  
      S17_sum_lhs = ZERO
      S18_sum_lhs = ZERO
      S20_sum_lhs = ZERO
!
      UsM_e = U_S(IJK,M)
      UsM_w = U_S(IMJK,M)
      VsM_n = V_S(IJK,M)
      VsM_s = V_S(IJMK,M)
      WsM_t = W_S(IJK,M)
      WsM_b = W_S(IJKM,M)
      UsM_p = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VsM_p = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M) )
      WsM_p = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M) )
! 
      D_PM = D_P(IJK,M) 
      M_PM = (Pi/6.d0)*D_PM**3 * RO_S(M)
      NU_PM_p = ROP_S(IJK,M)/M_PM
      NU_PM_E = ROP_S(IJKE,M)/M_PM
      NU_PM_W = ROP_S(IJKW,M)/M_PM
      NU_PM_N = ROP_S(IJKN,M)/M_PM
      NU_PM_S = ROP_S(IJKS,M)/M_PM
      NU_PM_T = ROP_S(IJKT,M)/M_PM
      NU_PM_B = ROP_S(IJKB,M)/M_PM
!
!
!     Production by shear: (S:grad(vi))
!         Pi_s*tr(Di)

      S10_lhs = P_S_C(IJK,M) * ZMAX(TRD_S_C(IJK,M)) 
      S10_rhs = P_S_C(IJK,M) * ZMAX(-TRD_S_C(IJK,M))
!
!
!     Production by shear: (S:grad(vi))  
!         Mu_s*tr(Di^2)
      S15_rhs = 2.d0*Mu_s_c(IJK,M)*TRD_S2(IJK,M)
!
!
!     Production by shear: (S:grad(vi))  
!         Lambda_s*tr(Di)^2
      S16_lhs = (TRD_S_C(IJK,M)**2)*ZMAX( -LAMBDA_s_c(IJK,M) )
      S16_rhs = (TRD_S_C(IJK,M)**2)*ZMAX(  LAMBDA_s_C(IJK,M) )
	 
!
!
      DO L = 1, MMAX
          D_PL = D_P(IJK,L) 
          M_PL = (Pi/6.d0)*D_PL**3 * RO_S(L)
          NU_PL_p = ROP_S(IJK,L)/M_PL
          NU_PL_E = ROP_S(IJKE,L)/M_PL
          NU_PL_W = ROP_S(IJKW,L)/M_PL
          NU_PL_N = ROP_S(IJKN,L)/M_PL
          NU_PL_S = ROP_S(IJKS,L)/M_PL   
          NU_PL_T = ROP_S(IJKT,L)/M_PL
          NU_PL_B = ROP_S(IJKB,L)/M_PL

          T_PL_p = Theta_m(IJK,L)
          T_PL_E = Theta_m(IJKE,L)
          T_PL_W = Theta_m(IJKW,L)
          T_PL_N = Theta_m(IJKN,L)
          T_PL_S = Theta_m(IJKS,L)
          T_PL_T = Theta_m(IJKT,L)
          T_PL_B = Theta_m(IJKB,L)

          UsL_e = U_S(IJK,L)
          UsL_w = U_S(IMJK,L)
          VsL_n = V_S(IJK,L)
          VsL_s = V_S(IJMK,L)
          WsL_t = W_S(IJK,L)
          WsL_b = W_S(IJKM,L)
          UsL_p = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I)
          VsL_p = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L) )
          WsL_p = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L) )
!
!
!         Energy dissipation by collisions: Sum(Nip)
!              SUM( EDT_s_ip )
          S20_sum_lhs = S20_sum_lhs + EDT_s_ip(IJK,M,L)
!
!
!         Energy dissipation by collisions: SUM(Nip)
!              SUM( EDvel_sL_ip* div(vp) ) !Modified by sof to include trace of V_s_L
          
	  S11a_sum_lhs = S11a_sum_lhs + ZMAX(-EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)
          
	  S11a_sum_rhs = S11a_sum_rhs + ZMAX( EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)
!
!
!         Energy dissipation by collisions: Sum(Nip)
!              SUM( EDvel_sM_ip* div(vi) ) !Modified by sof to include trace of V_s_M
          
	  S12a_sum_lhs = S12a_sum_lhs + ZMAX(-EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)
          
	  S12a_sum_rhs = S12a_sum_rhs + ZMAX( EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)
!
!                    
          IF (M .NE. L) THEN
	        LM = FUNLM(L,M)
!
!
!              Production by shear: (S:grad(vi))  
!                   SUM(2*Mu_sL_ip*tr(Dk*Di) )
               S17_sum_lhs = S17_sum_lhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( - TRD_s2_ip(IJK,M,L) )
               S17_sum_rhs = S17_sum_rhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( TRD_s2_ip(IJK,M,L) )
!
! These two terms can be treated explicitly here by uncommenting the following
! two lines. They are currently treated by PEA algorithm in solve_granular_energy
!
!               S13_sum_lhs = S13_sum_lhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,M)
!               S13_sum_rhs = S13_sum_rhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,L)
!
!
!              Production by shear: (S:grad(vi))  
!                   SUM( (Xi_sL_ip-(2/3)*Mu_sL_ip)*tr(Dk)tr(Di) )
               S18_sum_lhs = S18_sum_lhs + ZMAX(-1.d0* &
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )
               S18_sum_rhs = S18_sum_rhs + ZMAX(&
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )
!
!
!              Part of Heat Flux: div (q)
!                   Kth_sL_ip*[grad(Tp)]
!              Note for L=M S21 terms cancel with similar term arising from
!                   grad(Ti)
               Kth_sL_e = AVG_X_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKE,M,L),I)
               Kth_sL_w = AVG_X_S(Kth_sL_ip(IJKW,M,L),Kth_sL_ip(IJK,M,L), IM)
               Kth_sL_n = AVG_Y_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKN,M,L),J)
               Kth_sL_s = AVG_Y_S(Kth_sL_ip(IJKS,M,L),Kth_sL_ip(IJK,M,L), JM)
               Kth_sL_t = AVG_Z_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKT,M,L),K)
               Kth_sL_b = AVG_Z_S(Kth_sL_ip(IJKB,M,L),Kth_sL_ip(IJK,M,L), KM)

               S21a_sum_rhs = S21a_sum_rhs + ( (Kth_sL_e*(T_PL_E-T_PL_p) )*&
                    ODX_E(I)*AYZ(IJK) - (Kth_sL_w*(T_PL_p-T_PL_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S21b_sum_rhs = S21b_sum_rhs + ( (Kth_sL_n*(T_PL_N-T_PL_p) )*&
                    ODY_N(J)*AXZ(IJK) - (Kth_sL_s*(T_PL_p-T_PL_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S21c_sum_rhs = S21c_sum_rhs + ( (Kth_sL_t*(T_PL_T-T_PL_p) )*&
                    ODZ_T(K)*OX(I)*AXY(IJK) - (Kth_sL_b*(T_PL_p-T_PL_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )
!
!
!              Part of Heat Flux: div (q)
!                   Knu_s_ip*[ni*grad(np)-np*grad(ni)]
!              Note S14 terms should evaluate to zero for particles from the same phase
               Knu_sL_e = AVG_X_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKE,M,L),I)
               Knu_sM_e = AVG_X_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKE,M,L),I)
               Knu_sL_w = AVG_X_S(Knu_sL_ip(IJKW,M,L),Knu_sL_ip(IJK,M,L), IM)
               Knu_sM_w = AVG_X_S(Knu_sM_ip(IJKW,M,L),Knu_sM_ip(IJK,M,L), IM)
               Knu_sL_n = AVG_Y_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKN,M,L),J)
               Knu_sM_n = AVG_Y_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKN,M,L),J)
               Knu_sL_s = AVG_Y_S(Knu_sL_ip(IJKS,M,L),Knu_sL_ip(IJK,M,L), JM)
               Knu_sM_s = AVG_Y_S(Knu_sM_ip(IJKS,M,L),Knu_sM_ip(IJK,M,L), JM)
               Knu_sL_t = AVG_Z_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKT,M,L),K)
               Knu_sM_t = AVG_Z_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKT,M,L),K)
               Knu_sL_b = AVG_Z_S(Knu_sL_ip(IJKB,M,L),Knu_sL_ip(IJK,M,L), KM)
               Knu_sM_b = AVG_Z_S(Knu_sM_ip(IJKB,M,L),Knu_sM_ip(IJK,M,L), KM)

               S14a_sum_rhs = S14a_sum_rhs + ( (Knu_sL_e*(NU_PL_E-NU_PL_p) - &
                    Knu_sM_e*(NU_PM_E-NU_PM_p) )*ODX_E(I)*AYZ(IJK) - (Knu_sL_w*&
                    (NU_PL_p-NU_PL_W) - Knu_sM_w*(NU_PM_p-NU_PM_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S14b_sum_rhs = S14b_sum_rhs + ( (Knu_sL_n*(NU_PL_N-NU_PL_p) - &
                    Knu_sM_n*(NU_PM_N-NU_PM_p) )*ODY_N(J)*AXZ(IJK) - (Knu_sL_s*&
                    (NU_PL_p-NU_PL_S) - Knu_sM_s*(NU_PM_p-NU_PM_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S14c_sum_rhs = S14c_sum_rhs + ( (Knu_sL_t*(NU_PL_T-NU_PL_p) - &
                    Knu_sM_t*(NU_PM_T-NU_PM_p) )*ODZ_T(K)*OX(I)*AXY(IJK) - &
                    (Knu_sL_b*(NU_PL_p-NU_PL_B) - Knu_sM_b*(NU_PM_p-NU_PM_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )
!
!
!              Part of Heat Flux: div (q)
!                   Kvel_s_ip*[vi-vp]
!              Note S9 terms should evaluate to zero for particles from the same phase
               Kvel_s_e = AVG_X_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKE,M,L),I)
               Kvel_s_w = AVG_X_H(Kvel_s_ip(IJKW,M,L),Kvel_s_ip(IJK,M,L), IM)
               Kvel_s_n = AVG_Y_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKN,M,L),J)
               Kvel_s_s = AVG_Y_H(Kvel_s_ip(IJKS,M,L),Kvel_s_ip(IJK,M,L), JM)
               Kvel_s_t = AVG_Z_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKT,M,L),K)
               Kvel_s_b = AVG_Z_H(Kvel_s_ip(IJKB,M,L),Kvel_s_ip(IJK,M,L), KM)
 
               S9_sum_rhs = S9_sum_rhs + ( Kvel_s_e*(UsM_e-UsL_e)*AYZ(IJK) - &
                    Kvel_s_w*(UsM_w-UsL_w)*AYZ(IMJK) + Kvel_s_n*(VsM_n-VsL_n)*AXZ(IJK)-&
                    Kvel_s_s*(VsM_s-VsL_s)*AXZ(IJMK) + Kvel_s_t*(WsM_t-WsL_t)*AXY(IJK)-&
                    Kvel_s_b*(WsM_b-WsL_b)*AXY(IJKM) )
!
!
          ENDIF    ! (IF M.NE.L)
!
      ENDDO
!
!  WARNING: The terms due to granular temperature gradients S21 (a,b,c) have caused
!           some converegence issues, remove them from LHS and RHS for debugging (sof).
!
      SOURCELHS = ( (S11a_sum_lhs+S11b_sum_lhs+S11c_sum_lhs+S12a_sum_lhs+&
          S12b_sum_lhs+S12c_sum_lhs) + (S10_lhs+S16_lhs+S17_sum_lhs+&
          S18_sum_lhs-S20_sum_lhs + S13_sum_lhs)*VOL(IJK) + &
	  ZMAX(S21a_sum_rhs+S21b_sum_rhs+S21c_sum_rhs)+ &
	  ZMAX(S14a_sum_rhs+S14b_sum_rhs+S14c_sum_rhs)+ ZMAX(S9_sum_rhs) ) / &
          Theta_m(IJK,M)
!
!
      SOURCERHS = ( S10_rhs+ S15_rhs + S16_rhs + S17_sum_rhs+S18_sum_rhs  + S13_sum_rhs) * VOL(IJK) + &
          S11a_sum_rhs+S11b_sum_rhs+S11c_sum_rhs+S12a_sum_rhs+S12b_sum_rhs+S12c_sum_rhs + &
          ZMAX(- (S14a_sum_rhs+S14b_sum_rhs+S14c_sum_rhs) ) + ZMAX(-S9_sum_rhs) + &
	  ZMAX(- (S21a_sum_rhs+S21b_sum_rhs+S21c_sum_rhs) )
	  
!	  
!
      RETURN  
      END SUBROUTINE SOURCE_IA_NONEP_GRANULAR_ENERGY
!-----------------------------------------------  

