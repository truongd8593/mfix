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
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE indices
      USE constant
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
      INTEGER          IJK, I, J, K, M
!
!                      Source terms to be kept on rhs
      DOUBLE PRECISION sourcerhs
!
!                      Source terms to be kept on lhs
      DOUBLE PRECISION sourcelhs
!
!                      Slip velocity
      DOUBLE PRECISION VSLIP
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
      SOURCERHS = (ZMAX(LAMBDA_S_C(IJK,M))*TRD_S_C(IJK,M)**2+2.*MU_S_C(IJK,M)*&
         TRD_S2(IJK,M)+P_S_C(IJK,M)*ZMAX((-TRD_S_C(IJK,M))))*VOL(IJK) 
!
      VSLIP = (U_S(IJK,M)-U_G(IJK))**2 + (V_S(IJK,M)-V_G(IJK))**2 + (W_S(IJK,M)&
         -W_G(IJK))**2 
      VSLIP = DSQRT(VSLIP) 
!
      SOURCERHS = SOURCERHS + (SWITCH*81D0*EP_S(IJK,M)*(MU_G(IJK)*VSLIP)**2D0/(&
         G_0(IJK,M,M)*D_P(M)**3D0*RO_S(M)*(PI*THETA_M(IJK,M)+SMALL_NUMBER)**&
         0.5D0))*VOL(IJK) 
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
      SOURCELHS = (SWITCH*3.*F_GS(IJK,M)+(48./DSQRT(PI))*ETA*(1.-ETA)*ROP_S(IJK&
         ,M)*EP_S(IJK,M)*G_0(IJK,M,M)*DSQRT(THETA_M(IJK,M))/D_P(M)+(ROP_S(IJK,M&
         )*(1D0+4.*ETA*EP_S(IJK,M)*G_0(IJK,M,M)))*ZMAX(TRD_S_C(IJK,M))+ZMAX((-&
         LAMBDA_S_C(IJK,M)))*TRD_S_C(IJK,M)**2/DSQRT(THETA_M(IJK,M)+&
         SMALL_NUMBER))*VOL(IJK) 
!
      RETURN  
      END SUBROUTINE SOURCE_GRANULAR_ENERGY 
