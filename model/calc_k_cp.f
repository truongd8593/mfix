!TO DO:
!1. Kcp needs to be defined for each solids phase (?).
!2. The part of Kcp from P_star should be based on the sum of EP_s of
!     close-packed solids.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_K_cp(Kcp, IER)                                    C
!  Purpose: Calculate and store dPodEp_s                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_K_cp(Kcp, IER)
!
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE indices
      USE pscor
      USE geometry
      USE constant
      USE run
      USE visc_s
      USE trace
      USE compar
      USE sendrecv
      IMPLICIT NONE
 
      INTEGER          IJK, M
!
!                      dPodEP_s
      DOUBLE PRECISION Kcp(DIMENSION_3)
!
!                      error index
      INTEGER          IER
 
!              Radial distribution function
      DOUBLE PRECISION g0, G_0
 
!                      Other variables
      DOUBLE PRECISION Pc, DPcoDEPs, Mu, Mu_b, Mu_zeta, ZETA
      DOUBLE PRECISION F2, DF2oDEPs, DEPs2G_0oDEPs, Pf, Pfmax
      DOUBLE PRECISION  DZETAoDEPs, DG_0DNU
      double precision calc_ep_star
      
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'


!      DO 200 M = 1, MMAX
        M = Mcp

        IF(CLOSE_PACKED(M)) THEN

!$omp     parallel do  firstprivate(M) &
!$omp&    private( IJK, DEPs2G_0oDEPs, &
!$omp&             Pc, DPcoDEPS,  Mu, Mu_b, Mu_zeta, ZETA, &
!$omp&             F2, DF2oDEPs, Pf, Pfmax )
          DO 100 IJK = ijkstart3, ijkend3
            IF(.NOT.WALL_AT(IJK))THEN
 
! start anuj 4/20
               DEPs2G_0oDEPs= EP_s(IJK,M)*EP_s(IJK,M)*&
                              DG_0DNU(EP_s(IJK,M)) +&
                              2d0*EP_s(IJK,M)*G_0(IJK,M,M)
 
               IF (FRICTION) THEN
 
                  IF (EP_s(IJK,M).GT.EPS_f_min) THEN
 
	             IF (EP_s(IJK,M).GT.EPS_max) THEN
              	        Pc = 1d25*((EP_s(IJK,M) - EPS_max)&
                                                      **10d0)
                        DPcoDEPS =&
                             1d26*((EP_s(IJK,M) - EPS_max)**9d0)
 
	             ELSE
                        Pc = Fr*((EP_s(IJK,M) - EPS_f_min)**N_Pc)/&
                          ((EPS_max - EP_s(IJK,M) + SMALL_NUMBER)&
                           **D_Pc)
 
                        DPcoDEPs =&
                           Fr*((EP_s(IJK,M) - EPS_f_min)**(N_Pc - 1.))&
                           *(N_Pc*(EPS_max - EP_s(IJK,M)) -&
                             D_Pc*(EP_s(IJK,M) - EPS_f_min))&
      	                   / ((EPS_max - EP_s(IJK,M) + SMALL_NUMBER)**&
                              (D_Pc + 1.))
                     ENDIF
 
 
                     Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))&
                           *D_p(IJK,M)*RO_s(M))/96d0
 
 
                     Mu_b = (256d0*Mu*EP_s(IJK,M)*EP_s(IJK,M)&
                          *G_0(IJK,M,M))/(5d0*Pi)
 
                     IF (SAVAGE.EQ.1) THEN
            	      Mu_zeta =&
                           ((2d0+ALPHA)/3d0)*((Mu/(Eta*(2d0-Eta)*&
                           G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                           G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                           EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b*Eta))
 
 
                      ZETA = ((48d0*Eta*(1d0-Eta)*RO_s(M)*EP_s(IJK,M)*&
                            EP_s(IJK,M)*G_0(IJK,M,M)*&
                            (Theta_m(IJK,M)**1.5d0))/&
                            (SQRT_Pi*D_p(IJK,M)*2d0*Mu_zeta))**0.5d0
 
                     ELSEIF (SAVAGE.EQ.0) THEN
                      ZETA = (SMALL_NUMBER +&
                             trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0))**0.5d0
 
                     ELSE
                      ZETA = ((Theta_m(IJK,M)/(D_p(IJK,M)*D_p(IJK,M))) +&
                             (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                             trD_s_C(IJK,M))/3.d0)))**0.5d0
 
                     ENDIF
 
                     IF ((trD_s_Co(IJK,M)/(ZETA*N_Pf*&
                         DSQRT(2d0)*&
                         Sin_Phi)) .GT. 1d0) THEN
                      F2 = 0d0
                      DF2oDEPs = ZERO
 
                     ELSE
                      F2 = (1d0 - (trD_s_Co(IJK,M)/(ZETA*N_Pf*&
                          DSQRT(2d0)*Sin_Phi)))**(N_Pf-1d0)
 
                      IF (SAVAGE.EQ.1) THEN
 
                       DF2oDEPs = (N_Pf-1d0)*(F2**(N_Pf-2d0))*&
                         trD_s_Co(IJK,M)&
                         *DZETAoDEPs(EP_s(IJK,M), IJK, M)&
                         / (ZETA*ZETA*&
                         N_Pf*DSQRT(2d0)*Sin_Phi)
 
                      ELSE
                       DF2oDEPs=ZERO
                      ENDIF
 
                      Pf = Pc*F2
 
                      Pfmax = Pc*((N_Pf/(N_Pf-1d0))**(N_Pf-1d0))
 
                      IF (Pf> Pfmax) THEN
                         F2 = (N_Pf/(N_Pf-1d0))**(N_Pf-1d0)
                         DF2oDEPS = ZERO
                      ENDIF
 
                     ENDIF
 
!    Contributions to Kcp(IJK) from kinetic theory have been left out in the expressions below
!    as they cause convergence problems at low solids volume fraction.
 
                      Kcp(IJK) =&
                            F2*DPcoDEPS&
                          + Pc*DF2oDEPS
 
		  ELSE
 
			Kcp(IJK) = ZERO
 
                  ENDIF
 
 
               ELSE ! FRICTION = .FALSE.
!GERA*****************************
                 if (MMAX >= 2)EP_star = Calc_ep_star(IJK, IER)
!GERA_END************************* 
                 IF(EP_g(IJK) .LT. EP_star) THEN
                    Kcp(IJK) = dPodEP_s(EP_s(IJK, M))
 
		 ELSE
 		    Kcp(IJK) = ZERO	
 
                 ENDIF

	       ENDIF
 
!//SP
            ELSE
	       Kcp(IJK) = ZERO
 
            ENDIF
 
 100     CONTINUE
        ENDIF
!   end anuj 4/20
 
! 200   CONTINUE

      CALL send_recv(Kcp, 2)
      
      RETURN
      END
 
 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DZETAoDEPs (EPs,IJK,M)                                 C
!  Purpose: Calculate derivative of zeta                               C
!           w.r.t granular volume fraction                             C
!                                                                      C
!  Author: A. Srivastava                              Date: 8-JUNE-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION DZETAoDEPs(EPs, IJK, M)
!
      USE param 
      USE param1 
      USE fldvar
      USE physprop
      USE indices
      USE pscor
      USE geometry
      USE constant
      USE run
      USE visc_s
      USE trace
      USE compar
      IMPLICIT NONE
!
!                      solids volume fraction
      DOUBLE PRECISION, intent(in) :: EPs
!
!                      radial distribution function
      DOUBLE PRECISION g0
 
!                      Other variables
      DOUBLE PRECISION Mu, Mu_b,  DEPs2G_0oDEPs, F1, DF1oDEPs
 
      DOUBLE PRECISION DG_0DNU
      INTEGER , intent(in) :: IJK,M

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: G_0
!-----------------------------------------------
 
      INCLUDE 'ep_s1.inc'
      INCLUDE 's_pr1.inc'
      INCLUDE 'function.inc'
      INCLUDE 's_pr2.inc'
      INCLUDE 'ep_s2.inc'
 
 
      g0 = G_0(IJK, M, M)
 
      Mu = (5d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(IJK,M)*RO_s(M))/96d0
 
      Mu_b = (256d0*Mu*EPs*EPs*g0/(5d0*Pi))
 
      DEPs2G_0oDEPs = EPs*EPs*DG_0DNU(EPs) + 2d0*EPs*g0
 
      F1 = ((2d0+ALPHA)/3d0)*((2*Mu/(Eta*(2d0-Eta)*&
           g0))*(1d0+1.6d0*Eta*EPs*g0)*(1d0+1.6d0*Eta*(3d0*Eta-2d0)&
           *EPs*g0)+(1.2d0*Mu_b*Eta))
 
      DF1oDEPs = ((2d0+ALPHA)/3d0)*((2*Mu/(Eta*(2d0-Eta))*&
         ((-DG_0DNU(EPs)/(g0*g0)) + (1.6d0*Eta*(3d0*Eta-1d0))&
        + (64d0*Eta*Eta*(3d0*Eta-2d0)*DEPs2G_0oDEPs/25d0))) +&
        3.2d0*Eta*RO_s(M)*D_p(IJK,M)*((Theta_m(IJK,M)/Pi)**0.5d0)&
          *DEPs2G_0oDEPs)
 
      DZETAoDEPs = 0.5d0*((48d0*Eta*(1d0-Eta)*RO_s(M)*F1*&
              (Theta_m(IJK,M)**1.5d0)/&
           (SQRT_Pi*D_p(IJK,M)*EPs*EPs*g0))**0.5d0)*&
           (F1*DEPs2G_0oDEPs - EPs*Eps*g0*DF1oDEPs)/(F1*F1)
 
 
      RETURN
      END

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
