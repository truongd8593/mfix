!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  Module name: CALC_IA_NONEP_DRAG_ss(L, M, IER)
!  Purpose: This module computes the coefficient of drag between
!           solids phase m and solids phase l using Iddir & 
!           Arastoopour (2005) kinetic theory model
!
!  Literature/Document References:    
!	Iddir, Y.H., "Modeling of the multiphase mixture of particles 
!	  using the kinetic theory approach," PhD Thesis, Illinois
!	  Institute of Technology, Chicago, Illinois, 2004
!    Iddir, Y.H., & H. Arastoopour, "Modeling of Multitype particle 
!      flow using the kinetic theory approach," AIChE J., Vol 51,
!      no. 6, June 2005
!                                                         
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      SUBROUTINE CALC_IA_NONEP_DRAG_SS(L, M, IER) 
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE constant
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar 
      USE sendrecv 
      USE drag
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------   
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK
! 
!                      Index of solids phases 
      INTEGER          L, M 
! 
!                      Index for storing solids-solids drag coefficients 
!                      in the upper triangle of the matrix 
      INTEGER          LM 
!
!                      Particle diameters 
      DOUBLE PRECISION D_PM, D_PL 
! 
!                      Sum of particle diameters 
      DOUBLE PRECISION DPSUM 
! 
      DOUBLE PRECISION M_PM, M_PL, MPSUM, DPSUMo2, NU_PL, NU_PM
      DOUBLE PRECISION Ap_lm, Dp_lm, R0p_lm, R2p_lm, R3p_lm, R4p_lm, R10p_lm, Bp_lm
      DOUBLE PRECISION Fss_ip, Fnus_ip, FTsM_ip, FTsL_ip, F_common_term
!
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION , EXTERNAL :: G_0 
!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------
!
      DO IJK = ijkstart3, ijkend3
!         
          IF (.NOT.WALL_AT(IJK)) THEN 
!		  
               LM = FUNLM(L,M)
!
               IF (M == L) THEN
                    F_SS(IJK,LM) = ZERO
                    Fnu_s_ip(IJK,M,L) = ZERO
                    FT_sM_ip(IJK,M,L) = ZERO
                    FT_sL_ip(IJK,M,L) = ZERO
!
               ELSE 
                    D_PM = D_P(IJK,M) 
                    D_PL = D_P(IJK,L) 
                    DPSUM = D_PL + D_PM 

                    M_PM = (Pi/6.d0) * D_PM**3 *RO_S(M)
                    M_PL = (Pi/6.d0) * D_PL**3 *RO_S(L)
                    MPSUM = M_PM + M_PL
                    DPSUMo2 = DPSUM/2.d0
                    NU_PM = ROP_S(IJK,M)/M_PM
                    NU_PL = ROP_S(IJK,L)/M_PL
!
                 IF(Theta_m(IJK,M) > ZERO .AND. Theta_m(IJK,L) > ZERO) THEN
!
                    Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/&
                          2.d0
                    Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-Theta_m(IJK,M) ))/&
                         (2.d0*MPSUM)
                    Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+M_PL*Theta_m(IJK,L) ))/&
                         (2.d0*MPSUM*MPSUM)
                    
		    R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+ &
                              ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                              ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 * Dp_lm**4.5 ) )
                    
		    R2p_lm = ( 1.d0/( 2.d0*Ap_lm**1.5 * Dp_lm*Dp_lm ) )+&
                              ( (3.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**3 ) )+&
                              ( (15.d0*Bp_lm**4)/( 2.d0*Ap_lm**3.5 * Dp_lm**4 ) )
                    
		    R3p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) ) )+&
                              ( (21.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**2.5 * Dp_lm**4.5 ) )+&
                              ( (315.d0*Bp_lm**4)/( 8.d0 * Ap_lm**3.5 *Dp_lm**5.5 ) )
                    
		    R4p_lm = ( 3.d0/( Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                              ( (35.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**3.5 * Dp_lm**4.5 ) )+&
                              ( (441.d0*Bp_lm**4)/( 8.d0 * Ap_lm**4.5 * Dp_lm**5.5 ) )
                    
		    R10p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**2.5 ) )+&
                              ( (25.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**3.5 * Dp_lm**3.5 ) )+&
                              ( (1225.d0*Bp_lm**4)/( 24.d0* Ap_lm**4.5 * Dp_lm**4.5 ) )
!
!
                    F_common_term = (DPSUMo2*DPSUMo2/4.d0)*(M_PM*M_PL/MPSUM)*&
                         G_0(IJK,M,L)*(1.d0+C_E)*(M_PM*M_PL)**1.5
!
!
!                   Momentum source associated with relative velocity 
!                   between solids phase m and solid phase l
                    Fss_ip = F_common_term*NU_PM*NU_PL*DSQRT(PI)*R2p_lm*&
                           (Theta_m(IJK,M)*Theta_m(IJK,L))**2
!
!                   Momentum source associated with the difference in the gradients in
!                   number density of solids phase m and all other solids phases
                    Fnus_ip = F_common_term*(PI*DPSUMo2/12.d0)*R0p_lm*&
                           (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5
!
!                   Momentum source associated with the gradient in granular temperature of
!                   solid phase M
                    FTsM_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                              (Theta_m(IJK,M)**1.5 * Theta_m(IJK,L)**2.5) *&
                              (  (-1.5d0/12.d0*R0p_lm)+&
                              Theta_m(IJK,L)/16.d0*(  (-M_PM*R10p_lm) - &
                              ((5.d0*M_PL*M_PL*M_PM/(192.d0*MPSUM*MPSUM))*R3p_lm)+&
                              ((5.d0*M_PM*M_PL)/(96.d0*MPSUM)*R4p_lm*Bp_lm)  )  )
!
!                   Momentum source associated with the gradient in granular temperature of
!                   solid phase L ! no need to recompute (sof Aug 30 2006)
                    FTsL_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                             (Theta_m(IJK,L)**1.5 * Theta_m(IJK,M)**2.5) *&
                              (  (1.5d0/12.d0*R0p_lm)+&
                              Theta_m(IJK,M)/16.d0*(  (M_PL*R10p_lm)+&
                              (5.d0*M_PM*M_PM*M_PL/(192.d0*MPSUM*MPSUM)*R3p_lm)+&
                              (5.d0*M_PM*M_PL/(96.d0*MPSUM) *R4p_lm*Bp_lm)  )  )

		    F_SS(IJK,LM) = Fss_ip
               
                    Fnu_s_ip(IJK,M,L) = Fnus_ip 
!
!      WARNING: the following two terms have caused some convergence problems earlier
!               Set them to ZERO for debugging in case of converegence issues. (sof)
                    FT_sM_ip(IJK,M,L) = FTsM_ip  ! ZERO

                    FT_sL_ip(IJK,M,L) = FTsL_ip  ! ZERO
                 ELSE
		    F_SS(IJK,LM) = ZERO
                    Fnu_s_ip(IJK,M,L) = ZERO
                    FT_sM_ip(IJK,M,L) = ZERO
                    FT_sL_ip(IJK,M,L) = ZERO
                 ENDIF
               ENDIF
          ENDIF
      ENDDO 
           
      RETURN  
      END SUBROUTINE CALC_IA_NONEP_DRAG_SS 



!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3







