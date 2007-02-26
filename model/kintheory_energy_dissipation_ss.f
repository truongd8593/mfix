!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      SUBROUTINE CALC_IA_NONEP_ENERGY_DISSIPATION_SS(M, IER) 
!
!-----------------------------------------------
!     Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar 
      USE indices 
      USE physprop
      USE run
      USE constant
      USE toleranc
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------  
!     Error index
      INTEGER          IER
!
!     Index
      INTEGER          IJK, I, J, K
!     
!     Solids phase
      INTEGER          L,M
!
!                      Index for storing solids-solids drag coefficients 
!                      in the upper triangle of the matrix                  
      INTEGER          LM 
!
!     variables for IA equipartition model
      DOUBLE PRECISION ED_common_term
      DOUBLE PRECISION EDvel_sL, EDvel_sM
      DOUBLE PRECISION M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION Ap_lm, Dp_lm, R1p_lm, R10p_lm, R3p_lm, R4p_lm,&
                       R5p_lm, Bp_lm
!----------------------------------------------- 
!     Function subroutines
!----------------------------------------------- 
      DOUBLE PRECISION G_0
!-----------------------------------------------
!     Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------   
!
      DO IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)       
!     
          IF ( FLUID_AT(IJK) ) THEN
!
               D_PM = D_P(IJK,M)
               M_PM = (PI/6.d0)*(D_PM**3)*RO_S(M)
               NU_PM = ROP_S(IJK,M)/M_PM

               DO L = 1, MMAX
!
                    LM = FUNLM(L,M)
! 
                    D_PL = D_P(IJK,L)
                    M_PL = (PI/6.d0)*(D_PL**3)*RO_S(L)
                    MPSUM = M_PM + M_PL
                    DPSUMo2 = (D_PM+D_PL)/2.d0
                    NU_PL = ROP_S(IJK,L)/M_PL
!
!
                    ED_common_term = (3.d0/4.d0)*(DPSUMo2*DPSUMo2)*(1.d0+C_E)*&
                         G_0(IJK,M,L)*NU_PM*NU_PL*(M_PM*M_PL/MPSUM)*&
                         ((M_PM*M_PL)**1.5)
!
!                                  
                    IF (M .eq. L) THEN
!
                         Ap_lm = MPSUM/(2.d0)
                         Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
                         R1p_lm = 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) )
                         R3p_lm = 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) )
!
!
                         ! Dissipation associated with the difference in temperature
                         ! (e.g. interphase transfer term).  For the case case (L=M) 
                         ! the term ED_s * (TL-TM) cancels. Therefore, explicity set
                         ! the term to zero for this case.
                         ED_ss_ip(IJK,LM) = ZERO                        
!
!
                         ! Dissipation associated with temperature
                         EDT_s_ip(IJK,M,L) = -ED_common_term* (1.d0-C_E)*(M_PL/MPSUM)*&
                              (DSQRT(PI)/6.d0)*R1p_lm*(Theta_m(IJK,M)**1.5)
!
!
                         ! Dissipation associated with divergence of
                         ! velocity of solid phase L: do not explicity include
                         ! terms that cancel when EDvel_sL is summed with EDvel_sM
                         EDvel_sL = ED_common_term*((1.d0-C_E)*(M_PL/MPSUM)*&
                              (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm)
                         EDvel_sL_ip(IJK,M,L) = EDvel_sL*(Theta_m(IJK,L))
!
!
                         ! Dissipation associated with divergence of
                         ! velocity of solid phase M: do not explicity include
                         ! terms that cancel when EDvel_sL is summed with EDvel_sM
                        
			! commented by sof, no need to re-do computation
			! EDvel_sM = ED_common_term*((1.d0-C_E)*(M_PL/MPSUM)*&
                        !      (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm)
                        ! EDvel_sM_ip(IJK,M,L) = EDvel_sM*(Theta_m(IJK,M))
                         
			 EDvel_sM_ip(IJK,M,L) =  EDvel_sL_ip(IJK,M,L)
!
                    ELSE
!
                         Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/&
                               2.d0
                         Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-Theta_m(IJK,M) ))/&
                              (2.d0*MPSUM)
                         Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+M_PL*Theta_m(IJK,L) ))/&
                              (2.d0*MPSUM*MPSUM)
                         
			 R1p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) ) )+ &
                              ( (9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4 ) )+&
                              ( (30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 * Dp_lm**5 ) )
                         
			 R3p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) ) )+&
                              ( (21.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**2.5 * Dp_lm**4.5 ) )+&
                              ( (315.d0*Bp_lm**4)/( 8.d0 * Ap_lm**3.5 *Dp_lm**5.5 ) )
                         
			 R4p_lm = ( 3.d0/( Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                              ( (35.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**3.5 * Dp_lm**4.5 ) )+&
                              ( (441.d0*Bp_lm**4)/( 8.d0 * Ap_lm**4.5 * Dp_lm**5.5 ) )
                         
			 R5p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**3 ) )+ &
                              ( (5.d0*Bp_lm*Bp_lm)/( Ap_lm**3.5 * Dp_lm**4 ) )+&
                              ( (14.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5 ) )
                         
			 R10p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**2.5 ) )+&
                              ( (25.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**3.5 * Dp_lm**3.5 ) )+&
                              ( (1225.d0*Bp_lm**4)/( 24.d0* Ap_lm**4.5 * Dp_lm**4.5 ) )
!
!
                         ! Dissipation associated with the difference in temperature
                         ! (e.g. interphase transfer term). to solved using PEA
                         ED_ss_ip(IJK,LM) = ED_common_term*DSQRT(PI)*(M_PM*M_PL/&
                             (2.d0*MPSUM))*R5p_lm*( (Theta_M(IJK,M)*Theta_M(IJK,L))**3 )
!
!
                         ! Dissipation associated with temperature
                         EDT_s_ip(IJK,M,L) = -ED_common_term* (1.d0-C_E)*(M_PL/MPSUM)*&
                              (DSQRT(PI)/6.d0)*R1p_lm*( (Theta_m(IJK,M)*Theta_m(IJK,L))**3 )
!
!
                         ! Dissipation associated with divergence of
                         ! velocity of solid phase L  
                         EDvel_sL = ED_common_term*( ((3.d0*DPSUMo2*PI/40.d0)*M_PL*&
                              R10p_lm)+( (DPSUMo2*PI/4.d0)*(M_PM*M_PL/MPSUM)*Bp_lm*&
                              R4p_lm)-( (1.d0-C_E)*(M_PL/MPSUM)*(DPSUMo2*PI/16.d0)*&
                              M_PL*Bp_lm*R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*&
                              (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm))
                         
			 EDvel_sL_ip(IJK,M,L) = EDvel_sL*( Theta_m(IJK,M)**3.5 *&
                               Theta_m(IJK,L)**2.5 )

                         ! Dissipation associated with divergence of
                         ! velocity of solid phase M  
                         EDvel_sM = ED_common_term*( (-(3.d0*DPSUMo2*PI/40.d0)*M_PM*&
                              R10p_lm)+( (DPSUMo2*PI/4.d0)*(M_PM*M_PL/MPSUM)*Bp_lm*&
                              R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*(DPSUMo2*PI/16.d0)*&
                              M_PM*Bp_lm*R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*&
                              (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm))
                         
			 EDvel_sM_ip(IJK,M,L) = EDvel_sM*( Theta_m(IJK,M)**2.5 *&
                              Theta_m(IJK,L)**3.5 )
!
                    ENDIF

               ENDDO

          ENDIF
      ENDDO

      RETURN  
      END SUBROUTINE CALC_IA_NONEP_ENERGY_DISSIPATION_SS 
