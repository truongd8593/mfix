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
!      Yin, X, Sundaresan, S. (2008).  AIChE			       C
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
! --- Correction factors for implementing polydisperse drag model 
!     proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION FA_cor, FB_cor, FA, FB, tmp_sum, tmp_fac

! --- Correction factor for implementing HYS drag law
      DOUBLE PRECISION del_U_mix,tmp_del_U_mix 
!
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
!     
!     

!***********************************************************
!     Declaration of variables relevant to the HYS drag correlation
!***********************************************************

!     Index for particles of other species
      INTEGER j

!     Polydisperse correction factor for YS drag relation
      DOUBLE PRECISION a_YS

!     Cell center value of x-particle velocity (HYS drag relation)
      DOUBLE PRECISION USCM_HYS

!     Cell center value of y-particle velocity (HYS drag relation)
      DOUBLE PRECISION VSCM_HYS

!     Cell center value of z-particle velocity (HYS drag relation)
      DOUBLE PRECISION WSCM_HYS

!     Lubrication interaction prefactor in YS drag relation	
      DOUBLE PRECISION alpha_YS

!     Friction coefficient for a particle of type i (HYS drag relation)
      DOUBLE PRECISION beta_i_HYS

!     Friction coefficient for a particle of type j (HYS drag relation)
      DOUBLE PRECISION beta_j_HYS

!     Magnitude of gas-solids relative velocity for polydisperse HYS drag 
!     implementation
      DOUBLE PRECISION VREL_poly 

!     Stokes drag of a particle of type j
      DOUBLE PRECISION	FSTOKES_j

!     Diameter ratio used as temp variable
      DOUBLE PRECISION Y_i_temp

!     Variable for Beetstra et. al. drag relation
      DOUBLE PRECISION F_D_BVK

!     Variable for YS drag relation
      DOUBLE PRECISION F_YS

!     Minimum particle diameter in mixture
      DOUBLE PRECISION min_D_p

!      End of variable definition for HYS drag relation
!***********************************************************
!     Current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER) ! Tsuji drag
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
            IF (P_FLOW_AT(IJK)) THEN
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
!--------------------------Begin GIDASPOW_PCF --------------------
!    Additional modifications to apply the polydisperse correction 
!    factor (PCF) proposed by Beetstra et al. (2007)    09-26-07
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_PCF') THEN
!
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM)
	       ENDDO 
               D_p_av = ZERO
               tmp_sum = ZERO
	       tmp_fac = ZERO
	       DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN
                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)
                 ELSE
                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells
		 ENDIF
               ENDDO 
               D_p_av = ONE / tmp_sum

               Y_i = D_p(IJK,M)/D_p_av

               IF (Mu > ZERO) THEN
                    RE_G = D_p_av*VREL*ROP_G(IJK)/Mu
               ELSE
                    RE_G = LARGE_NUMBER 
               ENDIF
!
               IF(EP_g(IJK) .LE. 0.8D0) THEN
                    DgA = 150D0 * (ONE - EP_g(IJK)) * Mu &
                         / ( EP_g(IJK) * D_p_av**2 ) &
                         + 1.75D0 * RO_g(IJK) * VREL / D_p_av
               ELSE
                    IF(RE_G .LE. 1000D0)THEN
                         C_d = (24.D0/(RE_G+SMALL_NUMBER)) * (ONE + 0.15D0 * RE_G**0.687D0)
                    ELSE
                         C_d = 0.44D0
                    ENDIF
                    DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                         /D_p_av 
               ENDIF


               ! see the associated erratum by Beetstra et al. (2007) :
               FA_cor = Y_i 
               IF (M .EQ. 1) THEN
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2) 
               ELSE
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) 
               ENDIF

               FA = ONE/(Y_i*Y_i) * DgA * FA_cor 
               FB = ONE/(Y_i*Y_i) * DgA * FB_cor 

               IF (RE_G .EQ. ZERO) THEN
                    FA = ZERO
                    FB = ZERO
               ENDIF

!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    F_gstmp = FB * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    F_gstmp = FA * EP_s(IJK, M)
               ENDIF
!
!              
!--------------------------End GIDASPOW_PCF ----------------------
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
!
!-----------------------End Gidaspow_blend -----------------------
!     
!----------------------Begin Gidaspow_BLEND_PCF ------------------
!    Additional modifications to apply the polydisperse correction 
!    factor (PCF) proposed by Beetstra et al. (2007)    09-26-07
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND_PCF') THEN
!
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM)
	       ENDDO 
               D_p_av = ZERO
               tmp_sum = ZERO
	       tmp_fac = ZERO
	       DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN
                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)
                 ELSE
                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells
		 ENDIF
               ENDDO 
               D_p_av = ONE / tmp_sum

               Y_i = D_p(IJK,M)/D_p_av

               IF (Mu > ZERO) THEN
                    RE_G = D_p_av*VREL*ROP_G(IJK)/Mu
               ELSE
                    RE_G = LARGE_NUMBER
               ENDIF
!
!              Dense phase - EP_g < 0.8
               Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu &
                    / ( EP_g(IJK) * D_p_av**2 ) &
                    + 1.75D0 * RO_g(IJK) * VREL / D_p_av
!
!              Dilute phase - EP_g >= 0.8
               IF(RE_G .LE. 1000D0)THEN
                    C_d = (24.D0/(RE_G+SMALL_NUMBER)) * (ONE + 0.15D0 * RE_G**0.687D0)
               ELSE
                    C_d = 0.44D0
               ENDIF
               WenYu = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                    /D_p_av
!     
!              Switch function
               PHI_gs = ATAN(150D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
!
!              Blend the models
               DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
!               
               ! see the associated erratum by Beetstra et al. (2007) :
               FA_cor = Y_i 
               IF (M .EQ. 1) THEN
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2) 
               ELSE
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) 
               ENDIF

               FA = ONE/(Y_i*Y_i) * DgA * FA_cor
               FB = ONE/(Y_i*Y_i) * DgA * FB_cor 

               IF (RE_G .EQ. ZERO) THEN
                    FA = ZERO
                    FB = ZERO
               ENDIF

!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    F_gstmp = FB * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    F_gstmp = FA * EP_s(IJK, M)
               ENDIF
!
!               
!---------------------End Gidaspow_BLEND_PCF ---------------------
!     
!--------------------------Begin WEN_YU --------------------------
            ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU') then
                IF(Re_G .LE. 1000D0)THEN
                   C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
                ELSE
                   C_d = 0.44D0
                ENDIF
!
                DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                  /D_p(IJK,M)
               
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                  F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                  F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
               
!--------------------------End WEN_YU ----------------------------
! 
!--------------------------Begin WEN_YU_PCF ----------------------
!    Additional modifications to apply the polydisperse correction 
!    factor (PCF) proposed by Beetstra et al. (2007)    09-26-07
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU_PCF') then
!
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM)
	       ENDDO 
               D_p_av = ZERO
               tmp_sum = ZERO
	       tmp_fac = ZERO
	       DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN
                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)
                 ELSE
                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells
		 ENDIF
               ENDDO 
               D_p_av = ONE / tmp_sum

               Y_i = D_p(IJK,M)/D_p_av

               IF (Mu > ZERO) THEN
                    RE_G = D_p_av*VREL*ROP_G(IJK)/Mu
               ELSE
                    RE_G = LARGE_NUMBER
               ENDIF
!
               IF(RE_G .LE. 1000D0)THEN
                    C_d = (24.D0/(RE_G+SMALL_NUMBER)) * (ONE + 0.15D0 * RE_G**0.687D0)
               ELSE
                    C_d = 0.44D0
               ENDIF

               DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
                    /D_p_av

               ! see the associated erratum by Beetstra et al. (2007) :
               FA_cor = Y_i 
               IF (M .EQ. 1) THEN
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2) 
               ELSE
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) 
               ENDIF

               FA = ONE/(Y_i*Y_i) * DgA * FA_cor
               FB = ONE/(Y_i*Y_i) * DgA * FB_cor 

               IF (RE_G .EQ. ZERO) THEN
                    FA = ZERO
                    FB = ZERO
               ENDIF

!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    F_gstmp = FB * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    F_gstmp = FA * EP_s(IJK, M)
               ENDIF
!
!              
!------------------------End WEN_YU_PCF --------------------------
    
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
            ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') THEN
!     
               F_STOKES = 18D0*Mu*EP_g(IJK)*EP_g(IJK)/D_p(IJK,M)**2
	       
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM) ! this is slightly /= one-ep_g due to round-off
	       ENDDO 
               w = EXP(-10.0D0*(0.4D0-phis)/phis)
	   
               IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
                    F_0 = (1.0D0-w) *                                             &
                         (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis   &
                         *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*           &
                         phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *          &
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
                    F_2 = (1.0D0-w) *                                            &
                         (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis  &
                         *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*          &
                         phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *        &
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
	   
!              This is a check for phis (or eps_(ijk,m)) to be within physical range
               IF(phis <= ZERO .OR. phis > ONE) F = zero
!
               DgA = F * F_STOKES
!   
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
!               
!     
!---------------------End Koch & Hill (2001) ----------------------
!
!-------------------- Begin Koch_Hill_PCF -------------------------
!    Additional modifications to apply the polydisperse correction 
!    factor (PCF) proposed by Beetstra et al. (2007)    09-26-07
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL_PCF') THEN
!     
               F_STOKES = 18D0*Mu*EP_g(IJK)*EP_g(IJK)/D_p(IJK,M)**2
	       
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM)
	       ENDDO 
               D_p_av = ZERO
               tmp_sum = ZERO
	       tmp_fac = ZERO
	       DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN
                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)
                 ELSE
                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells
		 ENDIF
               ENDDO 
               D_p_av = ONE / tmp_sum

               Y_i = D_p(IJK,M)/D_p_av

               w = EXP(-10.0D0*(0.4D0-phis)/phis)
	   
               IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
                    F_0 = (1.0D0-w) *                                             &
                         (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis   &
                         *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*           &
                         phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *          &
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
                    F_2 = (1.0D0-w) *                                            &
                         (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis  &
                         *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*          &
                         phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *        &
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

               IF (Mu > ZERO) THEN
!                   Note Reynolds' number for Hill and Koch has an 
!                   additional factor of 1/2 & ep_g
                    RE_kh = 0.5D0 * D_p_av * VREL * ROP_G(IJK) / Mu
               ELSE
                    RE_kh = LARGE_NUMBER
               ENDIF
	   
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
	   
!              This is a check for phis (or eps_(ijk,m)) to be within physical range
               IF(phis <= ZERO .OR. phis > ONE) F = zero
!
               ! see the associated erratum by Beetstra et al. (2007) :
               FA_cor = Y_i 
               IF (M .EQ. 1) THEN
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2) 
               ELSE
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) 
               ENDIF

               FA = FA_cor * F
               FB = FB_cor * F 
	      
               IF(Re_kh == ZERO) THEN
                    FA = ZERO
                    FB = ZERO
               ENDIF

!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    DgA = FB * F_STOKES
                    F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    DgA = FA * F_STOKES
                    F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
!                
!     
!--------------------- End Koch_Hill_PCF --------------------------
!
!---- Begin Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62 (Jan 2007) -----
! 
!
            ELSE IF(TRIM(DRAG_TYPE).EQ.'BVK') THEN
!
!              eq(9) BVK J. fluid. Mech. 528, 2005 (this F_Stokes is /= of Koch_Hill by a factor of ep_g)
               F_STOKES = 18D0*Mu*EP_g(IJK)/D_p(IJK,M)**2 
	       
               phis = ZERO
	       DO IM = 1, MMAX
	         phis = phis + EP_S(IJK,IM)
	       ENDDO 
               D_p_av = ZERO
               tmp_sum = ZERO
	       tmp_fac = ZERO
	       DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN
                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)
                 ELSE
                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells
		 ENDIF
               ENDDO 
               D_p_av = ONE / tmp_sum

               Y_i = D_p(IJK,M)/D_p_av

               IF (Mu > ZERO) THEN	      
                    RE = D_p_av*VREL*ROP_G(IJK)/Mu
               ELSE
                    RE = LARGE_NUMBER
               ENDIF
	      
               F = 10d0 * phis / EP_g(IJK)**2 + EP_g(IJK)**2 * (ONE+1.5d0*DSQRT(phis))
               F = F + 0.413d0*Re/(24d0*EP_g(IJK)**2) * (ONE/EP_G(IJK) + 3d0*EP_G(IJK) &
                    *phis + 8.4d0/Re**0.343) / (ONE+10**(3d0*phis)/Re**(0.5+2*phis))

               ! see the associated erratum by Beetstra et al. (2007) :
               ! the correction factor differs for model A versus model B
               ! application of the correction factor for model A is found from
               ! the correction factor for model B and neglects the Y_i**3 term
               FA_cor = Y_i
               IF (M .EQ. 1) THEN
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2) 
               ELSE
                    FB_cor = (EP_g(IJK)*Y_i + phis*Y_i**2 + 0.064d0*EP_g(IJK)*Y_i**3) 
               ENDIF
            
               FA = FA_cor * F
               FB = FB_cor * F 
	      
               IF(Re == ZERO) THEN
                    FA = ZERO
                    FB = ZERO
               ENDIF
  
!              Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
               IF(Model_B)THEN
                    DgA = FB * F_STOKES
                    F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
               ELSE
                    DgA = FA * F_STOKES
                    F_gstmp = DgA * EP_s(IJK, M)
               ENDIF
!     
!---- End Beetstra, van der Hoef, Kuipers, Chem. Eng. Science 62 (Jan 2007) -----
! 

!-------------------------- Begin HYS drag relation -----------------------------

	ELSE IF(TRIM(DRAG_TYPE).EQ.'HYS') THEN
	
        	F_STOKES = 18D0*Mu*EP_s(IJK,M)*EP_g(IJK)/D_p(IJK,M)**2 
                F_gstmp = ZERO 
                phis = ZERO
	        
                DO IM = 1, MMAX

	          phis = phis + EP_S(IJK,IM)

	        ENDDO 		

		D_p_av = ZERO
                tmp_sum = ZERO
	        tmp_fac = ZERO

	        DO IM = 1, MMAX
                 IF (phis .GT. ZERO) THEN

                   tmp_fac = EP_S(IJK,Im)/phis
                   tmp_sum = tmp_sum + tmp_fac/D_p(IJK,Im)

                 ELSE

                   tmp_sum = tmp_sum + ONE/ D_p(IJK,Im) ! not important, but will avoid NaN's in empty cells

		 ENDIF
               ENDDO 

               D_p_av = ONE / tmp_sum
               Y_i = ZERO
               Y_i = D_p(IJK,M)/D_p_av
	       a_YS = 1d0 - 2.66d0*phis + 9.096d0*phis**2 - 11.338d0*phis**3 

!	       Calculate smallest particle diameter in the mixture
	       min_D_p= MIN(D_p(IJK,1),D_p(IJK,2))

!	       Calculate smallest diameter if number of particle types is greater than 2
	       IF (MMAX > 2) THEN
		  DO IM=3,MMAX
		     min_D_p = MIN(min_D_p,D_p(IJK,IM))
		  ENDDO
	       ENDIF

!		Calculate the prefactor of the off-diagonal friction coefficient
!		Use default value of lamdba if there are no particle asparities
	       IF (use_def_lam_HYS)THEN
                  IF (TRIM(UNITS).EQ.'CGS')THEN

		     lam_HYS = 0.0001d0
		     alpha_YS = 1.313d0*LOG10(min_D_p/lam_HYS) - 1.249d0

		  ELSEIF (TRIM(UNITS).EQ.'SI')THEN

		     lam_HYS = 0.000001d0
		     alpha_YS = 1.313d0*LOG10(min_D_p/lam_HYS) - 1.249d0

		  ENDIF
	        ELSE

		   alpha_YS = 1.313d0*LOG10(min_D_p/lam_HYS) - 1.249d0

		ENDIF

!	Calculate velocity components of each species
		USCM_HYS = ZERO
		VSCM_HYS = ZERO	
		WSCM_HYS = ZERO
		
		DO IM = 1, MMAX

	           USCM_HYS = USCM_HYS + EP_S(IJK,Im)*(UGC - AVG_X_E(U_S(IMJK,Im),U_S(IJK,Im),I))
		   VSCM_HYS = VSCM_HYS + EP_S(IJK,Im)*(VGC - AVG_Y_N(V_S(IJMK,Im),V_S(IJK,Im)))
		   WSCM_HYS = WSCM_HYS + EP_S(IJK,Im)*(WGC - AVG_Z_T(W_S(IJKM,Im),W_S(IJK,Im))) 
                   
	        ENDDO 
		
		USCM_HYS = USCM_HYS/phis
		VSCM_HYS = VSCM_HYS/phis
		WSCM_HYS = WSCM_HYS/phis
		
		VREL_poly = SQRT(USCM_HYS**2 + VSCM_HYS**2 + WSCM_HYS**2)
		
		IF (Mu > ZERO) THEN	      

                    RE = D_p_av*VREL_poly*ROP_G(IJK)/Mu

                ELSE

 		    RE = LARGE_NUMBER

                ENDIF

!	Beetstra correction for monodisperse drag
                F_D_BVK = ZERO

		F = 10d0 * phis / EP_g(IJK)**2 + EP_g(IJK)**2 * (ONE+1.5d0*DSQRT(phis))

        	F_D_BVK = F + 0.413d0*Re/(24d0*EP_g(IJK)**2) * (ONE/EP_G(IJK) + 3d0*EP_G(IJK) &
                	*phis + 8.4d0/Re**0.343d0) / (ONE+10**(3d0*phis)/Re**(0.5d0+2*phis))
	
!	YS correction for polydisperse drag
		beta_i_HYS = ZERO
                F_YS = ZERO

		F_YS = 1d0/EP_g(IJK) + (F_D_BVK - 1d0/EP_g(IJK))*(a_YS*Y_i+(1d0-a_YS)*Y_i**2)
		F_YS = F_YS*F_STOKES
		beta_i_HYS = F_YS
              
	
		DO j = 1,MMAX

	   	   IF (j /= M)THEN

		      beta_j_HYS = ZERO
		      Y_i_temp = ZERO
		      FSTOKES_j = ZERO
		      Y_i_temp = D_p(IJK,j)/D_p_av
		      
		      beta_j_HYS = 1d0/EP_g(IJK) + (F_D_BVK - 1d0/EP_g(IJK))*(a_YS*Y_i_temp+(1d0-a_YS)*Y_i_temp**2)
		      FSTOKES_j = 18D0*Mu*EP_s(IJK,j)*EP_g(IJK)/D_p(IJK,j)**2 
		      
		      beta_j_HYS = beta_j_HYS*FSTOKES_j
                      
!	Calculate off-diagonal friction coefficient
                      
		      beta_ij(IJK,M,j) = ZERO
                      
!	This if statement prevents NaN values from appearing for beta_ij
                      IF ((EP_S(IJK,M) == ZERO) .or. (EP_S(IJK,j)==ZERO))THEN
                      
                         beta_ij(IJK,M,j) = ZERO
                     
                      ELSE

		      beta_ij(IJK,M,j) = (2d0*alpha_YS*EP_S(IJK,M)*EP_S(IJK,j))/(EP_S(IJK,M)/beta_i_HYS + EP_S(IJK,j)/beta_j_HYS) 
         
                      ENDIF
                      
                      F_YS = F_YS + beta_ij(IJK,M,j)
	   
	            ENDIF
		ENDDO

		IF (Model_B)THEN

	   	   F_gstmp = F_YS/EP_g(IJK)

	   	ELSE

	           F_gstmp = F_YS

		ENDIF


!***************************END HYS drag relation**************************************


            ELSE
              CALL START_LOG 
              IF(.not.DMP_LOG)call open_pe_log(ier)
              if(mype == pe_io) WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
              CALL END_LOG 
              call mfix_exit(myPE)  
            ENDIF
      
            F_gs(IJK, M) = (ONE - UR_F_gs) * F_gs(IJK, M) + UR_F_gs * F_gstmp
	! Copy above line for beta_ij to keep drag correct for polydisperse case    
	    IF(TRIM(KT_TYPE) == 'GHD') THEN
	      IF(M==1) THEN
	        F_gs(IJK, MMAX) = F_gs(IJK, M)
	      ELSE
	        F_gs(IJK, MMAX) = F_gs(IJK, MMAX) + F_gs(IJK, M)
	      ENDIF
	    ENDIF
         
         ELSE 
            F_gs(IJK, M) = ZERO 
	    IF(TRIM(KT_TYPE) == 'GHD') F_gs(IJK, MMAX) = ZERO
         ENDIF 

      END DO
      
      RETURN  
      END SUBROUTINE DRAG_GS 

!//   Comments on the modifications for DMP version implementation      
!//   001 Include header file and common declarations for parallelization
!//   350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!     
